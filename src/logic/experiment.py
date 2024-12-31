
import json
import time
import logging
import requests

from uuid import uuid4
from dataclasses import dataclass, field, asdict
from typing import Dict, Any, List
from redis import Redis
from fastapi import WebSocket

from config import config
from static import REDIS_KEYS_PREFIX, ActionType
from utils import Singleton, timer, get_current_datetime_str
from logic.alphatool_sdk import AlphaTool, AlphaToolSdkException


logger = logging.getLogger(__name__)



RETRY_TIMES = 3
# polling run status interval, seconds
CHECK_INTERVAL = 1
# Exception machine status
EXCEPTION_STATUS = (
    'failed',
    'blocked-by-open-door'
)
# Mapping of exception types and corresponding handling methods
ISSUE_TYPE_ACTION_MAPPING = {
    'type1': 'redo',
    'type2': 'stop',
}
# redis key
EXPERIMENT_KEY_PREFIX = REDIS_KEYS_PREFIX['experiment']


class MachineStatus:

    IDLE_STATUS = ('idle', )
    PAUSED_STATUS = ('paused', )
    DOING_STATUS = ('running', 'finishing')
    EXCEPTION_STATUS = ('blocked-by-open-door', 'failed')
    STABLE_STATUS = ('stopped', 'succeeded')


class ExperimentException(Exception):
    pass


@dataclass
class Experiment(object):
    """experiment"""

    # conversationid
    conversation_uuid: str
    uuid: str = field(default_factory=lambda: uuid4().hex)
    name: str = 'New Experiment'
    experiment_type: str = 'alphatool'
    create_time: str = field(default_factory=get_current_datetime_str)
    update_time: str = field(default_factory=get_current_datetime_str)
    # alphatool Machine connection information, by default taken from the configuration file
    alphatool_baseurl: str = config['ALPHATOOL_BASEURL']
    # jetson Camera service url address
    jetson_baseurl: str = config['JETSON_BASEURL']

    # protocols
    protocols: List[List[str]] = field(default_factory=list)
    # Board information list
    panel_list: List[List[int]] = field(default_factory=list)
    # List of running processes
    runs: List[Dict] = field(default_factory=list)

    # current protocol
    current_protocol_index: int = 0
    # Currently running task information
    current_run: Dict = field(default_factory=dict)
    status: str = 'idle'

    @property
    def sdk(self) -> AlphaTool:
        return AlphaTool(self.alphatool_baseurl)

    def asdict(self) -> Dict[str, Any]:
        """asdict dataclass Object to Dictionary"""
        return asdict(self)

    def get_current_run(self) -> str:
        """Get the currently running program"""
        resp = self.sdk.get_current_run()
        current_run_id = resp['data']['current']['id']
        return current_run_id


class ExperimentManager(metaclass=Singleton):
    """
    Run the experiment
1. Check the machine status: If there is a run running, stop the current run
- Perform the pre-check task **The number of protocols checked is the same as the number of protocols in the experiment, both specified by current_protocol_index
2. Loop through the experimental steps of protocols
1. Create a Protocol and save the protocol information, protocol_id -- protocol_info
2. Create a Run and save the run information run_id -- run_info
1. While True: Run the action of the Run -- current_action_index
2. Poll to check the machine status:
Current action index: currentCommandIndex,
Current run status: status:
idle Created, not executed
running Executing
pause-requested Requesting to pause (not paused yet)
paused Paused
stop-requested Requesting to stop (not stopped yet)
stopped Stopped
failed Failed
finishing Finishing
succeeded Successful
blocked-by-open-door Pause (because the door is open)
3. After checking that the action is in the paused state, send the current state information to the detection service
4. Process according to the return comparison table of abnormal state and processing method
5. If it still fails after processing retry_times, exit and report an error, otherwise enter the next action
    """

    def __init__(self, redis_conn: Redis) -> None:
        self.redis_conn = redis_conn
        self.ws_connections: Dict[str, WebSocket] = {}
        self.conversations = {}
        self.load_experiments()

    def load_experiments(self):
        """Read all experimental data from redis"""
        experiment_uuids = self.redis_conn.keys(
            EXPERIMENT_KEY_PREFIX + '*'
        )  # type: Any
        for e_uuid in experiment_uuids:
            experiment_info = self.redis_conn.get(e_uuid)  # type: Any
            experiment = Experiment(**json.loads(experiment_info))
            self.conversations[experiment.conversation_uuid] = experiment

    def save_experiment(self, experiment: Experiment):
        """Save experiment to redis"""
        experiment.update_time = get_current_datetime_str()
        self.redis_conn.set(
            EXPERIMENT_KEY_PREFIX + experiment.uuid,
            json.dumps(experiment.asdict()),
        )

    def set_ws_conn(self, conversation_uuid: str, ws: WebSocket):
        """Set websocket connection with conversation id"""
        self.ws_connections[conversation_uuid] = ws

    def get_ws_conn(self, conversation_uuid: str):
        """Get websocket connection with conversation id"""
        return self.ws_connections.get(conversation_uuid)

    def load_experiment(self, experiment_uuid: str):
        """Read redis experimental data->experiment object"""
        experiment_info = self.redis_conn.get(
            EXPERIMENT_KEY_PREFIX + experiment_uuid
        )  # type: Any
        experiment = Experiment(**json.loads(experiment_info))
        return experiment

    def get_experiment_by_conversation(self, conversation_uuid: str):
        """Get experimental subjects based on conversation"""
        experiment = self.conversations.get(conversation_uuid)
        if not experiment:
            experiment = Experiment(conversation_uuid=conversation_uuid)
            self.save_experiment(experiment)
            self.conversations[conversation_uuid] = experiment
        return experiment

    def get_current_run_status(self, experiment: Experiment, run_id: str):
        """Get the current run status and save it"""
        run_info = experiment.sdk.get_run(run_id=run_id)['data']
        experiment.current_run = run_info
        experiment.status = run_info['status']
        self.save_experiment(experiment)
        return run_info

    def wait_until_done(
            self,
            experiment: Experiment,
            run_id: str | None = None,
            check_interval: float = CHECK_INTERVAL
    ):
        """Waiting for the action to complete"""
        if not run_id:
            run_id = experiment.get_current_run()
        status = self.get_current_run_status(experiment, run_id)['status']
        logger.info(f"[wait_until_done] current_status: {status}")
        while status in MachineStatus.DOING_STATUS:

            logger.debug(f"[wait_until_done] sleep for {CHECK_INTERVAL}")
            time.sleep(check_interval)

            status = self.get_current_run_status(
                experiment, run_id
            )['status']
        logger.debug(f"[wait_until_done] current_status: {status}")
        return status

    def try_stop_current_run(self, experiment: Experiment):
        """Try to stop the current task"""
        resp = experiment.sdk.get_current_run()
        current = resp['data'].get('current')
        if current:
            if current['status'] in MachineStatus.IDLE_STATUS + MachineStatus.STABLE_STATUS + MachineStatus.EXCEPTION_STATUS:
                return True
            result = experiment.sdk.run_action(current['id'], ActionType.stop)
            experiment.current_run = {}
            self.save_experiment(experiment)
            return result.get('success')
        experiment.current_run = {}
        self.save_experiment(experiment)
        return True

    async def send_ws_data(
            self,
            ws: WebSocket,
            code: int,
            msg: str,
            experiment: Experiment,
            data_type: str,
            generating: bool = False,
            command_options: List[str] = [],
            state: str = 'continue',
    ):
        """Send data to websocket: text or image data"""
        experiment_info = experiment.asdict()
        # protocols
        experiment_info.pop('protocols')
        resp = {
            'responses': {
                'code': code,
                'response': '',
                'data': {
                    'content': msg,
                    # text or image
                    'type': data_type,
                    'experiment_info': experiment_info,

                    'generating': generating,
                    'command_options': command_options,
                    'state': state,
                }
            }
        }
        # logger.info(f"[send_ws_data]: {resp}")
        await ws.send_json(resp)

    async def start_experiment(
            self,
            ws: WebSocket,
            experiment: Experiment,
            retry_times: int = RETRY_TIMES,
    ):
        """Start the experiment from the protocol_index of the experiment record"""

        logger.info(f"try_stop_current_run: {experiment.uuid}")
        self.try_stop_current_run(experiment)


        start_protocol_index = experiment.current_protocol_index
        await self.send_ws_data(
            ws,
            200,
            "[STEP] Start Experiment",
            experiment,
            data_type="text",
            generating=True,
        )

        logger.info(f"ProtocolsCounts: {len(experiment.protocols)}")
        for _i, protocols in enumerate(experiment.protocols):
            if _i < start_protocol_index:
                logger.info(f"The experiment is currently running: {start_protocol_index}, Skip the current {_i} protocol")
                continue
            logger.info(f"StartExperiment: {experiment.uuid}: {_i} protocol")
            experiment.current_protocol_index = _i
            self.save_experiment(experiment)
            for _j, protocol in enumerate(protocols):
                # TODO: The machine returns to the origin
                self.wait_until_done(experiment)
                # The machine returns the status "Stopped". The execution will still prompt "Running" so wait for 1 second.
                # time.sleep(1)
                experiment.sdk.move_position('8', 0, 0, 0)

                # TODO: Board layout identification
                await self.send_ws_data(
                    ws,
                    200,
                    "[STEP] Start board layout recognition Board Recognition",
                    experiment,
                    data_type="text",
                    generating=True,
                )
                result = self.check_machine_status('BOARD_RECOGNITION')
                if result.get('error'):
                    raise ExperimentException(

                    )

                # TODO: Pipette inspection
                await self.send_ws_data(
                    ws,
                    200,
                    "[STEP] Start pipette check Liquid Movement",
                    experiment,
                    data_type="text",
                    generating=True,
                )
                result = self.check_machine_status('LIQUID_MOVEMENT')
                if result.get('error'):
                    raise ExperimentException(
                        f"Pipette Check Failed: {result.get('message')}"
                    )

                # Run protocol
                t1 = time.time()
                try:
                    # Pause on error
                    status = await self.run_protocol(
                        ws=ws,
                        experiment=experiment,
                        protocol=protocol,
                        retry_times=retry_times,
                    )
                except ExperimentException as e:
                    await self.send_ws_data(
                        ws,
                        500,
                        f"[ERROR] Problems encountered in the experiment: {e}",
                        experiment,
                        data_type="text",
                        generating=False,
                        command_options=['retry', 'resume', 'stop']
                    )

                dt1 = time.time() - t1


                logger.info(f"RunExperimentProtocol: {experiment.uuid}-{_i}-{_j} in {dt1} seconds")
                await self.send_ws_data(
                    ws,
                    200,
                    f"[INFO] Running Protocol {_j + 1} cost time: {dt1}",
                    experiment,
                    data_type="text",
                    generating=False,
                    command_options=['start']
                )

            logger.info(f"RunExperimentProtocol: [INFO] Protocol {_i + 1} over")
            await self.send_ws_data(
                ws,
                200,
                f"[INFO] Protocol {_i + 1} over",
                experiment,
                data_type="text",
                generating=False,
                command_options=['start']
            )

        logger.info("RunExperiment: Experiment execution completed")
        await self.send_ws_data(
            ws,
            200,
            "[INFO] The experiment has ended",
            experiment,
            data_type="text",
            generating=False,
            state='stop',
        )

    @timer
    async def run_protocol(
        self,
        ws: WebSocket,
        experiment: Experiment,
        protocol: str,
        retry_times: int = 3
    ):
        """Run Protocol"""

        protocol_resp = experiment.sdk.create_protocol_from_str(
            json_string=protocol,
            debug=True,
        )
        # logger.debug(
        #     f"[run_protocol]: protocol_resp: {protocol_resp}"
        # )
        protocol_id = protocol_resp['data']['id']
        await self.send_ws_data(
            ws,
            200,
            f"[STEP] Create Protocol: {protocol_id}",
            experiment,
            data_type="text",
            generating=True,
        )

        run_resp = experiment.sdk.create_run(
            protocol_id=protocol_id,
        )
        logger.debug(
            f"[run_protocol]: run_resp: {run_resp}"
        )
        run_id = run_resp['data']['id']
        await self.send_ws_data(
            ws,
            200,
            f"[STEP] Create Run: {run_id}",
            experiment,
            data_type="text",
            generating=True,
        )

        status = await self.start_run(
            ws=ws,
            experiment=experiment,
            run_id=run_id,
            retry_times=retry_times,
        )
        return status

    @timer
    async def start_run(
        self,
        ws: WebSocket,
        experiment: Experiment,
        run_id: str,
        retry_times: int,
    ):
        """Run experiment steps"""
        # Initialize status
        run_info = self.get_current_run_status(
            experiment, run_id
        )
        status = run_info['status']
        # Loop status until stop/exception status
        while status not in MachineStatus.STABLE_STATUS + MachineStatus.EXCEPTION_STATUS:
            # Run next step
            experiment.sdk.run_action(
                run_id=run_id,
                action_type=ActionType.play,
            )

            run_info = self.get_current_run_status(
                experiment, run_id
            )
            status = run_info['status']

            logger.info(f"[start_run] current_step: {run_info['currentCommandIndex']}")
            logger.info(f"[start_run] current_status: {status}")
            await self.send_ws_data(
                ws,
                200,
                f"[INFO] Running RUN[{run_id}] Step: {run_info['currentCommandIndex']}, Status: {status}",
                experiment,
                data_type="text",
                generating=True,
            )
            status = self.wait_until_done(experiment, run_id)
            logger.info(f"[start_run] current_status: {status}")

            if status in MachineStatus.PAUSED_STATUS:
                # action.end
                # TODO: Wait for instruction each time paused state is reached
                # return status
                # TODO: Paused state: Upload current machine status and action type & handle machine exceptions
                check_result = self.check_machine_status(run_info)
                # TODO: Check if params.isMixEnd is true for current step # Determine if current or previous step
                # 1. Get params.labwareId
                # 2. Take photo - wait for completion
                # 3. Return params.labwareId
                logger.debug(
                    f"[start_run] check_machine_status: {check_result}"
                )
                await self.send_ws_data(
                    ws,
                    200,
                    f"[STEP] JETSON check RUN[{run_id}] status: {check_result}",
                    experiment,
                    data_type="text",  # TODO: maybe image
                    generating=True,
                )
                # TODO: Determine status requiring action
                if check_result.get('status') == 'error':
                    # Exception status
                    handle_result = False
                    for _i in range(retry_times):
                        logger.debug(
                            f"[start_run] check_machine_status for the {_i + 1} time: {check_result}"
                        )                        
                        handle_result = self.handle_machine_status(check_result)
                        # TODO: Handle exception
                        await self.send_ws_data(
                            ws,
                            200,
                            f"[STEP] Handle Run[{run_id}] exception result: {handle_result}",
                            experiment,
                            data_type="text",
                            generating=True,
                        )
                        if handle_result:
                            break
                    # Handle failed
                    if not handle_result:
                        raise ExperimentException(
                            f"Experiment failed with exception: {check_result}, and handle_result is {handle_result}"
                        )
            # action.start
            check_result = self.check_machine_status(run_info)
            # status = self.wait_until_done(experiment, run_id)

            # await self.send_ws_data(
            #     ws,
            #     200,
            #     f"[STEP] JETSON检查RUN[{run_id}]状态: {check_result}",
            #     experiment,
            #     data_type="text",
            #     generating=True,
            # )

        if status in MachineStatus.STABLE_STATUS:
            # status == 'succeeded' | 'stopped'
            logger.debug(f"[start_run] Current task reached terminal status: {run_id}: {status}")
            await self.send_ws_data(
                ws,
                200,
                f"[INFO] Current task executed successfully: {run_id}",
                experiment,
                data_type="text",
                generating=False,
            )

        if status in MachineStatus.EXCEPTION_STATUS:
            raise ExperimentException(
                f"[start_run] experiment failed with status: {status}")

        return status

    @timer
    async def redo_action(
        self,
        ws_conn: WebSocket,
        experiment: Experiment,
    ):
        """Redo current operation"""
        current_run_id = experiment.current_run['id']
        resp = experiment.sdk.run_redo(current_run_id)
        logger.debug(f"[redo_action] Redo current operation: {current_run_id}, {resp}")
        await self.send_ws_data(
            ws_conn,
            200,
            f"[INFO] Redo current operation: {current_run_id}",
            experiment,
            data_type="text",
            generating=True,
        )
        status = self.wait_until_done(experiment, current_run_id)
        await self.send_ws_data(
            ws_conn,
            200,
            f"[INFO] Redo current operation: {current_run_id} result: {status}",
            experiment,
            data_type="text",
            generating=False,
        )
        return status

    @timer
    async def resume_run(
        self,
        ws: WebSocket,
        experiment: Experiment,
        retry_times: int = 3,
    ):
        """Resume experiment"""
        current_run_id = experiment.current_run['id']
        await self.send_ws_data(
            ws,
            200,
            f"[INFO] Resume running task: {current_run_id}",
            experiment,
            data_type="text",
            generating=True,
        )
        status = await self.start_run(
            ws=ws,
            experiment=experiment,
            run_id=current_run_id,
            retry_times=retry_times
        )
        if status == MachineStatus.STABLE_STATUS[1]:
            # status == 'succeeded'
            logger.debug(f"[resume_run] Current task executed successfully: {current_run_id}")
            await self.send_ws_data(
                ws,
                200,
                f"[INFO] Current task executed successfully: {current_run_id}: {status}",
                experiment,
                data_type="text",
                generating=False,
                command_options=['start', 'stop']
            )
            # After current run succeeds, prepare to execute next protocol
            experiment.current_protocol_index += 1
            experiment.status = status
            await self.start_experiment(
                ws=ws,
                experiment=experiment,
            )
        else:
            logger.debug(f"[resume_run] Current task status: {current_run_id}: {status}")
            await self.send_ws_data(
                ws,
                200,
                f"[INFO] Current task did not succeed: {current_run_id}: {status}",
                experiment,
                data_type="text",
                generating=False,
                command_options=['retry', 'stop']
            )

    @timer
    def check_machine_status(self, type):
        """Call video detection service, return current exception status"""
        # TODO: Call interface
        
        return {}

    @timer
    def handle_machine_status(self, machine_status):
        """Take corresponding action based on camera return status"""
        # TODO: Corresponding logic
        # Successfully handled
        return True
