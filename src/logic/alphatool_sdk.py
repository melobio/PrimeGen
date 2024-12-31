import logging
from functools import wraps
from urllib.parse import urljoin
from typing import Optional, Callable, Dict

import requests
import requests.status_codes

from static import ActionType


logger = logging.getLogger(__name__)

TIMEOUT = 3


class AlphaToolSdkException(Exception):
    pass


def parse_response(func) -> Callable[..., Dict]:
    """response handle"""
    @wraps(func)
    def wrapper(*args, **kwargs):
        response: requests.Response = func(*args, **kwargs)

        if response.status_code not in range(200, 300):
            logger.debug(
                f"AlphaTool request Exception: \n"
                f"request: {response.request.__dict__} \n"
                f"response: {response.__dict__}"
            )
            raise AlphaToolSdkException(
                f"AlphaTool response status code error: \n"
                f"response: [{response.status_code}] {response.content}"
            )
        try:
            alphatool_response = response.json()
        except Exception:
            raise AlphaToolSdkException(
                f"AlphaTool response json serialize failed: {response.content}"
            )
        if not alphatool_response.get('success'):
            raise AlphaToolSdkException(
                f"AlphaTool api get an error: {alphatool_response}"
            )
        else:
            return alphatool_response

    return wrapper


class AlphaTool(object):

    def __init__(self, base_url: str, timeout: int = TIMEOUT) -> None:
        self.request_session = requests.Session()
        self.base_url = base_url
        self.timeout = timeout

    @parse_response
    def create_protocol(
        self,
        files,
        debug=True,
    ) -> requests.Response:
        """Create a protocol"""
        path = '/protocols'
        debug = 'true' if debug else 'false'
        req_body = {
            'debug': debug,
        }
        return self.request_session.post(
            url=urljoin(self.base_url, path),
            files=files,
            data=req_body,
        )

    @parse_response
    def create_protocol_from_str(
        self,
        json_string: str,
        debug=True,
    ) -> requests.Response:
        """Create a protocol based on the JSON string of the protocol"""
        from io import StringIO
        s_io = StringIO(json_string)
        path = '/protocols'
        debug = 'true' if debug else 'false'
        req_body = {
            'debug': debug,
        }
        return self.request_session.post(
            url=urljoin(self.base_url, path),
            files={
                'files': ('pre_check.json', s_io),
            },
            data=req_body,
        )

    @parse_response
    def create_run(
        self,
        protocol_id: str,
        server_project_id: Optional[str] = None,
        server_protocol_id: Optional[str] = None,
        server_url: Optional[str] = None,
        server_token: Optional[str] = None,
        server_device_id: Optional[str] = None,
    ) -> requests.Response:
        """Create an execution task
        reference: https://apifox.com/apidoc/shared-0dd65fe0-ecd3-4427-a19c-323b9503e1ac/api-192817161
        """
        path = '/runs'
        req_body = {
            'protocolId': protocol_id,
        }
        if server_project_id:
            req_body['serverProjectId'] = server_project_id
        if server_protocol_id:
            req_body['serverProtocolId'] = server_protocol_id
        if server_url:
            req_body['serverUrl'] = server_url
        if server_token:
            req_body['serverToken'] = server_token
        if server_device_id:
            req_body['serverDeviceId'] = server_device_id
        return self.request_session.post(
            urljoin(self.base_url, path),
            json={
                'data': req_body
            }
        )

    @parse_response
    def get_run(
        self,
        run_id: str,
    ) -> requests.Response:
        """Get current task information: status, steps, etc.
        reference: https://apifox.com/apidoc/shared-0dd65fe0-ecd3-4427-a19c-323b9503e1ac/api-192817163
        """
        path = f'/runs/{run_id}'
        return self.request_session.get(
            url=urljoin(self.base_url, path),
            timeout=self.timeout
        )

    @parse_response
    def get_current_run(self):
        """Get the current run"""
        path = '/runs/current'
        return self.request_session.get(
            url=urljoin(self.base_url, path),
            timeout=TIMEOUT
        )

    @parse_response
    def run_action(
        self,
        run_id: str,
        action_type: ActionType,
    ) -> requests.Response:
        """Control Tasks
        reference: https://apifox.com/apidoc/shared-0dd65fe0-ecd3-4427-a19c-323b9503e1ac/api-192817164
        """
        path = f'/runs/{run_id}/actions'
        req_body = {
            'data': {
                'actionType': action_type.value,
            }
        }
        return self.request_session.post(
            url=urljoin(self.base_url, path),
            json=req_body,
        )

    @parse_response
    def run_redo(
        self,
        run_id: str,
        count: int = 1,
    ) -> requests.Response:
        """Redo last step
        reference: https://apifox.com/apidoc/shared-0dd65fe0-ecd3-4427-a19c-323b9503e1ac/api-222332955
        """
        path = f"/runs/{run_id}/redo"
        return self.request_session.post(
            url=urljoin(self.base_url, path),
            json={
                'count': count
            }
        )

    @parse_response
    def move_position(self, channel: str, x: float, y: float, z: float):
        """
        Move the gun tip to a certain position
        reference: https://apifox.com/apidoc/shared-0dd65fe0-ecd3-4427-a19c-323b9503e1ac/api-192817154
        """
        path = f"/pipette/{channel}/move-xyz/{x}/{y}/{z}"
        return self.request_session.post(
            url=urljoin(self.base_url, path)
        )

    @parse_response
    def capture_labware(
        self,
        run_id: str,
        labware_id: str,
        to_location: tuple[int, int, int] = (0, 322, 40)
    ):
        """
        Use the gripper to grab the Labware to somewhere
        reference: https://apifox.com/apidoc/shared-0dd65fe0-ecd3-4427-a19c-323b9503e1ac/api-237932740
        """
        path = f"/runs/{run_id}/captureLabwareTo"
        return self.request_session.post(
            url=urljoin(self.base_url, path),
            headers={'ig-run': 'true'},
            json={
                'labwareId': labware_id,
                'toLocation': {
                    "x": to_location[0],
                    "y": to_location[1],
                    "z": to_location[2]
                }
            }
        )

    @parse_response
    def release_labware(self, run_id: str, labware_id: str):
        """
        Put the labware back in place
        reference: https://apifox.com/apidoc/shared-0dd65fe0-ecd3-4427-a19c-323b9503e1ac/api-237934499
        """
        path = f"/runs/{run_id}/releaseLabware"
        return self.request_session.post(
            url=urljoin(self.base_url, path),
            headers={'ig-run': 'true'},
            json={
                'labwareId': labware_id,
            }
        )
