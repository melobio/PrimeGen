

from enum import Enum


REDIS_KEYS_PREFIX = {
    'conversation': 'CONVERSATION_',
    'experiment': 'EXP_',
}


class ActionType(Enum):
    play = 'play'
    pause = 'pause'
    stop = 'stop'
