

import os
import time
import logging
from functools import wraps
from datetime import datetime, timezone, timedelta

import redis
import fastapi

from config import config


logger = logging.getLogger(__name__)

DT_FORMAT = '%Y-%m-%d %H:%M:%S'
TIMEZONE = timezone(timedelta(hours=8), 'Asia/Shanghai')
REDIS_POOL = None


class Singleton(type):
    _instances = {}

    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super(Singleton, cls).__call__(*args, **kwargs)
        return cls._instances[cls]


def get_current_datetime_str() -> str:
    """获取当前时间字符串"""
    return datetime.now().astimezone(TIMEZONE).isoformat(" ", "seconds")


def timer(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        t1 = time.time()
        ret = func(*args, **kwargs)
        logger.info(
            f"[{func.__name__}] execution_time: {time.time() - t1}"
        )
        return ret
    return wrapper


def init_logger():

    if not os.path.exists(config['LOG_PATH']):
        os.makedirs(config['LOG_PATH'])
    fileHandler = logging.FileHandler(
        os.path.join(
            config['LOG_PATH'],
            'app.log'
        )
    )
    streamHandler = logging.StreamHandler()
    streamHandler.setLevel(logging.DEBUG)
    log_formatter = logging.Formatter(
        "%(asctime)s [%(processName)s: %(process)d] [%(threadName)s: %(thread)d] [%(levelname)s] %(name)s: %(message)s"
    )
    streamHandler.setFormatter(log_formatter)
    fileHandler.setFormatter(log_formatter)

    logging.basicConfig(
        handlers=[
            fileHandler,
            streamHandler,
        ],
        level=config.get("LOG_LEVEL", logging.INFO),
        force=True
    )


def init_redis_cli() -> redis.ConnectionPool:
    redis_url = config.get('REDIS_URL', 'redis://localhost:6379')
    global REDIS_POOL
    if not REDIS_POOL:
        REDIS_POOL = redis.ConnectionPool.from_url(redis_url)
    return REDIS_POOL


def get_redis_conn():
    """获取redis连接"""
    global REDIS_POOL
    init_redis_cli()
    if REDIS_POOL:
        return redis.Redis(connection_pool=REDIS_POOL)
    raise ConnectionError("Cannot connect to Redis server")
