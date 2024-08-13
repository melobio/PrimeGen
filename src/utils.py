

import os

def init_logger():
    import logging
    from config import config

    if not os.path.exists(config['LOG_PATH']):
        os.makedirs(config['LOG_PATH'])
    fileHandler = logging.FileHandler(
        os.path.join(
            config['LOG_PATH'],
            'app.log'
        )
    )
    log_formatter = logging.Formatter(
        "%(asctime)s [%(processName)s: %(process)d] [%(threadName)s: %(thread)d] [%(levelname)s] %(name)s: %(message)s"
    )
    fileHandler.setFormatter(log_formatter)

    logging.basicConfig(handlers=[fileHandler], level=logging.INFO)


def init_openai():
    """初始化openai接口设置"""
    global openai
    import openai

    from config import config

    openai.api_key = config['OPENAI_API_KEY']
    openai.api_type = config['OPENAI_API_TYPE']
    openai.api_base = config['OPENAI_API_BASE']
    openai.api_version = config['OPENAI_API_VERSION']


def init_everything():
    init_logger()
    init_openai()
