

import os
import json

CONFIG_PATH = os.environ.get('CONFIG_PATH') or 'config.json'
CONFIG_EXAMPLE = {
    'LOG_PATH': '../logs',
    'OPENAI_API_TYPE': '',
    'OPENAI_API_BASE': '',
    'OPENAI_API_VERSION': '',
    'OPENAI_API_KEY': '',
}


def load_config():
    if os.path.exists(CONFIG_PATH):
        with open(CONFIG_PATH) as f:
            config = json.load(f)
    else:
        config = CONFIG_EXAMPLE.copy()
    # load variable from environment in case for docker/k8s usage 从环境变量中获取

    for key in config.keys():
        if key in os.environ:
            config[key] = os.environ.get(key)
    return config


config = load_config()
