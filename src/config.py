

import os
import json


DEFAULT_CONFIG_PATH = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    'config.json',
)

CONFIG_PATH = os.environ.get('CONFIG_PATH') or DEFAULT_CONFIG_PATH
CONFIG_EXAMPLE = {
    'LOG_LEVEL': 'INFO',
    'LOG_PATH': '../logs',
    "ALPHATOOL_BASEURL": "",
    "JETSON_BASEURL": "",
    "REDIS_URL": "redis://localhost:6379",
    "DEFAULT_LLM": "azure",
    "AZURE_API_KEY": "",
    "AZURE_BASE_URL": "",
    "AZURE_MODEL": "gpt-4o",
    "AZURE_API_VERSION": "2024-02-15-preview",
    "QWEN_API_KEY": "",
    "QWEN_BASE_URL": "https://dashscope.aliyuncs.com/compatible-mode/v1",
    "QWEN_MODEL": "qwen-plus",
    "ZHIPU_API_KEY": "",
    "ZHIPU_BASE_URL": "https://open.bigmodel.cn/api/paas/v4/",
    "ZHIPU_MODEL": "glm-4-plus"
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
            config[key] = os.environ.get(key, '')
    return config


config = load_config()
