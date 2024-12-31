
import openai
from openai.types.chat.chat_completion import ChatCompletion
from openai.types.chat.chat_completion_message_param import ChatCompletionMessageParam
from config import config
from openai.lib.azure import AzureOpenAI
from enum import EnumType
from typing import Any, Optional, List


class LLMType(EnumType):
    AZURE = 'AZURE'
    QWEN = 'QWEN'
    ZHIPU = 'ZHIPU'


def get_llm_client(
    api_key: str,
    base_url: str,
    llm_type: LLMType = LLMType.AZURE, # type: ignore
    **kwargs,
) -> Any:
    """
    Get llm client according to llm_type and initial client
    """
    if llm_type == LLMType.AZURE:
        client = AzureOpenAI(api_key=api_key, azure_endpoint=base_url, **kwargs)
    elif llm_type == LLMType.QWEN:
        # compatible with openai
        client = openai.Client(api_key=api_key, base_url=base_url, **kwargs)
        # or dashscope sdk
    elif llm_type == LLMType.ZHIPU:
        client = openai.Client(api_key=api_key, base_url=base_url, **kwargs)
        # or zhipuai sdk
        # import zhipuai
        # client = zhipuai.ZhipuAI(api_key=api_key, base_url=base_url, **kwargs)
    else:
        raise NotImplementedError(f"Unsupported LLM type: {llm_type}")
    return client


def get_llm_from_config(config: dict, llm_type: LLMType = LLMType.AZURE) -> Any: # type: ignore
    """

    Get llm client from config
    """
    llm_type = config.get('DEFAULT_LLM') or llm_type
    llm_type = str(llm_type)  # type: ignore
    client_kwargs = {}
    for k, v in config.items():
        #
        if k.startswith(f"{llm_type}_"):
            key = k.replace(f"{llm_type}_", '')
            #
            client_kwargs[key.lower()] = v
    client_kwargs.pop('model', None)
    return get_llm_client(llm_type=llm_type, **client_kwargs)


def get_llm_chat_completion(
        messages: List[ChatCompletionMessageParam],
        client: Optional[AzureOpenAI | openai.Client] = None,
        temperature=0.2,
        model=None,
        **kwargs,
) -> ChatCompletion:
    llm_type = config.get('DEFAULT_LLM') or LLMType.AZURE
    if client is None:
        client = get_llm_from_config(config, llm_type)
    if not client:
        raise ValueError("LLM client is None, please check")
    if not model:
        model = config[f"{llm_type}_MODEL"]
    print('****************model*********************',model)
    response: ChatCompletion = client.chat.completions.create(
        model=model,
        messages=messages,
        temperature=temperature,
        **kwargs,
    )
    return response
