

import logging
import traceback
import mgi_alphatool
from openai.types.chat import ChatCompletionSystemMessageParam, ChatCompletionUserMessageParam
from llm_utils.utils import get_llm_chat_completion


logger = logging.getLogger(__name__)


GET_PANEL_LIST_PROMPT = """
I will give you a text description of the experimental operation sequence. Please extract all the slot numbers or board numbers involved.
## Rules
1. Only the slot number or board number, no other content
2. The return format is: [1, 2, 3, 4, 5]
3. If there is no slot number or board number, return an empty list
4. Do not repeat the board number
"""


def tips_test(pos: int):
    ctx = mgi_alphatool.init('alphatool')
    tips200_1 = ctx.load_labware('mdk_96_filtertiprack_200ul', pos)
    # Pipettes
    right_pipette = ctx.load_pipette("p200_multi", "left")
    right_pipette.pick_up_tip(tips200_1[0])
    ctx.pause(2)
    right_pipette.drop_tip(tips200_1[0])
    right_pipette.move_to(tips200_1[0])
    json_data = ctx._export()
    return json_data


def get_panel_list_from_txt(txt_content: str) -> list[int]:
    """Use LLM to read the panel list from a txt file"""
    initial_prompt = [
        ChatCompletionSystemMessageParam(
            content=GET_PANEL_LIST_PROMPT,
            role='system',
            name='system',
        ),
        ChatCompletionUserMessageParam(
            content=txt_content,
            role='user',
            name='user',
        )
    ]
    llm_resp = get_llm_chat_completion(initial_prompt)
    if not llm_resp.choices:
        return []
    try:
        content = llm_resp.choices[0].message.content
        if not content:
            return []
        panel_list = eval(content)
    except Exception:
        logger.error(f"get_panel_list_from_txt error: {traceback.format_exc()}")
        return []
    else:
        return panel_list
