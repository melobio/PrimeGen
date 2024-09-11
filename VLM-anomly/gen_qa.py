from mimetypes import guess_type
import openai 
import json
import os
from tqdm import tqdm
import pandas as pd

#Similar to 3D-LLM instruction set construction, the user needs to input the OpenAI API and configure the engine, etc.
def GPT4o(prompts):

    openai.api_type = ""
    openai.api_base = ""
    openai.api_version = ""
    openai.api_key = ""

    response = openai.ChatCompletion.create(
        engine="gpt-4o", 
        response_format={"type": "json_object"},
        messages=[
            {"role": "system",
                 "content": "You're an editorial expert"},
            {"role": "user", 
                "content": f"{prompt}"}],
        temperature=0.7,
        top_p=0.95,
        max_tokens=4096
    )
    # 生成结果
    output = (response.choices[0].message.content)
    return output


output_path = "./Global_img/"
if not os.path.exists(output_path):
    os.makedirs(output_path)

with open('../Data_process/annotations_data.jsonl','r') as f:
    data = [json.loads(line) for line in f]

augment_data_dict = {"Q":[],"A":[],"img_path":[]}

#["Tip1","Tip8","Damaged","Open","Close","empty","full"]
for i in tqdm(data):
    label = i['labels']
    img_path =  i['img_path']
    print("true label: ",label)
    if "Tip8" in label:
        tip = 'The image shows 8 pipette tip being securely held by the OT-2 Gen-2 robotic arm gripper.[Holding/Released],Holding.'
        channel = "Eight pipette tips are being firmly grasped by OT-2's robotic arm gripper, OT-2 Gen-2.[One/Eight/None],Eight"

    if "Tip1" in label:
        tip = 'The image shows a pipette tip being securely held by the OT-2 Gen-2 robotic arm gripper.[Holding/Not Released],Holding.'
        channel = "OT-2's robotic arm gripper, OT-2 Gen-2, holds only one pipette tip, not eight.[One/Eight/None],One."


    if 'Damaged' in label:
        damage = 'The underside of this pipette has damaged. This damaged pipette may affect the accuracy of liquid transfer during experimental operation, so it needs to be replaced promptly.[Undamaged/Damaged/None], Damaged.'
    else:
        damage = 'This pipette tip is undamaged.[Undamaged/Damaged/None], Undamaged.'
         
    if "Tip8" not in label and "Tip1" not in label and "Damaged" not in label:
        tip = 'The OT-2 Gen-2 robotic arm gripper is holding any pipette tip.[Holding/Released],Released.'
        channel = "Since the OT-2 Gen-2 robotic arm gripper does not grab the pipette tip, it is not necessary to judge the pipette tip channel,[One/Eight/None],None."
        damage = "Since the OT-2 Gen-2 robotic arm gripper does not grab the pipette tip, it is not necessary to judge the pipette tip condition,[Undamaged/Damaged/None],None."
    
    if "Close" in label:
        pcr = "In the current image, the heated lid of the thermocycler module is in the closed state.[Open/Closed], Closed."
    
    if "Open" in label:
        pcr = "In the current image, the heated lid of the thermocycler module is open, so we can clearly see the 96-Well Plate inside.[Open/Closed], Open."
    
    if 'empty' in label:
        trash = "In the current picture, the trash can is empty.[Empty/Full], Empty."
    
    if 'full' in label:
        trash = "In the current image, the trash can is full and contains many discarded pipette tips. [Empty/Full], Full."
    if 'full' not in label and 'empty' not in label:
        trash = "In the current picture, the trash can is empty.[Empty/Full], Empty."

    prompt = f"""
In this image, we see a OT-2 experimental platform.
In the picture Description:
<1>{tip};
<2>{damage};
<3>{pcr};
<4>{channel};
<5>{trash}
##
Help me generate a question and answer based on the above picture description. 
The question require models that can be well trained in cognition, understanding, and detection.
Note that every question must be included all states, for example, if there is no pipette, it is Released, and the quantity is None.

##
Here is an example of the question and answer format：
"Q": "Can the OT-2 Gen-2 robotic arm gripper hold the pipette tips?What is the status of the pipette tips?How many pipette tips can the robotic arm gripper hold?Is the lid of the Thermalcycler module open or closed?Is the waste container empty or full?[Holding/Released , Undamaged/Damaged/None, Open/Closed, One/Eight/None, Empty/Full]", \
"A": "[Holding , Undamaged , Open , One, Empty]"

## 
The output format is json: {{"Q": " ", "A": " "}}
"""


    output = GPT4o(prompt)
    data_dict = json.loads(output)
    
    Q = data_dict["Q"]
    A = data_dict["A"]
    print('[generated Q]: ',Q)
    print('[generated Q]: ',A)

    csv_file_path = 'Global_VQA.csv'
    file_exists = os.path.exists(csv_file_path)
    if not file_exists:
        pd.DataFrame({"Q": [Q], "A": [A], "img_path": [img_path]}).to_csv(csv_file_path, mode='a', header=True, index=False)
    else:
        pd.DataFrame({"Q": [Q], "A": [A], "img_path": [img_path]}).to_csv(csv_file_path, mode='a', header=False, index=False)
   
