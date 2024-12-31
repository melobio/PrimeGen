# server.py
import asyncio
import websockets
import openai
import os
from fastapi import FastAPI,WebSocket
import uvicorn
from pydantic import BaseModel
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
import pandas as pd
import json
import traceback
from tools import get_protocol_param_dict,calculating_solution_volume,get_protocol_process
import requests
import logging
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)
from llm_utils.utils import get_llm_chat_completion


all_protocol_name = ["ATOPlex DNA Library Prep Kit",
"ATOPlex RNA Library Prep Kit","ATOPlex HIV Library Prep Kit","ATOPlex MPXV Library Prep Kit","MGIEasy Fast Enzymatic DNA Library Prep Kit","MGIEasy Fast PCR-FREE Enzymatic Library Prep Kit","MGIEasy Universal DNA Library Prep Kit","MGIEasy Enzymatic DNA Library Prep Kit","MGIEasy Cell-free DNA Library Prep Kit","MGIEasy Exome Enzymatic Library Prep Kit","MGICare Chromosomal Copy Number Variation Detection Kit","MGIEasy Fast RNA Library Prep Kit","MGIEasy Directional RNA Library Prep Kit","MGICare Single Cell Chromosomal Copy Number Variation Detection Kit"]

all_protocol_file_dict = {"ATOPlex DNA Library Prep Kit":'H-T-126ATOPlex_DNA建库试剂盒套装说明书中文受控版24_4_9.txt',
"ATOPlex RNA Library Prep Kit":'H-T-127ATOPlex_RNA建库试剂盒套装说明书中文受控版24_4_9.txt',
"ATOPlex HIV Library Prep Kit":'ATOPlexHIV-1建库试剂盒套装使用说明书_中文_RUO_SZ24_3_13.txt',                          
"ATOPlex MPXV Library Prep Kit":'ATOPlex_MPXV引物池V2_0使用说明书-非受控版24_3_13.txt',
"MGIEasy Fast Enzymatic DNA Library Prep Kit":'940-001193-00_940-001194-00_940-001196-00_MGIEasy_Fast酶切⽂库制备试剂套装使用说明书3_0.txt',
"MGIEasy Fast PCR-FREE Enzymatic Library Prep Kit":'940-000886-0_940-000884-00_940-000882-00_MGIEasy_Fast_PCR-FREE酶切文库制备试剂套装使用说明书1_0.txt',
"MGIEasy Enzymatic DNA Library Prep Kit":'1000006987_1000006988_1000017572-MGIEasy酶切DNA文库制备试剂套装使用说明书B6.txt',
"MGIEasy Cell-free DNA Library Prep Kit":'940-000184-0_940-000185-00_MGIEasy_游离DNA文库制备试剂套装使用说明书1_0.txt',
"MGIEasy Exome Enzymatic Library Prep Kit":'1000009658-MGIEasy外显子组酶切文库制备试剂套装使用说明书-B3.txt',
"MGICare Chromosomal Copy Number Variation Detection Kit":'',
"MGIEasy Fast RNA Library Prep Kit":'',
"MGIEasy Directional RNA Library Prep Kit":'',
"MGICare Single Cell Chromosomal Copy Number Variation Detection Kit":''}


def modify_code_file(code_file_path,param_dict):
    # Output string
    new_str = []
    with open(code_file_path,'r') as f:
        lines = f.readlines()
        for ll in lines:
            flag = False
            temp_p = ''
            for param in param_dict.keys():
                if param in ll:
                    flag = True
                    temp_p = param
            if flag:
                tt = str(param_dict[temp_p])
                new_str.append(ll.replace(temp_p, tt))
            else:
                new_str.append(ll)
        f.close()
    dir_path = '/reference_data/protocol_generate_code'
    file_name = os.path.basename(code_file_path).split('.')[0]
    new_name = file_name+'_modify.py'
    new_path = os.path.join(dir_path,new_name)
    with open(new_path, 'w') as f:
        for line in new_str:
            f.write(line)
        f.close()

    return new_path

# Get parameter variables from code
# Get parameter variables from code 
def get_param_name(code_name_list):
    template_code_path = '/reference_data/protocol_block_code_alphatool'
    param_name_list = []
    user_input_list = []
    for name in code_name_list:
        if '*' in name:
            temp_name = name.split('**')
            for temp_n in temp_name:
                temp_path = os.path.join(template_code_path,temp_n)
                with open(temp_path,'r') as f:
                    lines = f.readlines()
                    for line in lines:
                        if '=' in line and 'input_' in line and '/' not in line :
                            param_name_list.append(line.split('=')[-1].strip())
                        if  '=' in line and 'user_' in line:
                            user_input_list.append(line.split('=')[-1].strip())
                    f.close()
        else:
            temp_path = os.path.join(template_code_path,name)
            with open(temp_path,'r') as f:
                lines = f.readlines()
                for line in lines:
                    if '=' in line and 'input_' in line and '/' not in line:
                        param_name_list.append(line.split('=')[-1].strip())
                    if '=' in line and 'user_' in line:
                            user_input_list.append(line.split('=')[-1].strip())
                f.close()
    param_list = []
    for pa in param_name_list:
        if pa not in param_list:
            param_list.append(pa)
    return param_list,user_input_list


def template_protocol_design(instruction: dict, stage=1):

    # Initialize protocol design type
    protocol_design_type = 'template_protocol_design'

    all_protocol_param_dict={"ATOPlex DNA Library Prep Kit":{"sample_nums":96,
                        "spike_in_control_pcr_block":False,
                        "amplicon_multiplicity":50,
                        "panel_amp_length":300,
                         "barcode_type":'single_pcr_barcode'},      
                      "ATOPlex RNA Library Prep Kit":{"sample_nums":96,
                        "spike_in_control_pcr_block":False,
                        "amplicon_multiplicity":50,
                        "panel_amp_length":300,
                         "barcode_type":'single_pcr_barcode'},       
                      "ATOPlex HIV Library Prep Kit":{"sample_nums":96,},
                      "ATOPlex MPXV Library Prep Kit":{"sample_nums":96,},
                      "MGIEasy Exome Enzymatic Library Prep Kit":{"sample_nums":96,},
                      "MGIEasy Cell-free DNA Library Prep Kit":{"sample_nums":96,},
                      "MGICare Single Cell Chromosomal Copy Number Variation Detection Kit":{"sample_nums":96,},
                      "MGIEasy Universal DNA Library Prep Kit":{"sample_nums":96,},       
                      "MGIEasy Enzymatic DNA Library Prep Kit":{"sample_nums":96,},       
                      "MGICare Chromosomal Copy Number Variation Detection Kit":{"sample_nums":96,},       
                      "MGIEasy Fast Enzymatic DNA Library Prep Kit":{"sample_nums":96,},       
                      "MGIEasy Fast PCR-FREE Enzymatic Library Prep Kit":{"sample_nums":96,}, 
                      "MGIEasy Fast RNA Library Prep Kit":{"sample_nums":96,},      
                      "MGIEasy Directional RNA Library Prep Kit":{"sample_nums":96,},       
                             
     }
    describe_dict = {"sample_nums": "Sample quantity",
                     "spike_in_control_pcr_block": "Whether to add quality control amplification reagent",
                     "amplicon_multiplicity": "Panel amplicon multiplicity",
                     "panel_amp_length": "Panel amplicon length",
                     "barcode_type": "PCR barcode type"}
    if stage ==1:
        stage += 1
        response = f'Next is the protocol design stage. Please select the library construction kit you want to use from the following. If there is no kit that suits you, you can also manually upload your kit description file.'

        response_dict = {
            'response':response,
            'protocol_design_type': protocol_design_type,
            "operations":[
                {"title": "Protocol library","key": "selected_protocol","type": ["select","single"],"value":[] ,"options": all_protocol_name},
                {"title":"upload protocol file", "options":[], "type":["file"], "key":"protocol_file_path", "value":[]}
            ],
            'data':{},
            'stage':stage,
            'state':'continue'
        }
        return response_dict
    
    elif stage == 2:
        stage += 1
        protocol_design_type = ''
        if instruction.get("protocol_file_path"):
            protocol_design_type = 'auto_protocol_design'
        else:
            protocol_design_type = 'template_protocol_design'
        response_dict = {}
        if protocol_design_type == 'template_protocol_design':
            # Find corresponding execution steps based on selected template protocol
            # Let users provide all parameters directly
            template_protocol_name = instruction['selected_protocol'][0]
            protocol_file_name = all_protocol_file_dict[template_protocol_name]
            protocol_file_path = os.path.join('/reference_data/pdf_to_txt',protocol_file_name)
            txt_file = ''
            with open(protocol_file_path,'r')as f:
                lines = f.readlines()
                txt_file =' '.join(lines)# Also remove, keep only the main content, e.g. "First round PCR reaction" becomes 'PCR reaction'
                f.close()
            # Extract library construction process information and extract summary of each step content
            experiment_type_prompt = '请仔细分析用户提供的文本内容，这是一个文库构建的protocol。请提取目录中能够在OT-2上执行的步骤，并保留每个步骤的摘要描述（包括所有移液信息）。如果有单标签和双标签模块，请仅保留单标签步骤。如果目录中包含以下步骤，也请加入流程：杂交后PCR和杂交后PCR产物纯化，但不要重复加入。对于磁珠双选和磁珠单选，采用磁珠双选的方法。如果目录中的步骤名称有"第一轮"、"第二轮"等修饰词，请保留，但去掉括号内的内容。若目录中有"磁珠片段筛选"步骤，则根据前一步修改名称：如果前一步是酶切打断，则修改为"打断产物纯化"；如果前一步是PCR反应，则修改为"PCR产物纯化"。请跳过质检、文库混合、酶切消化、环化和试剂配制的步骤，并按先后顺序存放在JSON中，返回一个JSON。JSON格式如下：{'第一步'：{"名称"：步骤的名称，"描述"：该步骤的具体内容，包括所使用到的试剂类别和使用的体积，以及每次度溶液的操作步骤}，'第二步'：{"名称"：，"描述"：}}，注意：不能出现名称和描述都完全一样的两个步骤，具体例子如下：{'第一步'：{"名称"："RT-PCR 及多重扩增"，"描述"：首先使用25 μL RT-PCR Buffer ,2.5 μL RT-PCR Enzyme Mix, 2 μL ATOPlex HIV-1 Drug Resistance Primer Pool, 0.5 μL Nuclease-Free Water混合制备成RT-PCR 及多重扩增反应液 单个反应体积为30 μL。接着吸取20μL样本到新的样本管中，再吸取30μL反应液到样本管，最后进行反应，条件为50℃ 30 min，94℃ 3 min，1个循环；94℃ 30 s，58℃ 45 s，72℃ 2 min，42个循环；72℃ 5 min，12℃ Hold。1个循环；}'
            experiment_info_response = get_llm_chat_completion(
               messages=[
                    {"role": "system", "content": experiment_type_prompt},
                    {"role": "user", "content": txt_file}
                ],
                response_format={"type": "json_object"},
                temperature=0.2,
             )
            print('experiment_info_response', experiment_info_response.choices[0].message.content)
            experiment_process = experiment_info_response.choices[0].message.content
            print('experiment_process',type(experiment_process))

            total_tokens = experiment_info_response.usage.total_tokens
            print(f"************Total tokens used: {total_tokens}***************")
            logging.info(f'*************Total tokens used*******：\n{total_tokens}')
            # experiment_process = json.loads(experiment_process)
            process_list = get_protocol_process(template_protocol_name)
            # Let users provide all parameters at once
            option_list = []
            # Get parameter dictionary needed for this kit
            temp_parm_dict = all_protocol_param_dict[template_protocol_name]

            for param in temp_parm_dict.keys():
                option_list.append({"title":param,"options":[],"describe":describe_dict[param],"type":["single","input"],"key":param,"value":[temp_parm_dict[param]]})

            new_list = [x.split('.')[0] for x in process_list]
            

            response = f'Using {template_protocol_name} to construct library, it can be divided into {len(process_list)} steps as shown below: {experiment_process}. Next I will design the execution code for each step. According to the experimental requirements of {template_protocol_name}, you need to determine the following experimental parameters.'
            response_dict = {
                "response":response,
                'protocol_design_type': protocol_design_type,
                "operations":option_list,
                'data':{"protocol_param_dict":temp_parm_dict,"operation_steps":process_list,"selected_protocol":template_protocol_name},
                "stage":stage,
                "state":"continue"
            }
        elif protocol_design_type == 'auto_protocol_design':

            protocol_file = instruction['protocol_file_path'][0]
            txt_file = ''
            with open(protocol_file, 'r') as f:
                lines = f.readlines()
                txt_file = ' '.join(lines)  # Also remove, keep only the main content
                f.close()
            # Extract library construction process information and extract summary of each step content
            experiment_type_prompt = '请仔细分析用户提供的文本内容，这是一个文库构建的protocol。请提取目录中能够在OT-2上执行的步骤，并保留每个步骤的摘要描述（包括所有移液信息）。如果有单标签和双标签模块，请仅保留单标签步骤。如果目录中包含以下步骤，也请加入流程：杂交后PCR和杂交后PCR产物纯化，但不要重复加入。对于磁珠双选和磁珠单选，采用磁珠双选的方法。如果目录中的步骤名称有"第一轮"、"第二轮"等修饰词，请保留，但去掉括号内的内容。若目录中有"磁珠片段筛选"步骤，则根据前一步修改名称：如果前一步是酶切打断，则修改为"打断产物纯化"；如果前一步是PCR反应，则修改为"PCR产物纯化"。请跳过质检、文库混合、酶切消化、环化和试剂配制的步骤，并按先后顺序存放在JSON中，返回一个JSON。JSON格式如下：{'第一步'：{"名称"：步骤的名称，"描述"：该步骤的具体内容，包括所使用到的试剂类别和使用的体积，以及每次度溶液的操作步骤}，'第二步'：{"名称"：，"描述"：}}，注意：不能出现名称和描述都完全一样的两个步骤，具体例子如下：{'第一步'：{"名称"："RT-PCR 及多重扩增"，"描述"：首先使用25 μL RT-PCR Buffer ,2.5 μL RT-PCR Enzyme Mix, 2 μL ATOPlex HIV-1 Drug Resistance Primer Pool, 0.5 μL Nuclease-Free Water混合制备成RT-PCR 及多重扩增反应液 单个反应体积为30 μL。接着吸取20μL样本到新的样本管中，再吸取30μL反应液到样本管，最后进行反应，条件为50℃ 30 min，94℃ 3 min，1个循环；94℃ 30 s，58℃ 45 s，72℃ 2 min，42个循环；72℃ 5 min，12℃ Hold。1个循环；}'
            experiment_info_response = get_llm_chat_completion(
                messages=[
                    {"role": "system", "content": experiment_type_prompt},
                    {"role": "user", "content": txt_file}
                ],
                response_format={"type": "json_object"},
                temperature=0.2,
            )
            print('experiment_info_response', experiment_info_response.choices[0].message.content)
            experiment_response = experiment_info_response.choices[0].message.content

            total_tokens = experiment_info_response.usage.total_tokens
            print(f"************Total tokens used: {total_tokens}***************")
            logging.info(f'*************Total tokens used*******：\n{total_tokens}')

            experiment_process = json.loads(experiment_response)
            # Convert each step content into vectors using embedding and retrieve the best matching vector from vector database to return corresponding code block
            query_text_list = []
            name_list = []
            for process in experiment_process.keys():
                name = experiment_process[process]['名称']
                query_text = experiment_process[process]['描述']
                query_text_list.append(query_text)
                name_list.append(name)
            query_file_path = "/reference_data/query_text_code.csv"
            temp_df = pd.DataFrame({"text": query_text_list, "name": name_list})
            temp_df.to_csv(query_file_path, index=False)

            url = "http://10.49.60.23:8010/upload/"
            files = []
            files.append(('files', (os.path.basename(query_file_path), open(query_file_path, 'rb'))))
            # Send POST request
            response = requests.post(url, files=files)

            # Convert Python object to DataFrame
            match_df = pd.read_json(response.json())
            # Print returned results
            # print('response',response.json())
            # script_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'embedding_app.py')
            # cmd = f'python {script_path} --query_file_path {query_file_path}'
            # print('-' * 30, 'cmd: \n', cmd)
            # process = os.popen(cmd)  # return file
            # output = process.read()
            # process.close()
            # match_df = pd.read_csv(query_file_path)
            # Get final matched common code blocks
            code_list = list(match_df['code_block'])
            code_block_str = ','.join(code_list)
            code_block_path = '/reference_data/protocol_block_code_alphatool'
            code_path_list = [os.path.join(code_block_path, x) for x in code_list]
            # Extract parameters that need to be passed by users
            all_user_list = []
            for ii in range(len(code_list)):
                temp = [code_list[ii]]
                para_list, user_input_list = get_param_name(temp)
                if len(user_input_list) != 0:
                    for xx in user_input_list:
                        if xx not in all_user_list:
                            all_user_list.append(xx)


            response = f'The library construction steps in your provided text are as follows: {experiment_response}. We have matched corresponding code blocks from the code library for each specific operation, as shown below: {code_block_str}. Please confirm and provide the corresponding input parameters.'
            option_list = []
            for param in all_user_list:
                option_list.append(
                    {"title": param, "options": [], "type": ["single", "input"], "key": param, "value": []})

            response_dict = {'response': response,
                             'data':{'code_block': code_path_list,'protocol_file_path': protocol_file,
                                     "user_input_list": all_user_list,"operation_name":name_list },
                             "operations": option_list,
                             "protocol_design_type": protocol_design_type,
                             "stage": stage, "state": "continue"}
        return response_dict

    elif stage==3:
        stage+=1
        protocol_design_type = instruction["protocol_design_type"]
        if protocol_design_type == 'template_protocol_design':
            # Modify parameter dictionary based on user uploaded parameters
            temp_parm_dict = instruction['data']['protocol_param_dict']
            new_param_dict ={}
            for kk in temp_parm_dict.keys():
                new_param_dict[kk] = instruction[kk][0]
            template_protocol_name = instruction['data']['selected_protocol']
            process_list = instruction['data']['operation_steps']
            new_param_dict = calculating_solution_volume(new_param_dict,template_protocol_name)
            print('new_param_dict',new_param_dict)
            # Extract layout information for all steps
            Layout_info_dict =[]
            for temp_process in process_list:
                txt_filename = temp_process.replace('py','txt')
                txt_path = '/reference_data/block_txt'
                txt_file = os.path.join(txt_path, txt_filename)
                Layout_info_dict.append(txt_file)

            # Modify block code based on parameter dictionary
            code_path = '/reference_data/template_protocol_code'
            all_code_dict = []
            for temp_p in process_list:
                code_file = os.path.join(code_path, temp_p)
                new_code_path = modify_code_file(code_file, new_param_dict)
                all_code_dict.append(new_code_path)

            response = f"Below are the initial layout arrangements and corresponding execution code for each step of the experiment. Please place your experimental consumables and equipment according to the layout requirements, and start executing according to the code for each step after ensuring proper placement."
            response_dict = {
                "response":response,
                'protocol_design_type': protocol_design_type,
                'data':{"new_code":all_code_dict,"layout_info":Layout_info_dict,
                        "protocol_param_dict":new_param_dict},
                "stage":stage,
                'state':'stop'
            }
        elif protocol_design_type == 'auto_protocol_design':

            code_list = instruction['data']['new_code']
            print('Start parameter extraction and solution calculation')
            name_list = instruction['data']['operation_name']
            protocol_file = instruction['data']['protocol_file_path']
            user_input_list = instruction['data']['user_input_list']

            txt_file = ''
            with open(protocol_file, 'r') as f:
                lines = f.readlines()
                txt_file = ' '.join(lines)  #
                f.close()

            all_param_list = []
            for ii in range(len(code_list)):
                temp = [os.path.basename(code_list[ii])]
                temp_name = name_list[ii]
                para_list, _ = get_param_name(temp)

                print('Parameters to extract', para_list)
                experiment_type_prompt = (
                    f'请仔细阅读和分析文档中的内容，这是一个构建标准文库实验的protocol，目录中是文库构建的执行步骤，现在你需要聚焦于步骤{temp_name}中的所有内容，提取出该步骤中吸取（移液）溶剂的体积和样本体积。需要抽取的溶液参数列表如下：{para_list}\n'
                    f'User将会提供文档，你需要基于文档内容，找到当前步骤下参数列表中每种溶液吸取的体积。注意：参数input_Single_sample_volume 指的是上一个步骤样本管中剩余的溶液体积，input_final_transfer_volume指的是当前步骤最终移取上清液的体积，如果热循环中的温度是hold，则把他修改为60秒，如果参数是MIX的混合液，则只提取total的总体积，并以JSON格式返回。'
                    f'以下是一个具体的例子：'
                    f'User提供的文本为： 3.5.10将离心管瞬时离心，置于磁力架上，静置2-5 min至波体澄清，用移液器吸取19 μL上清液转移到新的0.2mL离心管中。| 3.6 PCR |  |  |  A 注意：PCR扩增步骤需要严格控制扩增循环数。循环数不足，会导致文库产出不足；循环数过多， | 又会影响后续数据性能表现。表12给出获得500ng和1ugPCR产物推荐的扩增循环数。当基因组DNA质量较差、主带较长时，需适当提高循环数以获取足量产物。 表12获得500 ng 和1μgPCR产物推荐的扩增循环数 |  |  |  | 打断后磁珠片段 对应产量所需循环数 |  |   基因组DNA（ng） | 筛选/纯化 | 500 ng 1 μg |  | 400 ng | 磁珠片段筛选 | 4-6 6-8 7-9 |  | 200 ng | 磁珠片段筛选 | 5-7 100 ng | 纯化 | 5-7 |  | 50 ng | 纯化 | 7-9 |  |13 | 7-9 9-11 |  | 3.6.1在冰上配制PCR反应液（见表13）。 |  |  |  |  |  | | 表13 PCR反应液的配制 |  |  |  |  |  | 组分 |  |  |  |  |  | PCR Enzyme Mix | 体积 25 μL |  |  |  |  | PCR Primer Mix | 6 μL |   | Total | 31μL |  |  |  | 3.6.2用移液器吸取31μL配制好的PCR反应液加入步骤3.5.10的PCR管中，涡旋震荡3次，每次3 S，瞬时离心将反应液收集至管底。 |  |  |  3.6.3 | 将步骤3.6.2所述PCR管置于PCR仪上，按照表14的条件进行PCR反应。 |  |  |  |  | PCR扩增反应条件 |  |  |  |  | 表14 | 时间 | 循环数 |  |  | 温度 | on |  |  |  | | 热盖 95°C | 3 min | 1循环 |  |  | 98°C |  |  |  |  | | 60°C |  | 20 s | 8循环 |  | 72°C |  | 15 s |  |   | 72°C |  | 30 s | 1循环 |    | 4°℃ |  | 10 min |  |  | A |  | Hold |  |  |  | | 注意：表中的循环数以200ng建库起始量为标准，不同的起始量请参考表12进行调整。 |  |  |  |   3.6.4 |  | 反应结束后，瞬时离心将反应液收集至管底，吸取全部反应液转移到新的1.5mL离心管中。 |  |  |  | 3.7PCR产物纯化 | A |  |  |  |  |  |  | 注意：操作前请仔细阅读附录A关于磁珠及纯化。 |  |  |  | 3.7.1 |  | 提前30 min取出 DNA Clean Beads置于室温，使用前充分震荡混匀。 |  |  |  | 3.7.2 | 吸取50 μL DNA Clean Beads至步骤3.6.4 的PCR产物中，充分震荡混匀至所有磁珠悬浮。 |  |  |  |  |3.7.3 | 室温孵育5 min。 |  |  |  |  | 3.7.4 |  | 将离心管瞬时离心，置于磁力架，静置2-5min至液体澄清，用移液器小心吸取上清并丢弃。 | 14 |  |  | 3.7.5 | 取上清并丢弃。 | 保持离心管置于磁力架上，加入200uL新鲜配制的80%乙醇漂洗磁珠及管壁，静置30s后小心吸 |  |  |  | 3.7.6 |  | 重复步骤3.7.5，尽量吸干管内液体，有少量残留在管壁时可将离心管瞬时离心，在磁力架上分离后， |  |  |  | | MGI |用小量程的移液器将管底液体吸干。 |  | 3.7.7 | 保持离心管固定于磁力架上，打开离心管管盖，室温干燥，直至磁珠表面无反光、无开裂。 将离心管从磁力架上取下，加入32μL TE Buffer进行DNA洗脱，用移液器轻轻吹打至少10次至 | 3.7.8 所有磁珠悬浮。 |  | 3.7.9 室温下孵育5 min。 | 3.7.10将离心管瞬时离心，置于磁力架，静置2-5min至液体澄清，将30 μL上清液转移到新的1.5 mL |'
                    f'需要计算的参数列表为["input_Single_sample_volume","input_PCR_MIX_volume","input_thermocycler_param","input_DNA_Beads_volume","input_TE_Buffer_volume","input_Ethanol_80_volume","input_final_transfer_volume"]。'
                    f'返回的JSON为：{json.dumps({"input_Single_sample_volume": 19, "input_PCR_MIX_volume": 31, "input_thermocycler_param": [{"temper_time": [(95, 180)], "repetition": 1}, {"temper_time": [(98, 20), (60, 15), (72, 30)], "repetition": 8}, {"temper_time": [(72, 600), (4, 60)], "repetition": 1}], "input_DNA_Beads_volume": 50, "input_TE_Buffer_volume": 32, "input_Ethanol_80_volume": 200, "input_final_transfer_volume": 30, }, indent=4)}'
                )
                experiment_info_response = get_llm_chat_completion(
                    messages=[
                        {"role": "system", "content": experiment_type_prompt},
                        {"role": "user", "content": txt_file}
                    ],
                    response_format={"type": "json_object"},
                    temperature=0.2,
                )
                print('experiment_type_prompt', experiment_type_prompt)
                print('experiment_info_response', experiment_info_response.choices[0].message.content)
                experiment_process = experiment_info_response.choices[0].message.content
                experiment_process = json.loads(experiment_process)
                all_param_list.append(experiment_process)
                # all_param_dict = all_param_dict + str(experiment_process)

            user_input_dict = {}
            for user_input in user_input_list:
                user_input_dict[user_input] = instruction[user_input][0]
            # Get the number of samples from user
            Sample_num = int(instruction['user_sample_nums'][0])
            new_para_list = []
            for ii in range(len(all_param_list)):
                temp_param = all_param_list[ii]
                for kk in temp_param.keys():
                    if 'Ethanol_80' in kk:
                        temp_param[kk] = temp_param[kk] * (Sample_num + 8) * 2
                    elif '_volume' in kk and 'transfer' not in kk and 'Single_' not in kk:
                        temp_param[kk] = temp_param[kk] * (Sample_num + 8)
                new_para_list.append(temp_param)

            print('new_para_list', new_para_list)
            option_list = []
            for ii in range(len(new_para_list)):
                temp_dict = new_para_list[ii]
                step_str = f'step{str(ii + 1)}-'
                for param in temp_dict.keys():
                    option_list.append(
                        {"title": step_str + param, "options": [], "type": ["single", "input"], "key": step_str + param,
                         "value": [temp_dict[param]]})
            response = f'Based on the sample parameters provided by the user, we have calculated the required liquid volume parameters for each code block as follows, please confirm.'
            response_dict = {"response": response,
                             'data':{"new_code": code_list, 'new_para_list': new_para_list,'user_input_dict': user_input_dict},
                             "operations": option_list,
                             'protocol_design_type': protocol_design_type,
                             "stage": stage, 'state': 'continue'}
        else:
            # log the protocol_design_type is incorrect
            response_dict = {'state': 'stop'}
        return response_dict

    elif stage == 4:
        stage += 1
        print('Modifying code according to parameter list')
        code_list = instruction['data']['new_code']
        new_modify_para_list = []
        new_para_list = instruction['data']['new_para_list']
        for ii in range(len(new_para_list)):
            new_dict = {}
            temp_dict = new_para_list[ii]
            step_str = f'step{str(ii + 1)}-'
            for key in temp_dict.keys():
                new_dict[key] = instruction[step_str + key][0]
            new_modify_para_list.append(new_dict)

        user_input_dict = instruction['data']['user_input_dict']
        template_code_path = '/reference_data/protocol_block_code_alphatool'
        new_code_path = []
        for ii in range(len(code_list)):
            code_path = os.path.join(template_code_path, code_list[ii])
            txt_file = ''
            with open(code_path, 'r') as f:
                lines = f.readlines()
                txt_file = ''.join(lines)
                f.close()
            temp_dict = new_modify_para_list[ii]
            temp_dict.update(user_input_dict)
            experiment_type_prompt = (
                f'You need to complete a keyword replacement task. The user will provide a code snippet and a replacement dictionary.\n'
                f'The task is to replace keywords in the code with corresponding values from the dictionary. Keep other parts of the code unchanged and return the modified code.\n'
                f'Note: Only return the modified code, no other explanations needed\n'
                f'Replacement dictionary：{temp_dict}\n'
                f'Here is an example：\n'
                '''
                Sample_nums = input_sample_nums
                MIX_volume = input_PCR_MIX_volume
                Single_sample_volume = input_Single_sample_volume
                thermocycler_param = input_thermocycler_param

                transfer_volume_1 = 10
                transfer_volume_2 = MIX_volume / (Sample_nums + 8)

                use_col_nums = int(Sample_nums / 8)
                if Sample_nums % 8 != 0:
                    use_col_nums += 1
                temp_module = protocol.load_module('temperature module', '3')
                temp_plate = temp_module.load_labware('biorad_96_wellplate_200ul_pcr')
                plate1 = protocol.load_labware('biorad_96_wellplate_200ul_pcr', '1')
                thermocycler_module = protocol.load_module('thermocycler', '7')
                thermocycler_plate = thermocycler_module.load_labware('biorad_96_wellplate_200ul_pcr')
                '''
                f'Replacement dictionary：{{"input_sample_nums": 96, "input_PCR_MIX_volume": 3224, "input_Single_sample_volume": 50, "input_thermocycler_param":[{{"temper_time": [[37, 1800], [65, 900], [4, 60]], "repetition": 1}}]}}\n'
                f'Result after replacement：\n'
                '''
                Sample_nums = 96
                MIX_volume = 3224
                Single_sample_volume = 50
                thermocycler_param = [{"temper_time": [[37, 1800], [65, 900], [4, 60]], "repetition": 1}]

                transfer_volume_1 = 10
                transfer_volume_2 = MIX_volume / (Sample_nums + 8)

                use_col_nums = int(Sample_nums / 8)
                if Sample_nums % 8 != 0:
                    use_col_nums += 1
                temp_module = protocol.load_module('temperature module', '3')
                temp_plate = temp_module.load_labware('biorad_96_wellplate_200ul_pcr')
                plate1 = protocol.load_labware('biorad_96_wellplate_200ul_pcr', '1')
                thermocycler_module = protocol.load_module('thermocycler', '7')
                thermocycler_plate = thermocycler_module.load_labware('biorad_96_wellplate_200ul_pcr')
                '''
            )
            experiment_info_response = get_llm_chat_completion(
                messages=[
                    {"role": "system", "content": experiment_type_prompt},
                    {"role": "user", "content": txt_file}
                ],
                temperature=0.2,
            )
            print('experiment_type_prompt', experiment_type_prompt)
            print('experiment_info_response', experiment_info_response.choices[0].message.content)
            experiment_process = experiment_info_response.choices[0].message.content or ''
            new_code_name = os.path.basename(code_list[ii]).split('.')[0]+'_modofy.py'
            new_path = os.path.join('/reference_data/protocol_generate_code',new_code_name)
            new_code_path.append(new_path)
            with open(new_path, 'w') as f:
                lines = experiment_process.split('\n')
                for ll in lines:
                    if '```' not in ll:
                        f.write(ll + '\n')
                f.close()

        # Use LLM to get layout information
        layout_info_list = []
        for ii in range(len(new_code_path)):
            code_path = new_code_path[ii]
            txt_file = ''
            with open(code_path, 'r') as f:
                lines = f.readlines()
                txt_file = ''.join(lines)
                f.close()

            experiment_type_prompt = (
                f'You need to extract the layout information from the experimental code provided by the user and return a layout information summary in natural language\n'
                f'Note: Please follow the template in the example\n'
                f'Here is an example：\n'
                f'Experimental code：\n'
                '''
                Sample_nums = user_sample_nums
                Beads_volume = input_DNA_Beads_volume
                TE_Buffer_volume = input_TE_Buffer_volume
                Ethanol_80_volume = input_Ethanol_80_volume
        
                beads_transfer_volume = Beads_volume/(Sample_nums+8)
                TE_Buffer_transfer_volume = TE_Buffer_volume/(Sample_nums+8)
                Ethanol_80_transfer_volume = Ethanol_80_volume/(Sample_nums+8)/2
                Single_sample_volume = 50
                use_col_nums = int(Sample_nums/8)
                if Sample_nums%8 !=0:
                    use_col_nums = use_col_nums+1
        
                # Load 200ul filter tip rack
                tips200_1 = ctx.load_labware('mdk_96_filtertiprack_200ul', 4)
                tips200_2 = ctx.load_labware('mdk_96_filtertiprack_200ul', 5)
                tips200_3 = ctx.load_labware('mdk_96_filtertiprack_200ul', 6)
                tips200_4 = ctx.load_labware('mdk_96_filtertiprack_200ul', 7)
                tips200_5 = ctx.load_labware('mdk_96_filtertiprack_200ul', 8)
                tips200_6 = ctx.load_labware('mdk_96_filtertiprack_200ul', 9)
                tips200_7 = ctx.load_labware('mdk_96_filtertiprack_200ul', 10)
                tips200_8 = ctx.load_labware('mdk_96_filtertiprack_200ul', 11)
        
                plate1 = ctx.load_labware('biorad_96_wellplate_200ul_pcr', 3)
                # Load magnetic module and PCR plate on it
                mag_module = ctx.load_module("magneticModuleV2", 1)
                mag_plate = mag_module.load_labware('biorad_96_wellplate_200ul_pcr')
                reservoir_plate = ctx.load_labware('usascientific_12_reservoir_22ml', 2)
        
                # Pipettes
                right_pipette = ctx.load_pipette("p200_multi", "left")
        
                # Execute instructions
                # Transfer
                transfer_times = use_col_nums
                transfer_info_right1 = []
                for ii in range(transfer_times):
                    transfer_info_right1.append((0,ii, beads_transfer_volume))
        
                transfer_all(right_pipette,reservoir_plate,mag_plate,tips200_1,transfer_info_right1,mix_flag=1,mix_value=[10,50])
                ctx.pause(300)
                mag_module.engage(9.0)
                ctx.pause(300)
        
                transfer_info_right2 = []
                for ii in range(transfer_times):
                    transfer_info_right2.append((ii,4,beads_transfer_volume+Single_sample_volume))
        
                transfer_all(right_pipette, mag_plate,reservoir_plate,tips200_1,transfer_info_right2)
        
                transfer_info_right3 = []
                for ii in range(transfer_times):
                    transfer_info_right3.append((2, ii, Ethanol_80_transfer_volume))
        
                transfer_all(right_pipette, reservoir_plate, mag_plate,tips200_1, transfer_info_right3)
                ctx.pause(30)
        
                transfer_info_right4 = []
                for ii in range(transfer_times):
                    transfer_info_right4.append((ii, 5, Ethanol_80_transfer_volume))
        
                transfer_all(right_pipette, mag_plate, reservoir_plate, tips200_1, transfer_info_right4)
                '''
                f'Return layout information：\n'
                '''
                Layout:
                    (1)Place 200ul mdk_96_filtertiprack_200ul tip racks in slots 4,5,6,7,8,9,10,11.
                    (2)Place magnetic module in slot 1, with a 200μL biorad_96_wellplate_200ul_pcr plate on it.
                    (3)Place a 200μL biorad_96_wellplate_200ul_pcr plate in slot 3.
                    (4)Place a 22ml usascientific_12_reservoir_22ml reservoir in slot 2.
                    (5)Use P200 multichannel pipette on the left side, using tip racks from slots 4,5,6,7,8,9,10,11.
                '''
            )
            experiment_info_response = get_llm_chat_completion(
                messages=[
                    {"role": "system", "content": experiment_type_prompt},
                    {"role": "user", "content": txt_file}
                ],
                temperature=0.2,
            ) 
            print('experiment_type_prompt', experiment_type_prompt)
            print('experiment_info_response', experiment_info_response.choices[0].message.content)
            experiment_process = experiment_info_response.choices[0].message.content or ''
            temp_txt_path = os.path.join(os.path.dirname(code_path),os.path.basename(code_path).split('.')[0]+'.txt')
            with open(temp_txt_path, 'w') as f:
                f.write(experiment_process)
                f.close()
            layout_info_list.append(temp_txt_path)
        # Convert python to json
        json_file_list =[]
        for ii in range(len(new_code_path)):
            new_code = new_code_path[ii]
            temp_path = os.path.join(os.path.dirname(new_code),os.path.basename(new_code).replace('.py','.json'))
            json_file_list.append(temp_path)
            cmd = f'mgi_alphatool -i {new_code} -o {temp_path}'
            process = os.popen(cmd)  # return file
            output = process.read()
            print(output)
            process.close()
        response = 'Below are the executable code and layout information generated according to user requirements. Users can import them into OT-2 machine to start the library construction experiment.'

        response_dict = {"response": response, 'data':{"new_code": new_code_path,'layout_info':layout_info_list,'json_file':json_file_list}, "stage": stage, 'state': 'stop'}

        return response_dict


app = FastAPI()

# Parse input data
class InputModel(BaseModel):
    instruction: str #user input
    conversation: list #context dialogue
    type: int
    stage: int #dialogue rounds
    upload_file_flag: bool
    species_identification_dict: dict
    experiment_select_strain: list
    non_target_select: list
    select_target_gene: list
    strain_select: list

# Heartbeat test to ensure server connection status
async def heartbeat(websocket: WebSocket):
    while True:
        await asyncio.sleep(30)
        try:
            print("HEARTBEAT")
            await websocket.send_text("heartbeat")
        except Exception as e:
            print("Heartbeat failed, retrying:", e)
            break

executor = ProcessPoolExecutor(max_workers=4)#ThreadPoolExecutor()
               
@app.websocket("/protocol_design")
async def websocket_endpoint(websocket: WebSocket):
    await websocket.accept()
    try:
        while True:
            data = await websocket.receive_text()
            print('input_data',data)
            input_dict = json.loads(data)
            print('data',input_dict)
            stage = input_dict["stage"]
            instruction = input_dict["instruction"]
            responses_dict = template_protocol_design(instruction, stage)
            print('responses_dict',responses_dict)
            await websocket.send_text(json.dumps({'content':responses_dict, 'state':'stop'},ensure_ascii=False)) 
            
    except Exception as e:
        traceback.print_exc()
    finally:
        try:
            await websocket.close()
        except RuntimeError:
            # Ignore this exception if connection is already closed
            pass
            
if __name__ == "__main__":
    openai.api_type = "azure"
    openai.api_base = ""
    openai.api_version = ""
    openai.api_key = ""

    os.system('rm log')
    # Start service
    uvicorn.run(app, host="0.0.0.0", port=8032,timeout_keep_alive=300)
