
import os

import json
from typing import List, Optional
import ast
import uvicorn
from Bio import Entrez
import asyncio
from fastapi import FastAPI,WebSocket
from pydantic import BaseModel
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor

from fastapi.responses import HTMLResponse, FileResponse
from fastapi import FastAPI, File, UploadFile, HTTPException, WebSocketDisconnect
# API请求有频率限制，TODO：多账号混合、时延机制
from protocol_agent import modify_code_file,template_protocol_design,get_param_name
from primer_design import snp_primer_design,redesign_primer
from SNP_Genotyping_search import SNP_Genotyping_search
from download_files import download_files
from cancer_search import cancer_search
from genetic_disorder_search import genetic_disorder_search
from whole_genome_search import whole_genome_search
from protein_mutation_search import protein_mutation_search
from species_identification_search import species_identification_search
from pathogen_drug_resistance_search import pathogen_drug_resistance_search
import traceback

import logging
from openai.types.chat.chat_completion import ChatCompletion
from openai.types.chat.chat_completion_message import ChatCompletionMessage
from logic.experiment import ExperimentManager
from llm_utils.utils import get_llm_chat_completion
from utils import init_logger, init_redis_cli, get_redis_conn


logger = logging.getLogger(__name__)


Entrez.email = "912837656@qq.com"
Entrez.api_key = "24e5b06635f7b32dbd9b89c178a3caf3f009"#"db4b31edb00e220c3dd378908403eddbc308"


def init_everything(app: FastAPI):
    init_logger()
    redis_conn = init_redis_cli()
    app.state.redis = redis_conn
    app.state.experiment_manager = ExperimentManager(get_redis_conn())


def parse_llm_response(response: ChatCompletion) -> ChatCompletionMessage:
    return response.choices[0].message


primer_functions= [
{
        "name": "download_files",  # 基因分型检索
        "description": " If you are downloading files now, this tool can help you download the relevant files.",
        "parameters": {
            "type": "object",
            "properties": {
                "name_list": {
                    "type": "string",
                    "description": "the name list."
                },

            },
            "required": ["name_list"]
        }
    },
    {
        "name": "snp_primer_design",
        "description": " ",
        "parameters": {
            "type": "object",
            "properties": {
                
                "instruction": {
                    "type": "string",
                    "description": ""
                },
            },
            "required": ["instruction"]
        }
    },
    {
        "name": "redesign_primer",
        "description": " redesign primer",
        "parameters": {
            "type": "object",
            "properties": {

                "instruction": {
                    "type": "string",
                    "description": ""
                },
            },
            "required": ["instruction"]
        }
    },
]
functions= [  
    {
        "name": "protein_mutation_search",#蛋白质突变检索
        "description": "If you want to use multiplex PCR methods for tNGS library construction and sequencing on various enzyme mutant libraries obtained through protein engineering, you can use this tool to obtain protein and gene information.",
        "parameters": {
            "type": "object",
            "properties":{
                "term":{
                    "type":"string",
                    "description":"The search query or terms. When the user has multiple information terms(like species name and protein name), “AND” should be used for splicing. For example, if you get “luciferase and Gaussia”, you need to extract “luciferase AND Gaussia”.",
                },
                
            },
            "required":["term"]
        }
    },
    {
        "name": "whole_genome_search",  # 全基因组检索
        "description": "If you want to use multiplex PCR methods for tNGS library construction and sequencing on whole genome detection of viruses, you can use this tool to obtain gene information.",
        "parameters": {
            "type": "object",
            "properties": {
                "term": {
                    "type": "string",
                    "description": "",
                },

            },
            "required": ["term"]
        }
    },
	{
        "name": "genetic_disorder_search",#遗传病检索
        "description": "If the ongoing experiment is a multiplex PCR experiment for genetic disorder analysis, this tool can be used to obtain target gene sequence and snp information. ",
        "parameters": {
            "type": "object",
            "properties": {
                "genetic_disease": {
                    "type": "string",
                    "description": "the genetic disease name."
                },
                
            },
            "required": ["genetic_disease"]
        }
    },
    {
        "name":"cancer_search", #肿瘤检索
        "description":"If the ongoing experiment is a multiplex PCR experiment for cancer analysis, this tool can be used to obtain target gene sequence and snp information. ",
        "parameters":{
            "type":"object",
            "properties":{
                "term":{
                    "type":"string",
                    "description":"Extract cancer-related terms from the user's input, which may be the name of the cancer.",
                },
            },
            "required":["term"]
        }
    },
    {
        "name":"pathogen_drug_resistance_search", #病原耐药性检索
        "description":"If the ongoing experiment is a multiplex PCR experiment for pathogen drug resistance analysis, this tool can be used to obtain target gene sequence and snp information. ",
        "parameters":{
            "type":"object",
            "properties": {
                "instruction": {
                    "type":"string",
                    "description":"the species name,When the user provides an abbreviation of the species name, you need to convert it to the full name of the species.",
                    },
                },
            "required":["instruction"]
            }
    },
    {
        "name":"species_identification_search", #物种鉴定检索 Species Identification
        "description":"If the ongoing experiment is a multiplex PCR experiment for Species Identification analysis, this tool can be used to obtain target gene sequence and snp information. ",
        "parameters": {
            "type": "object",
            "properties": {
                "instruction": {
                    "type": "string",
                    "description": "the specified species name.If the user provides an abbreviation, please help me change it to the full name."
                },
                
            },
            "required": ["instruction"]
            }
    },
    {
        "name": "SNP_Genotyping_search",  # 基因分型检索
        "description": "If the ongoing experiment is a multiplex PCR experiment for SNP typing, you can use this tool to obtain the target gene sequence and SNP information. ",
        "parameters": {
            "type": "object",
            "properties": {
                "specie_name": {
                    "type": "string",
                    "description": "the  species name."
                },

            },
            "required": ["specie_name"]
        }
    },
    {
        "name": "download_files",  # 基因分型检索
        "description": " If you are downloading files now, this tool can help you download the relevant files.",
        "parameters": {
            "type": "object",
            "properties": {
                "name_list": {
                    "type": "string",
                    "description": "the name list."
                },

            },
            "required": ["name_list"]
        }
    },
    {
        "name": "snp_primer_design",
        "description": " ",
        "parameters": {
            "type": "object",
            "properties": {
                "max_amp_len": {
                    "type": "number",
                    "description": "the Maximum amplicon length."
                },
                "min_amp_len": {
                    "type": "number",
                    "description": "the Minimum amplicon length."
                },
                "temperature": {
                    "type": "number",
                    "description": "the temperature."
                },
                "min_GC": {
                    "type": "number",
                    "description": "the Minimum GC content."
                },
                "max_GC": {
                    "type": "number",
                    "description": "the Maximum GC content."
                },
                "max_len": {
                    "type": "number",
                    "description": "the Maximum length of a primer. "
                },
            },
            "required": ["max_amp_len","min_amp_len","temperature","min_GC","max_GC","max_len"]
        }
    },
]


#定义X-chat的对话逻辑
def predict(input_data, stage=1, search_type=None, search_prompt=""):
    history_conversation = input_data["conversation"]
    window_history = history_conversation[-50:]
    try:
        if search_type==None:
            system_prompt = [{'role':'system','content':"""You are an artificial intelligence assistant designed to help users retrieve biometric information sequences. 
    When a user requests protein mutant testing, 
    you should call the protein_mutation_search function to obtain protein and gene information; when a user requests a multiplex PCR experiment for genetic disease analysis, 
    you should call the genetic_disorder_search function to obtain the target gene sequence and SNP information; When a user requests cancer analysis or wants to obtain cancer-related genes, 
    you should call the cancer_search function to obtain the target gene sequence and SNP information; 
    when the user requests a  antimicrobial resistance research When performing a multiplex PCR experiment, 
    you should call the pathogen_drug_resistance_search function to obtain the target gene sequence and SNP information; 
    When users want to do virus-related whole genome testing,you should call the whole_genome_search function to obtain the target gene sequence ;
    When a user wants to obtain or download a gene sequence or genome sequence, you should  call the download_files to meet the user's needs;
    when the user requests a multiplex PCR experiment for pathogen species identification analysis, you should call the species_identification_search function；In addition, the language in which you reply to the user should be consistent with the language used by the user."""
                            }]
            messages = system_prompt + window_history
            response = get_llm_chat_completion(
                messages=messages,
                temperature=0.7,
                top_p=0.95,
                max_tokens=4096,
                stream=False,
                functions=functions,
                function_call="auto",
            )
            search_response = parse_llm_response(response)
            if search_response.function_call:
                print('*******************Successful activate function*******************')
                # Call the function. The JSON response may not always be valid so make sure to handle errors
                function_name = search_response.function_call.name
                available_functions = {"protein_mutation_search": protein_mutation_search,
                                    "genetic_disorder_search": genetic_disorder_search,
                                    "cancer_search": cancer_search,
                                    "download_files": download_files,
                                    "pathogen_drug_resistance_search":pathogen_drug_resistance_search,
                                    "whole_genome_search":whole_genome_search,
                                    "species_identification_search":species_identification_search,
                                     "snp_primer_design":snp_primer_design,
                                    }
                function_to_call = available_functions[function_name]
                print("Function name: ",search_response.function_call.name)
                print("Function input: ",search_response.function_call.arguments)
                
                
                
                if search_response.function_call.name == "protein_mutation_search":
                
                    function_args = json.loads(search_response.function_call.arguments)
                
                    #instruction = function_args['term']
                        #user_input = window_history[-1]['content']#function_args['user_input']
                        #data_info_dict = function_to_call(stage=1,instruction=instruction,user_input=user_input)
                    responses = protein_mutation_search(input_data,stage=1)
                    #user_response = pathogen_drug_resistance_search(stage=1,instruction=species_name)
                    #responses = user_response
                    history_conversation.append({'role': 'assistant', 'content': responses})
                    responses_dict = {'responses': responses, 'history_conversation': history_conversation}
                    return responses_dict
                        #messages.append( # adding assistant response to messages
                        #        {
                        #            "role": search_response["role"],
                        #            "function_call": {
                       #                 "name": function_name,
                       #                 "arguments": search_response.function_call.arguments,
                       #                 },
                       #             "content": None
                       #             }
                       #         )
                       # messages.append( # adding assistant response to messages
                       #         {
                       #             "role": "function",
                       #             "name":function_name,
                       #             "content": data_info_dict["data_list"],
                       #             }
                       #         )
                        
                        #Modify user query in multiple rounds of self-reflection
                       # if data_info_dict["data_list"]==[]:
                           # i=1
                            #query_now = function_args["term"]
                            ##while True:
                                #print(f"""If no NCBI results are retrieved, it may be that the query is too long or too short. The current query is {query_now}. 
                                #Please reflect on the modification. You only need to return my new term without unnecessary description.""")
                                #messages.append({"role":'user','content':f"""If no NCBI results are retrieved, it may be that the query is too long or too short. The current query is {query_now}. 
                                #Please reflect on the modification. You only need to return my new term without unnecessary description."""})
                                #second_response = openai.ChatCompletion.create(
                                #                messages=messages,
                                #                deployment_id="XMGI-Chat1-GPT4"
                                #                )
                                #print(f'The modified term after the {i}-th round of reflection is:',second_response["choices"][0]["message"]['content'])
                                #data_info_dict = function_to_call(instruction=second_response["choices"][0]["message"]['content'], stage=1)
                                #if data_info_dict["data_list"]!=[]:
                                #    break
                                #i+=1
                                #if i >5:
                                #    break
                                # print(f'The search results of round {i} are:',data_info_dict["data_list"])
                                # query_now=second_response["choices"][0]["message"]['content']
    
    
                        # responses = 'Based on your question, I searched for the relevant sequences as follows, please select according to your needs:\n'+str(data_info_dict["options"])
    
                        #history_conversation.append({'role':'assistant','content':data_info_dict})
                        #responses_dict = {'responses':data_info_dict,'history_conversation':history_conversation}
                        #return responses_dict
                        # return data_info_dict
                elif search_response.function_call.name == "download_files":
                    function_args = json.loads(search_response.function_call.arguments)
                    name = function_args['name_list']
                    user_response = download_files(stage=1, instruction=name)
                    responses = user_response
                    history_conversation.append({'role': 'assistant', 'content': responses})
                    responses_dict = {'responses': responses, 'history_conversation': history_conversation}
                    return responses_dict

                elif search_response.function_call.name == "whole_genome_search":
                    function_args = json.loads(search_response.function_call.arguments)

                    user_response = whole_genome_search(stage=1, instruction=input_data)
                    responses = user_response
                    history_conversation.append({'role': 'assistant', 'content': responses})
                    responses_dict = {'responses': responses, 'history_conversation': history_conversation}
                    # responses_dict={'target_cds_path':[],'non_target_cds_path':[],'fna_file':[fna_path],'snp_file':[snp_path],'download_path':[],'responses':responses,'history_conversation':history_conversation}
                    return responses_dict

                elif search_response.function_call.name == "pathogen_drug_resistance_search":
                    function_args = json.loads(search_response.function_call.arguments)
                    #species_name = function_args['instruction']
                    responses = pathogen_drug_resistance_search(input_data,stage=1)
                    #user_response = pathogen_drug_resistance_search(stage=1,instruction=species_name)
                    #responses = user_response
                    history_conversation.append({'role': 'assistant', 'content': responses})
                    responses_dict = {'responses': responses, 'history_conversation': history_conversation}
                    return responses_dict
    
                elif search_response.function_call.name == "snp_primer_design":
                    print('成功调用primer design')
                    function_args = json.loads(search_response.function_call.arguments)
                    print('xxx',history_conversation[-1]['content'])
                    instruction = ast.literal_eval(history_conversation[-1]['content'])
                    # history_conversation = json.loads()
                    # instruction = history_conversation
                    print('xxxxinstruction',instruction)

                    responses = snp_primer_design(instruction, stage)
                    history_conversation.append({'role': 'assistant', 'content': responses})
                    responses_dict = {'responses': responses, 'history_conversation': history_conversation}
                    return responses_dict

                elif search_response.function_call.name == "genetic_disorder_search":
                    function_args = json.loads(search_response.function_call.arguments)
                    genetic_disease = function_args['genetic_disease']

                    user_response = genetic_disorder_search(input_data,stage=1)
                    responses = user_response
                    history_conversation.append({'role':'assistant','content':responses})
                    responses_dict = {'responses': responses, 'history_conversation': history_conversation}
                    return responses_dict

                elif search_response.function_call.name == "cancer_search":
                    function_args = json.loads(search_response.function_call.arguments)
                    #instruction = function_args['term']
                    
                    user_response = cancer_search(stage=1, instruction=input_data)
                    responses = user_response
                    history_conversation.append({'role': 'assistant', 'content': responses})
                    responses_dict = {'responses': responses, 'history_conversation': history_conversation}
                    # responses_dict={'target_cds_path':[],'non_target_cds_path':[],'fna_file':[fna_path],'snp_file':[snp_path],'download_path':[],'responses':responses,'history_conversation':history_conversation}
                    return responses_dict

                elif search_response.function_call.name =="species_identification_search":
                    function_args = json.loads(search_response.function_call.arguments)
                    instruction = history_conversation#function_args['instruction']++++++++++++++++++++++++++

                    responses = species_identification_search(input_data,stage=1)
                    history_conversation.append({'role':'assistant','content':responses})
                    responses_dict = {'responses': responses, 'history_conversation': history_conversation}
                    # responses_dict['history_conversation'] = history_conversation
                    return responses_dict

                else:
                    print("error funcation calling")

            else:
                responses = search_response.content
                history_conversation.append({'role':'assistant','content':responses})
                responses_dict={'responses': '没有调用函数', 'history_conversation': history_conversation}
                return responses_dict

        else:
            system_prompt = [{'role':'system','content':f"""You are an artificial intelligence assistant designed to help users retrieve biometric information sequences. {search_prompt}"""}]
            if search_type=="cancer_type":

                  #instruction = json.loads(history_conversation[-1]['content'])
                  #instruction = instruction['selected_options']
                  responses = cancer_search(input_data,stage)
                  #responses = cancer_search(instruction, stage)
                  history_conversation.append({'role': 'assistant', 'content': responses})
                  responses_dict = {'responses': responses, 'history_conversation': history_conversation}
                  return responses_dict

            elif search_type=="genetic_disorder_type":

                  responses = genetic_disorder_search(input_data, stage)
                  history_conversation.append({'role': 'assistant', 'content': responses})
                  responses_dict = {'responses': responses, 'history_conversation': history_conversation}
                  return responses_dict

            elif search_type=="whole_genome_type":

                  responses = whole_genome_search(input_data, stage)
                  history_conversation.append({'role': 'assistant', 'content': responses})
                  responses_dict = {'responses': responses, 'history_conversation': history_conversation}
                  return responses_dict
            elif search_type=="protein_mutation_type":
                #if stage==2:
                  #instruction = json.loads(history_conversation[-1]['content'])
                  #instruction = instruction['selected_options']
                responses = protein_mutation_search(input_data,stage)
                history_conversation.append({'role':'assistant','content':responses})
                responses_dict = {'responses': responses, 'history_conversation': history_conversation}
                return responses_dict

            elif search_type == "download_type":
                print('开始正式下载文件')
                #instruction = json.loads(history_conversation[-1]['content'])
                #instruction = instruction['selected_options']
                responses = download_files(input_data,stage)
                #responses = download_files(instruction, stage)
                history_conversation.append({'role': 'assistant', 'content': responses})
                responses_dict = {'responses': responses, 'history_conversation': history_conversation}
                return responses_dict



            elif search_type=="pathogen_drug_resistance_type":
                instruction = history_conversation
                responses = pathogen_drug_resistance_search(input_data,stage)
                history_conversation.append({'role':'assistant','content':responses})
                responses_dict = {'responses': responses, 'history_conversation': history_conversation}
                return responses_dict

            elif search_type =="species_identification_type":
                instruction = history_conversation
                responses = species_identification_search(input_data,stage)
                history_conversation.append({'role':'assistant','content':responses})
                responses_dict = {'responses': responses, 'history_conversation': history_conversation}
                return responses_dict

            else:
                print("error search type")
    except Exception as e:
        print("Exception:"+str(e))
        print("Traceback:")
        traceback.print_exc()  # 打印异常调用堆栈
        return "Exception:"+str(e)

# 定义X-chat的对话逻辑
def primer_predict(input_data, primer_type=None,stage=1, search_prompt=""):
    history_conversation = input_data["conversation"]
    window_history = history_conversation[-50:]
    try:
        if primer_type == None:
            system_prompt = [{'role': 'system', 'content': """You are an artificial intelligence assistant designed to help users design specific primers. 
    When the user asks for follow-up work or subsequent primer design , 
    you should call the snp_primer_design function to design primers;
    When the user replies "redesign Primer" in the history conversation, you should call the redesign_primer function to redesign primers;
    """ }]
            messages = system_prompt + window_history
            response = get_llm_chat_completion(messages=messages, functions=primer_functions)
            search_response = parse_llm_response(response)

            if search_response.function_call:
                print('*******************Successful activate function*******************')
                # Call the function. The JSON response may not always be valid so make sure to handle errors
                function_name = search_response.function_call.name
                available_functions = {
                                       "snp_primer_design": snp_primer_design,
                                        "redesign_primer":redesign_primer,
                                       }
                function_to_call = available_functions[function_name]
                print("Function name: ", search_response.function_call.name)
                print("Function input: ", search_response.function_call.arguments)
    
                if search_response.function_call.name == "snp_primer_design":
                    print('成功调用primer design')


                    responses = snp_primer_design(input_data,stage=stage)
                    history_conversation.append({'role': 'assistant', 'content': responses})
                    responses_dict = {'responses': responses, 'history_conversation': history_conversation}
                    return responses_dict
                elif search_response.function_call.name == "redesign_primer":
                    print('成功调用redesign')
                    responses = redesign_primer(input_data, stage=stage)
                    history_conversation.append({'role': 'assistant', 'content': responses})
                    responses_dict = {'responses': responses, 'history_conversation': history_conversation}
                    return responses_dict
                else:
                    print("error funcation calling")

            else:
                responses = search_response.content
                history_conversation.append({'role': 'assistant', 'content': responses})
                responses_dict = {'target_cds_path': [], 'non_target_cds_path': [], 'snp_file': [], 'fna_file': [],
                                  'download_path': [], 'responses': responses, 'history_conversation': history_conversation}
                return responses_dict

        elif primer_type == 'redesign_primer_type':
              print('进入函数')
              responses = redesign_primer(input_data,stage)
              history_conversation.append({'role': 'assistant', 'content': responses})
              responses_dict = {'responses': responses, 'history_conversation': history_conversation}
              return responses_dict
        else:
              print('进入函数')   
              responses = snp_primer_design(input_data,stage)
              history_conversation.append({'role': 'assistant', 'content': responses})
              responses_dict = {'responses': responses, 'history_conversation': history_conversation}
              return responses_dict

    except Exception as e:

        print("Exception:"+str(e))
        print("Traceback:")
        traceback.print_exc()  # 打印异常调用堆栈
        return "Exception:"+str(e)

app = FastAPI()


# 解析输入数据
class InputModel(BaseModel):
    instruction: str #用户输入
    conversation: list #上下文对话
    type: int
    stage: int #对话轮次
    upload_file_flag: bool
    species_identification_dict: dict
    experiment_select_strain: list
    non_target_select: list
    select_target_gene: list
    strain_select: list


#心跳测试，保证服务器连接状态
async def heartbeat(websocket: WebSocket):
    while True:
        await asyncio.sleep(30)
        try:
            print("HEARTBEAT")
            await websocket.send_text("heartbeat")
        except Exception as e:
            print("Heartbeat failed, retrying:", e)
            break

executor = ProcessPoolExecutor(max_workers=50)#ThreadPoolExecutor()


# WebSocket 路由
@app.websocket("/ncbi")
async def websocket_endpoint(websocket: WebSocket):
    await websocket.accept()

    asyncio.create_task(heartbeat(websocket))
    try:
        while True:
            data = await websocket.receive_text()
            #解析为JSON
            input_data = json.loads(data)
            logger.info(f"Endpoint: '/ncbi' 接收到数据：{input_data}")
            input_data["instruction"] = json.loads(input_data["instruction"])
            all_history_conversation = input_data["instruction"]["conversation"]
            stage = input_data["instruction"].get("stage") or 1


            #第一轮对话
            if stage == 1:
                loop = asyncio.get_event_loop()
                input_data["conversation"] = all_history_conversation
                responses_dict = await loop.run_in_executor(executor, predict, input_data, stage)

            #多轮对话
            else:

                input_data["conversation"] = all_history_conversation
                search_type = input_data["instruction"]["search_type"]

                temp_flag = False
                if search_type != "":
                    temp_flag = True

                loop = asyncio.get_event_loop()
                if temp_flag:
                    logger.info(f"Endpoint: '/ncbi' 进行第二次回复：{input_data}")
                    print('xxxxx')
                    responses_dict = await loop.run_in_executor(executor, predict, input_data, stage, search_type,'')
                else:
                    responses_dict = await loop.run_in_executor(executor, predict, input_data)

            await websocket.send_text(json.dumps({'content':responses_dict, 'state':'stop'},ensure_ascii=False))
                    
    except Exception as e:
        traceback.print_exc()
    finally:
        print('finished!')

@app.websocket("/primer")
async def websocket_endpoint(websocket: WebSocket):
    await websocket.accept()
    asyncio.create_task(heartbeat(websocket))
    #all_history_conversation = []
    try:
        while True:
            data = await websocket.receive_text()
            #print("接收到数据：",data)
            input_data = json.loads(data) #解析为JSON
            logger.info(f"Endpoint: '/primer' 接收到数据：{input_data}")
            input_data["instruction"] = json.loads(input_data["instruction"])
            all_history_conversation = input_data["instruction"]["conversation"]
            stage = input_data["instruction"].get("stage") or 1

            # 第一轮对话
            if stage == 1:
                all_history_conversation = input_data["instruction"]["conversation"]
                loop = asyncio.get_event_loop()
                input_data["conversation"] = all_history_conversation
                responses_dict = await loop.run_in_executor(executor, primer_predict, input_data)
            # 多轮对话
            else:
                print('进入第二轮对话')
                input_data["conversation"] = all_history_conversation
                primer_type = input_data["instruction"]["primer_type"]
                
                print('stage',stage)
                temp_flag = False
                if primer_type != "":
                    temp_flag = True

                loop = asyncio.get_event_loop()
                if temp_flag:
                    print('xxxxx进行第二次回复')
                    responses_dict = await loop.run_in_executor(executor, primer_predict, input_data, primer_type,stage)
                else:
                    responses_dict = await loop.run_in_executor(executor, primer_predict, input_data)
                #print("当前多轮对话用户输入：", input_)
                #history_conversation = [{'role': his_['role'], 'content': his_['content']} for his_ in history]
                # history_conversation.append({'role': 'user', 'content': f'{input_}'})
                #history_conversation.append({'role': 'user', 'content': f'{data}'})
                # print("当前用户输入：", input_)
                #print("当前多轮对话历史：", history_conversation)
                #loop = asyncio.get_event_loop()
                #responses_dict = await loop.run_in_executor(executor, primer_predict, history_conversation)

            await websocket.send_text(
                json.dumps({'content': responses_dict, 'state': 'stop'}, ensure_ascii=False))

    except Exception as e:
        traceback.print_exc()
    finally:
        print('finished!')


@app.websocket("/protocol_design")
async def websocket_endpoint(websocket: WebSocket):
    await websocket.accept()
    try:
        while True:
            data = await websocket.receive_text()
            #解析为JSON
            input_data = json.loads(data)
            logger.info(f"Endpoint: '/protocol_design' 接收到数据：{input_data}")
            input_data["instruction"] = json.loads(input_data["instruction"])
            # all_history_conversation = input_data["instruction"]["conversation"]
            responses_dict = template_protocol_design(
                input_data["instruction"],
                input_data["instruction"].get("stage", 1)
            )
            print('responses_dict', responses_dict)
            await websocket.send_text(json.dumps({'content': {'responses': responses_dict}, 'state': 'stop'}, ensure_ascii=False))
    except Exception as e:
        traceback.print_exc()
    finally:
        try:
            await websocket.close()
        except RuntimeError:
            # 如果连接已经关闭，忽略这个异常
            pass

#上传文件功能：
UPLOAD_FOLDER = '/reference_data/upload_file' # 设置上传文件的保存目录
ALLOWED_EXTENSIONS = {'txt', 'pdf', 'fasta', 'fa', 'csv'}  # 允许上传的文件类型


# # 检查文件类型是否允许上传
def allowed_file(filename):
    return '.' in filename and \
        filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

@app.post("/upload")
async def upload_file(files: List[UploadFile] = File(...)):
    uploaded_filenames = []
    for file in files:
        try:
            print('开始上传')
            contents = await file.read()
            file_path = os.path.join(UPLOAD_FOLDER, file.filename)
            with open(file_path, 'wb') as f:
                f.write(contents)
            uploaded_filenames.append(file_path)

        except Exception as e:
            raise HTTPException(status_code=500, detail=str(e))
    return {"filenames": uploaded_filenames, "success": 'True'}

#@app.get("/v1/files/xsearchdev/{file_path:path}")
@app.get("/reference_data/{file_path:path}")#"/download_reference_data")
async def download_reference_data(file_path:str):
    print("download: ", file_path)
    #return FileResponse(file_path)
    print('/reference_data/' + file_path)
    # pd_ = pd.read_csv('/reference_data/' + file_path)
    # print('已读取pd',pd_)
    # response = FileResponse('/reference_data/' + file_path)
    response = FileResponse(file_path)
    # response.headers["Content-Disposition"] = f'attachment; filename="{file_path.split("/")[-1]}"'

    from urllib.parse import quote
    filename = file_path.split("/")[-1]
    quoted_filename = quote(filename)
    response.headers["Content-Disposition"] = f'attachment; filename*=UTF-8\'\'{quoted_filename}'

    return response


class ConversationExperimentModel(BaseModel):
    conversation_uuid: str
    # protocol design 产生的结果, txt文件列表, 用LLM解析板位信息
    panel_list: List[str]
    protocol_path: List[str]


class ExperimentResponse(BaseModel):

    class ExperimentData(BaseModel):
        experiment_uuid: str

    data: ExperimentData | None
    code: int
    error_message: str


@app.post("/api/conversation/experiment")
async def get_experiment(
    conversation_experiment: ConversationExperimentModel
) -> ExperimentResponse:
    """根据对话id获取实验信息"""
    # create or get experiment
    from logic.alphatool_checks import tips_test, get_panel_list_from_txt

    def load_protocol_file(file_path: str):
        with open(file_path, 'r') as protocol_file:
            return protocol_file.read()

    logger.info(
        f"[get_experiment] {conversation_experiment}"
    )

    experiment_manager = app.state.experiment_manager  # type: ExperimentManager
    experiment = experiment_manager.get_experiment_by_conversation(
        conversation_experiment.conversation_uuid
    )

    logger.info(
        f"[get_experiment] saveExperiment: {experiment.uuid}"
    )
    protocols = [
        json.loads(load_protocol_file(protocol_file))
        for protocol_file in conversation_experiment.protocol_path
    ]
    # FIXME: 删除测试代码
    # protocols = [
    #     load_protocol_file(protocol_file)
    #     for protocol_file in [
    #         '/Users/gala/work/Genomics/code/pcr-agent/tests/protocol_debug.json',
    #     ]
    # ]

    # positions = []
    logger.info(
        f"[get_experiment]: Panels: {conversation_experiment.panel_list},"
        f"protocolNumbers: {len(protocols)}"
    )
    panel_number_list = []

    for panel_file in conversation_experiment.panel_list:
        with open(panel_file, 'r') as panel_file:
            panel_content = panel_file.read()
            logger.info(f"[get_experiment] panel_content: {panel_content}")
            panel_number_list.append(
                get_panel_list_from_txt(panel_content)
            )

    # FIXME: 删除测试代码
    # for panel_file in ['/Users/gala/work/Genomics/code/pcr-agent/tests/test_panel.txt']:
    #     with open(panel_file, 'r') as panel_file:
    #         panel_content = panel_file.read()
    #         logger.info(f"[get_experiment] panel_content: {panel_content}")
    #         panel_number_list.append(
    #             get_panel_list_from_txt(panel_content)
    #         )

    experiment.panel_list = panel_number_list
    experiment.protocols = []
    experiment.current_run = {}
    experiment.current_protocol_index = 0
    experiment_manager.save_experiment(experiment)
    # protocol 列表
    for position_list, protocol in zip(experiment.panel_list, protocols):
        # 检查各个板位吸头的protocol
        pre_check_protocols = [
            tips_test(position)
            for position in position_list
        ]
        # 该protocol加上检查的protocols 得到总的要运行的protocol列表
        one_complete_protocol = pre_check_protocols + [protocol]
        experiment.protocols.append(one_complete_protocol)
        experiment_manager.save_experiment(experiment)

    return ExperimentResponse(
        data=ExperimentResponse.ExperimentData(experiment_uuid=experiment.uuid),
        code=200,
        error_message=''
    )


@app.websocket('/ws/experiment')
async def ws_experiment(
    *,
    ws_conn: WebSocket,
    conversation_uuid: str,
):
    """
    Experiment code execution websocket
    request msg format:
    {
        "command": "start" | "resume" | "stop",
    }
    response msg format:
    {
        "responses": {
            "code": 200,
            "response": "success",
            "data": {
                "content": "response content",
                "type": "text" | "image",
                "experiment_info": experiment_info,
                "generating": true | false
                "command_options": ["start", "resume", "stop"] | [],
                "state": "running" | "succeeded" | "failed" | "paused",
            }
        }
    }
    """
    await ws_conn.accept()

    experiment_manager = app.state.experiment_manager
    experiment_manager.set_ws_conn(
        conversation_uuid=conversation_uuid,
        ws=ws_conn,
    )
    experiment = experiment_manager.get_experiment_by_conversation(
        conversation_uuid
    )
    try:
        async for message in ws_conn.iter_json():
            logger.info(f"[ws_experiment] receive: {message}")
            command = message.get('command')
            if command not in ('start', 'resume', 'stop'):
                await experiment_manager.send_ws_data(
                    ws=ws_conn,
                    code=400,
                    msg=f"[ERROR] Invalid command: {command}",
                    experiment=experiment,
                    data_type='text',
                    generating=False,
                )
            elif command == 'start':
                await experiment_manager.start_experiment(ws_conn, experiment)
            elif command == 'resume':
                await experiment_manager.resume_run(ws_conn, experiment)
            elif command == 'redo':
                await experiment_manager.redo_action(ws_conn, experiment)
            elif command == 'stop':
                try:
                    stop_result = experiment_manager.try_stop_current_run(experiment)
                except Exception as e:
                    await experiment_manager.send_ws_data(
                        ws=ws_conn,
                        code=500,
                        msg=f"[ERROR] Failed to stop the experiment: {e}",
                        experiment=experiment,
                        data_type='text',
                        generating=False,
                    )
                else:
                    await experiment_manager.send_ws_data(
                        ws=ws_conn,
                        code=200 if stop_result else 400,
                        msg=f'[INFO] Stop experiment result: {stop_result}',
                        experiment=experiment,
                        data_type='text',
                        generating=False,
                    )
    except WebSocketDisconnect as e:
        logger.error(f"ConversationWSDisconnected: {conversation_uuid}, reason: {e}")
        await experiment_manager.send_ws_data(
            ws=ws_conn,
            code=500,
            msg=f"[ERROR] Websocket Disconnected {e}",
            experiment=experiment,
            data_type='text',
            generating=False,
        )
        await ws_conn.close()
    except Exception as e:
        logger.error(f"ConversationWSError: {conversation_uuid}, reason: {e}")
        logger.error(traceback.format_exc())
        await experiment_manager.send_ws_data(
            ws=ws_conn,
            code=500,
            msg=f"[ERROR] Experiment executed failed: {e}",
            experiment=experiment,
            data_type='text',
            generating=False,
            command_options=['start', 'stop'],
        )
        await ws_conn.close()


if __name__ == "__main__":
    # 开启服务
    init_everything(app)
    uvicorn.run(app, host="0.0.0.0", port=8081, timeout_keep_alive=300)
