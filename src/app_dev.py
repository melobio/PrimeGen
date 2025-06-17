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
from direct_protein_gene_search import direct_protein_gene_search
import traceback

import logging
from openai.types.chat.chat_completion import ChatCompletion
from openai.types.chat.chat_completion_message import ChatCompletionMessage
from logic.experiment import ExperimentManager
from llm_utils.utils import get_llm_chat_completion
from utils import init_logger, init_redis_cli, get_redis_conn


logger = logging.getLogger(__name__)


Entrez.email = ""
Entrez.api_key = ""


def init_everything(app: FastAPI):
    init_logger()
    redis_conn = init_redis_cli()
    app.state.redis = redis_conn
    app.state.experiment_manager = ExperimentManager(get_redis_conn())


def parse_llm_response(response: ChatCompletion) -> ChatCompletionMessage:
    return response.choices[0].message


primer_functions= [
    {
        "type": "function",
        "function":{
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
    }},
    {
        "type": "function",
        "function":{
        "name": "snp_primer_design",
        "description": "when user want to start primer design,please choose this function",
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
    }},
    {
        "type": "function",
        "function":{
        "name": "redesign_primer",
        "description": "when user want to redesign primer ,please choose this function",
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
    }},
]

tools = [
    {
        "type": "function",
        "function":{
        "name": "direct_protein_gene_search",
        "description": "If you want to search for protein and gene information directly without performing specific experiments like protein mutation search, whole genome search, genetic disorder search, or cancer search, you can use this tool. This is a general-purpose search tool suitable for basic protein and gene information queries.",
        "parameters": {
            "type": "object",
            "properties": {
                "instruction": {
                    "type": "string",
                    "description": "The search query or terms. When the user has multiple information terms(like species name and protein name), \"AND\" should be used for splicing. For example, if you get \"luciferase and Gaussia\", you need to extract \"luciferase AND Gaussia\"."
                }
            },
            "required": ["instruction"]
        }
    }},

    {
        "type": "function",
        "function":{
        "name": "protein_mutation_search", 
        "description": "If you want to use multiplex PCR methods for tNGS library construction and sequencing on various enzyme mutant libraries obtained through protein engineering, you can use this tool to obtain protein and gene information.",
        "parameters": {
            "type": "object",
            "properties": {
                "term": {
                    "type": "string",
                    "description": "The search query or terms. When the user has multiple information terms(like species name and protein name), \"AND\" should be used for splicing. For example, if you get \"luciferase and Gaussia\", you need to extract \"luciferase AND Gaussia\".",
                },

            },
            "required": ["term"]
        }
    }},

    {
        "type": "function",
        "function":{
        "name": "whole_genome_search", 
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
    }},

    {
        "type": "function",
        "function":{
        "name": "genetic_disorder_search",
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
    }},
    {
        "type": "function",
        "function":{
        "name": "cancer_search", 
        "description": "If the ongoing experiment is a multiplex PCR experiment for cancer analysis, this tool can be used to obtain target gene sequence and snp information. ",
        "parameters": {
            "type": "object",
            "properties": {
                "term": {
                    "type": "string",
                    "description": "Extract cancer-related terms from the user's input, which may be the name of the cancer.",
                },
            },
            "required": ["term"]
        }
    }},
        {
            "type": "function",
            "function":{
        "name": "pathogen_drug_resistance_search",
        "description": "If the ongoing experiment is a multiplex PCR experiment for pathogen drug resistance analysis, this tool can be used to obtain target gene sequence and snp information. ",
        "parameters": {
            "type": "object",
            "properties": {
                "instruction": {
                    "type": "string",
                    "description": "the species name,When the user provides an abbreviation of the species name, you need to convert it to the full name of the species.",
                },
            },
            "required": ["instruction"]
        }
    }},
    {
        "type": "function",
        "function":{
        "name": "species_identification_search",
        "description": "If the ongoing experiment is a multiplex PCR experiment for Species Identification analysis, this tool can be used to obtain target gene sequence and snp information. ",
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
    }},
    {
        "type": "function",
        "function":{
        "name": "SNP_Genotyping_search",
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
    }},
    {
        "type": "function",
        "function":{
        "name": "download_files", 
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
    }},
    {
        "type": "function",
        "function":{
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
            "required": ["max_amp_len", "min_amp_len", "temperature", "min_GC", "max_GC", "max_len"]
        }
    }},
]

def predict(input_data, stage=1, search_type=None, search_prompt=""):
    logging.info(f'-----------start stage: {stage}')
    logging.info(f'-----------input_data: {input_data}')
    history_conversation = input_data["conversation"]
    logging.info(f'-----------history_conversation: {history_conversation}')
    window_history = history_conversation[-50:]
    logging.info(f'-----------window_history: {window_history}')
    try:
        if search_type==None:
            system_prompt = [{'role':'system','content':"""You are an artificial intelligence assistant designed to help users retrieve biometric information sequences. 
    When a user requests protein and gene information directly, you should call the direct_protein_gene_search function to obtain protein and gene information;
    When a user requests protein mutant testing, you should call the protein_mutation_search function to obtain protein and gene information; 
    when a user requests a multiplex PCR experiment for genetic disease analysis, you should call the genetic_disorder_search function to obtain the target gene sequence and SNP information; 
    When a user requests cancer analysis or wants to obtain cancer-related genes, you should call the cancer_search function to obtain the target gene sequence and SNP information; 
    when the user requests a  antimicrobial resistance research When performing a multiplex PCR experiment, you should call the pathogen_drug_resistance_search function to obtain the target gene sequence and SNP information; 
    When users want to do virus-related whole genome testing,you should call the whole_genome_search function to obtain the target gene sequence ;
    When a user wants to obtain or download a gene sequence or genome sequence, you should  call the download_files to meet the user's needs;
    when the user requests a multiplex PCR experiment for pathogen species identification analysis, you should call the species_identification_search function；In addition, the language in which you reply to the user should be consistent with the language used by the user."""
                            }]
            # window_history = [{'role': 'user', 'content': 'Hi!'},] +window_history
            user_input = []
            for his in window_history:
                if his['role'] == 'user':
                    user_input.append(his)
            # messages = system_prompt + window_history
            messages = user_input
            response = get_llm_chat_completion(
                messages=messages,
                temperature=0.7,
                top_p=0.95,
                max_tokens=4096,
                stream=False,
                tools=tools,
                # functions=functions,
                tool_choice="required",
                # function_call="auto",
            )
            search_response = parse_llm_response(response)
            logger.info(f'*******************search_response:{search_response}*******************')
            # if search_response.function_call:
            if search_response.tool_calls:
                logger.info('*******************Successful activate function*******************')
                # Call the function. The JSON response may not always be valid so make sure to handle errors
                # function_name = search_response.function_call.name
                function_name = search_response.tool_calls[0].function.name
                available_functions = {"direct_protein_gene_search": direct_protein_gene_search,
                                    "protein_mutation_search": protein_mutation_search,
                                    "genetic_disorder_search": genetic_disorder_search,
                                    "cancer_search": cancer_search,
                                    "download_files": download_files,
                                    "pathogen_drug_resistance_search":pathogen_drug_resistance_search,
                                    "whole_genome_search":whole_genome_search,
                                    "species_identification_search":species_identification_search,
                                    "snp_primer_design":snp_primer_design,
                                    }
                logger.info(f'*******************function_name {function_name}*******************')
                # print("Function name: ",search_response.function_call.name)
                # print("Function input: ",search_response.function_call.arguments)
                
                
                
                if function_name == "protein_mutation_search":
                    # function_args = json.loads(search_response.function_call.arguments)
                    responses = protein_mutation_search(input_data,stage=1)
                    history_conversation.append({'role': 'assistant', 'content': responses})
                    responses_dict = {'responses': responses, 'history_conversation': history_conversation}
                    return responses_dict
                
                elif function_name == "direct_protein_gene_search":
                    # function_args = json.loads(search_response.function_call.arguments)
                    responses = direct_protein_gene_search(input_data,stage=1)
                    history_conversation.append({'role': 'assistant', 'content': responses})
                    responses_dict = {'responses': responses, 'history_conversation': history_conversation}
                    return responses_dict

                elif function_name == "whole_genome_search":
                    # function_args = json.loads(search_response.function_call.arguments)
                    user_response = whole_genome_search(stage=1, instruction=input_data)
                    responses = user_response
                    history_conversation.append({'role': 'assistant', 'content': responses})
                    responses_dict = {'responses': responses, 'history_conversation': history_conversation}
                    return responses_dict

                elif function_name == "pathogen_drug_resistance_search":
                    # function_args = json.loads(search_response.function_call.arguments)
                    responses = pathogen_drug_resistance_search(input_data,stage=1)
                    history_conversation.append({'role': 'assistant', 'content': responses})
                    responses_dict = {'responses': responses, 'history_conversation': history_conversation}
                    return responses_dict
    
                elif function_name == "snp_primer_design":
                    print('successful calling primer design')
                    # function_args = json.loads(search_response.function_call.arguments)
                    print('xxx',history_conversation[-1]['content'])
                    instruction = ast.literal_eval(history_conversation[-1]['content'])
                    print('xxxxinstruction',instruction)

                    responses = snp_primer_design(instruction, stage)
                    history_conversation.append({'role': 'assistant', 'content': responses})
                    responses_dict = {'responses': responses, 'history_conversation': history_conversation}
                    return responses_dict

                elif function_name == "genetic_disorder_search":
                    # function_args = json.loads(search_response.function_call.arguments)
                    # genetic_disease = function_args['genetic_disease']

                    user_response = genetic_disorder_search(input_data,stage=1)
                    responses = user_response
                    history_conversation.append({'role':'assistant','content':responses})
                    responses_dict = {'responses': responses, 'history_conversation': history_conversation}
                    return responses_dict

                elif function_name == "cancer_search":
                    # function_args = json.loads(search_response.function_call.arguments)
                    #instruction = function_args['term']
                    
                    user_response = cancer_search(stage=1, instruction=input_data)
                    responses = user_response
                    history_conversation.append({'role': 'assistant', 'content': responses})
                    responses_dict = {'responses': responses, 'history_conversation': history_conversation}
                    # responses_dict={'target_cds_path':[],'non_target_cds_path':[],'fna_file':[fna_path],'snp_file':[snp_path],'download_path':[],'responses':responses,'history_conversation':history_conversation}
                    return responses_dict

                elif function_name =="species_identification_search":
                    # function_args = json.loads(search_response.function_call.arguments)
                    # instruction = history_conversation#function_args['instruction']++++++++++++++++++++++++++

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
                  responses = cancer_search(input_data,stage)
                  history_conversation.append({'role': 'assistant', 'content': responses})
                  responses_dict = {'responses': responses, 'history_conversation': history_conversation}
                  return responses_dict
            
            elif search_type=="direct_protein_gene_search_type":
                responses = direct_protein_gene_search(input_data,stage)
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
                responses = protein_mutation_search(input_data,stage)
                history_conversation.append({'role':'assistant','content':responses})
                responses_dict = {'responses': responses, 'history_conversation': history_conversation}
                return responses_dict

            elif search_type == "download_type":
                print('开始正式下载文件')
                responses = download_files(input_data,stage)
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
        traceback.print_exc()
        return "Exception:"+str(e)

def primer_predict(input_data, primer_type=None,stage=1, search_prompt=""):
    history_conversation = input_data["conversation"]
    window_history = history_conversation[-50:]

    logging.info(f'-------- primer window_history: {window_history}')
    try:
        if primer_type == None:
            system_prompt = [{'role': 'system', 'content': """You are an artificial intelligence assistant designed to help users design specific primers. 
    When the user asks for follow-up work or start primer design in the history conversation , 
    you should call the snp_primer_design function to design primers;
    When the user replies "redesign Primer" in the history conversation, you should call the redesign_primer function to redesign primers;
    """ }]

            # messages = system_prompt + window_history
            user_input = []
            for his in window_history:
                if his['role'] == 'user':
                    user_input.append(his)

            messages = user_input
            response = get_llm_chat_completion(messages=messages,
                                               tools=primer_functions,
                                               tool_choice="required",)

            search_response = parse_llm_response(response)
            logger.info(f'*************search_response:{search_response}*************************')
            # if search_response.function_call:
            if search_response.tool_calls:
                logger.info('*******************Successful activate function*******************')
                # Call the function. The JSON response may not always be valid so make sure to handle errors
                # function_name = search_response.function_call.name
                function_name = search_response.tool_calls[0].function.name
                available_functions = {
                                       "snp_primer_design": snp_primer_design,
                                        "redesign_primer":redesign_primer,
                                       }
                # function_to_call = available_functions[function_name]
                logger.info(f"********************Function name:{function_name}************")
                # print("Function input: ", search_response.function_call.arguments)
    
                if function_name == "snp_primer_design":
                    logging.info('成功调用primer design')
                    responses = snp_primer_design(input_data,stage=stage)
                    history_conversation.append({'role': 'assistant', 'content': responses})
                    responses_dict = {'responses': responses, 'history_conversation': history_conversation}
                    return responses_dict
                elif function_name == "redesign_primer":
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
        traceback.print_exc() 
        return "Exception:"+str(e)

app = FastAPI()


class InputModel(BaseModel):
    instruction: str 
    conversation: list 
    type: int
    stage: int 
    upload_file_flag: bool
    species_identification_dict: dict
    experiment_select_strain: list
    non_target_select: list
    select_target_gene: list
    strain_select: list

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



@app.websocket("/ncbi")
async def websocket_endpoint(websocket: WebSocket):
    await websocket.accept()

    asyncio.create_task(heartbeat(websocket))
    try:
        while True:
            data = await websocket.receive_text()
            # parse to JSON
            input_data = json.loads(data)
            logger.info(f"Endpoint: '/ncbi' 接收到数据：{input_data}")
            input_data["instruction"] = json.loads(input_data["instruction"])
            all_history_conversation = input_data["instruction"]["conversation"]
            stage = input_data["instruction"].get("stage") or 1
            
            if "search_type"in input_data["instruction"]:
                search_type = input_data["instruction"]["search_type"]
            else:
                search_type = ""


            # 1st stage
            if stage == 1 and search_type == "":
                loop = asyncio.get_event_loop()
                input_data["conversation"] = all_history_conversation
                responses_dict = await loop.run_in_executor(executor, predict, input_data, stage)

            # multi-stage
            else:

                input_data["conversation"] = all_history_conversation
                search_type = input_data["instruction"]["search_type"]

                temp_flag = False
                if search_type != "":
                    temp_flag = True

                loop = asyncio.get_event_loop()
                if temp_flag:
                    logger.info(f"Endpoint: '/ncbi' 进行第二次回复：{input_data}")
                    
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
            logger.info(f"'/primer' 接收到数据：{data}")
            
            input_data = json.loads(data) 
            

            input_data["instruction"] = json.loads(input_data["instruction"])
            all_history_conversation = input_data["instruction"]["conversation"]
            stage = input_data["instruction"].get("stage") or 1

            
            if stage == 1:
                all_history_conversation = input_data["instruction"]["conversation"]
                loop = asyncio.get_event_loop()
                input_data["conversation"] = all_history_conversation
                logger.info(f"Endpoint: '/primer' input_data：{input_data}")
                logger.info("*************************************************")
                responses_dict = await loop.run_in_executor(executor, primer_predict, input_data)
            
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
            
            input_data = json.loads(data)
            logger.info(f"Endpoint: '/protocol_design' 接收到数据：{input_data}")
            logger.info("****************************************************")
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
            
            pass


UPLOAD_FOLDER = '/reference_data/upload_file' 
ALLOWED_EXTENSIONS = {'txt', 'pdf', 'fasta', 'fa', 'csv'}  


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


class Query(BaseModel):
    query: str


@app.post("/api/test_llm")
async def test_llm(query: Query):
    """Test LLM model"""
    from llm_utils.utils import get_llm_chat_completion
    response = get_llm_chat_completion(
        messages=[
            {"role": "user", "content": query.query},
        ],
        temperature=0.7,
        top_p=0.95,
        max_tokens=4096,
        stream=False,
    )
    return response


class ConversationExperimentModel(BaseModel):
    conversation_uuid: str
    # protocol design result, txt files list, use LLM to extract panel list
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
    """get experiment by conversation"""
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
    # protocols = [
    #     json.loads(load_protocol_file(protocol_file))
    #     for protocol_file in conversation_experiment.protocol_path
    # ]
    # FIXME: delete test code
    protocols = [
        load_protocol_file(protocol_file)
        for protocol_file in [
            '/reference_data/test_protocol.json',
            # '/Users/gala/work/Genomics/code/pcr-agent/tests/test_protocol.json'
        ]
    ]

    # positions = []
    logger.info(
        f"[get_experiment]: Panels: {conversation_experiment.panel_list},"
        f"protocolNumbers: {len(protocols)}"
    )
    panel_number_list = []


    experiment.panel_list = []
    experiment.protocols = [protocols, ]
    experiment.current_run = {}
    experiment.current_protocol_index = 0
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
            if command not in ('start', 'resume', 'stop', 'redo'):
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
            command_options=['resume', 'stop'],
        )
        await ws_conn.close()


if __name__ == "__main__":
    # start server
    init_everything(app)
    uvicorn.run(app, host="0.0.0.0", port=8081, timeout_keep_alive=300)
