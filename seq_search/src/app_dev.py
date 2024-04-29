#import glob2 as glob
import shutil
import glob
import time
import zipfile
import os
import time as time_
import pickle
import openai
from tenacity import retry,wait_fixed
import traceback
import json
import requests
from typing import List
import ast
from bs4 import BeautifulSoup
import uvicorn
from Bio import Entrez
import asyncio
from fastapi import FastAPI,WebSocket
from pydantic import BaseModel
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
import pandas as pd
import socket
import re
# from flask import Flask, request, redirect, url_for, flash
# from werkzeug.utils import secure_filename
import threading
from fastapi.responses import HTMLResponse
from fastapi import FastAPI, File, UploadFile, HTTPException
# API请求有频率限制，TODO：多账号混合、时延机制

from primer_design import snp_primer_design,primer3_design
from SNP_Genotyping_search import SNP_Genotyping_search
from download_files import download_files
from cancer_search import cancer_search
from genetic_disorder_search import genetic_disorder_search
from protein_mutation_search import protein_mutation_search
from species_identification_search import species_identification_search
from pathogen_drug_resistance_search import pathogen_drug_resistance_search

Entrez.email = "912837656@qq.com"
Entrez.api_key = "24e5b06635f7b32dbd9b89c178a3caf3f009"#"db4b31edb00e220c3dd378908403eddbc308"

#free gpt user 20 request/min
#@retry(wait=wait_fixed(20))
def get_chatgpt_response(messages, functions=[], stream_flag=False):#input,a_):
    response = openai.ChatCompletion.create(
            temperature=0.6,
            engine='XMGI-Chat1-GPT4',
            messages=messages,
            stream=stream_flag,
            functions = functions,
            function_call="auto",
            )
    
    if stream_flag:
        return  response
    elif functions!=[]:
        return response['choices'][0]['message']
    else:
        return response['choices'][0]['message']['content']

primer_functions= [
{
        "name": "download_files",  # 基因分型检索
        "desctription": " If you are downloading files now, this tool can help you download the relevant files.",
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
        "desctription": " ",
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
        "desctription":"If the ongoing experiment is a multiplex PCR experiment for cancer analysis, this tool can be used to obtain target gene sequence and snp information. ",
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
        "desctription":"If the ongoing experiment is a multiplex PCR experiment for pathogen drug resistance analysis, this tool can be used to obtain target gene sequence and snp information. ",
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
        "desctription":"If the ongoing experiment is a multiplex PCR experiment for Species Identification analysis, this tool can be used to obtain target gene sequence and snp information. ",
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
        "desctription": "If the ongoing experiment is a multiplex PCR experiment for SNP typing, you can use this tool to obtain the target gene sequence and SNP information. ",
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
        "desctription": " If you are downloading files now, this tool can help you download the relevant files.",
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
        "desctription": " ",
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







#Protein first-level search: continue searching according to user needs, returning protein name, species introduction and protein length
# def search_uniprot_one(term):
#     get_info_url = f"http://rest.uniprot.org/uniprotkb/search?query=%28protein_name%3A{term}%29+AND+%28reviewed%3Atrue%29&format=json"
#     # get_info_url = f"http://rest.uniprot.org/uniprotkb/search?query={term}&format=json"
#     response = requests.get(get_info_url)
#     search_info_list = []
#     data_list = []
#     if response.status_code != 200:
#         print("Error:", response.status_code)
#         return search_info_list, data_list
#         exit()
#
#     # parse JSON
#     data = response.json()
#
#
#     i = 0
#
#     for data_info in data['results']:
#         i += 1
#         data_dict = {}
#         search_info = {}
#         data_dict["Function information"] = {}
#
#         if 'comments' in list(data_info.keys()):
#             for comments_ in data_info["comments"]:
#                 try:
#                     if 'FUNCTION' in comments_["commentType"]:
#                         data_dict["Function information"][comments_['commentType']] = comments_['texts'][0]['value']
#                     elif 'CATALYTIC ACTIVITY' in comments_["commentType"]:
#                         data_dict["Function information"][comments_['commentType']] = comments_['reaction']['name']
#                     elif 'COFACTOR' in comments_["commentType"]:
#                         data_dict["Function information"][comments_['commentType']] = comments_['cofactors'][0]['name']
#                     elif 'ACTIVITY REGULATION' in comments_["commentType"]:
#                         data_dict["Function information"][comments_['commentType']] = comments_['texts'][0]['value']
#                 except KeyError as e:
#                     print(f"An error occurred: {e}")
#                     traceback.print_exc()
#         else:
#             data_dict["Function information"] = "No further information provided"
#
#         data_dict["Protein Sequence"] = data_info["sequence"]["value"]
#         data_dict["Protein ID"] = data_info["primaryAccession"]
#         data_dict["Protein Name fasta"] = data_info["uniProtkbId"]
#
#         # A protein may have multiple DNA sequences or be inserted into some large gene sequence
#         reference_dna_id = []
#         for refs_ in data_info["uniProtKBCrossReferences"]:
#             if "EMBL" == refs_["database"]:
#                 reference_dna_id.append(refs_["id"])
#         data_dict["DNA List"] = reference_dna_id
#         if "recommendedName" not in data_info["proteinDescription"]:
#             protein_name = data_info["proteinDescription"]["submissionNames"][0]['fullName']["value"]
#         else:
#             protein_name = data_info["proteinDescription"]["recommendedName"]["fullName"]["value"]
#         organism_name = data_info["organism"]["scientificName"]
#         protein_len = data_info["sequence"]["length"]
#         if data_dict["Function information"] == "No further information provided":
#             print(data_dict["Function information"])
#
#             search_info['Protein_name'] = protein_name
#             search_info['Organism'] = organism_name
#             search_info['Protein_length'] = protein_len
#             search_info['Function_information'] = data_dict["Function information"]
#         else:
#             Function_content = '\n'.join(
#                 f'    {na_}:  {content}' for na_, content in data_dict["Function information"].items())
#             search_info['Protein_name'] = protein_name
#             search_info['Organism'] = organism_name
#             search_info['Protein_length'] = protein_len
#             search_info['Function_information'] = Function_content
#
#         search_info_list.append(search_info)
#         data_dict["Protein Name"] = protein_name
#         data_dict["Organism"] = organism_name
#         data_dict["Length"] = protein_len
#         data_list.append(data_dict)
#     print("Extracted functional information", len(data_list))
#     return search_info_list, data_list





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
    when the user requests a pathogen resistance analysis When performing a multiplex PCR experiment, 
    you should call the pathogen_drug_resistance_search function to obtain the target gene sequence and SNP information; 
    When a user wants to obtain or download a gene sequence or genome sequence, you should  call the download_files to meet the user's needs;
    when the user requests a multiplex PCR experiment for pathogen species identification analysis, you should call the species_identification_search function；In addition, the language in which you reply to the user should be consistent with the language used by the user."""
                            }]
            messages = system_prompt + window_history
            search_response = get_chatgpt_response(messages,functions)
            if search_response.get("function_call"):
                print('*******************Successful activate function*******************')
                # Call the function. The JSON response may not always be valid so make sure to handle errors
                function_name = search_response["function_call"]["name"]
                available_functions = {"protein_mutation_search": protein_mutation_search,
                                    "genetic_disorder_search": genetic_disorder_search,
                                    "cancer_search": cancer_search,
                                    "download_files": download_files,
                                    "pathogen_drug_resistance_search":pathogen_drug_resistance_search,
                                    "species_identification_search":species_identification_search,
                                     "snp_primer_design":snp_primer_design,
                                    }
                function_to_call = available_functions[function_name]
                print("Function name: ",search_response["function_call"]["name"])
                print("Function input: ",search_response["function_call"]["arguments"])
                
                
                
                if search_response["function_call"]["name"] == "protein_mutation_search":
                
                    function_args = json.loads(search_response["function_call"]["arguments"])
                
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
                       #                 "arguments": search_response["function_call"]["arguments"],
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
                elif search_response["function_call"]["name"] == "download_files":
                    function_args = json.loads(search_response["function_call"]["arguments"])
                    name = function_args['name_list']
                    user_response = download_files(stage=1, instruction=name)
                    responses = user_response
                    history_conversation.append({'role': 'assistant', 'content': responses})
                    responses_dict = {'responses': responses, 'history_conversation': history_conversation}
                    return responses_dict
    
    
                elif search_response["function_call"]["name"] == "pathogen_drug_resistance_search":
                    function_args = json.loads(search_response["function_call"]["arguments"])
                    #species_name = function_args['instruction']
                    responses = pathogen_drug_resistance_search(input_data,stage=1)
                    #user_response = pathogen_drug_resistance_search(stage=1,instruction=species_name)
                    #responses = user_response
                    history_conversation.append({'role': 'assistant', 'content': responses})
                    responses_dict = {'responses': responses, 'history_conversation': history_conversation}
                    return responses_dict
    
                elif search_response["function_call"]["name"] == "snp_primer_design":
                    print('成功调用primer design')
                    function_args = json.loads(search_response["function_call"]["arguments"])
                    print('xxx',history_conversation[-1]['content'])
                    instruction = ast.literal_eval(history_conversation[-1]['content'])
                    # history_conversation = json.loads()
                    # instruction = history_conversation
                    print('xxxxinstruction',instruction)
    
                    responses = snp_primer_design(instruction)
                    history_conversation.append({'role': 'assistant', 'content': responses})
                    responses_dict = {'responses': responses, 'history_conversation': history_conversation}
                    return responses_dict
    
                elif search_response["function_call"]["name"] == "genetic_disorder_search":
                    function_args = json.loads(search_response["function_call"]["arguments"])
                    genetic_disease = function_args['genetic_disease']
                    user_response = genetic_disorder_search(stage=1,instruction=genetic_disease)
                    responses = user_response
                    history_conversation.append({'role':'assistant','content':responses})
                    responses_dict = {'responses': responses, 'history_conversation': history_conversation}
                    return responses_dict
                
                elif search_response["function_call"]["name"] == "cancer_search":
                    function_args = json.loads(search_response["function_call"]["arguments"])
                    #instruction = function_args['term']
                    
                    user_response = cancer_search(stage=1, instruction=input_data)
                    responses = user_response
                    history_conversation.append({'role': 'assistant', 'content': responses})
                    responses_dict = {'responses': responses, 'history_conversation': history_conversation}
                    # responses_dict={'target_cds_path':[],'non_target_cds_path':[],'fna_file':[fna_path],'snp_file':[snp_path],'download_path':[],'responses':responses,'history_conversation':history_conversation}
                    return responses_dict
    
                elif search_response["function_call"]["name"] =="species_identification_search":
                    function_args = json.loads(search_response["function_call"]["arguments"])
                    instruction = history_conversation#function_args['instruction']++++++++++++++++++++++++++
    
                    responses = species_identification_search(input_data,stage=1)
                    history_conversation.append({'role':'assistant','content':responses})
                    responses_dict = {'responses': responses, 'history_conversation': history_conversation}
                    # responses_dict['history_conversation'] = history_conversation
                    return responses_dict
    
                else:
                    print("error funcation calling")
    
            else:
                responses = search_response['content']
                history_conversation.append({'role':'assistant','content':responses})
                responses_dict={'没有调用函数'}
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
                
                    # print('xxxxxxxxxxxxxhistory_conversation', history_conversation)
                  #instruction = json.loads(history_conversation[-1]['content'])
                  #instruction = instruction['selected_options']
                  #print('instruction',instruction)
                  responses = genetic_disorder_search(input_data, stage)
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
            
                #instruction = json.loads(history_conversation[-1]['content'])
                #instruction = instruction['selected_options']
                #responses = pathogen_drug_resistance_search(instruction, stage)
                #history_conversation.append({'role': 'assistant', 'content': responses})
                #responses_dict = {'responses': responses, 'history_conversation': history_conversation}
                #return responses_dict
    
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
        return "Exception:"+str(e)

# 定义X-chat的对话逻辑
def primer_predict(input_data, primer_type=None,stage=1, search_prompt=""):
    history_conversation = input_data["conversation"]
    window_history = history_conversation[-50:]
    #try:
    if primer_type == None:
        system_prompt = [{'role': 'system', 'content': """You are an artificial intelligence assistant designed to help users design specific primers. 
When the user asks for follow-up work or subsequent primer design , 
you should call the snp_primer_design function to design primers;"""
                          }]
        messages = system_prompt + window_history
        search_response = get_chatgpt_response(messages, primer_functions)
        if search_response.get("function_call"):
            print('*******************Successful activate function*******************')
            # Call the function. The JSON response may not always be valid so make sure to handle errors
            function_name = search_response["function_call"]["name"]
            available_functions = {
                                   "snp_primer_design": snp_primer_design,
                                   }
            function_to_call = available_functions[function_name]
            print("Function name: ", search_response["function_call"]["name"])
            print("Function input: ", search_response["function_call"]["arguments"])


            if search_response["function_call"]["name"] == "snp_primer_design":
                print('成功调用primer design')


                responses = snp_primer_design(input_data,stage=stage)
                history_conversation.append({'role': 'assistant', 'content': responses})
                responses_dict = {'responses': responses, 'history_conversation': history_conversation}
                return responses_dict

            else:
                print("error funcation calling")

        else:
            responses = search_response['content']
            history_conversation.append({'role': 'assistant', 'content': responses})
            responses_dict = {'target_cds_path': [], 'non_target_cds_path': [], 'snp_file': [], 'fna_file': [],
                              'download_path': [], 'responses': responses, 'history_conversation': history_conversation}
            return responses_dict

    else:
          print('进入函数')   
          responses = snp_primer_design(input_data,stage)
          history_conversation.append({'role': 'assistant', 'content': responses})
          responses_dict = {'responses': responses, 'history_conversation': history_conversation}
          return responses_dict
    #except Exception as e:
     #   print("Exception:"+str(e))
     #   return "Exception:"+str(e)
        
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

executor = ProcessPoolExecutor(max_workers=4)#ThreadPoolExecutor()




# WebSocket 路由
@app.websocket("/ncbi")
async def websocket_endpoint(websocket: WebSocket):
    await websocket.accept()

    asyncio.create_task(heartbeat(websocket))
    try:
        while True:
            data = await websocket.receive_text()
            #print("接收到数据：",data)
            input_data = json.loads(data) #解析为JSON
            print("接收到数据：",input_data)
            
            input_data["instruction"] = json.loads(input_data["instruction"])
            all_history_conversation = input_data["instruction"]["conversation"]
            stage = input_data["instruction"]["stage"]


            #第一轮对话
            if stage == 1:
                loop = asyncio.get_event_loop()
                input_data["conversation"] = all_history_conversation
                responses_dict = await loop.run_in_executor(executor, predict, input_data, stage)

            #多轮对话
            else:
                #print("当前对话历史：",history_conversation)
                input_data["conversation"] = all_history_conversation
                search_type = input_data["instruction"]["search_type"]

                temp_flag = False
                if search_type != "":
                    temp_flag = True

                loop = asyncio.get_event_loop()
                if temp_flag:
                    print('xxxxx进行第二次回复')
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
            print("接收到数据：",input_data)
            
            input_data["instruction"] = json.loads(input_data["instruction"])
            all_history_conversation = input_data["instruction"]["conversation"]
            stage = input_data["instruction"]["stage"]
            

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
                primer_type = 'xxx'#input_data["instruction"]["primer_type"]
                
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





if __name__=="__main__":
    #openai api    
    openai.api_type = "azure"
    openai.api_base = "https://xmgi-chat1.openai.azure.com/"
    openai.api_version = "2023-07-01-preview"
    openai.api_key = "42cb08b9d67747d9a60c3d4722168108"

    os.system('rm log')


    #开启服务
    uvicorn.run(app, host="0.0.0.0", port=8081,timeout_keep_alive=300)


