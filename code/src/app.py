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
from concurrent.futures import ThreadPoolExecutor
import pandas as pd
import socket
# from flask import Flask, request, redirect, url_for, flash
# from werkzeug.utils import secure_filename
import threading
from fastapi.responses import HTMLResponse
from fastapi import FastAPI, File, UploadFile, HTTPException
# API请求有频率限制，TODO：多账号混合、时延机制
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
                "user_input":{
                    "type":"string",
                    "description":"All user input information"
                },
            },
            "required":["term","user_input"]
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


def snp_primer_design(instruction,parameters=[]):


    snp_primer_design_prompt = """
here is the designed primer file.
                    """
    out_path = '/reference_data/primer_result'
    dna_file_path = '/reference_data/dna_file'
    snp_file_path = '/reference_data/snp_file'
    # if search_type == "genetic_disorder_type":
    search_type = instruction['search_type']

    print('search_type',search_type)
    # if stage == 1:
    #
    # #获取引物设计所必需的参数
    # parameters = instruction['parameters']
    # parameters = []
    if search_type == 'genetic_disorder_type':
        print(f'当前结果来自：{search_type}')
        temp_dict = instruction['search_responses']['target_gene_info']

        res = pd.DataFrame.from_dict(temp_dict)

        print('res',res)
        gene_names = list(res['gene'])
        end_list = list(res['end'])
        start_list = list(res['start'])
        chr_list = list(res['chrom'])
        # 定义文件路径
        fasta_file = "/reference_data/GCF_000001405.40_GRCh38.p14_genomic_filter.fna"
        result = []
        # 打开文件并读入内容
        with open(fasta_file, 'r') as file:
            file_content = file.read().splitlines()
            file.close()
        for ii in range(len(gene_names)):
            na = gene_names[ii]
            temp_seq = get_dna_seq(chr_list[ii][3:], start_list[ii], end_list[ii],file_content)
            temp_path = os.path.join(dna_file_path,f'{na}.fasta')
            with open(temp_path,'w') as f:
                f.write(f'>{na}\n')
                f.write(temp_seq)
                f.close()

            fasta_path = temp_path
            snp_path = ''
            if parameters ==[]:
                cmd = 'python /reference_data/plasmid.py --fasta_path {}  --out_path {} --gene_name {}'.format(fasta_path, out_path,na)
            else:
                cmd = 'python /reference_data/plasmid.py --fasta_path {}  --out_path {} --gene_name {} --max_amp_len {} --min_amp_len {} --temperature {} --min_GC {} --max_GC {} --max_len {}'.format(
                    fasta_path, out_path, na,parameters[0],parameters[1],parameters[2],parameters[3],parameters[4],parameters[5])
            print('-' * 30, 'cmd: \n', cmd)
            process = os.popen(cmd)  # return file
            output = process.read()
            print(output)
            process.close()
            for ff in os.listdir(out_path):
                if ff.endswith('.csv') and na in ff:
                    result.append(os.path.join(out_path,ff))


    elif search_type == 'protein_mutation_type':
        print(f'当前结果来自：{search_type} ')
        path_list_select = instruction['search_responses']['protein_gene_seq_path']
        path_list = instruction['upload_file']
        if len(path_list_select) != len(path_list):
            print('上传文件数量与原文件数量不相等。')
            print('path_list_select',path_list_select)
            print('path_list', path_list)
        name_list = []
        for pa in path_list:
            temp_name = os.path.basename(pa)
            temp_name = temp_name.split('.')[0]
            temp_name = temp_name.replace('(','').replace(')','')
            name_list.append(temp_name)
            # cmd = 'python /reference_data/plasmid.py --fasta_path {}  --out_path {} --gene_name {}'.format(pa,
            #                                                                                               out_path, temp_name)
            if parameters ==[]:
                cmd = 'python /reference_data/plasmid.py --fasta_path {}  --out_path {} --gene_name {}'.format(pa, out_path,temp_name)
            else:
                cmd = 'python /reference_data/plasmid.py --fasta_path {}  --out_path {} --gene_name {} --max_amp_len {} --min_amp_len {} --temperature {} --min_GC {} --max_GC {} --max_len {}'.format(
                    pa, out_path, temp_name,parameters[0],parameters[1],parameters[2],parameters[3],parameters[4],parameters[5])
            print('-' * 30, 'cmd: \n', cmd)
            process = os.popen(cmd)  # return file
            output = process.read()
            print(output)
            process.close()
            time.sleep(3)
        result =[]
        for temp_name in name_list:
            for file in os.listdir(out_path):
                if temp_name in file and '.csv' in file:
                    result.append(os.path.join(out_path,file))
    elif search_type == 'pathogen_drug_resistance_type':
        print(f'当前结果来自：{search_type} ')
        res_path = instruction['search_responses']['filter_result_path']
        res_df = pd.read_csv(res_path)

        gene_list = list(set(list(res_df['Gene'])))
        print('gene_list',gene_list)
        for gen in gene_list:
            temp = res_df[res_df['Gene']==gen]
            snp_list = list(temp['Position'])
            seq_temp = list(temp['Sequence'])[0]
            snp_path = os.path.join(snp_file_path,gen+'_var.csv')
            fa_path = os.path.join(dna_file_path,gen+'_seq.fasta')
            freq = [1 for x in range(len(snp_list))]
            snp_df = pd.DataFrame({'START':snp_list,'STOP':snp_list,'FREQ':freq})
            snp_df = snp_df.drop_duplicates()
            snp_df.to_csv(snp_path,index=False)
            with open(fa_path,'w') as f:
                f.write('>'+gen+'\n')
                f.write(seq_temp)
                f.close()
            # cmd = 'python /reference_data/plasmid.py --fasta_path {}  --snp_path {} --out_path {} --gene_name {}'.format(fa_path,snp_path,out_path,gen)
            if parameters ==[]:
                cmd = 'python /reference_data/plasmid.py --fasta_path {} --snp_path {} --out_path {} --gene_name {}'.format(fa_path,snp_path, out_path,gen)
            else:
                cmd = 'python /reference_data/plasmid.py --fasta_path {} --snp_path {} --out_path {} --gene_name {} --max_amp_len {} --min_amp_len {} --temperature {} --min_GC {} --max_GC {} --max_len {}'.format(
                    fa_path,snp_path, out_path,gen,parameters[0],parameters[1],parameters[2],parameters[3],parameters[4],parameters[5])
            print('-' * 30, 'cmd: \n', cmd)
            process = os.popen(cmd)  # return file
            output = process.read()
            print(output)
            process.close()
            time.sleep(3)
        result = []
        for temp_name in gene_list:
            for file in os.listdir(out_path):
                if temp_name in file and '.csv' in file and 'primer_design' in file:
                    result.append(os.path.join(out_path, file))

    elif search_type == 'cancer_type':
        print(f'当前结果来自：{search_type}')
        res_dict = instruction['search_responses']['target_gene_info']
        res_df = pd.DataFrame.from_dict(res_dict)
        gene_list = list(set(list(res_df['Gene name'])))
        print('gene_list', gene_list)
        # 定义文件路径
        fasta_file = "/reference_data/GCF_000001405.40_GRCh38.p14_genomic_filter.fna"

        # 打开文件并读入内容
        with open(fasta_file, 'r') as file:
            file_content = file.read().splitlines()
            file.close()
        for gen in gene_list:
            temp_df = res_df[res_df['Gene name']==gen]
            snp_path = os.path.join(snp_file_path, gen + '_var.csv')
            snp_list = list(temp_df['Start position in gene'])
            freq = [1 for x in range(len(snp_list))]
            snp_df = pd.DataFrame({'START': snp_list, 'STOP': snp_list, 'FREQ': freq})
            snp_df = snp_df.drop_duplicates()
            snp_df.to_csv(snp_path, index=False)

            temp_seq = get_dna_seq(list(temp_df['Chromosome'])[0], list(temp_df['Gene start'])[0], list(temp_df['Gene end'])[0], file_content)
            temp_path = os.path.join(dna_file_path,gen+'_seq.fa')
            with open(temp_path, 'w') as f:
                f.write(f'>{gen}\n')
                f.write(temp_seq)
                f.close()

            fasta_path = temp_path
            print(f'{gen}基因序列提取完成。开始设计引物')
            # snp_path = ''
            # cmd = 'python /reference_data/plasmid.py --fasta_path {} --snp_path {} --out_path {} --gene_name {}'.format(fasta_path,snp_path,
            #                                                                                                out_path, gen)
            if parameters ==[]:
                cmd = 'python /reference_data/plasmid.py --fasta_path {} --snp_path {} --out_path {} --gene_name {}'.format(fasta_path,snp_path,out_path,gen)
            else:
                cmd = 'python /reference_data/plasmid.py --fasta_path {} --snp_path {} --out_path {} --gene_name {} --max_amp_len {} --min_amp_len {} --temperature {} --min_GC {} --max_GC {} --max_len {}'.format(
                    fasta_path,snp_path, out_path,gen,parameters[0],parameters[1],parameters[2],parameters[3],parameters[4],parameters[5])
            print('-' * 30, 'cmd: \n', cmd)
            process = os.popen(cmd)  # return file
            output = process.read()
            print(output)
            process.close()
        result = []
        for temp_name in gene_list:
            for file in os.listdir(out_path):
                if temp_name in file and '.csv' in file and 'primer_design' in file:
                    result.append(os.path.join(out_path, file))
    elif search_type == 'species_identification_type':

        print(f'当前结果来自：{search_type}')
        if type == '':
            target_path = instruction['search_responses']['target_gene_info'][0]
            nontarget_path  = instruction['search_responses']['target_gene_info'][1]
            cmd = f'python /reference_data/run_primer_design.py --target_path {target_path} --nontarget_path {nontarget_path}'
            print('-' * 30, 'cmd: \n', cmd)
            process = os.popen(cmd)  # return file
            output = process.read()
            print(output)
            process.close()
        #直接上传单个基因序列文件进行引物设计。
        elif type =='':
            files_list= instruction['search_responses']['upload']
            output_path = out_path
            result = primer3_design(files_list,output_path)


    if len(result)==0:
        responses = {'response': 'Sorry, primer design failed.', 'primer': result, 'state': 'stop',
                     "primer_type": "snp_primer_design_type",
                     'primer_design_prompt': ''}
    else:
        responses = {'response':snp_primer_design_prompt,'primer': result, 'state': 'stop', "primer_type": "snp_primer_design_type",
                     'primer_design_prompt': snp_primer_design_prompt}

    return responses

def SNP_Genotyping_search(species_name):
    #获得参考基因组文件
    if species_name == '':
        snp_info = pd.read()
    #读取分型snp数据表，确定你要做那些亚型
    print()
    #组合成文件，传给引物设计Agent


def download_files(instruction,stage):
    download_system_prompt = ''
    if stage == 1:
        stage += 1
        protein_name = instruction
        print('name', protein_name)
        protein_info, data_list = search_uniprot_one(protein_name)
        if data_list ==[]:
            response = f'No {instruction} relevant sequence was found on the NCBI website. Sorry for the inconvenience, you can change the name to get the information you want.'
            return {"response": response, "options": [], "data_list": data_list, "stage": stage,
                    "search_type": "download_type", "search_system_prompt": download_system_prompt,
                    "state": 'stop'}
        else:
            response = f'I will provide the protein information related to {instruction}'
            return {"response": response, "options": protein_info[:30], "data_list": data_list, "stage": stage,
                    "search_type": "download_type", "search_system_prompt": download_system_prompt,
                    "state": 'continue'}
    else:
        data_list = instruction  # user select
        # user_select_protien_ncbi_download=""
        aa_dna_content, download_path = user_select_protein_ncbi_download(data_list)
        response = f"Protein gene selection has been completed as needed, the details are as follows: {aa_dna_content} "
        print('结果',download_path)
        return {"response": response, "protein_gene_seq_path": download_path,
                "search_type": "download_type", "search_system_prompt": download_system_prompt,
                "state": 'stop'}




#cancer search tool
def cancer_search(instruction, stage):
    #description cancer search tool
    cancer_search_system_prompt="""The current Search Agent has entered the cancer retrieval phase and 
needs to interact with the user to obtain the genomic and corresponding drug resistance locus information for the cancer. 
The main task here is to interact with the user, guiding them to select the Gene name, Primary site, Site subtype 1, 
Primary histology, and Histology subtype 1. Based on the selected content and filtering through the Cosmic database, 
it acquires the gene name and corresponding SNP information, and also retrieves the genomic sequence information from NCBI using the gene name."""
    #cancer name --> gene name list
    if stage==1:
        stage+=1
        cancer_name = instruction
        #加载一个目录数据集，让用户选择
        print('文件是否存在',os.path.exists("/reference_data/cancer_directory.json"))
        with open("/reference_data/cancer_directory.json","r") as f:
            Tumor_catalog = json.load(f)
        f.close()
        # cancer_name_catalog = Tumor_catalog[cancer_name]
        response = "I will provide a catalog of tumors, please confirm the cancer name, Primary site, Site subtype 1, Primary histology, Histology subtype 1."
        return {"response":response,"cancer_list":Tumor_catalog,"stage":stage,"search_type":"cancer_type","search_system_prompt":cancer_search_system_prompt,"state":'continue'}
    
    #user select --> gene sequence and snp
    else:
        stage+=1
        #user_select_list = instruction
        print('instruction',instruction)
        data = pd.read_csv("/reference_data/Cosmicfilter.tsv", header=0, sep='\t')
        #文件名中的空格用下划线代替
        primary_site = instruction[0].lower().strip().replace(' ','_')
        site_subtype_1 = instruction[1].lower().strip().replace(' ','_')
        primary_histology = instruction[2].lower().strip().replace(' ','_')
        histology_subtype_1 = instruction[3].lower().strip().replace(' ','_')

        #catalog list
        filtered_data = data[
        (data['Primary site'].str.lower() == primary_site) &
        (data['Site subtype 1'].str.lower() == site_subtype_1) &
        (data['Primary histology'].str.lower() == primary_histology) &
        (data['Histology subtype 1'].str.lower() == histology_subtype_1)
        ]
        name_mutation_csv = filtered_data.loc[:,['Gene name', 'Chromosome', 'Gene start',
                                                 'Gene end', 'Start position in gene', 'End position in gene']]
        result_data = name_mutation_csv.to_dict()
        # file_path_list = []
        # file_path_mutation_dict={}
        # for gene_name, mutation_cds in name_mutation_csv:
            # fna_path, user_response = search_ncbi("genome", gene_name, f"{primary_site}_{site_subtype_1}_{primary_histology}_{histology_subtype_1}_cancer")
            
            # file_path_list.append(fna_path)
        if len(name_mutation_csv) ==0:
            print('没有找到相关的基因信息')
            stage = stage -1
            with open("/reference_data/cancer_directory.json", "r") as f:
                Tumor_catalog = json.load(f)
            f.close()
            # cancer_name_catalog = Tumor_catalog[cancer_name]
            response = "I will provide a catalog of tumors, please confirm the cancer name, Primary site, Site subtype 1, Primary histology, Histology subtype 1."
            return {"response": response, "cancer_list": Tumor_catalog, "stage": stage, "search_type": "cancer_type",
                    "search_system_prompt": cancer_search_system_prompt, "state": 'continue'}

        else:
            response ='I will provide the snp information'
                #save snp
                # target_gene_dict[gene_seq]=SNP
            return{"response":response,"target_gene_info":result_data,"search_type":"cancer_type","search_system_prompt":cancer_search_system_prompt,"state":'stop'}


#genetic disorder search tool
def genetic_disorder_search(instruction, stage):
    #description genetic disorder search tool
    genetic_disorder_system_prompt="""
    Current search agents are involved in search projects related to human genetic diseases.First, it requires user interaction to extract the entered human genetic disease name and retrieve the gene names related to the genetic disease on Omim.
And feed back the retrieved information, such as genetic disease phenotypes, related genes and other information.
Further user interaction will select the genetic disease and associated genes required for study. Finally, according to the user’s choice,
The gene sequence of the corresponding gene and the position information on the genome were downloaded from NCBI.
    """

    #genetic disorder name --> gene name list
    if stage==1:
        response = f'I will provide the name of the genetic disease and the corresponding gene name related to {instruction}'
        stage+=1
        # cancer_name = instruction
        #OMIM dataset
        #根据 disorder name来抽取gene names
        gene_result = get_disease_gene(instruction)
        if isinstance(gene_result,str):
            return {"response": gene_result, "options": gene_result, "stage": stage,
                    "search_type": "genetic_disorder_type", "search_system_prompt": genetic_disorder_system_prompt,
                    "state": 'stop'}
        else:
            #返回gene names 、 stage、type
            return {"response":response,"options":gene_result ,"stage":stage,"search_type":"genetic_disorder_type","search_system_prompt":genetic_disorder_system_prompt,"state":'continue'}
    #user select --> gene name list --> gene sequence and snp
    else:
    # elif stage==2:
        stage+=1
        target_gene_names_list = instruction
        #打开 ClinVar dataset
        gene_info_df = pd.read_csv('/reference_data/gene_info.tsv',sep='\t')
        print('target_gene_names_list',target_gene_names_list)
        # target_gene_dict = {}
        df_list =[]
        for gene_name in target_gene_names_list:
            temp_df = gene_info_df[gene_info_df['gene']==gene_name]
            df_list.append(temp_df)
        #返回目标基因在染色体上的位置信息
        target_gene_dict = pd.concat(df_list).to_dict()
        response = f'I will provide the location information of the target gene on the whole genome'
        return{"response":response,"target_gene_info":target_gene_dict,"search_type":None,"search_system_prompt":genetic_disorder_system_prompt,"state":'stop'}
    
    #user select --> snp/region --> continue primer design
    # else:
    #     select_gene_dict = instruction
    #     snp_gene_path_dict = {}
    #     for snp_, seq_ in select_gene_dict.items():
    #         #local save snp files and gene sequence fasta
    #         snp_path=""
    #         gene_path=""
    #         snp_gene_path_dict[gene_path]=snp_path
    #
    #     response=f"Genetic disorder gene and SNP selection has been completed as needed, the details are as follows: 1) 2)"
    #
    #     return{"response":response,"snp_gene_path":snp_path,"type":"genetic_disorder_type","search_system_prompt":genetic_disorder_system_prompt,"state":'stop'}

#protein engineering search tool
def protein_mutation_search(instruction, stage,user_input):
    #description protein engineering search tool
    protein_engineering_system_prompt="""The current search agent is involved in search projects related to various enzyme mutant libraries. 
First it requires user interaction to extract the entered protein name and retrieve the relevant sequence on the UniProt. 
And feed back the retrieved sequence information, such as function, species, length and other information. 
Further user interaction will select the desired protein. Finally, based on the user's selection, 
the corresponding protein's gene sequence is downloaded from NCBI. 
This sequence is then provided to the user for use in multiplex PCR methods, facilitating tNGS library construction and sequencing."""

    #protein name --> protein information
    if stage==1:
        stage+=1
        protein_name = instruction
        #判断用户是否需要上传
        # messages = []
        # messages.append({"role": 'user',
        #                  'content': f"""You need to determine the meaning of the sentence {user_input} entered by the user. If the user means the sequence file he already has, or he wants to upload the file himself, then you return the character True, otherwise you return False, and there is no need to return other information. """})
        # # mess_all = messages +instruction
        # # print('message:',mess_all)
        # second_response = openai.ChatCompletion.create(
        #     messages=messages,
        #     deployment_id="XMGI-Chat1-GPT4"
        # )
        # print(f'The modified term  reflection is:',
        #       second_response["choices"][0]["message"]['content'])
        # flag = second_response["choices"][0]["message"]['content']
        #
        # if 'True' in flag:
        #     print('可以上传')
        #     response = ''
        #     return {"response": response,  "stage": stage,"data_list": ['xxx'],
        #             "upload": "True", "upload_title": "If you already have a sequence file, you can click the button below to upload it. If not, please ignore it. ","search_type": "protein_mutation_type",
        #             "search_system_prompt": protein_engineering_system_prompt, "state": 'stop'}
        print('protein_name',protein_name)
        protein_info, data_list = search_uniprot_one(protein_name)
        print('搜索结果数量',len(data_list))
        if len(data_list) ==0:
            response = f'Sorry, the current search results for the protein {instruction} on NCBI are empty. Please adjust the name and search again.'
            return {"response": response, "options": protein_info[:20], "data_list": data_list, "stage": stage-1,
                    "upload": "True", "upload_title": " ","search_type": "protein_mutation_type",
                    "search_system_prompt": protein_engineering_system_prompt, "state": 'continue'}

        else:
            # print('xxxxxprotein_info',len(protein_info))
            response = f'I will provide the protein information related to {instruction}'
            return {"response":response,"options":protein_info[:20],"upload_title": " ","data_list":data_list,"stage":stage,"upload":"True","search_type":"protein_mutation_type","search_system_prompt":protein_engineering_system_prompt,"state":'continue'}
    
    #user select proetin number --> gene sequence --> continue primer design
    elif stage==2:
        data_list = instruction  #user select
        #user_select_protien_ncbi_download=""
        try:
            aa_dna_content, download_path = user_select_protein_ncbi_download(data_list)
        except Exception as e:
            raise Exception('文件下载失败')
        response=f"Protein gene selection has been completed as needed, the details are as follows: {aa_dna_content} "
        if download_path==[]:
            return {"response": 'Sorry, the gene sequence corresponding to the protein you selected was not found.', "protein_gene_seq_path": download_path,
                    "search_type": "protein_mutation_type", "upload":"True","upload_title": "You can manually upload the file or select a different protein.","search_system_prompt": protein_engineering_system_prompt,
                    "state": 'stop'}
        else:
            return {"response":response,"protein_gene_seq_path":download_path,
                    "search_type":"protein_mutation_type","upload":"True","upload_file":[],"upload_title": "If you want to modify the downloaded gene sequence so that it can be better used for primer design, please click the upload button below to upload your modified sequence file, and keep the number of file uploads and file downloads equal. If you do not need Modification, please ignore the upload button.",
                    "search_system_prompt":protein_engineering_system_prompt,"state":'stop'}
    else:
        new_file_list = instruction
        return {"response": response, "protein_gene_seq_path": new_file_list, "upload":"True","upload_file":[],"search_type": "protein_mutation_type",
                "search_system_prompt": protein_engineering_system_prompt, "state": 'stop'}

#species identification search tool
# def species_identification_search(genus_name):
#     #genus_name =
#     #species_identification search tool
#     species_identification_system_prompt="""The current Search Agent has entered into the pathogen species identification search.
# It requires interaction with the user to extract the inputted pathogen genus name.
# By using all species names under the genus name, it completes the download of whole-genome files for all species,
# to be used for subsequent species identification."""
#
#     #询问是group还是unqiue
#     # if stage ==1:
#     #     messages =[]
#     #     messages.append({"role": 'user',
#     #                      'content': f""" You first need to determine whether the experimental target provided by the user is a species-level name or a strain-level name. If it is a species-level name, you need to inform the user that what he provided is a species name and ask him to provide a more detailed strain name.
#     #                       and you also need to use json to return 10 common strain names of this species. """})
#     #     # mess_all = messages +instruction
#     #     # print('message:',mess_all)
#     #     second_response = openai.ChatCompletion.create(
#     #         messages=messages,
#     #         deployment_id="XMGI-Chat1-GPT4"
#     #     )
#     #     print(f'The modified term  reflection is:',
#     #           second_response["choices"][0]["message"]['content'])
#     #     user_response = second_response["choices"][0]["message"]['content']
#     #     # data_info_dict = function_to_call(instruction=second_response["choices"][0]["message"]['content'], stage=1)
#     #     # user_response = []
#     #
#     #     return {"response":user_response,"options":['Group-specific','Microbe-specific'],"choice_type":'Single',"state":'continue',"search_type":"species_identification_type"}
#     # elif stage ==2:
#     #     select_type = instruction['selected']
#     #     user_response = 'Do you want to do Group-specific or Microbe-specific primer design?'
#     #     if select_type =='Microbe-specific':
#     #         #llm判断是物种还是菌株
#     #         print('xxxx')
#     #
#     #     else:
#     #         pass
#
#     target_flie_path, non_target_path, user_response = get_specie_cds_file(genus_name)
#     response=f"species identification searching has been completed as needed"
#
#     return {"response":user_response,"target_flie_path":target_flie_path,"non_target_path":non_target_path,"search_type":"species_identification_type","search_system_prompt":species_identification_system_prompt,"state":'stop'}


# species identification search tool
def species_identification_search(instruction, stage):
    # 第一轮对话，实体提取，看用户是否说明物种、菌株和Group或Microbe实验方向，提取信息在后续search中使用
    if stage == 1:
        species_identification_dict = {
            "Experiment Direction": "species identification",
            "Experiment Type": "",
            "Genus": "",
            "Species Name": "",
            "Strain": "",
            "Gene Name": "",
            "Taget Need": "",
        }

        """
        #Examples1：Mycobacterium tuberculosis
        new_species_identification_dict = {
        "Experiment Direction": "species identification",
        "Experiment Type": "Group",
        "Genus":"Mycobacterium",
        "Species Name":"Mycobacterium tuberculosis",
        "Strain":"H37Rv",
        "Gene Name":"katG",
        "Taget Need":"",
        }
        #Examples2：Influenza A virus
        new_species_identification_dict = {
        "Experiment Direction": "species identification",
        "Experiment Type": "Group",
        "Genus":"Influenza",
        "Species Name":"Influenza A virus",
        "Strain":"infA_H1N1",
        "Gene Name":"Matrix protein 1",
        "Taget Need":"",
        }
        """

        conversation = instruction['conversation']

        # 根据上下文对话，提取实验类型、实验类型是否提供、物种名称、菌株名称
        experiment_type_prompt = """Based on the user's contextual dialogue, you need to perform information extraction. The information to be extracted includes:  
1.Experiment type (Microbe or Group)  
2.Whether the experiment type is provided (True or False, where False is only when the experiment type is not mentioned)  
3.Species Name (which the user may or may not mention; if not mentioned, it is None)  
4.Strain Name (which the user may or may not mention; if not mentioned, it is None)  
The results should be returned to me in JSON format.  JSON example is shown below:
{  
  "Experiment type": "Microbe",  
  "Experiment type provided": true,  
  "Species Name": "E. coli",  
  "Strain Name": null  
}
"""
        conversation_str = "\n\n".join([i["role"] + ':\n' + i['content'] for i in conversation])
        experiment_info_response = openai.ChatCompletion.create(
            messages=[{"role": "system", "content": experiment_type_prompt},
                      {"role": "user", "content": conversation_str}
                      ],
            deployment_id="XMGI-Chat1-GPT4"
        )
        experiment_info_json_data = extract_json_response(experiment_info_response)
        species_identification_dict["Species Name"] = experiment_info_json_data["Species Name"]
        species_identification_dict["Strain"] = experiment_info_json_data["Strain Name"]
        species_identification_dict["Experiment Type"] = experiment_info_json_data["Experiment Type"]

        if species_identification_dict["Strain"].lower() is None:
            recommended_system_prompt = """Based on the contextual content, 
the current user only describes the species content. 
The experiment requires strain-level information. Please recommend strains to the user. 
The results are returned in JSON format. If the user only said "Staphylococcus aureus", you need to provide the following JSON file: {
"recommended_strains": [
"Staphylococcus aureus subsp. aureus Rosenbach",
"Staphylococcus aureus subsp. aureus MRSA252",
"Staphylococcus aureus subsp. aureus N315",
"Staphylococcus aureus subsp. aureus Mu50",
"Staphylococcus aureus subsp. aureus USA300",
"Staphylococcus aureus subsp. aureus MSSA476",
"Staphylococcus aureus subsp. aureus Newman",
"Staphylococcus aureus subsp. aureus COL",
"Staphylococcus aureus subsp. aureus NCTC 8325"]}"""

            conversation_str = "\n\n".join([i["role"] + ':\n' + i['content'] for i in conversation])
            recommended_strains_response = openai.ChatCompletion.create(
                messages=[{"role": 'system', 'content': recommended_system_prompt},
                          {"role": "user", 'content': conversation_str}],
                deployment_id="XMGI-Chat1-GPT4"
            )
            recommended_strains_json_data = extract_json_response(recommended_strains_response)

        # 用户已明确实验方向
        if 'microbe' in species_identification_dict["Experiment Type"].lower():
            stage += 1
            response = f"You've chosen to conduct an experiment focusing on {species_identification_dict['Experiment Type']}, specifically targeting the {species_identification_dict['Species Name']} species."
            if species_identification_dict["Strain"] is None:
                response += "You are currently set to do a Microbe experiment. Please continue to select a strain from the recommendations below, or you can directly tell me what you need, and I can help you obtain it."
                response += """When you confirm the research strain, 
                we will recommend potential target genes you might need. 
                You can choose to use these target genes to enhance primer design, 
                or you can upload your own target gene files.
                Alternatively, you may opt not to use target genes, 
                in which case our AI will automatically identify target genes and design primers for you."""

                return {"response": response, "options": recommended_strains_json_data, "upload_file": True,
                        "choice_type": 'Single', "stage": stage,
                        "state": 'continue', "search_type": "species_identification_type",
                        "species_identification_dict": species_identification_dict}
            else:
                response += f"You intend to study the strain {species_identification_dict['Strain']}. "
                response += """we will recommend potential target genes you might need. 
                You can choose to use these target genes to enhance primer design, 
                or you can upload your own target gene files.
                Alternatively, you may opt not to use target genes, 
                in which case our AI will automatically identify target genes and design primers for you."""

                return {"response": response, "choice_type": 'Single', "stage": stage,
                        "state": 'continue', "search_type": "species_identification_type",
                        "species_identification_dict": species_identification_dict}

        elif 'group' in species_identification_dict["Experiment Type"].lower():
            stage += 1
            response = f"You've chosen to conduct an experiment focusing on {species_identification_dict['Experiment Type']}, specifically targeting the {species_identification_dict['Species Name']} species."

            if species_identification_dict["Strain"] is None:
                response += "You are currently set to do a Group experiment. Please continue to select multiple strains from the recommendations below, or you can directly tell me what you need, and I can help you obtain it."
                response += """When you confirm the research strain, 
                we will recommend potential target genes you might need. 
                You can choose to use these target genes to enhance primer design, 
                or you can upload your own target gene files.
                Alternatively, you may opt not to use target genes, 
                in which case our AI will automatically identify target genes and design primers for you."""

                return {"response": response, "options": recommended_strains_json_data, "upload_file": True,
                        "choice_type": 'Single', "stage": stage,
                        "state": 'continue', "search_type": "species_identification_type",
                        "species_identification_dict": species_identification_dict}
            else:
                response += f"You intend to study the strain {species_identification_dict['Strain']}. "
                response += """we will recommend potential target genes you might need. 
                You can choose to use these target genes to enhance primer design, 
                or you can upload your own target gene files.
                Alternatively, you may opt not to use target genes, 
                in which case our AI will automatically identify target genes and design primers for you."""

                return {"response": response, "choice_type": 'Single', "stage": stage,
                        "state": 'continue', "search_type": "species_identification_type",
                        "species_identification_dict": species_identification_dict}


        # 用户未明确实验方向
        else:
            user_response = 'Do you want to do group or microbe primer design?'
            # 如果用户没有说明物种名，就中断Search，state是stop，controller要进一步提问咨询。
            if species_identification_dict["Species Name"] is None:
                user_response += "And the species type you want to do is not found."
                return {"response": user_response, "options": ['group', 'microbe'], "choice_type": 'Single',
                        "stage": stage,
                        "state": 'stop', "search_type": "controller_type",
                        "species_identification_dict": species_identification_dict}
            else:
                stage += 1
                user_response += "Please continue to select multiple strains from the recommendations below, or you can directly tell me what you need, and I can help you obtain it."
                user_response += """When you confirm the type of experiment (Group/Microbe) and the research strain, 
                we will recommend potential target genes you might need. 
                You can choose to use these target genes to enhance primer design, 
                or you can upload your own target gene files.
                Alternatively, you may opt not to use target genes, 
                in which case our AI will automatically identify target genes and design primers for you."""

                return {"response": user_response, "options": ['group', 'microbe', recommended_strains_json_data],
                        "upload_file": True, "choice_type": 'Single', "stage": stage,
                        "state": 'continue', "search_type": "species_identification_type",
                        "species_identification_dict": species_identification_dict}


    # 判断用户是否需要进行靶基因确认
    elif stage == 2:
        species_identification_dict = instruction['species_identification_dict']
        # 用户未提供实验方向和starin时，上轮的选择
        if species_identification_dict["Experiment Type"] == 'null':
            species_identification_dict["Experiment Type"] = instruction['experiment_select_strain']
        if species_identification_dict["Strain"] == 'null':
            species_identification_dict["Strain"] = instruction['strain_select_strain']

        conversation = instruction['conversation']

        target_gene_recommended_system_prompt = """Context-based dialogue, please provide target genes for customer research content and return it to me in JSON format.
The JSON format is as follows:
{
    "recommended_target_genes": [
        "lacZ (beta-galactosidase gene)",
        "araC (arabinose operon regulatory protein gene)",
        "recA (recombinase A gene)",
        "rpoS (RNA polymerase sigma-S factor gene)",
        "gyrA (DNA gyrase subunit A gene)",
        "lon (ATP-dependent protease La gene)",
        "ompC (outer membrane porin C gene)",
        "phoA (alkaline phosphatase gene)",
        "uvrA (UvrABC system protein A gene)"
    ]
}"""
        conversation_str = "\n\n".join([i["role"] + ':\n' + i['content'] for i in conversation])
        target_gene_recommended_response = openai.ChatCompletion.create(
            messages=[{"role": 'system', 'content': target_gene_recommended_system_prompt},
                      {"role": "user", 'content': conversation_str}],
            deployment_id="XMGI-Chat1-GPT4"
        )
        target_gene_recommended_json_data = extract_json_response(target_gene_recommended_response)

        # 无论如何都要下载参考基因组
        strain_info = species_identification_dict["Strain"]
        for strain_ in strain_info:
            strain_search_info, strain_file_path = ncbi_search(strain_,
                                                               'genome')  # 直接基于菌株下载参考序列 #TODO：设计ncbi search genome
        stage += 1
        strain_response = f"""The current sequence search requirements have been completed. 
The following is the strain reference genome file. 
Please confirm whether it is correct. 
If correct, the primer design process will be carried out.
search results:{strain_search_info}"""

        if selespecies_identification_dict["Experiment Type"].lower() == 'group':
            user_response = strain_response + """For the current Group experiment, 
you can provide common genes of multiple strains as targets genes to improve the success rate of primer design. 
You can refer to the target gene content I recommend, or you can specify a specific target gene and provide its relevant information, 
and I will search it on NCBI. If you do not provide the target gene, I will design the primers myself. The target gene names I recommend are as follows:"""

        elif species_identification_dict["Experiment Type"].lower() == 'microbe':
            user_response = strain_response + """For current Microbe experiments, you can provide the target gene of the strain, 
which will increase the success rate of primer design. You can refer to the target gene content I recommend, 
or you can specify a specific target gene and provide its relevant information, and I will search it on NCBI. 
If you do not provide the target gene, I will design the primers myself. The target gene names I recommend are as follows:"""

        user_response += """If you choose or upload target genes, we will search for and download information based on those target genes and proceed with primer design based on them. 
        If you choose not to use target genes, we will complete the primer design based on the reference genome information of the strain you have obtained."""

        return {"response": user_response, "options": target_gene_recommended_json_data,
                "strain_path": strain_file_path, "upload_file": True, "choice_type": 'Single', "stage": stage,
                "state": 'continue', "search_type": "species_identification_type",
                "species_identification_dict": species_identification_dict}

        # 基于用户靶基因的抉择，判断后续任务：
    elif stage == 3:
        upload_file_flag = instruction['upload_file_flag']  # 用户上传则为True，我不需要走search流程

        select = instruction['select_target_gene']
        species_identification_dict = instruction['species_identification_dict']

        species_identification_dict["Genus"] = select
        conversation = instruction['conversation']

        flag_system_prompt = """From the current context dialogue, determine whether the user needs to further confirm the target gene information to design primers. You just need to reply yes or no."""
        conversation_str = "\n\n".join([i["role"] + ':\n' + i['content'] for i in conversation])
        flag_response = openai.ChatCompletion.create(
            messages=[{"role": 'system', 'content': flag_system_prompt},
                      {"role": "user", 'content': conversation_str}],
            deployment_id="XMGI-Chat1-GPT4"
        )
        flag_ = flag_response["choices"][0]["message"]['content']

        # non-target推荐与选择
        non_target_gene_recommended_system_prompt = """Based on our conversation and to ensure the efficiency of the primer design, 
I will recommend some non-target species for you. This will help you design primers with higher specificity. 
Below, I'll provide the recommendations in JSON format.Only need to return JSON recommendations, no redundant descriptions are needed. 
The JSON format is as follows:  
{    
  "recommendations": [  
      "Salmonella enterica",    
      "Shigella flexneri",    
      "Enterobacter cloacae",    
      "Klebsiella pneumoniae",    
      "Proteus mirabilis"    
    ]    
} """
        # conversation_str = "\n\n".join([i["role"]+':\n'+i['content']for i in conversation])
        non_target_gene_recommended_response = openai.ChatCompletion.create(
            messages=[{"role": 'system', 'content': non_target_gene_recommended_system_prompt},
                      {"role": "user", 'content': conversation_str}],
            deployment_id="XMGI-Chat1-GPT4"
        )
        non_target_gene_recommended_json_data = extract_json_response(non_target_gene_recommended_response)
        non_target_gene_response = """We have confirmed the target sequence information. 
        To enhance the specificity of primer design, you can select from my recommended non-target information or directly upload a non-target file. 
        If you choose not to use non-target information, we will proceed with the primer design based solely on the target sequence. 
        Please inform me of your decision or requirements. My recommendations for non-target information are as follows:"""

        if upload_file_flag == True:
            stage += 1
            user_response = "The target gene file uploaded by the user has been identified. Next, I will design primers based on the reference genome file of the strain."
            user_response += non_target_gene_response

            return {"response": user_response, "options": non_target_gene_recommended_json_data, "upload_file": True,
                    "choice_type": 'Multi', "stage": stage,
                    "state": 'continue', "search_type": "species_identification_type",
                    "species_identification_dict": species_identification_dict}

            # 用户不需要靶基因,继续问用户non-target事情
        if 'no' in flag_.lower():
            user_response = f"""Through our communication, it has been established that you do not need target genes for primer design. Next, we will need you to decide whether you require non-target genome information to optimize the specificity of the primers."""
            user_response += non_target_gene_response
            return {"response": user_response, "options": non_target_gene_recommended_json_data, "upload_file": True,
                    "stage": stage, "state": 'continue',
                    "search_type": "species_identification_type",
                    "species_identification_dict": species_identification_dict}

            # 用户需要靶基因，进一步给出检索结果，让用户可看到，继续问用户non-target事情
        else:
            Genus_name = species_identification_dict["Genus"]
            target_gene_search_info_list = []
            target_gene_search_path_list = []
            # for循环包含单个的microbe情况和多个的Group情况
            for sel in Genus_name:
                search_info, search_path = ncbi_search(sel,
                                                       'gene')  # 直接基于菌株下载参考序列 #TODO：设计ncbi search genome #这要做好target的过滤操作，根据description过滤 "Mycobacterium tuberculosis H37Rv" AND katG [Gene Name]  #这里要注意Gene Name的正确性，要让GPT4check，检索不到要反馈用户错误。
                target_gene_search_info_list.append(search_info)
                target_gene_search_path_list.append(search_path)
            stage += 1
            user_response = f"""Based on our conversation, 
we've looked up the target gene information for {Genus_name} and found the following sequence results for you:{target_gene_search_info_list}"""

            user_response += "Through our communication, it has been established that you do not need target genes for primer design. Next, we will need you to decide whether you require non-target genome information to optimize the specificity of the primers."
            user_response += non_target_gene_response

            return {"response": user_response, "target_gene_search_results": target_gene_search_path_list,
                    "optinons": non_target_gene_recommended_json_data, "choice_type": 'Multi', "stage": stage,
                    "state": 'continue', "search_type": "species_identification_type",
                    "species_identification_dict": species_identification_dict}

            # 判断non-target是否需要，并给出推荐
    elif stage == 4:
        conversation = instruction['conversation']
        species_identification_dict = instruction['species_identification_dict']
        user_select = instruction['non_target_select']

        # 用户没有选择
        if user_select == []:
            non_target_extract_system_prompt = """ Based on our conversation and the context, 
I will extract the non-target information User interested in. If User haven't specified any non-targets you want to focus on, 
I will return 'no'. If User have clearly mentioned User desired non-targets, I will provide them to you in JSON format.   
"""
            conversation_str = "\n\n".join([i["role"] + ':\n' + i['content'] for i in conversation])
            non_target_extract_response = openai.ChatCompletion.create(
                messages=[{"role": 'system', 'content': non_target_extract_system_prompt},
                          {"role": "user", 'content': conversation_str}],
                deployment_id="XMGI-Chat1-GPT4"
            )

            non_target_extract_json_data = extract_json_response(non_target_extract_response)
            stage += 1

            non_target = non_target_extract_json_data["non-target information"]

            if non_target == 'no':
                user_response = """You haven't provided any non-target information, so we will proceed with primer design based solely on the target information."""
                user_response += """Thus far, through our communication, we have confirmed the information regarding target . 
                Next, I will inform the Controller to communicate with you to determine the next steps related to primer design."""
                return {"response": user_response, "choice_type": 'Single', "stage": stage, "state": 'stop',
                        "search_type": "species_identification_type",
                        "species_identification_dict": species_identification_dict}

            else:
                non_target_search_info_list = []
                for non_tar_ in non_target:
                    search_info = search_ncbi(non_tar_, "genome")  # 检索物种参考基因组
                    non_target_search_info_list.append(search_info)

                user_response = """Based on the non-target information you provided, I have downloaded the following reference genomes."""
                user_response += """Thus far, through our communication, we have confirmed the information regarding both target and non-target. 
                Next, I will inform the Controller to communicate with you to determine the next steps related to primer design."""

                return {"response": user_response, "non_target_results": non_target_search_info_list,
                        "choice_type": 'Multi', "stage": stage,
                        "state": 'top', "search_type": "species_identification_type",
                        "species_identification_dict": species_identification_dict}

                # 用户在选择栏目中输入no
        elif user_select[0] == 'no':
            stage += 1
            user_response = """You haven't provided any non-target information, so we will proceed with primer design based solely on the target information."""
            user_response += """Thus far, through our communication, we have confirmed the information regarding target . 
            Next, I will inform the Controller to communicate with you to determine the next steps related to primer design."""
            return {"response": user_response, "choice_type": 'Single', "stage": stage, "state": 'stop',
                    "search_type": "species_identification_type",
                    "species_identification_dict": species_identification_dict}

            # 用户选择non-target
        else:
            stage += 1
            non_target_search_info_list = []
            non_target_path_list = []
            for non_tar_ in user_select:
                search_info, non_target_path = search_ncbi(non_tar_, "genome")  # 检索物种参考基因组
                non_target_search_info_list.append(search_info)
                non_target_path_list.append(non_target_path)
            user_response = f"""Based on the non-target information you selected, I have downloaded the following reference genomes:{non_target_search_info_list} so we will proceed with primer design based solely on the target information."""
            user_response += """Thus far, through our communication, we have confirmed the information regarding both target and non-target. 
                Next, I will inform the Controller to communicate with you to determine the next steps related to primer design."""

            return {"response": user_response, "non_target_results": non_target_path_list, "choice_type": 'Single',
                    "stage": stage, "state": 'stop', "search_type": "species_identification_type",
                    "species_identification_dict": species_identification_dict}

#增加药品名称，增加可选择SNP
def pathogen_drug_resistance_search(instruction,stage):

    pathogen_drug_resistance_prompt = """ 
    Currently, Search Agent has entered pathogen resistance search. It requires interaction with the user to obtain the name of the pathogen species, and by querying the pathogen library, it can find the SNP and sequence information of the corresponding species for subsequent design of primers for drug resistance analysis.
    """
    if stage==1:
        stage+=1
        response = f'I will provide resistance information and site information related to species {instruction}'
        species_name = instruction.lower()
        dfpatho = pd.read_csv('/reference_data/pathodrugdb.csv', header=0, sep=',')
        patholist = list(set(list(dfpatho['Species'])))
        patholist = [x.lower() for x in patholist]
        if species_name not in patholist:
            result = f"Warning: Pathogen name '{species_name}' is not in the predefined list of pathogens.Currently we only support drug resistance analysis of the following species {patholist}"
            return {"response":result,'options': [],
                "search_type": "pathogen_drug_resistance_type", "stage":stage,"search_system_prompt": pathogen_drug_resistance_prompt,
                "state": 'stop'}

        dfresult = dfpatho[dfpatho['Species'].str.lower() == species_name]
        print('dfresult',dfresult)
        df_unique = dfresult.drop_duplicates(subset=['Drug'])

        gen_list = list(df_unique['Gene'])
        pos_list = list(df_unique['Position'])
        Drug_list = list(df_unique['Drug'])

        result = []
        for ii in range(len(gen_list)):
            temp_dict ={}
            temp_dict['gene'] = gen_list[ii]
            temp_dict['position'] = pos_list[ii]
            temp_dict['Drug'] = Drug_list[ii]
            result.append(temp_dict)
        return {"response":response,"options": result, "species_name":species_name,"search_type": "pathogen_drug_resistance_type", "stage": stage,
                "search_system_prompt": pathogen_drug_resistance_prompt,"state": 'continue'}
    else:
        species_name = instruction[-1]
        drug_name = instruction[:-1]
        response = 'Based on your final choice, I give the following DNA sequence and snp information'
        dfpatho = pd.read_csv('/reference_data/pathodrugdb.csv', header=0, sep=',')
        dfresult = dfpatho[dfpatho['Species'].str.lower() == species_name]
        print('drug_name',drug_name)
        dfresult = dfresult[dfresult['Drug'].isin(drug_name)]
        print('xxxx',dfresult)
        save_path = '/reference_data/drug_filter_result.csv'
        dfresult.to_csv(save_path, index=False)
        if species_name == 'Mycobacterium tuberculosis':
            MTBseq = '/reference_data/MTB.fa'
            return {"response": response,"filter_result_path":save_path,'sequences':MTBseq,
                "search_type": "pathogen_drug_resistance_type", "stage":stage,"search_system_prompt": pathogen_drug_resistance_prompt,
                "state": 'stop'}
        else:
            return {"response": response,"filter_result_path":save_path,'sequences':'',"search_type": "pathogen_drug_resistance_type",
                    "stage":stage,"search_system_prompt": pathogen_drug_resistance_prompt,"state": 'stop'}

#针对单个基因序列文件设计引物
def primer3_design(fa_path_list,output_path):
    result = []
    for fa_path in fa_path_list:
        seqs = ''
        with open(fa_path,'r') as f:
            lines = f.readlines()
            for ll in lines:
                if '>' in ll:
                    continue
                else:
                    seqs =seqs + ll
        i_length = len(seqs)
        filename = os.path.basename(fa_path).split('.')[0]
        primers = primer3.bindings.design_primers(
                seq_args={
                    "SEQUENCE_ID": filename,
                    "SEQUENCE_TEMPLATE": seqs,
                    "SEQUENCE_INCLUDED_REGION": [0, int(i_length)],
                },
                global_args={
                    "PRIMER_OPT_SIZE": 20,
                    "PRIMER_MIN_SIZE": 18,
                    "PRIMER_MAX_SIZE": 22,
                    "PRIMER_OPT_TM": 60.0,
                    "PRIMER_MIN_TM": 58.0,
                    "PRIMER_MAX_TM": 63.0,
                    "PRIMER_PAIR_MAX_DIFF_TM": 2.0,
                    "PRIMER_MIN_GC": 40.0,
                    "PRIMER_MAX_GC": 60.0,
                    "PRIMER_MAX_POLY_X": 3,
                    "PRIMER_PICK_RIGHT_PRIMER": 1,
                    "PRIMER_PICK_LEFT_PRIMER": 1,
                    "PRIMER_PICK_INTERNAL_OLIGO": 0,
                    "PRIMER_PRODUCT_SIZE_RANGE": "75-150",
                    "PRIMER_GC_CLAMP": 1,
                    "PRIMER_NUM_RETURN": 4,
                },
            )
                    # Create a name specific to the bacterium and CDS name
        filenameTSV = filename + ".tsv"
        # Save the dictionary as a .tsv file
        res = os.path.join(output_path, filenameTSV)
        with open(res, "w") as f:
            for key in primers.keys():
                f.write("%s\t %s\n" % (key, primers[key]))
        result.append(res)
    return result

def get_nucleotide_info(name, aa_seq):
    # 搜索nucleotide数据库
    query = f"{name}"
    handle = Entrez.esearch(db="nucleotide", term=query, retmax=100000)
    record = Entrez.read(handle)
    handle.close()
    nucleotide_ids = record["IdList"]
    print('nucleotide_ids',nucleotide_ids)

    # 根据uniport提供的序列和名称来判断，取出真正的DNA序列
    for nucleotide_id in nucleotide_ids:
        handle = Entrez.efetch(db="nucleotide", id=nucleotide_id,  retmode="xml")
        records = Entrez.read(handle)
        handle.close()
        flag=False
        for record_info in records[0]["GBSeq_feature-table"]:
            if "CDS" == record_info["GBFeature_key"]:
                product_found = False
                for quals in record_info["GBFeature_quals"]:
                    if quals["GBQualifier_name"] == "translation":
                        pro_seq_d = quals["GBQualifier_value"]
                        if pro_seq_d == aa_seq:
                            loc_split = record_info["GBFeature_location"].split('(')[-1].split(')')[0].split(',')
                            if len(loc_split)==1:
                                print('GBFeature_location',record_info["GBFeature_location"])
                                loc = [int(x.replace('<','').replace('>',''))-1 for x in record_info["GBFeature_location"].split('(')[-1].split(')')[0].split('..')]
                                print(loc)
                                dna_seq = records[0]["GBSeq_sequence"][loc[0]:loc[1]]
                            else:
                                dna_seq=''
                                print(loc_split)
                                loc = [x.split('..') for x in loc_split]
                                print(loc)
                                for loc_ in loc:
                                    print(loc_)
                                    dna_seq += records[0]["GBSeq_sequence"][int(loc_[0]):int(loc_[1])]
                            dna_define = records[0]["GBSeq_definition"]
                            dna_name = records[0]["GBSeq_accession-version"]
                            flag=True
                            break
                        else:
                            dna_define = 'None'
                            dna_seq = 'None'
                            dna_name = 'None'
                            flag=False
                    
                if flag==True:
                    break
            else:
                dna_define = 'None'
                dna_seq = 'None'
                dna_name = 'None'
                flag=False


    return dna_seq, dna_define, dna_name

def get_disease_gene(disease_name):
    path = '/reference_data/pheno2gene.csv'
    gene_df = pd.read_csv(path)
    disease_name = disease_name.replace('"', "'")
    disease_gene_list = []
    gene_list = list(gene_df['HGNC'])
    phenotype_list = list(gene_df['Phenotype'])
    for ii in range(len(gene_list)):
        disease_genedict = {}
        if disease_name.lower() in phenotype_list[ii].lower():
            temp = phenotype_list[ii].replace('?', '').replace(':', '').replace('{', '').replace('}', '')
            disease_genedict['disease'] = temp
            disease_genedict['gene'] = gene_list[ii]
            disease_gene_list.append(disease_genedict)
    print('disease_gene_list',disease_gene_list)
    if len(disease_gene_list) == 0:
        return f'sorry, No causative gene related to {disease_name} has been found,You can select other genetic diseases for analysis.'
    else:
        return disease_gene_list

# def search_ncbi( term, work,db="genome",):
#     organism_name=term
#     print('User Input:',organism_name)
#     query = f"{organism_name}[Organism]"
#     handle = Entrez.esearch(db="genome", term=query, retmax=100000)
#     record = Entrez.read(handle)
#     handle.close()
#     genome_ids = record["IdList"]
#     name = []
#     for id in genome_ids:
#         esummary_handle = Entrez.esummary(db="genome", id=id)
#         esummary_record = Entrez.read(esummary_handle)
#         esummary_handle.close()
#         # 打印每个记录的摘要信息
#         assembly_name = esummary_record[0]["Assembly_Name"]
#         if assembly_name!='':
#             name.append(assembly_name)
#
#     for query in name:
#         handle = Entrez.esearch(db="assembly", term=query, retmax=100000)
#         record = Entrez.read(handle)
#         handle.close()
#         assembly_ids = record["IdList"]
#
#     os.makedirs(f"/reference_data/{work}}", exist_ok=True)
#     for assembly_id in assembly_ids:
#         handle = Entrez.efetch(db="assembly", id=assembly_id, rettype="docsum", retmode="xml")
#         records = Entrez.read(handle)
#         handle.close()
#         assembly_accessions = [doc['AssemblyAccession'] for doc in records['DocumentSummarySet']['DocumentSummary']]
#
#         #download
#         path = f"/reference_data/snp_experiment/{term.replace(' ','_')}"
#         os.makedirs(path, exist_ok=True)
#         print(f"datasets download genome accession {assembly_accessions[0]} --include gff3,rna,cds,protein,genome,seq-report --filename {path}/{assembly_accessions[0]}.zip")
#         url = f"datasets download genome accession {assembly_accessions[0]} --include gff3,rna,cds,protein,genome,seq-report --filename {path}/{assembly_accessions[0]}.zip"
#         try:
#             os.system(url)
#         except:
#             print('*****************download error, download again*****************')
#             os.system(f'rm {path}/{assembly_accessions[0]}.zip')
#             os.system(url)
#
#         file_path = f"{path}/{assembly_accessions[0]}.zip"
#         with zipfile.ZipFile(file_path, 'r') as zip_ref:
#             zip_ref.extractall(f"{path}/")
#
#         target_file = glob.glob(f"{path}/ncbi_dataset/data/G*")[0]
#         target_file = target_file+'/cds_from_genomic.fna'
#         file_to_copy = f"/reference_data/{work}/{term.replace(' ','_')}.fna"
#         print(f"Source file: {target_file}")
#         print(f"File exists: {os.path.exists(target_file)}")
#         print(f"Destination file: {file_to_copy}")
#         shutil.copy(target_file, file_to_copy)
#         print(f"finish{target_file} --> {file_to_copy}")
#
#     fna_path = file_to_copy
#     #snp_path = "/reference_data/snp_file/snp_file.xlsx"
#     user_response = f"已下载{assembly_accessions[0]}参考序列至/Reference_data文件夹下"
#     return fna_path, user_response, #snp_path,
def scrape_specie_name(specie_name):
    # 指定搜索的URL
    url = f"https://lpsn.dsmz.de/search?word={specie_name}"
    name_list = []
    # 发送请求
    response = requests.get(url)
    if response.status_code == 200:
        # 解析HTML内容
        soup = BeautifulSoup(response.text, 'html.parser')

        # 找到包含结果的HTML元素
        results = soup.find_all("a", class_="color-species")

        # 提取结果并保存到列表中
        species_results = [result.text for result in results]

        # 打印结果
        for species in species_results:
            if species.strip() == 'Species':
                pass
            elif '"' in species:
                name_list.append(species[1:-1].lower())
            else:
                name_list.append(species.strip().lower())
#             print(species)
    else:
        print("请求失败")
    name_list = list(set(name_list))
    return name_list


def get_specie_cds_file(genus_name):
    name_list = scrape_specie_name(genus_name)[:3]
    os.makedirs('/reference_data/target', exist_ok=True)
    os.makedirs('/reference_data/target_zip_file', exist_ok=True)
    for name in name_list:
        path = f'/reference_data/target_zip_file/{name.replace(" ","_")}'
        os.makedirs(path, exist_ok=True)
        search_ncbi_cds('genome',name, path)

    fail_name =[]
    cds_path = "/reference_data/target/"
    temp_name = os.listdir(cds_path)
    for na in name_list:
        if f"{na.replace(' ','_')}_cds_from_genomic.fna" not in temp_name:
            fail_name.append(na)
    print('fail_name',fail_name)
    for second in fail_name:
        temp = second.strip().split(' ')
        if len(temp) > 2:
            continue
        new_name = temp[0]+' tuberculosis variant '+temp[1]
        search_ncbi_cds('genome',new_name, path)
    
    non_target_path = "/reference_data/nontarget/"
    target_list = os.listdir(cds_path)
    target_flie_path =[]
    for tar in target_list:
        target_flie_path.append(os.path.join(cds_path,tar))

    user_response = f"There are {len(name_list)} species in the genus {genus_name} that have been found. The genus {genus_name} that failed to download has {len(fail_name)} species, namely {', '.join(fail_name)}"

    return target_flie_path, non_target_path, user_response

def search_ncbi_cds(db, term, path):
    print('User input：',term)
    query = f"{term}[Organism]"
    handle = Entrez.esearch(db="genome", term=query, retmax=100000)
    record = Entrez.read(handle)
    handle.close()
    genome_ids = record["IdList"]
    name = []

    for id in genome_ids:
        esummary_handle = Entrez.esummary(db="genome", id=id)
        esummary_record = Entrez.read(esummary_handle)
        esummary_handle.close()
        # 打印每个记录的摘要信息
        assembly_name = esummary_record[0]["Assembly_Name"]

        if assembly_name!='':
            name.append(assembly_name)

    if name==[]:
        return 

    for query in name:
        handle = Entrez.esearch(db="assembly", term=term + ' ' + query, retmax=100000)
        record = Entrez.read(handle)
        handle.close()
        assembly_ids = record["IdList"]

    for assembly_id in assembly_ids:

        handle = Entrez.efetch(db="assembly", id=assembly_id, rettype="docsum", retmode="xml")
        records = Entrez.read(handle)
        #当多个版本的数据时，选择全基因组的版本
        if len(assembly_ids)>1:
            if 'Complete Genome' not in records['DocumentSummarySet']['DocumentSummary'][0]['AssemblyStatus']:
                continue

        handle.close()
        assembly_accessions = [doc['AssemblyAccession'] for doc in records['DocumentSummarySet']['DocumentSummary']]
        #print('----------------------',assembly_accessions)
        print(f"datasets download genome accession {assembly_accessions[0]} --include gff3,rna,cds,protein,genome,seq-report --filename {path}/{term.replace(' ','_')}.zip")
        url = f"datasets download genome accession {assembly_accessions[0]} --include gff3,rna,cds,protein,genome,seq-report --filename {path}/{term.replace(' ','_')}.zip"
        try:
            os.system(url)
        except:
            print('*****************download error, download again*****************')
            os.system(f'rm {path}/{term.replace(" ","_")}.zip')
            os.system(url)

        file_path = f"{path}/{term.replace(' ','_')}.zip"
        with zipfile.ZipFile(file_path, 'r') as zip_ref:
            zip_ref.extractall(f"{path}/")

        target_file_path = glob.glob(f"{path}/ncbi_dataset/data/G*")[0]
        target_file = target_file_path+'/cds_from_genomic.fna'

        #file_to_copy = f"./non_target_cds/{term.replace(' ','_')}_cds_from_genomic.fna"
        file_to_copy = f"/reference_data/target/{term.replace(' ','_')}_cds_from_genomic.fna"
        print(f"Source file: {target_file}")
        print(f"File exists: {os.path.exists(target_file)}")
        print(f"Destination file: {file_to_copy}")
        shutil.copy(target_file, file_to_copy)
        print(f"finish{target_file} --> {file_to_copy}")


def get_dna_seq(chr_num, start, end,file_content):
    chr_name = f'chromosome {chr_num}'
    # fasta文件中的染色体名称
    all_chr_name = 'Homo sapiens ' + chr_name + ', GRCh38.p14 Primary Assembly'
    print(all_chr_name)
    # 用于存储提取的内容
    extracted_content = ""

    # 是否找到了chr
    found_chr = False

    # 在内存中搜索包含"all_chr_name"的行并提取内容
    for i, line in enumerate(file_content):
        if not found_chr and all_chr_name in line:
            found_chr = True
        elif found_chr:
            if line.startswith('>'):
                # 如果是新的标题行，则停止提取
                break
            extracted_content += line
    temp_seq = extracted_content[int(start):int(end)].upper()
    return temp_seq

#Protein first-level search: continue searching according to user needs, returning protein name, species introduction and protein length
def search_uniprot_one(term):
    get_info_url = f"http://rest.uniprot.org/uniprotkb/search?query=%28protein_name%3A{term}%29+AND+%28reviewed%3Atrue%29&format=json"
    # get_info_url = f"http://rest.uniprot.org/uniprotkb/search?query={term}&format=json"
    response = requests.get(get_info_url)
    search_info_list = []
    data_list = []
    if response.status_code != 200:
        print("Error:", response.status_code)
        return search_info_list, data_list
        exit()

    # parse JSON
    data = response.json()


    i = 0

    for data_info in data['results']:
        i += 1
        data_dict = {}
        search_info = {}
        data_dict["Function information"] = {}

        if 'comments' in list(data_info.keys()):
            for comments_ in data_info["comments"]:
                try:
                    if 'FUNCTION' in comments_["commentType"]:
                        data_dict["Function information"][comments_['commentType']] = comments_['texts'][0]['value']
                    elif 'CATALYTIC ACTIVITY' in comments_["commentType"]:
                        data_dict["Function information"][comments_['commentType']] = comments_['reaction']['name']
                    elif 'COFACTOR' in comments_["commentType"]:
                        data_dict["Function information"][comments_['commentType']] = comments_['cofactors'][0]['name']
                    elif 'ACTIVITY REGULATION' in comments_["commentType"]:
                        data_dict["Function information"][comments_['commentType']] = comments_['texts'][0]['value']
                except KeyError as e:
                    print(f"An error occurred: {e}")
                    traceback.print_exc()
        else:
            data_dict["Function information"] = "No further information provided"

        data_dict["Protein Sequence"] = data_info["sequence"]["value"]
        data_dict["Protein ID"] = data_info["primaryAccession"]
        data_dict["Protein Name fasta"] = data_info["uniProtkbId"]

        # A protein may have multiple DNA sequences or be inserted into some large gene sequence
        reference_dna_id = []
        for refs_ in data_info["uniProtKBCrossReferences"]:
            if "EMBL" == refs_["database"]:
                reference_dna_id.append(refs_["id"])
        data_dict["DNA List"] = reference_dna_id
        if "recommendedName" not in data_info["proteinDescription"]:
            protein_name = data_info["proteinDescription"]["submissionNames"][0]['fullName']["value"]
        else:
            protein_name = data_info["proteinDescription"]["recommendedName"]["fullName"]["value"]
        organism_name = data_info["organism"]["scientificName"]
        protein_len = data_info["sequence"]["length"]
        if data_dict["Function information"] == "No further information provided":
            print(data_dict["Function information"])

            search_info['Protein_name'] = protein_name
            search_info['Organism'] = organism_name
            search_info['Protein_length'] = protein_len
            search_info['Function_information'] = data_dict["Function information"]
        else:
            Function_content = '\n'.join(
                f'    {na_}:  {content}' for na_, content in data_dict["Function information"].items())
            search_info['Protein_name'] = protein_name
            search_info['Organism'] = organism_name
            search_info['Protein_length'] = protein_len
            search_info['Function_information'] = Function_content

        search_info_list.append(search_info)
        data_dict["Protein Name"] = protein_name
        data_dict["Organism"] = organism_name
        data_dict["Length"] = protein_len
        data_list.append(data_dict)
    print("Extracted functional information", len(data_list))
    return search_info_list, data_list


def ncbi_search(term, db_name):
    # 参考基因组检索
    if db_name == "genome":
        handle = Entrez.esearch(db=db_name, term=term, retmax=10)
        record = Entrez.read(handle)
        print(record)
        ids = record["IdList"]
        for id in ids:
            esummary_handle = Entrez.esummary(db="genome", id=id)
            esummary_record = Entrez.read(esummary_handle)
            esummary_handle.close()

            # 打印每个记录的摘要信息
            assembly_name = esummary_record[0]["Assembly_Name"]
            handle = Entrez.esearch(db="assembly", term=assembly_name, retmax=100000)
            record = Entrez.read(handle)
            handle.close()
            assembly_ids = record["IdList"]
            handle = Entrez.efetch(db="assembly", id=assembly_ids, rettype="docsum", retmode="xml")
            records = Entrez.read(handle)
            handle.close()
            assembly_accessions = [doc['AssemblyAccession'] for doc in records['DocumentSummarySet']['DocumentSummary']]
            print(assembly_accessions)
            # exit()
            url = f"datasets download genome accession {assembly_accessions[0]} --include gff3,rna,cds,protein,genome,seq-report --filename "  # {path}/{term.replace(' ','_')}.zip"
            print(url)
            exit()
            try:
                os.system(url)
            except:
                print('*****************download error, download again*****************')
                os.system(f'rm {path}/{term.replace(" ", "_")}.zip')
                os.system(url)

            file_path = f"{path}/{term.replace(' ', '_')}.zip"
            with zipfile.ZipFile(file_path, 'r') as zip_ref:
                zip_ref.extractall(f"{path}/")

            target_file_path = glob.glob(f"{path}/ncbi_dataset/data/G*")[0]
            target_file = target_file_path + '/cds_from_genomic.fna'

            # file_to_copy = f"./non_target_cds/{term.replace(' ','_')}_cds_from_genomic.fna"
            file_to_copy = f"/reference_data/target/{term.replace(' ', '_')}_cds_from_genomic.fna"
            print(f"Source file: {target_file}")
            print(f"File exists: {os.path.exists(target_file)}")
            print(f"Destination file: {file_to_copy}")
            shutil.copy(target_file, file_to_copy)
            print(f"finish{target_file} --> {file_to_copy}")

        return file_to_copy  # 返回下载好的数据位置

    # 靶基因检索
    elif db_name == "gene":
        term_str = f"""({term[0]}) AND {term[1]} [Text Word]"""

        # 先下载物种参考基因组
        handle = Entrez.esearch(db="genome", term=term[0], retmax=10)
        record = Entrez.read(handle)
        ids = record["IdList"]
        for id in ids:
            esummary_handle = Entrez.esummary(db="genome", id=id)
            esummary_record = Entrez.read(esummary_handle)
            esummary_handle.close()

            # 打印每个记录的摘要信息
            assembly_name = esummary_record[0]["Assembly_Name"]
            handle = Entrez.esearch(db="assembly", term=assembly_name, retmax=100000)
            record = Entrez.read(handle)
            handle.close()
            assembly_ids = record["IdList"]
            handle = Entrez.efetch(db="assembly", id=assembly_ids, rettype="docsum", retmode="xml")
            records = Entrez.read(handle)
            handle.close()
            assembly_accessions = [doc['AssemblyAccession'] for doc in records['DocumentSummarySet']['DocumentSummary']]
            print(assembly_accessions)
            exit()
            url = f"datasets download genome accession {assembly_accessions[0]} --include gff3,rna,cds,protein,genome,seq-report --filename {path}/{term.replace(' ', '_')}.zip"
            try:
                os.system(url)
            except:
                print('*****************download error, download again*****************')
                os.system(f'rm {path}/{term.replace(" ", "_")}.zip')
                os.system(url)

            file_path = f"{path}/{term.replace(' ', '_')}.zip"
            with zipfile.ZipFile(file_path, 'r') as zip_ref:
                zip_ref.extractall(f"{path}/")

            target_file_path = glob.glob(f"{path}/ncbi_dataset/data/G*")[0]
            target_file = target_file_path + '/cds_from_genomic.fna'

            # file_to_copy = f"./non_target_cds/{term.replace(' ','_')}_cds_from_genomic.fna"
            file_to_copy = f"/reference_data/target/{term.replace(' ', '_')}_cds_from_genomic.fna"
            print(f"Source file: {target_file}")
            print(f"File exists: {os.path.exists(target_file)}")
            print(f"Destination file: {file_to_copy}")
            shutil.copy(target_file, file_to_copy)
            print(f"finish{target_file} --> {file_to_copy}")

        # file_to_copy

        # kaishi1
        handle = Entrez.esearch(db=db_name, term=term_str, retmax=10)
        record = Entrez.read(handle)
        print(record)
        ids = record["IdList"]
        print(len(ids))
        for taxonomy_id in ids:
            esummary_handle = Entrez.esummary(db="gene", id=taxonomy_id)
            esummary_record = Entrez.read(esummary_handle)
            esummary_handle.close()
            # print(esummary_record)
            Name = esummary_record['DocumentSummarySet']['DocumentSummary'][0]['Name']
            if term[1] not in Name:  # 确保用户需要的靶基因在其中
                continue
            print(esummary_record)
            Description = esummary_record['DocumentSummarySet']['DocumentSummary'][0]['Description']
            Assembly_accessions = esummary_record['DocumentSummarySet']['DocumentSummary'][0]['GenomicInfo'][0][
                'ChrAccVer']
            Start_ = esummary_record['DocumentSummarySet']['DocumentSummary'][0]['GenomicInfo'][0]['ChrStart']
            End_ = esummary_record['DocumentSummarySet']['DocumentSummary'][0]['GenomicInfo'][0]['ChrStop']
            print(Name)
            print(Description)
            print(Assembly_accessions)
            print(Start_)
            print(End_)
            continue
        exit()
        pass


def extract_json_response(response):
    # 从response中提取JSON内容
    experiment_info_dict = response["choices"][0]["message"]['content']
    match = re.search(r'```json\n(.*?)\n```', experiment_info_dict, re.DOTALL)
    if match:
        experiment_info_json_str = match.group(1)
        experiment_info_json_data = json.loads(experiment_info_json_str)

    return experiment_info_json_data


# Protein secondary search: User specifies the protein, and NCBI downloads its reference gene sequence
def user_select_protein_ncbi_download(data_list):
    # print(human_decision)
    aa_dna_list = []
    for  idx in range(len(data_list)):
        pro_seq = data_list[idx]["Protein Sequence"]
        pro_id = data_list[idx]["Protein ID"]
        pro_name = data_list[idx]["Protein Name fasta"]
        for dna_id in data_list[idx]["DNA List"]:
            dna_seq, dna_define, dna_name = get_nucleotide_info(dna_id, pro_seq)
            if dna_seq != 'None':
                aa_dna_list.append((data_list[idx]["Protein Name"], data_list[idx]["Organism"], dna_name,
                                    dna_define, dna_seq, pro_id, pro_name, pro_seq))
                break

    # 一方面要前端展示，一方面要后端保存供给引物设计
    aa_dna_content = 'Based on your final choice, I give the following DNA sequence and amino acid sequence information:\n'
    download_path = []
    for idx, aa_dna_info in enumerate(aa_dna_list):
        # aa_dna_content += f"[{idx + 1}] {aa_dna_info[0]}protein({aa_dna_info[1]})：\n"
        # aa_dna_content += f"\>sp|{aa_dna_info[-3]}|{aa_dna_info[-2]}\n{aa_dna_info[-1]}\n\n"
        # aa_dna_content += f"\>{aa_dna_info[2]} {aa_dna_info[3]}\n{aa_dna_info[4]}\n\n"
        to_path = '/reference_data/dna_file'
        ss = aa_dna_info[0].strip().replace(' ', '_').replace('(', '').replace(')', '').replace('/','')
        print('ss',ss)
        temp_path = os.path.join(to_path, ss+'.fasta')
        download_path.append(temp_path)
        with open(temp_path, 'w') as f:
            f.write(f">{aa_dna_info[2]} {aa_dna_info[3]}\n{aa_dna_info[4]}")
        f.close()

    return aa_dna_content, download_path
# #Protein secondary search: User specifies the protein, and NCBI downloads its reference gene sequence
# def user_select_protien_ncbi_download(human_decision, data_list):
#     print(human_decision)
#     aa_dna_list=[]
#     for i, idx in enumerate(human_decision):
#         pro_seq = data_list[idx-1]["Protein Sequence"]
#         pro_id = data_list[idx-1]["Protein ID"]
#         pro_name = data_list[idx-1]["Protein Name fasta"]
#         for dna_id in data_list[idx-1]["DNA List"]:
#             dna_seq, dna_define, dna_name = get_nucleotide_info(dna_id, pro_seq)
#             if dna_seq != 'None':
#                 aa_dna_list.append((data_list[idx-1]["Protein Name"], data_list[idx-1]["Organism"],dna_name, dna_define, dna_seq, pro_id, pro_name,pro_seq))
#                 break
#
#     #一方面要前端展示，一方面要后端保存供给引物设计
#     aa_dna_content = 'Based on your final choice, I give the following DNA sequence and amino acid sequence information:\n'
#     download_path=[]
#     for idx, aa_dna_info in enumerate(aa_dna_list):
#         aa_dna_content += f"[{idx+1}] {aa_dna_info[0]}protein({aa_dna_info[1]})：\n"
#         aa_dna_content += f"\>sp|{aa_dna_info[-3]}|{aa_dna_info[-2]}\n{aa_dna_info[-1]}\n\n"
#         aa_dna_content += f"\>{aa_dna_info[2]} {aa_dna_info[3]}\n{aa_dna_info[4]}\n\n"
#         to_path='/reference_data/'
#         download_path.append(to_path+'/'+aa_dna_info[0]+'.fa')
#         with open(to_path+'/'+aa_dna_info[0]+'.fa','w') as f:
#             f.write(f"{aa_dna_info[2]} {aa_dna_info[3]}\n{aa_dna_info[4]}")
#         f.close()
#
#     return aa_dna_content, download_path


#定义X-chat的对话逻辑
def predict(history_conversation, stage=1, search_type=None, search_prompt=""):
    window_history = history_conversation[-50:]
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
                if stage==1:
                    instruction = function_args['term']
                    user_input = window_history[-1]['content']#function_args['user_input']
                    data_info_dict = function_to_call(stage=1,instruction=instruction,user_input=user_input)
                    messages.append( # adding assistant response to messages
                            {
                                "role": search_response["role"],
                                "function_call": {
                                    "name": function_name,
                                    "arguments": search_response["function_call"]["arguments"],
                                    },
                                "content": None
                                }
                            )
                    messages.append( # adding assistant response to messages
                            {
                                "role": "function",
                                "name":function_name,
                                "content": data_info_dict["data_list"],
                                }
                            )
                    
                    #Modify user query in multiple rounds of self-reflection
                    if data_info_dict["data_list"]==[]:
                        i=1
                        query_now = function_args["term"]
                        while True:
                            print(f"""If no NCBI results are retrieved, it may be that the query is too long or too short. The current query is {query_now}. 
                            Please reflect on the modification. You only need to return my new term without unnecessary description.""")
                            messages.append({"role":'user','content':f"""If no NCBI results are retrieved, it may be that the query is too long or too short. The current query is {query_now}. 
                            Please reflect on the modification. You only need to return my new term without unnecessary description."""})
                            second_response = openai.ChatCompletion.create(
                                            messages=messages,
                                            deployment_id="XMGI-Chat1-GPT4"
                                            )
                            print(f'The modified term after the {i}-th round of reflection is:',second_response["choices"][0]["message"]['content'])
                            data_info_dict = function_to_call(instruction=second_response["choices"][0]["message"]['content'], stage=1)
                            if data_info_dict["data_list"]!=[]:
                                break
                            i+=1
                            if i >5:
                                break
                            # print(f'The search results of round {i} are:',data_info_dict["data_list"])
                            # query_now=second_response["choices"][0]["message"]['content']


                    # responses = 'Based on your question, I searched for the relevant sequences as follows, please select according to your needs:\n'+str(data_info_dict["options"])

                    history_conversation.append({'role':'assistant','content':data_info_dict})
                    responses_dict = {'responses':data_info_dict,'history_conversation':history_conversation}
                    return responses_dict
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
                species_name = function_args['instruction']
                user_response = pathogen_drug_resistance_search(stage=1,instruction=species_name)
                responses = user_response
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
                instruction = function_args['term']
                user_response = cancer_search(stage=1, instruction=instruction)
                responses = user_response
                history_conversation.append({'role': 'assistant', 'content': responses})
                responses_dict = {'responses': responses, 'history_conversation': history_conversation}
                # responses_dict={'target_cds_path':[],'non_target_cds_path':[],'fna_file':[fna_path],'snp_file':[snp_path],'download_path':[],'responses':responses,'history_conversation':history_conversation}
                return responses_dict

            elif search_response["function_call"]["name"] =="species_identification_search":
                # function_args = json.loads(search_response["function_call"]["arguments"])
                instruction = function_args['instruction']
                responses = species_identification_search(instruction,stage=1)
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
            if stage == 2:
                instruction = json.loads(history_conversation[-1]['content'])
                instruction = instruction['selected_options']

                responses = cancer_search(instruction, stage)
                history_conversation.append({'role': 'assistant', 'content': responses})
                responses_dict = {'responses': responses, 'history_conversation': history_conversation}
                return responses_dict

        elif search_type=="genetic_disorder_type":
            if stage == 2:
                # print('xxxxxxxxxxxxxhistory_conversation', history_conversation)
                instruction = json.loads(history_conversation[-1]['content'])
                instruction = instruction['selected_options']
                print('instruction',instruction)
                responses = genetic_disorder_search(instruction, stage)
                history_conversation.append({'role': 'assistant', 'content': responses})
                responses_dict = {'responses': responses, 'history_conversation': history_conversation}
                return responses_dict
        elif search_type=="protein_mutation_type":
            if stage==2:
                instruction = json.loads(history_conversation[-1]['content'])
                instruction = instruction['selected_options']
                responses = protein_mutation_search(instruction, stage)
                history_conversation.append({'role': 'assistant', 'content': responses})
                responses_dict = {'responses': responses, 'history_conversation': history_conversation}
                return responses_dict
                # return responses
            elif stage==3:
                pass
        elif search_type == "download_type":
            if stage == 2:
                print('开始正式下载文件')
                instruction = json.loads(history_conversation[-1]['content'])
                instruction = instruction['selected_options']
                responses = download_files(instruction, stage)
                history_conversation.append({'role': 'assistant', 'content': responses})
                responses_dict = {'responses': responses, 'history_conversation': history_conversation}
                return responses_dict
                # return responses
            elif stage == 3:
                pass
            pass
        elif search_type=="pathogen_drug_resistance_type":
            if stage==2:
                instruction = json.loads(history_conversation[-1]['content'])
                instruction = instruction['selected_options']
                responses = pathogen_drug_resistance_search(instruction, stage)
                history_conversation.append({'role': 'assistant', 'content': responses})
                responses_dict = {'responses': responses, 'history_conversation': history_conversation}
                return responses_dict

        else:
            print("error search type")


# 定义X-chat的对话逻辑
def primer_predict(history_conversation, stage=1, search_type=None, search_prompt=""):
    window_history = history_conversation[-50:]
    if search_type == None:
        system_prompt = [{'role': 'system', 'content': """You are an artificial intelligence assistant designed to help users design specific primers. 
When the user asks for follow-up work or subsequent primer design and provides parameters for primer design, 
you should call the snp_primer_design function to design primers;
when the user requests a multiplex PCR experiment for pathogen species identification analysis, you should call the species_identification_search function."""
                          }]
        messages = system_prompt + window_history
        search_response = get_chatgpt_response(messages, primer_functions)
        if search_response.get("function_call"):
            print('*******************Successful activate function*******************')
            # Call the function. The JSON response may not always be valid so make sure to handle errors
            function_name = search_response["function_call"]["name"]
            available_functions = {
                                   "species_identification_search": species_identification_search,
                                   "snp_primer_design": snp_primer_design,
                                   }
            function_to_call = available_functions[function_name]
            print("Function name: ", search_response["function_call"]["name"])
            print("Function input: ", search_response["function_call"]["arguments"])


            if search_response["function_call"]["name"] == "snp_primer_design":
                print('成功调用primer design')

                function_args = json.loads(search_response["function_call"]["arguments"])
                print('xxx', history_conversation[-1]['content'])
                instruction = ast.literal_eval(history_conversation[-1]['content'])
                max_amp_len = function_args['max_amp_len']
                min_amp_len = function_args['min_amp_len']
                temperature = function_args['temperature']
                min_GC = function_args['min_GC']
                max_GC = function_args['max_GC']
                max_len = function_args['max_len']
                parameter = [max_amp_len, min_amp_len, temperature, min_GC, max_GC, max_len]
                print('parameter',parameter)
                # print('xxxxinstruction', instruction)

                responses = snp_primer_design(instruction,parameter)
                history_conversation.append({'role': 'assistant', 'content': responses})
                responses_dict = {'responses': responses, 'history_conversation': history_conversation}
                return responses_dict


            elif search_response["function_call"]["name"] == "species_identification_search":
                function_args = json.loads(search_response["function_call"]["arguments"])
                responses = species_identification_search(**function_args)
                history_conversation.append({'role': 'assistant', 'content': responses})
                responses_dict = {'responses': responses, 'history_conversation': history_conversation}
                # responses_dict['history_conversation'] = history_conversation
                return responses_dict

            else:
                print("error funcation calling")

        else:
            responses = search_response['content']
            history_conversation.append({'role': 'assistant', 'content': responses})
            responses_dict = {'target_cds_path': [], 'non_target_cds_path': [], 'snp_file': [], 'fna_file': [],
                              'download_path': [], 'responses': responses, 'history_conversation': history_conversation}
            return responses_dict



app = FastAPI()


# 解析输入数据
class InputModel(BaseModel):
    instruction: str
    history: list
    type: int

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

executor = ThreadPoolExecutor()




# WebSocket 路由
@app.websocket("/ncbi")
async def websocket_endpoint(websocket: WebSocket):
    await websocket.accept()

    asyncio.create_task(heartbeat(websocket))
    try:
        while True:
            data = await websocket.receive_json()
            print("接收到数据：",data)
            input_model = InputModel(**data)
            method_type = input_model.type
            input_ = input_model.instruction
            history = input_model.history
            
            #第一轮对话
            if history == []:
                history_conversation = []
                input_ = input_.replace('\n','').replace('   ','')
                history_conversation.append({'role':'user','content':f'{input_}'})
                print("当前用户输入：",input_)
                loop = asyncio.get_event_loop()
                responses_dict = await loop.run_in_executor(executor, predict, history_conversation)
                # history_conversation = responses_dict['history_conversation']
            #多轮对话
            else:
            # elif 'primer_type' not in input_:
                history_conversation = [{'role':his_['role'],'content':his_['content']}for his_ in history]
                history_conversation.append({'role':'user','content':f'{input_}'})
                print("当前用户输入：",input_)
                print("当前对话历史：",history_conversation)
                temp_flag = False
                if "search_type" in input_:
                    temp_flag = True

                loop = asyncio.get_event_loop()
                if temp_flag:
                    print('xxxxx进行第二次回复')
                    search_type = json.loads(input_)["search_type"]
                    responses_dict = await loop.run_in_executor(executor, predict, history_conversation,2, search_type,'')
                else:
                    responses_dict = await loop.run_in_executor(executor, predict, history_conversation)
                # history_conversation = responses_dict['history_conversation']
            # else:
            #     history_conversation = [{'role': his_['role'], 'content': his_['content']} for his_ in history]
            #     history_conversation.append({'role': 'user', 'content': f'{input_}'})
            #     print("当前用户输入：", input_)
            #     print("当前对话历史：", history_conversation)
            #     search_type = json.loads(input_)["search_type"]
            #     instruction = json.loads(history_conversation[-1]['content'])['selected_options']
            #     # if search_type == "genetic_disorder_type":
            #
            #
            #     responses_dict = await loop.run_in_executor(executor, snp_primer_design, search_type,instruction)

            # print("输出:", responses_dict['responses'])
            #循环？
            await websocket.send_text(json.dumps({'content':responses_dict,'type':method_type,'state':'stop'},ensure_ascii=False))
                    
    except Exception as e:
        traceback.print_exc()
    finally:
        print('finished!')

@app.websocket("/primer")
async def websocket_endpoint(websocket: WebSocket):
    await websocket.accept()

    asyncio.create_task(heartbeat(websocket))
    try:
        while True:
            data = await websocket.receive_json()
            print("primer接收到数据：", data)
            input_model = InputModel(**data)
            method_type = input_model.type
            input_ = input_model.instruction
            history = input_model.history

            # 第一轮对话
            if history == []:
                history_conversation = []
                history_conversation.append({'role': 'user', 'content': f'{data}'})
                print("当前用户输入：", input_)
                loop = asyncio.get_event_loop()
                responses_dict = await loop.run_in_executor(executor, primer_predict, history_conversation)
                # history_conversation = responses_dict['history_conversation']
            # 多轮对话
            else:
                print("当前多轮对话用户输入：", input_)
                history_conversation = [{'role': his_['role'], 'content': his_['content']} for his_ in history]
                # history_conversation.append({'role': 'user', 'content': f'{input_}'})
                history_conversation.append({'role': 'user', 'content': f'{data}'})
                # print("当前用户输入：", input_)
                print("当前多轮对话历史：", history_conversation)
                loop = asyncio.get_event_loop()
                responses_dict = await loop.run_in_executor(executor, primer_predict, history_conversation)

            await websocket.send_text(
                json.dumps({'content': responses_dict, 'type': method_type, 'state': 'stop'}, ensure_ascii=False))

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
    uvicorn.run(app, host="0.0.0.0", port=8082,timeout_keep_alive=300)


