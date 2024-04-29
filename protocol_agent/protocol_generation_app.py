from Bio import Entrez
import openai
import requests
import re
from bs4 import BeautifulSoup
import json
import os
import zipfile
import glob
import shutil
import Levenshtein
import math

def pre_calculate_two_step(experiment_stage, extract_information_json,pre_calculate_parameter_json={}):
    
    #human input
    dna_concentration = 
    sample_size_1st = 
    #用户上传Qubit产出的 DNA_ concentration
    
    plate_well = [f"{chr(65 + i)}{j}" for j in range(1, 12 + 1) for i in range(8)]
    #1st PCR formula
    if experiment_stage==1:
        #
        tip_rack_number = ['2','3','4','5','9'][:(sample_size_1st*extract_information_json["First_PCR_Amplification_Reaction_Mixture_Components_Numbers"]/96)]
        
        #17、	1st_tcmod_block_max_volume
        #1st_PCR Enzyme Mix + 1st_Barcode1_volume + DNA Template Volume + 1st_PCR Clean Enzyme
        tcmod_1st_block_max_volume = extract_information_json["PCR_1st_Enzyme_Mix"] + extract_information_json["Barcode1_1st_volume"] + extract_information_json["DNA_Template_Volume"] + extract_information_json["PCR_1st_Clean_Enzyme"]

        #19、	Barcode 1 well
        #Sample size = Barcode1 well 注意！按列移动溶液，所以是A1~H1，A2~H2
        Barcode_1_well = plate_well[:sample_size_1st]
        
        #21、	1st_PCR_MIX_Volume_fault_tolerance_volume
        #(1st_PCR Enzyme Mix + 1st_PCR Clean Enzyme + (1st_PCR Enzyme Mix + 1st_PCR Clean Enzyme)*0.1)*sample size
        PCR_1st_MIX_Volume_fault_tolerance_volume = (extract_information_json["PCR_1st_Enzyme_Mix"] + + extract_information_json["PCR_1st_Clean_Enzyme"] + (extract_information_json["DNA_Template_Volume"] + extract_information_json["Barcode1_1st_volume"])*0.1)*sample_size

        #22、	1st_PCR_MIX_Volume
        #1st_PCR Enzyme Mix + 1st_PCR Clean Enzyme
        PCR_1st_MIX_Volume = extract_information_json["PCR_1st_Enzyme_Mix"] + extract_information_json["PCR_1st_Clean_Enzyme"]

        #25、	1st_PCR_Mix_well计算
        #[1st_PCR_MIX_Volume] * sample size / 96_well_volume #注意！按列移动溶液，所以是A1~H1，A2~H2,…An~Hn
        PCR_1st_Mix_well = PCR_1st_MIX_Volume*sample_size_1st/extract_information_json["well96_volume"]

        #27、	Primer_Pool_well_Volume计算
        #Primer Pool Volume * sample size + (Primer Pool Volume * sample size) *0.1
        Primer_Pool_well_Volume = extract_information_json["Primer_Pool_Volume"] * sample_size_1st + (extract_information_json["Primer_Pool_Volume"]*sample_size_1st)*0.1

        #1st_DNA_TEMPLATE_WELL 计算
        #1st_DNA_TEMPLATE_WELL = sample size 注意！按列移动溶液，所以是A1~H1，A2~H2
        DNA_1st_TEMPLATE_WELL = plate_well[:sample_size_1st]
        
        extract_information_json["sample_size_1st"] = sample_size_1st
        pre_calculate_parameter_json={
                                      "tip_rack_number":tip_rack_number,
                                      "tcmod_1st_block_max_volume":tcmod_1st_block_max_volume,
                                      "Barcode_1_well":Barcode_1_well,
                                      "PCR_1st_MIX_Volume_fault_tolerance_volume":PCR_1st_MIX_Volume_fault_tolerance_volume,
                                      "PCR_1st_MIX_Volume":PCR_1st_MIX_Volume,
                                      "PCR_1st_Mix_well":PCR_1st_Mix_well,
                                      "Primer_Pool_well_Volume":Primer_Pool_well_Volume,
                                      "DNA_1st_TEMPLATE_WELL":DNA_1st_TEMPLATE_WELL
                                      }

        return extract_information_json, pre_calculate_parameter_json

    #1st purification formula
    if experiment_stage==2:

        #1st_purification_TE_Buffer_tolerance = 1st_purification_TE_Buffer_elute_DNA*2
        purification_1st_TE_Buffer_tolerance = extract_information_json["1st_purification_TE_Buffer_elute_DNA"]
        pre_calculate_parameter_json["purification_1st_TE_Buffer_tolerance"] = purification_1st_TE_Buffer_tolerance
        return extract_information_json, pre_calculate_parameter_json


    #2nd PCR formula
    if experiment_stage==3:
        #28、	U-R_well_Volume 计算
        #U-R Volume + (U-R Volume)*1
        U_R_well_Volume = extract_information_json["U_R_well_Volume"] + extract_information_json["U_R_well_Volume"]*1

        #18、	2nd_tcmod_block_max_volume
        #2nd_PCR Enzyme Mix + 2nd_PCR Clean Enzyme + PCR Addictive + TE Buffer+2nd_DNA_volume+2nd_Barcode2_volume+ U-R_well_Volume
        tcmod_2nd_block_max_volume = extract_information_json["PCR_2nd_Enzyme_Mix"] + extract_information_json["PCR_2nd_Clean_Enzyme"] + extract_information_json["PCR_Addictive"] + extract_information_json["TE_Buffer"] + extract_information_json["DNA_2nd_volume"] + extract_information_json["Barcode2_2nd_volume"] + extract_information_json["Barcode2_2nd_volume"] + U_R_well_Volume

        #Barcode 2 well计算
        #Sample size /96 = 向上取整 Barcode 2 well 
        #len(Barcode 2 well) --> sample size --> tip_rack_number = sample size * reagent_type_2nd_number
        #Barcode_1_well = plate_well[:sample_size]
        Barcode_2_well = plate_well[:math.ceil(extract_information_json["sample_size_1st"]/96)]
        sample_size_2nd = len(Barcode_2_well)
        tip_rack_number_2nd = sample_size_2nd * extract_information_json["Second_PCR_Amplification_Reaction_Mixture_Components_Numbers"]

        extract_information_json["sample_size_2nd"] = sample_size_2nd
        
        #23、	2nd_PCR_MIX_Volume_fault_tolerance_volume
        #2nd_PCR Enzyme Mix + 2nd_PCR Clean Enzyme + PCR Addictive + TE Buffer + (2nd_PCR Enzyme Mix + 2nd_PCR Clean Enzyme + PCR Addictive + TE Buffer)*0.1
        PCR_2nd_MIX_Volume_fault_tolerance_volume = extract_information_json["PCR_2nd_Enzyme_Mix"] + extract_information_json["PCR_2nd_Clean_Enzyme"] + extract_information_json["PCR_Addictive"] + extract_information_json["TE_Buffer"] +(extract_information_json["PCR_2nd_Enzyme_Mix"] + extract_information_json["PCR_2nd_Clean_Enzyme"] + extract_information_json["PCR_Addictive"] + extract_information_json["TE_Buffer"])*0.1


        #24、	2nd_PCR_MIX_Volume
        #2nd_PCR Enzyme Mix + 2nd_PCR Clean Enzyme + PCR Addictive + TE Buffer
        PCR_2nd_MIX_Volume = extract_information_json["PCR_2nd_Enzyme_Mix"] + extract_information_json["PCR_2nd_Clean_Enzyme"] + extract_information_json["PCR_Addictive"] + extract_information_json["TE_Buffer"]



        #26、	2nd_PCR_Mix_well计算
        #[2nd_PCR_MIX_Volume] * sample size / 96_well_volume
        PCR_2nd_Mix_well = PCR_2nd_MIX_Volume*sample_size_2nd/extract_information_json["well96_volume"]


        #2nd_DNA_TEMPLATE_WELL
        #2nd_DNA_TEMPLATE_WELL = sample size 注意！按列移动溶液，所以是A1~H1，A2~H2
        DNA_2nd_TEMPLATE_WELL = plate_well[:sample_size_2nd]

        pre_calculate_parameter_json["DNA_1st_TEMPLATE_WELL"] = DNA_2nd_TEMPLATE_WELL
        pre_calculate_parameter_json["PCR_2nd_Mix_well"] = PCR_2nd_Mix_well
        pre_calculate_parameter_json["Barcode_2_well"] = Barcode_2_well
        pre_calculate_parameter_json["tip_rack_number_2nd"] = tip_rack_number_2nd
        pre_calculate_parameter_json["PCR_2nd_MIX_Volume_fault_tolerance_volume"] = PCR_2nd_MIX_Volume_fault_tolerance_volume
        pre_calculate_parameter_json["PCR_2nd_MIX_Volume"] = PCR_2nd_MIX_Volume
        pre_calculate_parameter_json["tcmod_2nd_block_max_volume"] = tcmod_2nd_block_max_volume

        return extract_information_json, pre_calculate_parameter_json

    

    #2nd purification formula
    if experiment_stage==4:

        #2nd_purification_TE_Buffer_tolerance = 2nd_purification_TE_Buffer_elute_DNA*2
        purification_2nd_TE_Buffer_tolerance = extract_information_json["2nd_purification_TE_Buffer_elute_DNA"]*2
        #3、2nd_purification_supernatant_volumen 获取
        #RAG protocol中的提取
        purification_2nd_supernatant_volumen = extract_information_json["2nd_purification_supernatant_volumen"]
        
        pre_calculate_parameter_json["purification_2nd_TE_Buffer_tolerance"] = purification_2nd_TE_Buffer_tolerance
        pre_calculate_parameter_json["purification_2nd_supernatant_volumen"] = purification_2nd_supernatant_volumen

        return extract_information_json, pre_calculate_parameter_json
    
    #Make DNB
    if experiment_stage==5:
        # 2、	dsDNA_library_input 获取
        # RAG protocol的Make DNB中有提及
        # 30ng/ DNA_ concentration = dsDNA_library_input
        dsDNA_library_input = 30/dna_concentration
        
        # 3、	dsDNA_library_volume 计算
        # dsDNA_library_volume = dsDNA_library_input / DNA_ concentration
        dsDNA_library_volume = dsDNA_library_input/dna_concentration

        # 4、	dsDNA_library_volume_tolerance计算
        # dsDNA_library_volume_tolerance = dsDNA_library_volume*2
        dsDNA_library_volume_tolerance = dsDNA_library_volume*2

        # 5、	Low_TE_buffer_Volume计算
        # Low_TE_buffer_Volume = 20 - dsDNA_library_volume
        Low_TE_buffer_Volume = 20 - dsDNA_library_volume

        # 6、	Low_TE_buffer_Volume_tolerance计算
        # Low_TE_buffer_Volume_tolerance = Low_TE_buffer_Volume*2
        Low_TE_buffer_Volume_tolerance = Low_TE_buffer_Volume*2

        # 8、	 Make_DNB_Buffer_Volume_tolerance 计算
        # Make_DNB_Buffer_Volume_tolerance = Make_DNB_Buffer_Volume*2
        Make_DNB_Buffer_Volume_tolerance = extract_information_json["Make_DNB_Buffer_Volume"] * 2

        # 11、	Make_DNB_Enzyme_Mix_I_tolerance 计算
        # Make_DNB_Enzyme_Mix_I_tolerance = Make_DNB_Enzyme_Mix_I*2
        Make_DNB_Enzyme_Mix_I_tolerance = extract_information_json["Make_DNB_Enzyme_Mix_I"] * 2

        # 12、	Make_DNB_Enzyme_Mix_II_tolerance计算
        # Make_DNB_Enzyme_Mix_II_tolerance = Make_DNB_Enzyme_Mix_II*2
        Make_DNB_Enzyme_Mix_II_tolerance = extract_information_json["Make_DNB_Enzyme_Mix_II"] * 2


        # 14、	Make_DNB_Stop_Buffer_tolerance计算
        # Make_DNB_Stop_Buffer_tolerance = Make_DNB_Stop_Buffer*2
        Make_DNB_Stop_Buffer_tolerance = extract_information_json["Make_DNB_Stop_Buffer"] * 2
        
        pre_calculate_parameter_json["dsDNA_library_input"] = dsDNA_library_input
        pre_calculate_parameter_json["dsDNA_library_volume"] = dsDNA_library_volume
        pre_calculate_parameter_json["dsDNA_library_volume_tolerance"] = dsDNA_library_volume_tolerance
        pre_calculate_parameter_json["Low_TE_buffer_Volume"] = Low_TE_buffer_Volume
        pre_calculate_parameter_json["Low_TE_buffer_Volume_tolerance"] = Low_TE_buffer_Volume_tolerance
        pre_calculate_parameter_json["Make_DNB_Buffer_Volume_tolerance"] = Make_DNB_Buffer_Volume_tolerance
        pre_calculate_parameter_json["Make_DNB_Enzyme_Mix_I_tolerance"] = Make_DNB_Enzyme_Mix_I_tolerance
        pre_calculate_parameter_json["Make_DNB_Enzyme_Mix_II_tolerance"] = Make_DNB_Enzyme_Mix_II_tolerance
        pre_calculate_parameter_json["Make_DNB_Stop_Buffer_tolerance"] = Make_DNB_Stop_Buffer_tolerance
        
        return extract_information_json, pre_calculate_parameter_json


def pre_calculate_free_step():
    pass


def get_retrieval(retrieval_response):
    url = 'http://127.0.0.1:8000/search/'
    data = {"text": retrieval_response,
            "top_k": "10",
            "top_r": "1",
            "score_threshold": "0.5"}
    response = requests.post(url, json=data)

    if response.status_code == 200:
        retrieval = response.json()
        print(retrieval)  # 打印结果
    else:
        print("Failed to fetch data:", response.status_code, response.text)
    return retrieval


def extract_json_response(response):
    # 从response中提取JSON内容
    experiment_info_dict = response["choices"][0]["message"]['content']
    match = re.search(r'```json\n(.*?)\n```', experiment_info_dict, re.DOTALL)
    if match:
        experiment_info_json_str = match.group(1)
        experiment_info_json_data = json.loads(experiment_info_json_str)
    else:
        print("GPT构造JSON错误，情况如下所示：\n",experiment_info_dict)
        exit()
    return experiment_info_json_data

def extract_python_response(response):
    # 从response中提取JSON内容
    experiment_info_dict = response["choices"][0]["message"]['content']
    match = re.search(r'```python\n(.*?)\n```', experiment_info_dict, re.DOTALL)
    if match:
        experiment_info_json_str = match.group(1)
        experiment_info_json_data = json.loads(experiment_info_json_str)
    else:
        print("GPT构造JSON错误，情况如下所示：\n",experiment_info_dict)
        exit()
    return experiment_info_json_data


def protocol_generation(instruction, stage):
    if stage == 1:
        #stage1 用户可以选择试验方案SOP，也可以通过引物设计数量来确定（step free or two step）
        #暂时不考虑用户上传自己的protocol
        primer_num = instruction["primer_number"]

        if primer_num <=100:
            with open("","r") as f:
                step_free_protocol_doc = f.read()
            f.close()
            response = "It has been confirmed that you have completed the sequence retrieval and primer synthesis tasks, and you would like to begin designing Opentrons scripts based on an experimental protocol. Given your primer design results, with the number of primers being 100 or fewer, I recommend using a step-free protocol. I also offer the following options for you to freely choose the experimental protocol method. If you do not make a selection, I will proceed to design the Opentrons script using the default step-free protocol mode." 
            stage += 1
            return {
                "response": response, 
                "operations":[{"title":"recommend protocol","options":["two", "free"],"type":["single","input"],"key":"protocol_select","value":[]}],
                "stage": stage,
                "state": 'continue',
                "protocol_doc": step_free_protocol_doc}
        else:
            with open("","r") as f:
                two_step_protocol_doc = f.read()
            f.close()
            response = "It has been confirmed that you have completed the sequence retrieval and primer synthesis tasks, and you would like to begin designing Opentrons scripts based on an experimental protocol. Given the results of your primer design, with the number of primers exceeding 100, I recommend using a two-step protocol. I also offer the following options for you to freely choose the experimental protocol method. If you do not make a selection, I will default to designing the Opentrons script using the two-step protocol mode."
            stage += 1
            return {
                "response": response, 
                "operations":[
                        {"title":"recommend protocol",
                          "options":["two", "free"],
                          "type":["single","input"],
                          "key":"protocol_select","value":[]
                        },
                        {
                         "title":"upload protocol", 
                         "options":[], 
                         "type":["file"], 
                         "key":"protocol_path", 
                         "value":[]
                        }],
                "stage": stage,
                "state": 'continue', 
                "protocol_doc": two_step_protocol_doc}
    
    if stage == 2:
        #stage2 用户有可能上传Protocol时，那就RAG检索相似protocol；
        #stage2 用户直接选择protocol时，直接进一步提示用户需要提取的信息，用于后续的pre-calculate；

        conversation = instruction['conversation']
        primer_num = instruction["primer_number"]
        protocol_select = instruction["protocol_select"]
        
        
        #用户选择protocol实验
        if protocol_select != "":
            if 'two' in protocol_select:
                read_file_text = "****" #TODO设定多种文件的文本获取函数,extract(protocol_upload)
                            #extract information form retrieval doc
                retrieval_extraction_prompt="""I will provide the specific content of the article below. Please read the article carefully, and then extract the values for JSON from it. Return the JSON to me. Note! You do not need to include units; you only need to extract the numerical values.
'''json
{
    "1st_purification_TE_Buffer_elute_DNA": "",
    "2nd_purification_TE_Buffer_elute_DNA": "",
    "dsDNA_library_input": "",
    "Make_DNB_Buffer_Volume": "",
    "Make_DNB_Enzyme_Mix_I": "",
    "Make_DNB_Enzyme_Mix_II": "",
    "Make_DNB_Stop_Buffer": "",
    "Primer_hybridization_Reaction_step1_temperature": "",
    "Primer_hybridization_Reaction_step2_temperature": "",
    "Primer_hybridization_Reaction_step3_temperature": "",
    "Rolling_circle_amplification_step1_temperature": "",
    "Rolling_circle_amplification_step2_temperature": "",
    "Primer_Pool_Volume": "",
    "PCR_1st_Enzyme_Mix": "",
    "PCR_1st_Clean_Enzyme": "",
    "PCR_2nd_Enzyme_Mix": "",
    "PCR_2nd_Clean_Enzyme": "",
    "96_Well_Plate": "",
    "U_R_well_Volume": "",
    "PCR_Addictive": "",
    "TE_Buffer": "",
    "DNA_Template_Volume": "",
    "Barcode1_1st_volume": "",
    "Barcode2_2nd_volume": "",
    "DNA_2nd_volume": "",
    "PCR_1st_Cycles": "",
    "PCR_2nd_Cycles": "",
    "First_PCR_Amplification_Reaction_Mixture_Components_Numbers":"",
    "Second_PCR_Amplification_Reaction_Mixture_Components_Numbers":"",
    "2nd_purification_supernatant_volumen":""    
}
'''"""
                extraction_protocol_description = openai.ChatCompletion.create(
                    messages=[{"role": "system", "content": retrieval_extraction_prompt},
                            {"role": "user", "content": retrieval_extraction_prompt+f"\n##article:'''{relevance_doc}'''"+"\n##Notice! Don't lose objects, calculate carefully."}],
                    deployment_id="XMGI-Chat1-GPT4"
                )
                extraction_protocol_parameter = extract_python_response(extraction_protocol_description)
                print("-------------extraction protocol parameter information:\n", extraction_protocol_parameter)
                response = "You need to provide the number of samples, and subsequently assist me with calculating all the variables required for the experiment."  #下一步提取
                stage += 1
                return {
                "response": response, 
                "operations":[],
                "stage": stage,
                "state": 'continue',
                "protocol_doc": step_free_protocol_doc,
                "extraction_protocol_parameter":extraction_protocol_parameter}
            
            else:
                read_file_text = "****" #TODO设定多种文件的文本获取函数,extract(protocol_upload)
                #extract information form retrieval doc
                retrieval_extraction_prompt="""I will provide the specific content of the article below. Please read the article carefully, and then extract the values for JSON from it. Return the JSON to me. Note! You do not need to include units; you only need to extract the numerical values.
'''json
{
    "1st_purification_TE_Buffer_elute_DNA": "",
    "2nd_purification_TE_Buffer_elute_DNA": "",
    "dsDNA_library_input": "",
    "Make_DNB_Buffer_Volume": "",
    "Make_DNB_Enzyme_Mix_I": "",
    "Make_DNB_Enzyme_Mix_II": "",
    "Make_DNB_Stop_Buffer": "",
    "Primer_hybridization_Reaction_step1_temperature": "",
    "Primer_hybridization_Reaction_step2_temperature": "",
    "Primer_hybridization_Reaction_step3_temperature": "",
    "Rolling_circle_amplification_step1_temperature": "",
    "Rolling_circle_amplification_step2_temperature": "",
    "Primer_Pool_Volume": "",
    "PCR_1st_Enzyme_Mix": "",
    "PCR_1st_Clean_Enzyme": "",
    "PCR_2nd_Enzyme_Mix": "",
    "PCR_2nd_Clean_Enzyme": "",
    "96_Well_Plate": "",
    "U_R_well_Volume": "",
    "PCR_Addictive": "",
    "TE_Buffer": "",
    "DNA_Template_Volume": "",
    "Barcode1_1st_volume": "",
    "Barcode2_2nd_volume": "",
    "DNA_2nd_volume": "",
    "PCR_1st_Cycles": "",
    "PCR_2nd_Cycles": "",
    "First_PCR_Amplification_Reaction_Mixture_Components_Numbers":"",
    "Second_PCR_Amplification_Reaction_Mixture_Components_Numbers":"",
    "2nd_purification_supernatant_volumen":""    
}
'''"""
                extraction_protocol_description = openai.ChatCompletion.create(
                    messages=[{"role": "system", "content": retrieval_extraction_prompt},
                            {"role": "user", "content": retrieval_extraction_prompt+f"\n##article:'''{relevance_doc}'''"+"\n##Notice! Don't lose objects, calculate carefully."}],
                    deployment_id="XMGI-Chat1-GPT4"
                )
                extraction_protocol_parameter = extract_python_response(extraction_protocol_description)
                print("-------------extraction protocol parameter information:\n", extraction_protocol_parameter)

                response = "You need to provide the number of samples, and subsequently assist me with calculating all the variables required for the experiment."  #下一步提取
                stage += 1
                return {
                "response": response, 
                "operations":[],
                "stage": stage,
                "state": 'continue',
                "protocol_doc": step_free_protocol_doc,
                "extraction_protocol_parameter":extraction_protocol_parameter}
            
        
        
        #用户选择上传protocol
        else:
            protocol_upload = instruction["protocol_path"]
            read_file_text = "****" #TODO设定多种文件的文本获取函数,extract(protocol_upload)
            #get protocol summary to retrieval tempalte code
            experiment_type_prompt = """Based on the context provided and the protocol file uploaded by the user, please summarize the content of the uploaded protocol file. 
Generate a description of up to 512 words that captures the core content of the protocol.
"""

            conversation_str = "\n\n".join([i["role"] + ':\n' + i['content'] for i in conversation])
            upload_protocol_description = openai.ChatCompletion.create(
                messages=[{"role": "system", "content": experiment_type_prompt},
                            {"role": "user", "content": "Based on the context provided and the protocol file uploaded by the user, please summarize the content of the uploaded protocol file. Generate a description of up to 512 words that captures the core content of the protocol.\n"+"#conversation:"+conversation_str+f"#upload_file_text:{read_file_text}"}
                            ],
                deployment_id="XMGI-Chat1-GPT4"
            )
            print("-------------extract experiment information:\n", upload_protocol_description)
            
            relevance_doc = get_retrieval(upload_protocol_description)
            

            #extract information form retrieval doc
            retrieval_extraction_prompt="""I will provide the specific content of the article below. Please read the article carefully, and then extract the values for JSON from it. Return the JSON to me. Note! You do not need to include units; you only need to extract the numerical values.
'''json
{
    "1st_purification_TE_Buffer_elute_DNA": "",
    "2nd_purification_TE_Buffer_elute_DNA": "",
    "dsDNA_library_input": "",
    "Make_DNB_Buffer_Volume": "",
    "Make_DNB_Enzyme_Mix_I": "",
    "Make_DNB_Enzyme_Mix_II": "",
    "Make_DNB_Stop_Buffer": "",
    "Primer_hybridization_Reaction_step1_temperature": "",
    "Primer_hybridization_Reaction_step2_temperature": "",
    "Primer_hybridization_Reaction_step3_temperature": "",
    "Rolling_circle_amplification_step1_temperature": "",
    "Rolling_circle_amplification_step2_temperature": "",
    "Primer_Pool_Volume": "",
    "PCR_1st_Enzyme_Mix": "",
    "PCR_1st_Clean_Enzyme": "",
    "PCR_2nd_Enzyme_Mix": "",
    "PCR_2nd_Clean_Enzyme": "",
    "96_Well_Plate": "",
    "U_R_well_Volume": "",
    "PCR_Addictive": "",
    "TE_Buffer": "",
    "DNA_Template_Volume": "",
    "Barcode1_1st_volume": "",
    "Barcode2_2nd_volume": "",
    "DNA_2nd_volume": "",
    "PCR_1st_Cycles": "",
    "PCR_2nd_Cycles": "",
    "First_PCR_Amplification_Reaction_Mixture_Components_Numbers":"",
    "Second_PCR_Amplification_Reaction_Mixture_Components_Numbers":"",
    "2nd_purification_supernatant_volumen":""    
}
'''"""
            extraction_protocol_description = openai.ChatCompletion.create(
                messages=[{"role": "system", "content": retrieval_extraction_prompt},
                          {"role": "user", "content": retrieval_extraction_prompt+f"\n##article:'''{relevance_doc}'''"+"\n##Notice! Don't lose objects, calculate carefully."}],
                deployment_id="XMGI-Chat1-GPT4"
            )
            extraction_protocol_parameter = extract_python_response(extraction_protocol_description)
            print("-------------extraction protocol parameter information:\n", extraction_protocol_parameter)
            

            

            if 'two' in relevance_doc:
                response = f"We found that the file you uploaded is similar to the local protocol. We recommend{relevance_doc}. You need to provide the number of samples, and subsequently assist me with calculating all the variables required for the experiment."  #下一步提取
                stage += 1
                return {
                "response": response, 
                "operations":[],
                "stage": stage,
                "state": 'continue',
                "protocol_doc": relevance_doc,
                "extraction_protocol_parameter":extraction_protocol_parameter}
            
            else:
                response = f"We found that the file you uploaded is similar to the local protocol. We recommend{relevance_doc}. You need to provide the number of samples, and subsequently assist me with calculating all the variables required for the experiment."  #下一步提取
                stage += 1
                return {
                "response": response, 
                "operations":[],
                "stage": stage,
                "state": 'continue',
                "protocol_doc": relevance_doc,
                "extraction_protocol_parameter":extraction_protocol_parameter}
            

    if stage == 3:
        #stage3 提取用户的描述信息，然后激活pre-calculate计算流程，获得JSON文件。当用户没有描述完所有信息，将会进一步循环补充信息。
        
        conversation = instruction['conversation']
        protocol_select = instruction["protocol_doc"]
        extraction_protocol_parameter = instruction["extraction_protocol_parameter"]
        
        
        if "two" in protocol_select:
            two_step_pre_calculate_json=json.loads(extraction_protocol_parameter)
            two_step_pre_calculate_json["sample_size_1st"] = ""
            json_data = json.dumps(two_step_pre_calculate_json)
            pre_calculate_prompt = f"""Based on the context of the communication and the user's description, obtain the information required for the following JSON.'''json\n{json_data}\nNote that you only need to return the refined and extracted JSON file, without any additional descriptions. Ensure that no information is lost in the extraction; if the user does not provide specific data, that field should be null.'''"""

        else:
            step_free_pre_calculate_json=json.loads(extraction_protocol_parameter)
            step_free_pre_calculate_json["sample_size_1st"] = ""
            json_data = json.dumps(step_free_pre_calculate_json)
            pre_calculate_prompt = f"""Based on the context of the communication and the user's description, obtain the information required for the following JSON.'''json\n{json_data}\nNote that you only need to return the refined and extracted JSON file, without any additional descriptions. Ensure that no information is lost in the extraction; if the user does not provide specific data, that field should be null.'''"""
        

        conversation_str = "\n\n".join([i["role"] + ':\n' + i['content'] for i in conversation])
        pre_calculate_parameter_description = openai.ChatCompletion.create(
            messages=[{"role": "system", "content": pre_calculate_prompt},
                        {"role": "user", "content": f"""{pre_calculate_prompt}"""+"#conversation:"+conversation_str+f"""{pre_calculate_prompt}"""}
                        ],
            deployment_id="XMGI-Chat1-GPT4"
        )
        pre_calculate_parameter_data = extract_json_response(pre_calculate_parameter_description)
        print("-------------pre-calculate data information:\n", pre_calculate_parameter_data)

        flag_pre_calculate_parameter = [x for x,y in pre_calculate_parameter_data.items() if y==None]

        if flag_pre_calculate_parameter != []:
            response = "你缺少提供的参数如下***"
            return {
                "response": response, 
                "operations":[],
                "stage": stage,
                "state": 'continue',
                "protocol_doc": relevance_doc,
                "pre_calculate_parameter_data": pre_calculate_parameter_data,
                "extraction_protocol_parameter":extraction_protocol_parameter}

        else:
            response = "我已后续实验所需要"
            if "two" in protocol_select:
                extraction_protocol_parameter, pre_calculate_parameter_json = pre_calculate_two_step(pre_calculate_parameter_data,extraction_protocol_parameter)
            else:
                extraction_protocol_parameter, pre_calculate_parameter_json = pre_calculate_free_step(pre_calculate_parameter_data,extraction_protocol_parameter)
            stage += 1
            return {
                "response": response, 
                "operations":[],
                "stage": stage,
                "state": 'continue',
                "protocol_doc": relevance_doc,
                "pre_calculate_parameter_data": pre_calculate_parameter_json,
                "extraction_protocol_parameter":extraction_protocol_parameter}


    if stage == 4:
        #stage4 拆解Protocol多个Task，并且让用户确认
        conversation = instruction['conversation']
        protocol_select = instruction["protocol_doc"]
                
        decomposed_prompt = f"""Based on the previous communication and the content of the protocol article provided below, please break down the protocol article into different experimental steps. 
For example, the steps for a Two-step protocol could be as follows in JSON: {{"experiment step":{{"step1":"Amplification of target fragment", "step2":"purification", "step3":"Amplification and adapter ligation", "step4":"purification", "step5":"make dnb"}}}}. 
You only need to return the JSON file without any additional descriptions!"""
                
        conversation_str = "\n\n".join([i["role"] + ':\n' + i['content'] for i in conversation])
        decomposed_description = openai.ChatCompletion.create(
            messages=[{"role": "system", "content": decomposed_prompt},
                        {"role": "user", "content": f"""{decomposed_prompt}"""+"#conversation:"+conversation_str+f"#protocol:"+protocol_select+f"""{decomposed_prompt}"""}
                        ],
            deployment_id="XMGI-Chat1-GPT4"
        )
        decomposed_description_data = extract_json_response(decomposed_description)
        print("-------------decomposed description information:\n", decomposed_description_data)
        task_flow_response = '\n'.join([str(x)+": "+str(y) for x,y in decomposed_description_data.items()])
        response = f"The process for the task is as follows. Please confirm whether it is correct. If it is correct, proceed with the subsequent script generation tasks according to the workflow. If it is incorrect or information is missing, please provide the necessary additions.
"
        response += task_flow_response
        stage += 1

        return {
                "response": response, 
                "operations":[],
                "stage": stage,
                "state": 'continue',
                "protocol_doc": relevance_doc,
                "pre_calculate_parameter_data": experiment_parameter,
                "decomposed_description_data": decomposed_description_data,
                "extraction_protocol_parameter":extraction_protocol_parameter}


    if stage == 5:
        #stage5 Based on user feedback, check whether new content needs to be added. Simultaneously, begin the breakdown of the SOP task flow, retrieve template code, decompose template code blocks, and retrieve and recall corresponding block code.

        conversation = instruction['conversation']
        protocol_select = instruction["protocol_doc"]
        decomposed_description = instruction["decomposed_description_data"]
                
        judge_prompt = f"""Based on the context and the user's current feedback, 
        determine if there is any need for modifications or additional information. 
        If necessary, make the additions in the task flow JSON file provided below. 
        You only need to return "no" (if the user indicates no need for changes or additions) 
        or the updated task flow JSON file. Do not include any additional descriptions.
"""
                
        conversation_str = "\n\n".join([i["role"] + ':\n' + i['content'] for i in conversation])
        judge_description = openai.ChatCompletion.create(
            messages=[{"role": "system", "content": judge_prompt},
                        {"role": "user", "content": f"""{judge_prompt}"""+"#conversation:"+conversation_str+f"#task flow JSON:"+decomposed_description+f"""\n{judge_prompt}"""}
                        ],
            deployment_id="XMGI-Chat1-GPT4"
        )
        if 'no' in judge_description:
            #开始进行代码块检索设计
            task_list = [(x,y) for x,y in decomposed_description.tiems()]
            current_task = task_list[0]
            current_task_descrption = f"The current experimental step is the {current_task[0]}, and the task is {current_task[1]}."
            summary_prompt = f"""Based on the context, the user's current feedback, and the current task description, 
            please extract a summary that can determine which type of code script is needed. 
            Currently, there are two types of code scripts: a free-step script and a two-step script. 
            Keep the description under 250 words."""
            conversation_str = "\n\n".join([i["role"] + ':\n' + i['content'] for i in conversation])
            code_retrieval_response = openai.ChatCompletion.create(
            messages=[{"role": "system", "content": summary_prompt},
                    {"role": "user", "content": f"""{summary_prompt}"""+"#conversation:"+conversation_str+f"#current task description:"+current_task_descrption+f"""\n{summary_prompt}"""}
                    ],
                deployment_id="XMGI-Chat1-GPT4"
            )

            code_retrieval = get_retrieval(code_retrieval_response)
            
            #Decomposing template code
            code_decomposed_prompt = """The tasks will be broken down step-by-step in the provided Opentrons OT-2 script below. The decomposed tasks will be returned to me in the JSON format example for the first round of a two-step PCR reaction as shown below. Note that the comments in the OT-2 script clearly specify the code task steps; do not include any additional descriptions.
```json
[
  {"step": 1, "description": "Definition of models and labware placement"},
  {"step": 2, "description": "Reagent definition, which includes correct entry of round and accurate definition of volumes"},
  {"step": 3, "description": "Reagent loading, specifying the description of liquid loading, volumes loaded, and target wells"},
  {"step": 4, "description": "Pipetting, based on the loaded contents, specifies the pipetting actions, settings of temperature control modules, specific pipetting actions, and speed of aspiration and dispensing"},
  {"step": 5, "description": "PCR reaction configuration, which starts the thermal cycler after the above tasks are completed, initiating the PCR reaction. This requires input of the corresponding temperatures, cycle times, number of cycles, and thermal cycler well capacity volumes"}
]```"""
            code_decomposed_response = openai.ChatCompletion.create(
                        messages=[{"role": "system", "content": code_decomposed_prompt},
                                  {"role": "user", "content": f"""{code_decomposed_prompt}"""+f"#Opentrons OT-2 script:"+code_retrieval+f"""\n{code_decomposed_prompt}"""}
                                ],
                deployment_id="XMGI-Chat1-GPT4"
            )
            code_decomposed_data = extract_json_response(code_decomposed_response)
            print("-------------code decomposed information:\n", code_decomposed_data)

            
        else:
            judge_description_data = extract_json_response(judge_description)
            print("-------------djudge_description information:\n", judge_description_data)
            #retrieval task code
            decomposed_description = judge_description_data
            task_list = [(x,y) for x,y in decomposed_description.tiems()]
            current_task = task_list[0]
            current_task_descrption = f"The current experimental step is the {current_task[0]}, and the task is {current_task[1]}."
            summary_prompt = f"""Based on the context, the user's current feedback, and the current task description, 
            please extract a summary that can determine which type of code script is needed. 
            Currently, there are two types of code scripts: a free-step script and a two-step script. 
            Keep the description under 250 words."""
            conversation_str = "\n\n".join([i["role"] + ':\n' + i['content'] for i in conversation])
            code_retrieval_response = openai.ChatCompletion.create(
            messages=[{"role": "system", "content": summary_prompt},
                    {"role": "user", "content": f"""{summary_prompt}"""+"#conversation:"+conversation_str+f"#current task description:"+current_task_descrption+f"""\n{summary_prompt}"""}
                    ],
                deployment_id="XMGI-Chat1-GPT4"
            )

            code_retrieval = get_retrieval(code_retrieval_response)
            
            #Decomposing template code
            code_decomposed_prompt = """The tasks will be broken down step-by-step in the provided Opentrons OT-2 script below. The decomposed tasks will be returned to me in the JSON format example for the first round of a two-step PCR reaction as shown below. Note that the comments in the OT-2 script clearly specify the code task steps; do not include any additional descriptions.
```json
[
  {"step": 1, "description": "Definition of models and labware placement"},
  {"step": 2, "description": "Reagent definition, which includes correct entry of round and accurate definition of volumes"},
  {"step": 3, "description": "Reagent loading, specifying the description of liquid loading, volumes loaded, and target wells"},
  {"step": 4, "description": "Pipetting, based on the loaded contents, specifies the pipetting actions, settings of temperature control modules, specific pipetting actions, and speed of aspiration and dispensing"},
  {"step": 5, "description": "PCR reaction configuration, which starts the thermal cycler after the above tasks are completed, initiating the PCR reaction. This requires input of the corresponding temperatures, cycle times, number of cycles, and thermal cycler well capacity volumes"}
]```"""
            code_decomposed_response = openai.ChatCompletion.create(
                        messages=[{"role": "system", "content": code_decomposed_prompt},
                                  {"role": "user", "content": f"""{code_decomposed_prompt}"""+f"#Opentrons OT-2 script:"+code_retrieval+f"""\n{code_decomposed_prompt}"""}
                                ],
                deployment_id="XMGI-Chat1-GPT4"
            )
            code_decomposed_data = extract_json_response(code_decomposed_response)
            print("-------------code decomposed information:\n", code_decomposed_data)
        
        response = f"Based on our conversation, I have identified the content of the experiment you need to conduct. Next, I will proceed with the script block retrieval and synthesis for the first step of the task. The breakdown of the first step's script content is as follows::{code_decomposed_data}"        
        response += "\nWe are generating the corresponding Opentrons OT-2 execution script. Please be patient. Once the script is generated, I will display the code to you for confirmation."
        stage += 1
        return {
                "response": response, 
                "operations":[],
                "stage": stage,
                "state": 'continue',
                "protocol_doc": relevance_doc,
                "pre_calculate_parameter_data": experiment_parameter,
                "decomposed_description_data": decomposed_description_data,
                "code_decomposed_data": code_decomposed_data,
                "task_list":task_list,
                "extraction_protocol_parameter":extraction_protocol_parameter}


    if stage == 6:
        #stage6 Remind the user that the RAG script generation is slow and to patiently wait for the code to be produced.
        conversation = instruction['conversation']
        protocol_select = instruction["protocol_doc"]
        decomposed_description = instruction["decomposed_description_data"]
        pre_calculate_parameter_data = instruction["pre_calculate_parameter_data"]
        
        #rertieval decomposed code block
        code_retrieval_list = []
        code_rewrite_list = []
        for idx,x in enumerate(code_decomposed_data):
            code_decomposed_description = f"The task breakdown is for step {idx+1}, and its description is {x["description"]}."
            code_retrieval = get_retrieval(code_decomposed_description)
            code_retrieval_list.append(code_retrieval)

            #Decomposing template code
            code_rewrite_prompt = """You will receive the template code and a parameters information JSON. 
            The template code contains some variables that need to be modified. 
            These variables are detailed in the parameters information JSON. 
            You will need to extract the values for the corresponding variable fields from the JSON and rewrite the template code. 
            Note! You can only modify fields in the template code that appear in the JSON; variable names not present in the JSON must not be changed. 
            You only need to return the modified code to me, without any additional description!"""

            code_rewrite_response = openai.ChatCompletion.create(
                        messages=[{"role": "system", "content": code_rewrite_prompt},
                                  {"role": "user", "content": f"""{code_rewrite_prompt}"""+f"##Opentrons OT-2 script:"+code_retrieval+f"##parameters information JSON:\n{pre_calculate_parameter_data}"}
                                ],
                deployment_id="XMGI-Chat1-GPT4"
            )
            code_rewrite_data = extract_python_response(code_rewrite_response)
            print("-------------code rewrite information:\n", code_rewrite_data)
            code_rewrite_list.append(code_rewrite_data)
        
        all_rewrite_block_code = '\n'.join([f"###rewrite block code {str(idx+1)}:\n{x}" for idx,x in enumerate(code_rewrite_list)])
        merging_prompt = """You will receive several rewritten executable sub-blocks of code for the Opentrons OT-2, referred to as rewrite block code. 
        You are not to modify the variables within these sub-blocks. Your task is to sequentially assemble these sub-blocks into a runnable complete code. 
        You only need to return the complete code to me, without any additional description."""

        code_merging_response = openai.ChatCompletion.create(
                        messages=[{"role": "system", "content": merging_prompt},
                                  {"role": "user", "content": f"""{merging_prompt}"""+f"##Opentrons OT-2 rewrite_block_code:{all_rewrite_block_code}"}
                                ],
                deployment_id="XMGI-Chat1-GPT4"
            )
        code_merging_data = extract_python_response(code_merging_response)
        print("-------------code merging information:\n", code_merging_data)
        
        stage += 1        
        experiment_description = f"""Currently conducting step-{task_list[0]["step"]} of the experiment. Upon completion, you must inform the user that it has been completed. 
        If an anomaly occurs, alert the user and provide updates on the progress of subsequent experiments as they unfold."""

        task_list = task_list[1:]
        response = f"""Based on our communication, I have generated the following script code. Please confirm whether the code is correct. 
        If it is, we will submit it to the Experiment Agent to initiate the Opentrons experiment. 
        Note! Before starting the experiment, please ensure that materials, equipment, and other instruments are correctly positioned on the OT-2 device!
        OT-2 running script:\n{code_merging_data}"""

        experiment_description = f"""Currently conducting step-{} of the experiment. Upon completion, you must inform the user that it has been completed. 
        If an anomaly occurs, alert the user and provide updates on the progress of subsequent experiments as they unfold."""

        #protocol agent 和 code excute agent communcation term：
        return {
                "response": response, 
                "operations":[{"title":'experiment',"finished": "false","description":experiment_description, "task_list":task_list}],
                "stage": stage,
                "state": 'continue',
                "protocol_doc": relevance_doc,
                "pre_calculate_parameter_data": experiment_parameter,
                "decomposed_description_data": decomposed_description_data,
                "code_decomposed_data": code_decomposed_data,
                "task_list":task_list,
                "code_merging_data":code_merging_data,
                "extraction_protocol_parameter":extraction_protocol_parameter}



    if stage == 7:
        #stage7
        conversation = instruction['conversation']
        protocol_select = instruction["protocol_doc"]
        decomposed_description = instruction["decomposed_description_data"]
        pre_calculate_parameter_data = instruction["pre_calculate_parameter_data"]
        

        
        return {
                "response": response, 
                "operations":[],
                "stage": stage,
                "state": 'continue',
                "protocol_doc": relevance_doc,
                "pre_calculate_parameter_data": experiment_parameter,
                "decomposed_description_data": decomposed_description_data,
                "code_decomposed_data": code_decomposed_data,
                "task_list":task_list,
                "extraction_protocol_parameter":extraction_protocol_parameter}


        
        




