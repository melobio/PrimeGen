import time
from Bio import Entrez, SeqIO, SeqRecord
from Bio.Seq import Seq
import pandas as pd
import subprocess
import sys
import openai
import requests
import re
import json
import os
import zipfile
import glob
import shutil
import traceback


class CustomError(Exception):
    pass

# species identification search tool
def species_identification_search(instruction, stage):
    # 第一步确定用户是否已提供相关信息，根据提供的信息量选择不同的第二步
    if stage == 1:
         
        if 'search_type' not in instruction['instruction']:
            conversation = instruction['conversation']
            #初始化物种鉴定需要的Search内容
            species_identification_dict = {
                    "Experiment_Direction": "species identification",
                    "Experiment_Type": None, #General or Group or Microbe
                    "organism_names": "",
                    "genome_path":"",
                    "target_cds_path": "", #ALL
                    "non_target_cds_path": "", #ALL
                    "target_gene_path": "", #General or Microbe
                    "target_number":"" #single or multi
                    }

            # 尝试提取用户语言中的信息：实验类型、物种名称、菌株名称
            experiment_type_prompt = """Based on the user's contextual dialogue, you need to perform information extraction. The information to be extracted includes:  
1.Experiment type: Extract if 'precise' mentioned by the user("Microbe" or "Group" or "General"), return null if not 'precise' mentioned. The "General" primers described here are designed primarily for the reference genome of the species, such as Escherichia coli, which has been sequenced in detail and used as a baseline for research."Group" here is defined as a collection of several microbes that may come from the same genus or different genera.3) "Microbe" applies to a single colony whose entire genome has been sequenced, such as E. coli MC1061, a substrain of E. coli.

2.Organism Names: Extract all biological names mentioned by the user, including species and strain names. If the user provides an abbreviation, please complete it; return null if no biological name is mentioned. 
The results should be returned to me in JSON format. Example is shown below:
'''json
{  
  "Experiment_type": "Group",  
  "Organism_Names": ["Escherichia coli", "E. coli MC1061",...]
}
'''
Note:
1)The name of the organism must be complete and professional. Please fill in the complete and professional vocabulary. If it is a variant, strain, isolate, etc., there is no need to modify it and just extract it directly. Remember not to use brackets or other additional descriptions.
2)Notice! If the user hasn't mentioned it, fill in null, and don't speculate on the user's intentions!
"""
            conversation_str = "\n\n".join([i["role"] + ':\n' + i['content'] for i in conversation])
        
            #JSON MODE输出
            experiment_info_response = openai.ChatCompletion.create(
                messages=[{"role": "system", "content": experiment_type_prompt},
                          {"role": "user", "content": "##user's contextual dialogue:\n"+"'''\n"+conversation_str+"\n'''"}
                          ],
                response_format={"type": "json_object"},
                deployment_id="gpt-4o",
            )
            #结构化提取
            experiment_info_json_data = extract_json_response(experiment_info_response)
            #填充提取出的数据
            species_identification_dict["organism_names"] = experiment_info_json_data["Organism_Names"]
            #如果有物种的名称，引导用户选择实验类型
            if len(species_identification_dict["organism_names"])!=0:
                failed_list = []
                right_list = []
                for s_n in species_identification_dict["organism_names"]:
                    print(s_n)
                    if type(s_n)==list:
                        s_n=s_n[0]
                    search_result = get_embedding_recommended(s_n, type = "determine_species")
                    if type(search_result[0]) == str:
                        right_list.append(search_result)
                    else:
                        failed_list.append((s_n,search_result))
                
                if len(failed_list) != 0:
                    operations_list = []
                    species_identification_dict["organism_names"] = right_list
                    for i,f_l in enumerate(failed_list):
                        operations_list.append({"title":f"The {f_l[0]} is unclear, which may not be accurate. Please confirm again.(required)", "options":f_l[-1][0], "type":["single","required"],"key":f"species_select_{i+1}", "value":[]})
                    response = """We would like to inform you that based on the biological information search, the species name is unclear. Please complete the confirmation below, and we will conduct further experimental search content."""
            
                    return {"response": response, "operations": operations_list,"stage": stage,"state": 'continue',
                            "search_type": "species_identification_type", "species_identification_dict": species_identification_dict, 'type':'check_organism'}
                else:
                    species_identification_dict["organism_names"] = right_list
            else:
                # If the user does not specify the species name, the search will be interrupted, the state is stop, and the controller will ask further questions for consultation.
                user_response = "I was unable to find the species or strain information you provided. You need to provide this information so that I can proceed with gene sequence retrieval. At a minimum, you need to tell me the species name, and I will guide you through the subsequent steps."
                return {"response": user_response,
                        "operations":[],
                        "stage": stage+2,
                        "state": 'stop', "search_type": "controller_type",
                        "species_identification_dict": species_identification_dict}
        else:

            species_identification_dict = instruction['instruction']['species_identification_dict']
            species_selected_values = [value for key, value in instruction['instruction'].items() if key.startswith("species_select_")]
            if len(species_selected_values) >= 1:
                search_result = get_embedding_recommended(species_selected_values[0][0], type = "determine_species")
                species_identification_dict["organism_names"].append(search_result)
                
        # 物种不含有菌株，个别任务不能做，如microbe
        if species_identification_dict["organism_names"][-1] == "False" or "False" in species_identification_dict["organism_names"][-1]:
            response = """Currently, the species you mentioned does not have strains or lower levels, so only General and species-level Group experiments are supported. The “General” primer designs described here are mainly targeted at single species or discriminatory primer designs for multiple species."""
            return {"response": response, "operations":[],#{"title":"Select your experiment type (required)", "options":["general", "group"], "type":["single","required"],"key":"experiment_select", "value":[]},],
                    "stage": stage+1,"state": 'continue',
                    "search_type": "species_identification_type",
                    "species_identification_dict": species_identification_dict}
                
        #物种含有菌株，即所有任务可做
        else:
            response = """Would you like to design primers for type of "Gneral" or "Microbe"? Your choice will affect the subsequent retrieval and design task workflows.
1) The "General" primers described here are designed primarily for the reference genome of the species, such as Escherichia coli, which has been sequenced in detail and used as a baseline for research.
2) "Microbe" applies to a single colony whose entire genome has been sequenced, such as E. coli MC1061, a substrain of E. coli."""
            return {"response": response, "operations":[{"title":"Select your experiment type (required)", "options":["general", "microbe"], "type":["single","required"],"key":"experiment_select", "value":[]},],
                    "stage": stage+1,"state": 'continue',
                    "search_type": "species_identification_type",
                    "species_identification_dict": species_identification_dict}



    # 第二步根据用户选择的实验类型，选择不同的search路径
    if stage == 2:
        # Stage1已完成物种名称和实验类型的确定，stage2要实现1）基于name2class，判断物种类型；2）基于类型不同，给出不同search路线和推荐
        species_identification_dict = instruction['instruction']['species_identification_dict']
        experiment_type = instruction['instruction']['experiment_select'][0]
        species_identification_dict["Experiment_Type"] = experiment_type
        Species_name_list = species_identification_dict["organism_names"]        
        
        # 先判断物种类型
        species_type_list = []
        taxid_list = []
        name_list = []
        print(Species_name_list)
        for s_n, s_t, taxid, isolate in Species_name_list: #物种名、物种类型、isolate推荐
            species_type_list.append(s_t)
            taxid_list.append(taxid)
            name_list.append(s_n)
        

        # 用户选择General，在物种层面设计引物
        # TODO：让用户可参与上传rna、cds、genome文件
        if experiment_type == 'general':
            # 首先判断是单物种
            if len(species_identification_dict["organism_names"])==1:
                # 如果是 fungi或protozoa 只需要在下载Genome
                # TODO:待测试——7-11
                if species_type_list[0] in ["fungi", "protozoa"]:
                    species_path_list = []
                    failed_list = []
                    species = name_list[0]
                    taxid = taxid_list[0]
                    print(f'*********************{species} searching***************************')
                    #原生生物的search 直接考虑gene，没有则直接genome
                    species_file_path = ncbi_search(species, 'genome', taxid=taxid, target_type="target", search_type="rna_species_gene")
                    species_path_list = species_file_path
                    strain_response = "The current sequence search requirements have been completed. The following is the Species reference genome file."
                    species_identification_dict["genome_path"] = species_path_list
                    species_identification_dict["target_number"] = 'single'         
                    #core_cds_path = get_core_gene( ,type="Eukaryotes") #这个只有基因组
                else:
                    cds_result = get_embedding_recommended(taxid_list[0], type = "determine_cds", species_class=species_type_list[0])
                    # 判断是否存在cds文件
                    type_ = set(cds_result)
                    species_path_list = []
                    failed_list = []
                    species = species_identification_dict["organism_names"][0]
                     
                    print(f'*********************{species} species level searching***************************')
                    if "cds" in type_:
                        species_file_path = ncbi_search(species[0], 'genome', taxid=species[-2], target_type="target", search_type="species_cds")
                    elif "cds" not in type_ and "gtf" in type_:
                        species_file_path = ncbi_search(species[0], 'genome', taxid=species[-2], target_type="target", search_type="species_gene")
                    else:
                        species_file_path = ncbi_search(species[0], 'genome', taxid=species[-2], target_type="target", search_type="species_genome")
                    
                    species_path_list.append(species_file_path)
                    species_path_list = [x[0] for x in species_path_list]
                    strain_response = "The current sequence search requirements have been completed. The following is the Species reference genome file."
                    species_identification_dict["target_cds_path"] = species_path_list
                    species_identification_dict["target_number"] = 'single'
                    core_cds_path = get_core_gene( species_file_path,taxid=[taxid], type_cds="have_cds")
               
                species_identification_dict["target_cds_path"] = core_cds_path
                #user_response = strain_response + """For current General experiments (single species case), which has completed the acquisition of sequence information at the species level."""
                user_response = "Thus far, through our communication, we have confirmed the information regarding both cds gene. Next, I will inform the Controller to communicate with you to determine the next steps related to primer design."
                species_identification_dict["target_number"] = 'single'
                print('species_identification_dict:\n', species_identification_dict)
                #TODO:multi
                return {"response": user_response, "stage": stage+2,"operations":[],
                        "core_cds_path":[core_cds_path], "state": 'stop', 
                        "search_type": "species_identification_type", 
                        "species_identification_dict": species_identification_dict
                        }
            else:
                species_path_list = []
                core_cds_path_list=[]
                taxid_list = []
                for s_n, s_t, taxid, _  in species_identification_dict["organism_names"]:
                    # TODO：待测试——7-11
                    if s_t == "fungi" or s_t == "protozoa":
                        print(f'*********************{s_n} searching***************************')
                        species_file_path = ncbi_search(s_n, 'genome', taxid=taxid, target_type="target", search_type="rna_species_gene")
                        species_path_list.append(species_file_path)
                            
                    
                    else:
                        #search_result = get_embedding_recommended(species_identification_dict["Organism_Names"][0], type = "dna")
                        cds_result = get_embedding_recommended(taxid, type = "determine_cds", species_class=s_t)
                        type_ = set(cds_result)
                        print(f'*********************{s_n} searching***************************')
                        if "cds" in type_:
                            species_file_path = ncbi_search(s_n, 'genome', taxid=taxid, target_type="target", search_type="species_cds")
                        elif "cds" not in type_ and "gtf" in type_:
                            species_file_path = ncbi_search(s_n, 'genome', taxid=taxid, target_type="target", search_type="species_gene")
                        else:
                            species_file_path = ncbi_search(s_n, 'genome', taxid=taxid, target_type="target", search_type="species_genome")
                        species_path_list.append(species_file_path)
                        taxid_list.append(taxid)

                species_path_list = [x[0] for x in species_path_list]
                
                for species_path,taxid  in zip(species_path_list,taxid_list):

                    core_cds_path_list.append(get_core_gene(file_path=[species_path],taxid=[taxid], type_cds="have_cds", type_experiment = experiment_type))

                strain_response = "The current sequence search requirements have been completed. The following is the Species reference genome file."
                species_identification_dict["target_number"] = 'multi'
                species_identification_dict["target_cds_path"] = core_cds_path_list

                response_str = """Currently it is a General implementation, which has completed the acquisition of sequence information at the species level.Next,  I will inform the Controller to communicate with you to determine the next steps related to primer design."""
                user_response = strain_response + response_str
                
                return {"response": user_response, "stage": stage+2,"operations":[],
                    "core_cds_path": core_cds_path_list, "state": 'stop', 
                    "search_type": "species_identification_type","species_identification_dict": species_identification_dict
                    }

        
        # 用户选择Microbe模式，设计特异性引物
        else:
            # 用户选择，以此判断是单菌株还是多菌株情况
            # TODO:如果多个microbe情况下，直接下载cds，不需要推荐
            stage += 1
            response = f"You've chosen to conduct an experiment focusing on {species_identification_dict['Experiment_Type']}, specifically targeting the {species_identification_dict['Species_Name'][0][0]} species."
            
            # 用户未提及菌株
            # TODO: Microbe 也要区分真菌、原虫和其他类别
            if species_identification_dict["Strain_Name"] is None:
                recommended_strains_json_data = {"recommended_strains":species_identification_dict["organism_names"][0][-1]}
                response += """You are currently preparing to conduct a microbe experiment. I will recommend several candidate microbial strains for your selection. 
Alternatively, you may upload the .fna file of the strain you intend to investigate. 
[Note: The strain file should be in .fna format]. 
Upon confirmation of the selected strain, we will suggest potential target genes that may be required for your study. 
You have the option to utilize these target genes to enhance the design of your primers, 
or you may upload your own target gene file. [Note: The target gene file should also be in .fna format]. 
Alternatively, if you opt not to use specific target genes, our AI algorithm will automatically identify relevant target genes and design the necessary primers for your experiment."""
                
                return {"response": response, "operations":[{"title":"Recommend Strain: Choose 'Recommend Strain' or upload your own below","options":recommended_strains_json_data["recommended_strains"], "type":["multi","input","required"],"key":"strain_select","value":[]}, 
                       {"title":"If you have Strain cds.fna files, please upload them", "options":[], "type":["file"], "key":"strain_path", "value":[]}],
                        "stage": stage,
                        "state": 'continue', 
                        "search_type": "species_identification_type",
                        "species_identification_dict": species_identification_dict}
            
            # 用户已提及菌株
            # TODO：check bug
            else:
                response += f"You intend to study the strain {species_identification_dict['Strain_Name']}. "
                response += """we will suggest potential target genes that may be required for your study. You have the option to utilize these target genes to enhance the design of your primers, or you may upload your own target gene file. [Note: The target gene file should also be in .fna format]. Alternatively, if you opt not to use specific target genes, our AI algorithm will automatically identify relevant target genes and design the necessary primers for your experiment."""

                return {"response": response, "operations":[],
                        "stage": stage,"state": 'continue', "search_type": "species_identification_type",
                        "species_identification_dict": species_identification_dict}
    elif stage == 3:

        conversation = instruction['instruction']['conversation']

        if instruction['instruction']!=[]:
            if len(instruction['instruction']['operations'])==3:
                operation_temp = instruction['instruction']['operations'][1]['options'] 
            elif len(instruction['instruction']['operations'])==2:
                operation_temp = instruction['instruction']['operations'][0]['options']

        species_identification_dict = instruction['instruction']['species_identification_dict']
        experiment_type = species_identification_dict['Experiment_Type']
        
            
        #Group或Microbe 下载菌株的cds 
        upload_file_flag = instruction['instruction']['upload_file_flag']
        if upload_file_flag == True:
            strain_list = []
            strain_path_list = instruction['instruction']['strain_path']
            for strain_path in strain_path_list:
                strain_list.append(strain_path.split('/')[0].split('_cds')[0].replace("_", " "))
            species_identification_dict["Strain_Name"] = strain_list

            strain_file_path = instruction['instruction']['strain_path']
            strain_response = ""
            stage += 1

        else:

            select_strain_list = instruction['instruction']['strain_select'] 
            stage += 1
            cds_list=[]
            idx_list=[]
            species_class = species_identification_dict["Species_Name"][0][1]

            for strain_ in select_strain_list:
                idx = get_embedding_recommended(strain_, type = "get_isolate2id", species_class = species_class)
                idx_list.append(idx)
                cds_result = get_embedding_recommended(idx, type = "determine_cds")
                cds_list.append(cds_result)
            
            strain_path_list = []
            failed_list = []
            recommend_target_list = []
            #TODO: 这里的strain_name 包含名称、taxid和down_class
            for strain_, taxid, cds_flag in zip(select_strain_list, idx_list, cds_list):
                print(f'********************{strain_} searching***************************')
                if cds_flag:
                    recommend_target, strain_file_path = ncbi_search(strain_, 'genome', taxid=taxid, target_type="target", search_type="species_cds") #cds/gene/genome
                    strain_path_list.append(strain_file_path)
                    recommend_target_list.append(recommend_target)
                else:
                    recommend_target, strain_file_path = ncbi_search(strain_, 'genome', taxid=taxid, target_type="target", search_type="strain_genome")  #gene/genome
                    strain_path_list.append(strain_file_path)
                    recommend_target_list.append(recommend_target)
                
            strain_path_list = [x[0] for x in strain_path_list]
            strain_response = "The current sequence search requirements have been completed. The following is the strain reference genome file."
            species_identification_dict["target_cds_path"] = strain_path_list
        

        # Group、多物种General和多菌株Microbe，直接跳过target阶段
        if (species_identification_dict["Experiment_Type"].lower() == 'microbe' or species_identification_dict["Experiment_Type"].lower() == 'group') and len(species_identification_dict["target_cds_path"]) > 1:
            # non-target recommend
            non_target_gene_recommended_system_prompt = f"""You are conducting an {species_identification_dict["Experiment_Type"]} experiment on the {species_identification_dict["Species_Name"]} species. 
            To ensure the efficiency of the primer design, I will Recommend the top-20 bacteria and viruses that inhabit the same parasitic sites as {species_identification_dict["Species_Name"]}. 
            This will help you design primers with higher specificity.Below, I will provide recommendations in JSON format. 
            Only the JSON recommendations are required, no additional description is needed.The JSON format example is as follows:  
            '''json{{
             "recommendations non-target": [  
             "**",    
             "**", 
             "**"    
             ]    
             }}'''
             Note, recommendations for non-target species should be tailored based on the user's experiment type and target species, and should not be generic.
             """
            
            non_target_gene_recommended_response = openai.ChatCompletion.create(
                  messages=[{"role": 'system', 'content': non_target_gene_recommended_system_prompt},
                            {"role": "user", 'content': non_target_gene_recommended_system_prompt}],
                  response_format={"type": "json_object"},
                  deployment_id="gpt-4o",
                )
            
            if species_identification_dict["Experiment_Type"].lower() == 'group':
               response_str = """Currently it is a Group implementation, which has completed the acquisition of sequence information at the strain/isolute level. 
               Next, non-target information will be recommended, and you can make choices to improve the effect of the final primer design."""
            elif species_identification_dict["Experiment_Type"].lower() == 'microbe':
               response_str = """Currently it is a Micorbe implementation, which has completed the acquisition of sequence information at the strain/isolute level. 
               Next, non-target information will be recommended, and you can make choices to improve the effect of the final primer design."""
             
            else:
               response_str = """Currently it is a General implementation, which has completed the acquisition of sequence information at the species level.
                              Next, non-target information will be recommended, and you can make choices to improve the effect of the final primer design."""
                   
            non_target_gene_recommended_json_data = extract_json_response(non_target_gene_recommended_response)
            user_response = strain_response + response_str

            return {"response": user_response, "stage": stage,"operations":[
                {"title":"recommendations non-target genome/species","options":non_target_gene_recommended_json_data["recommendations non-target"],"type":["multi","input"],"key":"non_target_genome/species_select","value":[]},
                {"title":"Upload your non-target CDS .fna files if you have them.(Optional)", "options":[], "type":["file"], "key":"non_target_gene_path", "value":[]}],
                 "strain_path": strain_path_list, "state": 'continue', "search_type": "species_identification_type","species_identification_dict": species_identification_dict
                }
            
            
        # 单菌株Microbe 或 单物种General情况
        elif species_identification_dict["Experiment_Type"].lower() == 'microbe':
            #TODO:优化修改
            non_target_gene_recommended_system_prompt = f"""You are conducting an {species_identification_dict["Experiment_Type"]} experiment on the {species_identification_dict["Species_Name"]} species.
            To ensure the efficiency of the primer design, I will Recommend the top-20 bacteria and viruses that inhabit the same parasitic sites as {species_identification_dict["Species_Name"]}.
            This will help you design primers with higher specificity.Below, I will provide recommendations in JSON format.
            Only the JSON recommendations are required, no additional description is needed.The JSON format example is as follows:
            '''json{{
            "recommendations non-target": [
            "**",
            "**",
            "**"
            ]
            }}'''
            Note, recommendations for non-target species should be tailored based on the user's experiment type and target species, and should not be generic."""
            
            non_target_gene_recommended_response = openai.ChatCompletion.create(
                    messages=[{"role": 'system', 'content': non_target_gene_recommended_system_prompt},
                        {"role": "user", 'content': non_target_gene_recommended_system_prompt}],
                    response_format={"type": "json_object"},
                    deployment_id="gpt-4o",
                    )
            non_target_gene_recommended_json_data = extract_json_response(non_target_gene_recommended_response)
            user_response = strain_response + """For current Microbe experiments, you can provide the target gene of the strain, 
which will increase the success rate of primer design.
You can refer to the target gene content I recommend, or you can specify a specific target gene and provide its relevant information(gene name, like "lacZ"), and 
I will search it on NCBI. If you do not provide the target gene, non-target information must be provided next because specific primers need to be designed."""
            species_identification_dict["target_number"] = 'single'
            #TODO:multi
            return {"response": user_response, "operations":[{"title":"recommend target genes","options":recommend_target_list[0][:100],"type":["multi","input"],"key":"target_gene_select","value":[]},
                    {"title":"upload target genes", "options":[], "type":["file"], "key":"target_gene_path", "value":[]},
                    {"title":"recommendations non-target genome/species","options":non_target_gene_recommended_json_data["recommendations non-target"],"type":["multi","input"],"key":"non_target_genome/species_select","value":[]},
                    {"title":"Upload your non-target CDS .fna files if you have them.(Optional)", "options":[], "type":["file"], "key":"non_target_gene_path", "value":[]}],
                    "strain_path": strain_path_list, "stage": stage,
                    "state": 'continue', "search_type": "species_identification_type",
                    "species_identification_dict": species_identification_dict}
    # 基于Non-target的选择，进行下载沟通
    elif stage == 4:
        if instruction['instruction']!=[]:
            operation_temp = instruction['instruction']['operations'][0]['options'] #back
            
        conversation = instruction["instruction"]['conversation']
        upload_file_flag = instruction["instruction"]['upload_file_flag']  # If the user uploads it, it is True. I don’t need to go through the search process.
        species_identification_dict = instruction["instruction"]['species_identification_dict']
        target_path = species_identification_dict["target_cds_path"]
        #if 
        user_select = instruction["instruction"]['non_target_genome/species_select']
        
        if "target_gene_select" in instruction["instruction"]:
            if instruction["instruction"]["target_gene_select"]!=[]:#species_identification_dict["Experiment_type"]=="general":
                print("开始提取gene数据")
                _, target_gene_path = ncbi_search(instruction["instruction"]["target_gene_select"][0], 'gene', taxid='', cds_path=target_path[0] , target_type="target", search_type="strain_cds")
                species_identification_dict["target_gene_path"] = target_gene_path
        
        # User uploads non-target
        # 用户上传non-target
        if instruction["instruction"]['non_target_gene_path']!=[] and instruction["instruction"]['non_target_gene_path']!=['']:
            user_response = """Thus far, through our communication, we have confirmed the information regarding both target and non-target. 
            Next, I will inform the Controller to communicate with you to determine the next steps related to primer design."""
            species_identification_dict["non_target_cds_path"] = instruction["instruction"]['non_target_gene_path']
            print("finish work:",species_identification_dict)
            return {"response": user_response, "non_target_path": instruction["instruction"]['non_target_gene_path'], "operations":[],
                    "stage": stage, "state": 'stop', "search_type": "species_identification_type",
                    "species_identification_dict": species_identification_dict}

        # User selects non-target
        else:
            print("User selects non-target")
            stage += 1
            non_target_search_info_list = []
            non_target_path_list = []
            non_target_failed_list = []
            for non_tar_ in user_select:
                try:
                    print('++++++++++用户选择的non-target：', non_tar_)
                    non_target_name, non_target_path = ncbi_search(non_tar_, "genome",target_type="non_target", search_type="non-target_genome")
                    non_target_search_info_list.append(non_target_name)
                    non_target_path_list.append(non_target_path)
                except:
                    non_target_failed_list.append(non_tar_)
                    traceback.print_exc()
                    continue

                if len(non_target_search_info_list) == 0:
                    if len(non_target_failed_list) == 1:
                        response = f"Sorry, the target gene sequence({non_target_failed_list[0]})) was not found. You can manually upload the target gene cds files or select other target."
                    else:
                        response = f"Sorry, the target gene sequence({','.join(non_target_failed_list)})) was not found. You can manually upload the target gene cds files or select other target."
                    #TODO:
                    return {"response":response,
                            "operations":[{"title":"recommend non-target information","options":operation_temp,"type":["multi","input","required"],"key":"non_target_genome/species_select","value":[]},
                                          {"title":"upload non-target gene", "options":[], "type":["file"], "key":"non_target_gene_path", "value":[]}],
                            "stage": stage-1,"state": 'continue', "search_type": "species_identification_type",
                            "species_identification_dict": species_identification_dict}
                
            user_response = f"""Based on the non-target information you selected, I have downloaded the following reference genomes:{non_target_path_list} 
            so we will proceed with primer design based solely on the target information."""
            user_response += """Thus far, through our communication, we have confirmed the information regarding both target and non-target. 
            Next, I will inform the Controller to communicate with you to determine the next steps related to primer design."""
            non_target_path_list = [x[0] for x in non_target_path_list]
            species_identification_dict["non_target_cds_path"] = non_target_path_list
            print("finish work:",species_identification_dict)
            return {"response": user_response, "non_target_path": non_target_path_list, "target_path": species_identification_dict["target_cds_path"], "operations":[],
                    "stage": stage, "state": 'stop', "search_type": "species_identification_type",
                    "species_identification_dict": species_identification_dict}


def ncbi_search(term, db_name, taxid=None, target_type=None, cds_path=None, search_type=None):
    #任务是下载全基因组和CDS文件
    if db_name == "genome":
        #taget gene场景
        if target_type=="target":
            temp_path = "target_zip_file"
        
        #non-target gene场景
        elif target_type=="non_target":
            temp_path = "non_target_zip_file"
        # taxid
        id = taxid
        datasets_limitation = ["--assembly-source RefSeq --reference", "--assembly-source RefSeq", "--reference'",""]
        failed_list = []
        #cds download
        if search_type == "species_cds":
            term = term.replace('/','_').replace("(","_").replace(")","_").replace(' ','_')
            make_file_str = f"/reference_data/{temp_path}/{term}"
            print('-------------------------------:',term)
            for idx in range(len(datasets_limitation)):
                #记录错误序列情况
                if len(failed_list) != 0:
                    print(f'*****************{failed_list} download error, try change limitation *****************')
                    #全部检索都无，则失败
                if idx == len(datasets_limitation):
                    return f'{term} gene search failed.',''
                #选择检索模式
                datasets_limi_str_ = datasets_limitation[idx]
                #当前序列已经下载过
                if os.path.exists(f"{make_file_str}/{term}.zip"):
                    file_path = f"{make_file_str}/{term}.zip"
                #当前序列未下载，开始进行下载
                else:
                    print("开始构建文件夹")
                    os.system(f"mkdir {make_file_str}")
                    #print(url)
                    url = f"datasets download genome taxon {id} {datasets_limi_str_} --include cds --filename {make_file_str}/{term}.zip"
                    print(url)
                    os.system(url)
                    if not os.path.exists(f"{make_file_str}/{term}.zip"):
                        failed_list += datasets_limitation[idx]
                        continue
                    else:
                        print('--------完成genome下载')
                        file_path = f"{make_file_str}/{term}.zip"
                        break

        
        if search_type == "dna_species_gene" or search_type == "rna_species_gene":
            term = term.replace('/','_').replace("(","_").replace(")","_").replace(' ','_')
            make_file_str = f"/reference_data/{temp_path}/{term}"
            if search_type == "dna_species_gene":
                url = f"datasets download gene taxon {id} --include gene --filename {make_file_str}/{term}.zip"
            else:
                url = f"datasets download genome taxon {id} --filename {make_file_str}/{term}.zip"
            
            os.system(url)
            if not os.path.exists(f"{make_file_str}/{term}.zip"):
                for idx in range(len(datasets_limitation)):
                    #记录错误序列情况
                    if len(failed_list) != 0:
                        print(f'*****************{failed_list} download error, try change limitation *****************')
                    #全部检索都无，则失败
                    if idx == len(datasets_limitation):
                        return f'{term} gene search failed.',''
                    #选择检索模式
                    datasets_limi_str_ = datasets_limitation[idx]
                    #当前序列已经下载过
                    if os.path.exists(f"{make_file_str}/{term}.zip"):
                        file_path = f"{make_file_str}/{term}.zip"
                    else:
                        os.system(f"mkdir {make_file_str}")
                        url = f"datasets download genome taxon {id} {datasets_limi_str_} --filename {make_file_str}/{term}.zip"
                        os.system(url)
                        if not os.path.exists(f"{make_file_str}/{term}.zip"):
                            failed_list += datasets_limitation[idx]
                            continue
                        #下载成功，跳出循环下载模式
                        else:
                            print('----------genome 下载完成')
                            file_path = f"{make_file_str}/{term}.zip"
                            break
            else:
                print('----------gene下载完成')
        else:
            print(f"Search type出现问题：{search_type}")
            pass
        #解压缩包到指定的"make file str"路径内        
        with zipfile.ZipFile(file_path, 'r') as zip_ref:
            zip_ref.extractall(make_file_str)
        #通过解读report，定位最佳基因组
        best_gene_info = get_best_genome(f"{make_file_str}/ncbi_dataset/data/assembly_data_report.jsonl")
        #获取最佳基因组
        genomic_source_path = glob.glob(f"{make_file_str}/ncbi_dataset/data/{best_gene_info[0]}/*")[0]        
        #基因组要拷贝的位置
        genomic_copy_path = f"/reference_data/{temp_path}/{term.replace(' ', '_')}.fna"

        #完成基因组文件拷贝任务，用于后续
        shutil.copy(genomic_source_path, genomic_copy_path)
        #如果用户只要物种层面内容，此时即可返回全基因组文件
        #TODO：需要进一步细化任务类型
        if search_type in ["species_genome", "dna_species_gene", "rna_species_gene"]:
            return [genomic_copy_path]
        #非物种层面，还需要挖掘cds gene信息
        make_cds_file_str = f"/reference_data/{temp_path}/{term}/cds"
        #确定cds文件是否下载，没通过上一步获取的accession number，下载cds文件
        if not os.path.exists(make_cds_file_str):
            os.system(f"mkdir {make_cds_file_str}")

        target_gene_list = []
        #文件转储
        if search_type=="species_cds":
            #cds文件路径
            cds_path = genomic_source_path #f"/reference_data/{temp_path}/{term}/cds/ncbi_dataset/data/{best_gene_info[0]}/cds_from_genomic.fna"
            
            #转储新的cds文件路径
            if not os.path.exists(f"/reference_data/{temp_path}/{term}/cds/ncbi_dataset/data/{best_gene_info[0]}/"):
                os.makedirs(f"/reference_data/{temp_path}/{term}/cds/ncbi_dataset/data/{best_gene_info[0]}/")
            new_cds_path = f"/reference_data/{temp_path}/{term}/cds/ncbi_dataset/data/{best_gene_info[0]}/{term}_cds_from_genomic.fna"
        
        else:
            if not os.path.exists(make_cds_file_str):
                os.system(f"mkdir {make_cds_file_str}")
                get_cds_url = f"datasets download genome accession {best_gene_info[0]} --include cds --filename {make_cds_file_str}/{term}_cds.zip"
                os.system(get_cds_url)
            #解压cds文件
            with zipfile.ZipFile(f"{make_cds_file_str}/{term}_cds.zip", 'r') as zip_ref:
                zip_ref.extractall(make_cds_file_str)
            
            #cds文件路径
            cds_path = f"/reference_data/{temp_path}/{term}/cds/ncbi_dataset/data/{best_gene_info[0]}/cds_from_genomic.fna"
            #转储新的cds文件路径
            new_cds_path = f"/reference_data/{temp_path}/{term}/cds/ncbi_dataset/data/{best_gene_info[0]}/{term}_cds_from_genomic.fna"
        #文件重命名，已名称作为区分
        print(f"路径转储: cp {cds_path} {new_cds_path}")
        os.system(f"cp {cds_path} {new_cds_path}")
        return [new_cds_path]
def extract_json_response(response):
    # 从response中提取JSON内容
    experiment_info_dict = response["choices"][0]["message"]['content']
    try:
        #if json.loads(experiment_info_dict):
        return json.loads(experiment_info_dict)
    except:
        match = re.search(r'```json\n(.*?)\n```', experiment_info_dict, re.DOTALL)
        if match:
            experiment_info_json_str = match.group(1)
            experiment_info_json_data = json.loads(experiment_info_json_str)
        else:
            print("GPT构造JSON错误，情况如下所示：\n",experiment_info_dict)
            #exit()
        return experiment_info_json_data

def get_gca_gene(data,data_info_list,genome_dict):
    for data_info in data:
        if "contigN50" not in data_info["assemblyStats"]:
            continue
        if "assemblyLevel" not in data_info["assemblyInfo"]:
            continue

        accession = data_info["accession"]
        assemblyLevel = data_info["assemblyInfo"]["assemblyLevel"]
        contigN50 = data_info["assemblyStats"]["contigN50"]
        data_info_list.append((accession, assemblyLevel, contigN50))

        if assemblyLevel == "Complete Genome":
            genome_dict['complete'][accession] = (assemblyLevel, scaffoldN50)
            complete_flag = True
        else:
            genome_dict['noncomplete'][accession] = (assemblyLevel, scaffoldN50)
    return data_info_list,complete_flag,genome_dict
    

def get_best_genome(path):
    #GCF考虑scaffoldN50，assemblyLevel，前者最大，后者需要Complete Genome
    #无GCF下，考虑GCA，contigN50，assemblyLevel，有无complete都可，需要最大contigN50
    with open(path, 'r') as f:
        lines = f.readlines()
        data = [json.loads(line) for line in lines]

    data_info_list = []
    genome_dict = {'complete':{}, 'noncomplete':{}}
    complete_flag = False
    print("start to find GCF")
    for data_info in data:
        if "scaffoldN50" not in data_info["assemblyStats"]:
            continue
        if "assemblyLevel" not in data_info["assemblyInfo"]:
            continue
        accession = data_info["accession"]
        assemblyLevel = data_info["assemblyInfo"]["assemblyLevel"]
        scaffoldN50 = data_info["assemblyStats"]["scaffoldN50"]
        data_info_list.append((accession, assemblyLevel, scaffoldN50))
        if assemblyLevel == "Complete Genome":
            genome_dict['complete'][accession] = (assemblyLevel, scaffoldN50)
            complete_flag = True
        else:
            genome_dict['noncomplete'][accession] = (assemblyLevel, scaffoldN50)
    
    if data_info_list == []:
        print("NO GCF, refind GCA")
        data_info_list,complete_flag,genome_dict = get_gca_gene(data,data_info_list,genome_dict)

    if complete_flag:
        best_complete_genome = sorted(genome_dict['complete'].items(), key=lambda x: x[1][1], reverse=True)[0] if genome_dict['complete'] else None
        return best_complete_genome

    else:
        best_noncomplete_genome = sorted(genome_dict['noncomplete'].items(), key=lambda x: x[1][1], reverse=True)[0] if genome_dict['noncomplete'] else None
        return best_noncomplete_genome

def extract_gene(fna_content, gene_name):
    pattern = r">(.*?)\n([ATGC\n]+)"
    matches = re.findall(pattern, fna_content, re.MULTILINE)
    
    for header, sequence in matches:
        if f"[gene={gene_name}]" in header:
            seq = sequence.replace('\n','')
            return f">{header}\n{seq}"
    return ""

def get_embedding_recommended(text, type, species_class=""):
    #species_identification_dict["Species_Name"], type = "determine_isolate"
    # 定义请求数据
    data = {
            "text": text,
            "mode": type,
            "species_class": species_class,
            "top_k": 50,
            "top_r": 3,
            "score_threshold": 0.2
            }
    
    print("Sending data:", json.dumps(data, indent=4))

    # 发送POST请求
    response = requests.post(
                url='http://127.0.0.1:8119/search/',
                headers={'Content-Type': 'application/json', 'accept': 'application/json'},
                json=data
                )
    # 检查响应状态码
    if response.status_code == 200:
        # 打印响应内容
        if type == "determine_species":
            print(response.json())
            print('+++++++++++++++')
            search_results = response.json()["species_results"]
            print("recommended search results sampling: ",search_results)
            return search_results

        elif type == "determine_cds":
            search_results = response.json()["cds_results"]
            #print("cds yes or no:", search_results)
            return search_results

        elif type == "get_isolate2id":
            search_results = response.json()["isolate2id_results"]
            return search_results

    else:
        print(f"Request failed with status code {response.status_code}")

def getCMD(kmer_file, klen, db, ofile, qcov_hsp_perc=50, ident=15):
    if klen <= 50:
        cmd = 'time blastn -task blastn-short'
    elif klen <= 500:
        cmd = 'time blastn -task blastn'
    else:
        cmd = 'time blastn -task megablast'
    cmd += f' -query {kmer_file} -db {db} -num_threads 100 -perc_identity {ident} -qcov_hsp_perc {qcov_hsp_perc} -outfmt "6 qseqid sseqid qstart qend sstart send length pident mismatch gapopen evalue bitscore staxids sscinames qlen"  > {ofile}'
    print(cmd)
    return cmd

def genKmer(seq,  kmer_len):
    kmer = set()
    for i in range(0, len(seq)):
        if i + kmer_len > len(seq):     ## 最后一个窗口可能达不到kmer长度，跳过
            break
        kmer.add(seq[i:(i + kmer_len)])
    return kmer

def getkmer(infile, outfile):
    all_kmers = set()
    try:
        seq_list = []
        seq_len_list = []
        with open(outfile, 'w') as out_fasta:
            average_len = 150
            window = 0
            for r in SeqIO.parse(infile, 'fasta'):     ## 读取基因组文件，因为有的基因组组装不完整，会有很多条序列，所以需要迭代
                seq = str(r.seq)
                seq_list.append((seq, r.id))

            # 根据seq长度进行排序，选出前三个长度最长的元素
            top_3_longest = sorted(seq_list, key=lambda x: len(x[0]),reverse=True)[:3]

            for seq in top_3_longest:
                print(len(seq[0]))
                kmers = genKmer(seq[0], average_len)
                name = seq[1]
                for idx, kmer in enumerate(kmers):
                    if 'lcl' in name:
                        record = SeqRecord.SeqRecord(Seq(kmer), id=f"{name}_cds_kmer_{idx}", description="")
                    else:
                        record = SeqRecord.SeqRecord(Seq(kmer), id=f"lcl|{name}_cds_kmer_{idx}", description="")

                    SeqIO.write(record, out_fasta, 'fasta')
    except:
        traceback.print_exc()
        sys.exit(1)

#taxonkit来获取种菌分支taxid
def get_taxon_subtree(taxid):
    # 调用 taxonkit list 命令
    result = subprocess.run(["taxonkit", "list", "--ids", str(taxid)],
            capture_output=True,
            text=True,
            check=True)
    # 处理输出
    output = result.stdout.split('\n')
    subtree = [int(i.replace(' ','')) for i in output if i!='']
    return subtree

#物种和菌株不同
def filter_cds_blast_output(file_path, subtree, identity_threshold=90, type_='general_speices_blast'):
    #解析Blast结果
    column_names = ["qseqid", "sseqid", "qstart", "qend", "sstart", "send", "length", "pident", "mismatch", "gapopen", "evalue", "bitscore", "staxids", "sscinames", "qlen"]
    blast_output = pd.read_csv(file_path, sep="\t", header=None, names=column_names)
    #指定identity过滤值
    identity_filter_blast_output = blast_output[blast_output['pident'] >= identity_threshold ]
    #任务不同后续操作不同
    blast_dict = {}
    #做区分
    if type_ == "general_speices_blast":
        filter_blast_list = []
        #过滤同物种
        filtered_subtree_blast = blast_output[~blast_output["staxids"].isin(subtree)]
        #获取唯一query name
        filter_blast_list = filtered_subtree_blast["qseqid"].unique().tolist()
        print('ok')

    elif type_ == "Group_species_blast":
        filtered_subtree_blast = blast_output[~blast_output["staxids"].isin(subtree)]
        for idx, line_ in filtered_subtree_blast.iterrows():
            pass
    else:
        pass

    return filter_blast_list

def get_cds_blast_result(cds_name_list, inputfile, outputfile):
    with open(inputfile, 'r') as input_handle, open(outputfile, 'w') as output_handle:
        for record in SeqIO.parse(input_handle, 'fasta'):
            name = record.id
            if name not in cds_name_list:
                SeqIO.write(record, output_handle, 'fasta')

def get_best_core_gene(blastoutput, core_cds_path, output_path):
    column_names = ["qseqid", "sseqid", "qstart", "qend", "sstart", "send", "length", "pident", "mismatch", "gapopen", "evalue", "bitscore", "staxids", "sscinames", "qlen"]
    blast_output = pd.read_csv(blastoutput, sep="\t", header=None, names=column_names)
    name_list = []
    info_list = []
    print(core_cds_path)
    for record in SeqIO.parse(core_cds_path, 'fasta'):
        name = record.id.strip()
        print(name)
        name_list.append(name)
        info_list.append(record)
    print('----------blast:\n',blast_output)
    filter_blast_output = blast_output[blast_output['qseqid'].isin(name_list)]
    print('----------filter:\n',filter_blast_output)
    sorted_filter_blast_output = filter_blast_output.sort_values(by='pident')
    print('----------sorted:\n',sorted_filter_blast_output)
    result = sorted_filter_blast_output.head(1)
    index = name_list.index(result['qseqid'].values[0])
    with open(output_path, 'w') as output_handle:
        SeqIO.write(info_list[index], output_handle, 'fasta')

def get_core_gene(file_path, taxid=[],type_num="single", type_cds="",type_experiment=""):
    print('***************开始获取core gene***************')
    if type_cds == "have_cds":
        kmer_file = file_path[0]
    else:
        #当无cds或gtf文件时，要做kmer处理
        infile = file_path[0]
        kmer_file = ""
        getkmer(infile, outfile)
    #kmer长度
    start_time = time.time()
    len_list=[]
    for r in SeqIO.parse(kmer_file, 'fasta'):
        seqlen = len(str(r.seq))
        len_list.append(seqlen)
    klen = int(sum(len_list)/len(len_list))
    print('片段平均长度:', klen)
        
    #blast数据库
    db = '/home/YJhou/workspace/pcr_paper/search_base/blast_db/refseq'
    #blast输出
    path_lst = kmer_file.split('/')
    output_path = '/'.join(path_lst[:-1])
    name = path_lst[-1]
    blast_name = name.split('.fna')[0]+'_blast_output'
    blast_ofile = output_path+'/'+blast_name
    print('kmer_file：',kmer_file)
    cmd = getCMD(kmer_file, int(klen), db, blast_ofile, 80, 0)
    print('blast指令：', cmd)
    os.system(cmd)
    print('blast耗时：', time.time()-start_time)
        
    start_time = time.time()
    subtree = get_taxon_subtree(taxid[0]) 
    print("获取同物种分支耗时：",time.time()-start_time)
    start_time = time.time()
    filter_blast_list = filter_cds_blast_output(blast_ofile, subtree, 15)
    print("按照物种分支和identity过滤耗时：",time.time()-start_time)
        
    start_time = time.time()
    #过滤cds
    outputfile = output_path+'/'+blast_name+'_core_cds.fna'
    get_cds_blast_result(filter_blast_list, kmer_file, outputfile)
        
    #再次比对，找出一致性最小的
    blast_identity_ofile = output_path+'/'+blast_name+'_compare_identity'
    cmd = getCMD(outputfile, int(klen), db, blast_identity_ofile, 50, 0)
    print('blast指令：', cmd) 
    os.system(cmd)
    best_core_gene_path = "/reference_data/target/"+name.split('.fna')[0] + '_core_cds.fna'
    get_best_core_gene(blast_identity_ofile, outputfile, best_core_gene_path)
    return best_core_gene_path


