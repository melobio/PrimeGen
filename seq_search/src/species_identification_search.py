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

# species identification search tool
def species_identification_search(instruction, stage):
    # 第一轮对话，实体提取，看用户是否说明物种、菌株和Group或Microbe实验方向，提取信息在后续search中使用
    if stage == 1:
        species_identification_dict = {
            "Experiment_Direction": "species identification",
            "Experiment_Type": None,
            "Species_Name": "",
            "Strain": "",
            "Gene_Name": "",
            "Taget_Need": "",
        }

        conversation = instruction['conversation']

        # 根据上下文对话，提取实验类型、实验类型是否提供、物种名称、菌株名称
        experiment_type_prompt = """Based on the user's contextual dialogue, you need to perform information extraction. The information to be extracted includes:  
1.Experiment type: Identify whether it's "Microbe" or "Group" (return null if the user does not specify) 
2.Experiment type provided: Return a Boolean value; false only when the experiment type is not mentioned.
3.Species Name: Extract if mentioned by the user; return null if not mentioned. 
4.Strain Name: Extract if mentioned by the user; return null if not mentioned. 
The results should be returned to me in JSON format. Example is shown below:
'''json
{  
  "Experiment_type": null,  
  "Experiment type provided": false,  
  "Species_Name": "E. coli",  
  "Strain_Name": null  
}
'''
Note:
1)The “group” as defined here refers to several microorganisms that may be from the same genus or from different genera, while the “Microbe” refers to a well-defined species, such as a strain or substrain.
2)The species name and strain name should be complete. If the user provides an abbreviation, please complete it.
"""
        conversation_str = "\n\n".join([i["role"] + ':\n' + i['content'] for i in conversation])
        experiment_info_response = openai.ChatCompletion.create(
            messages=[{"role": "system", "content": experiment_type_prompt},
                      {"role": "user", "content": conversation_str}
                      ],
            deployment_id="XMGI-Chat1-GPT4"
        )
        experiment_info_json_data = extract_json_response(experiment_info_response)
        print("-------------extract experiment information:\n", experiment_info_json_data)

        species_identification_dict["Species_Name"] = experiment_info_json_data["Species_Name"]
        species_identification_dict["Strain"] = experiment_info_json_data["Strain_Name"]
        # species_identification_dict["Experiment Type"] = experiment_info_json_data["Experiment type"]

        if species_identification_dict["Strain"] is None:
            recommended_system_prompt = f"""Based on conversation, the user has provided species name {species_identification_dict["Species_Name"]}, but detailed strain-level information is required for our experiment. Please provide recommendations for specific strains of the mentioned species. The results should be formatted in JSON. For example, if the user mentions 'Staphylococcus aureus', you should return a JSON file with recommended strains like this:
'''json
{{
"recommended_strains": [
"Staphylococcus aureus subsp. aureus Rosenbach",
"Staphylococcus aureus subsp. aureus MRSA252",
"Staphylococcus aureus subsp. aureus N315",
"Staphylococcus aureus subsp. aureus Mu50",
"Staphylococcus aureus subsp. aureus USA300",
"Staphylococcus aureus subsp. aureus MSSA476",
"Staphylococcus aureus subsp. aureus Newman",
"Staphylococcus aureus subsp. aureus COL",
]}}
'''

if the user mentions 'Escherichia coli', you should return a JSON file with recommended strains like this:
'''json
{{"recommended_strains": [
"Escherichia coli O157:H7 str. Sakai",
"Escherichia coli O157:H7 str. EDL933",
"Escherichia coli O157:H7 str. TW14359"]}}
'''Please ensure that the strains listed are relevant to the species mentioned by the user.At the same time, the names of examples should be accurate and concise."""

            conversation_str = "\n\n".join([i["role"] + ':\n' + i['content'] for i in conversation])
            recommended_strains_response = openai.ChatCompletion.create(
                messages=[{"role": 'system', 'content': recommended_system_prompt},
                          {"role": "user", 'content': conversation_str}],
                deployment_id="XMGI-Chat1-GPT4"
            )
            recommended_strains_json_data = extract_json_response(recommended_strains_response)
            print("-------------recommended_strains:\n", recommended_strains_json_data)

        # 用户已明确实验方向
        if species_identification_dict["Experiment_Type"] is not None:
            if 'microbe' in species_identification_dict["Experiment_Type"].lower():
                stage += 1
                response = f"You've chosen to conduct an experiment focusing on {species_identification_dict['Experiment_Type']}, specifically targeting the {species_identification_dict['Species_Name']} species."
                if species_identification_dict["Strain"] is None:
                    response += "You are currently preparing to conduct a microbe experiment. I will recommend several candidate microbial strains for your selection. Alternatively, you may upload the .fna file of the strain you intend to investigate. [Note: The strain file should be in .fna format]. Upon confirmation of the selected strain, we will suggest potential target genes that may be required for your study. You have the option to utilize these target genes to enhance the design of your primers, or you may upload your own target gene file. [Note: The target gene file should also be in .fna format]. Alternatively, if you opt not to use specific target genes, our AI algorithm will automatically identify relevant target genes and design the necessary primers for your experiment."

                    return {"response": response, "operations":[{"title":"recommend strains","options":recommended_strains_json_data["recommended_strains"],"type":["single","input"],"key":"strain_select","value":[]},
                        {"title":"experiment type", "options":["group", "microbe"], "type":["single","input"],"key":"experiment_select", "value":[]},
                        {"title":"upload strains", "options":[], "type":["file"], "key":"strain_path", "value":[]}],"stage": stage,
                            "state": 'continue', "search_type": "species_identification_type",
                            "species_identification_dict": species_identification_dict}
                else:
                    response += f"You intend to study the strain {species_identification_dict['Strain']}. "
                    response += """we will suggest potential target genes that may be required for your study. You have the option to utilize these target genes to enhance the design of your primers, or you may upload your own target gene file. [Note: The target gene file should also be in .fna format]. Alternatively, if you opt not to use specific target genes, our AI algorithm will automatically identify relevant target genes and design the necessary primers for your experiment."""

                    return {"response": response, "operations":[{"title":"experiment type", "options":["group", "microbe"], "type":["single","input"],"key":"experiment_select", "value":[]}],
                            "stage": stage,"state": 'continue', "search_type": "species_identification_type",
                            "species_identification_dict": species_identification_dict}


            elif 'group' in species_identification_dict["Experiment_Type"].lower():
                stage += 1
                response = f"You've chosen to conduct an experiment focusing on {species_identification_dict['Experiment_Type']}, specifically targeting the {species_identification_dict['Species_Name']} species."

                if species_identification_dict["Strain"] is None:
                    response += """You are currently preparing to conduct a group experiment. To facilitate subsequent primer design, please provide information according to the options and fields below. You can:
1) I will recommend several candidate strains; please select multiple strains from these recommendations.
2) You may also directly upload the strain file you intend to work with.
3) Alternatively, you can directly inform me of the strains you need in the provided section; please ensure the names are correct, and I will retrieve and download the files for you.
Once you confirm the strains for your study, our AI will automatically identify target genes within the group and design primers for you.
"""
                    return {"response": response, "operations":[{"title":"recommend strains","options":recommended_strains_json_data["recommended_strains"], "type":["multi","input"],"key":"strain_select","value":[]},
                        {"title":"experiment type", "options":["group", "microbe"], "type":["single","input"],"key":"experiment_select", "value":[]},
                        {"title":"upload strains", "options":[], "type":["file"], "key":"strain_path", "value":[]}],
                            "stage": stage,"state": 'continue', "search_type": "species_identification_type",
                            "species_identification_dict": species_identification_dict}
                else:
                    response += f"You intend to study the strain {species_identification_dict['Strain']}. "
                    response += """our AI will automatically identify target genes within the group and design primers for you."""
                    #TODO:这里不需要靶基因推荐
                    return {"response": response, "operations":[{"title":"experiment type", "options":["group", "microbe"], "type":["single","input"],"key":"experiment_select", "value":[]}],
                            "stage": stage,"state": 'continue', "search_type": "species_identification_type",
                            "species_identification_dict": species_identification_dict}


        # 用户未明确实验方向
        else:
            # 如果用户没有说明物种名，就中断Search，state是stop，controller要进一步提问咨询。
            if species_identification_dict["Species_Name"] is None and species_identification_dict["Strain"] is None:
                user_response = "I was unable to find the species or strain information you provided. You need to provide this information so that I can proceed with gene sequence retrieval. At a minimum, you need to tell me the species name, and I will guide you through the subsequent steps."
                return {"response": user_response, 
                        "operations":[],
                        "stage": stage,
                        "state": 'stop', "search_type": "controller_type",
                        "species_identification_dict": species_identification_dict}
            else:
                stage += 1
                user_response = f"""You have indicated your interest in conducting experiments on the {species_identification_dict["Species_Name"]} species. """

                user_response += """Would you like to design primers for a "Group" or "Microbe" type? Your choice will affect the subsequent retrieval and design task workflows. 
1) "Group" here is defined as a collection of several microbes that may come from the same genus or different genera.
2) "Microbe" refers to a clearly defined species, such as a strain or sub-strain.

Additionally, I will recommend strains or sub-strains for your reference and selection. You can also directly enter the strain name in the section provided, and I will assist with the download. Alternatively, you can upload a gene file directly. [Note, the gene file should be in .fna format].

Further, once you confirm the experiment type (Group/Microbe) and the species/subspecies under study, we will recommend potential target genes you might need (no recommendations are provided for Group; in this case, AI must either design autonomously or the user must upload specific genes for the Group).
"""

                return {"response": user_response,
                        "operations":[{"title":"experiment type", "options":["group", "microbe"], "type":["single","input"],"key":"experiment_select", "value":[]},
                        {"title":"recommend strains","options":recommended_strains_json_data["recommended_strains"],"type":["multi","input"],"key":"strain_select","value":[]},
                        {"title":"upload strains", "options":[], "type":["file"], "key":"strain_path", "value":[]}
                            ],
                        "stage": stage,"state": 'continue', "search_type": "species_identification_type",
                        "species_identification_dict": species_identification_dict}


    # 判断用户是否需要进行靶基因确认
    elif stage == 2:
        species_identification_dict = instruction['instruction']['species_identification_dict']
        
        upload_file_flag = instruction['instruction']['upload_file_flag']

        # 用户上传则为True，提取菌株名称用于后续target gene推荐
        if upload_file_flag == True:
            strain_list = []
            strain_path_list = instruction['instruction']['strain_path']
            for strain_path in strain_path_list:
                strain_list.append(strain_path.split('/')[0].split('_cds')[0].replace("_", " "))
            species_identification_dict["Strain"] = strain_list

        elif species_identification_dict["Strain"] is None:
            species_identification_dict["Strain"] = instruction['instruction']['strain_select']

        # 用户未提供实验方向和starin时，上轮的选择
        if species_identification_dict["Experiment_Type"] is None:
            species_identification_dict["Experiment_Type"] = instruction['instruction']['experiment_select'][0]

        conversation = instruction['instruction']['conversation']
        # if species_identification_dict["Strain"] !="null" and species_identification_dict["Species_Name"] !="null":

        strain_ = str(species_identification_dict['Strain'])
        species_name = str(species_identification_dict["Species_Name"])

        target_gene_recommended_system_prompt = f"""Based on the conversational context, it is known that User wants to study {strain_}. Please provide the target gene of the customer's research content and return it to me in JSON format. The JSON format is as follows:
'''json
{{
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
}}'''
Note, the recommended target genes should be related to the strain being studied, as some target genes belong to other strains within the same species.
"""
        conversation_str = "\n\n".join([i["role"] + ':\n' + i['content'] for i in conversation])
        print("target_gene_recommended:",conversation_str)
        target_gene_recommended_response = openai.ChatCompletion.create(
            messages=[{"role": 'system', 'content': target_gene_recommended_system_prompt},
                      {"role": "user", 'content': conversation_str}],
            deployment_id="XMGI-Chat1-GPT4"
        )
        print("target_gene_recommended_response:",target_gene_recommended_response)
        target_gene_recommended_json_data = extract_json_response(target_gene_recommended_response)

        
        if upload_file_flag == True:
            strain_file_path = instruction['instruction']['strain_path']
            strain_response = ""
            stage += 1
        
        # 无论如何都要下载参考基因组
        else:
            strain_info = species_identification_dict["Strain"]
            strain_info_list = []
            for strain_ in strain_info:
                print(f'*********************{strain_} searching***************************')
                strain_search_info, strain_file_path = ncbi_search(strain_, 'genome', target_type="target")  # 直接基于菌株下载参考序列 #TODO：设计ncbi search genome
                print(strain_search_info)
                strain_info_list.append(strain_file_path)
            stage += 1
            strain_info_list = [x[0] for x in strain_info_list]
            strain_response = "The current sequence search requirements have been completed. The following is the strain reference genome file.  "

        if species_identification_dict["Experiment_Type"].lower() == 'group':
            user_response = strain_response + """For the current Group experiment, multiple common genes shared among the strains can be provided as target genes to enhance the success rate of primer design. You may directly upload potential common genes as target genes. If you do not provide target genes, I will design the primers myself."""
            return {"response": user_response, "operations":[{"title":"upload target genes", "options":[], "type":["file"], "key":"target_gene_path", "value":[]}],
                "strain_path": strain_info_list, "stage": stage,
                "state": 'continue', "search_type": "species_identification_type",
                "species_identification_dict": species_identification_dict}
        
        elif species_identification_dict["Experiment_Type"].lower() == 'microbe':
            user_response = strain_response + """For current Microbe experiments, you can provide the target gene of the strain, which will increase the success rate of primer design. You can refer to the target gene content I recommend, or you can specify a specific target gene and provide its relevant information(gene name, like "lacZ"), and I will search it on NCBI. If you do not provide the target gene, I will design the primers myself. The target gene names I recommend are as follows:"""
            return {"response": user_response, "operations":[{"title":"recommend target genes","options":target_gene_recommended_json_data["recommended_target_genes"],"type":["multi","input"],"key":"target_gene_select","value":[]},
                {"title":"upload target genes", "options":[], "type":["file"], "key":"target_gene_path", "value":[]}],
                "strain_path": strain_file_path, "stage": stage,
                "state": 'continue', "search_type": "species_identification_type",
                "species_identification_dict": species_identification_dict}

    # 基于用户靶基因的抉择，判断后续任务：
    elif stage == 3:
        upload_file_flag = instruction["instruction"]['upload_file_flag']
        species_identification_dict = instruction["instruction"]['species_identification_dict']

        species_identification_dict["Gene_Name"] = instruction["instruction"]["target_gene_select"]#'select_target_gene']
        conversation = instruction["instruction"]['conversation']

        # 用户上传则为True，提取菌株名称用于后续target gene推荐
        if upload_file_flag == True:
            gene_list = []
            target_gene_path_list = instruction["instruction"]['target_gene_path']
            for gene_path in target_gene_path_list:
                gene_list.append(gene_path.split('_gene')[0].replace("_", " "))
            species_identification_dict["Gene_Name"] = gene_list

        elif species_identification_dict["Gene_Name"] == 'null':
            species_identification_dict["Gene_Name"] = instruction["instruction"]['target_gene_select']

        flag_system_prompt = """From the current context dialogue, determine whether the user needs to further confirm the target gene information to design primers. You just need to reply yes or no."""
        conversation_str = "\n\n".join([i["role"] + ':\n' + i['content'] for i in conversation])
        flag_response = openai.ChatCompletion.create(
            messages=[{"role": 'system', 'content': flag_system_prompt},
                      {"role": "user", 'content': conversation_str}],
            deployment_id="XMGI-Chat1-GPT4"
        )
        flag_ = flag_response["choices"][0]["message"]['content']

        # non-target推荐与选择
        non_target_gene_recommended_system_prompt = f"""Based on our previous conversation, you are conducting an {species_identification_dict["Experiment_Type"]} experiment on the {species_identification_dict["Species_Name"]} species. 
To ensure the efficiency of the primer design, I will recommend some non-target species. This will help you design primers with higher specificity.
Below, I will provide recommendations in JSON format. Only the JSON recommendations are required, no additional description is needed.
The JSON format is as follows:  
'''json{{    
  "recommendations non-target": [  
      "Salmonella enterica",    
      "Shigella flexneri",    
      "Enterobacter cloacae",    
      "Klebsiella pneumoniae",    
      "Proteus mirabilis"    
    ]    
}}'''
Note, recommendations for non-target species should be tailored based on the user's experiment type and target species, and should not be generic.
"""

        non_target_gene_recommended_response = openai.ChatCompletion.create(
            messages=[{"role": 'system', 'content': non_target_gene_recommended_system_prompt},
                      {"role": "user", 'content': conversation_str}],
            deployment_id="XMGI-Chat1-GPT4"
        )
        non_target_gene_recommended_json_data = extract_json_response(non_target_gene_recommended_response)
        non_target_gene_response = """To enhance the specificity of your primer design, you may select from the non-target information I recommend (note that these recommendations are broad and may not necessarily be feasible in your experimental environment). We suggest that you directly upload non-target files (in .fna format) based on your experimental context. If you choose not to use the non-target information, we will design primers based solely on the target sequence. Please inform me of your decision or requirements. My suggestions for non-target information are as follows:"""
        if upload_file_flag == True:
            stage += 1
            user_response = """The target gene file uploaded by the user has been identified. Next, I will design primers based on the reference genome file of the strain."""
            user_response += non_target_gene_response

            return {"response": user_response, 
                    "stage": stage,
                    "operations":[
                        {
                            "title":"recommend non-target gene",
                            "options":non_target_gene_recommended_json_data["recommendations non-target"],
                            "type":["multi","input"],
                            "key":"non_target_gene_select",
                            "value":[]
                        },
                        {
                         "title":"upload non-target gene", 
                         "options":[], 
                         "type":["file"], 
                         "key":"non_target_gene_path", 
                         "value":[]
                         }
                        ],
                    "state": 'continue', 
                    "search_type": "species_identification_type",
                    "species_identification_dict": species_identification_dict
                    }

        # 用户不需要靶基因,继续问用户non-target事情
        if 'llm' in flag_.lower():
            user_response = f"""Through our communication, it has been established that you have no ideas on target genes for primer design. Next, we will need you to decide whether you require non-target genome information to optimize the specificity of the primers."""
            user_response += non_target_gene_response
            return {"response": user_response, "operations":[{"title":"recommend non-target gene","options":non_target_gene_recommended_json_data["recommendations non-target"],"type":["multi","input"],"key":"non_target_gene_select","value":[]},
                    {"title":"upload non-target gene", "options":[], "type":["file"], "key":"non_target_gene_path", "value":[]}],
                    "stage": stage, "state": 'continue',
                    "search_type": "species_identification_type",
                    "species_identification_dict": species_identification_dict}

        # 用户需要靶基因，进一步给出检索结果，让用户可看到，继续问用户non-target事情
        else:
            Genus_name = species_identification_dict["Gene_Name"]
            target_gene_search_info_list = []
            target_gene_search_path_list = []
            # for循环包含单个的microbe情况和多个的Group情况
            for sel in Genus_name:
                print([species_identification_dict["Species_Name"]+' '+species_identification_dict["Strain"][0], sel.split("(")[0]])
                if Levenshtein.distance(species_identification_dict["Species_Name"], species_identification_dict["Strain"][0])<len(species_identification_dict["Species_Name"]):
                    search_info, search_path = ncbi_search([species_identification_dict["Species_Name"]+' '+species_identification_dict["Strain"][0], sel.split("(")[0]], 'gene')  # 直接基于菌株下载参考序列 #TODO：设计ncbi search genome #这要做好target的过滤操作，根据description过滤 "Mycobacterium tuberculosis H37Rv" AND katG [Gene Name]  #这里要注意Gene Name的正确性，要让GPT4check，检索不到要反馈用户错误。
                else:
                     search_info, search_path = ncbi_search([species_identification_dict["Strain"][0], sel.split("(")[0]], 'gene')
                target_gene_search_info_list.append(search_info)
                target_gene_search_path_list.append(search_path)
            target_gene_search_path_list = [x[0] for x in target_gene_search_path_list]
            stage += 1
            user_response = f"""Based on our conversation, we've looked up the target gene information for {Genus_name} and found the following sequence results for you:{target_gene_search_info_list}"""

            user_response += non_target_gene_response

            return {"response": user_response, "target_gene_path": target_gene_search_path_list,
                    "operations":[{"title":"recommend non-target gene","options":non_target_gene_recommended_json_data["recommendations non-target"],"type":["multi","input"],"key":"non_target_gene_select","value":[]},
                    {"title":"upload non-target gene", "options":[], "type":["file"], "key":"non_target_gene_path", "value":[]}],
                    "stage": stage,"state": 'continue', "search_type": "species_identification_type",
                    "species_identification_dict": species_identification_dict}

    # 判断non-target是否需要，并给出推荐
    elif stage == 4:
        conversation = instruction["instruction"]['conversation']
        upload_file_flag = instruction["instruction"]['upload_file_flag']  # 用户上传则为True，我不需要走search流程
        species_identification_dict = instruction["instruction"]['species_identification_dict']
        user_select = instruction["instruction"]['non_target_gene_select']

        # 用户没有选择
        if user_select == []:
            non_target_extract_system_prompt = """ Based on our conversation and the context, I will extract the non-target information User interested in. If User haven't specified any non-targets you want to focus on, I will return 'no'. If User have clearly mentioned User desired non-targets, I will provide them to you in JSON format.   
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
                user_response = """You have confirmed that you do not need to use non-target information to assist in designing primers. Therefore, we will proceed with the primer design based solely on the target information. Next, I will notify the Controller to communicate with you to determine the subsequent steps related to primer design."""

                return {"response": user_response, "operations":[], "stage": stage, "state": 'stop',
                        "search_type": "species_identification_type",
                        "species_identification_dict": species_identification_dict}

            else:
                non_target_search_info_list = []
                for non_tar_ in non_target:
                    search_info = ncbi_search(non_tar_, "genome", target_type="non_target")  # 检索物种参考基因组
                    non_target_search_info_list.append(search_info)

                user_response = """Based on the non-target information you provided, I have downloaded the following reference genomes."""
                user_response += """Thus far, through our communication, we have confirmed the information regarding both target and non-target. Next, I will inform the Controller to communicate with you to determine the next steps related to primer design."""

                return {"response": user_response, "operations":[], "non_target_results": non_target_search_info_list,
                        "choice_type": 'Multi', "stage": stage,
                        "state": 'top', "search_type": "species_identification_type",
                        "species_identification_dict": species_identification_dict}

        # 用户在选择栏目中输入no
        elif user_select[0].lower() == 'no':
            stage += 1
            user_response = """You haven't provided any non-target information, so we will proceed with primer design based solely on the target information."""
            user_response += """Thus far, through our communication, we have confirmed the information regarding target. Next, I will inform the Controller to communicate with you to determine the next steps related to primer design."""
            return {"response": user_response, "operations":[], "stage": stage, "state": 'stop',
                    "search_type": "species_identification_type",
                    "species_identification_dict": species_identification_dict}
        
        #用户上传non-target
        elif instruction["instruction"]['non_target_gene_path']!=[]:
            user_response = """Thus far, through our communication, we have confirmed the information regarding both target and non-target. Next, I will inform the Controller to communicate with you to determine the next steps related to primer design."""
            return {"response": user_response, "non_target_path": instruction["instruction"]['non_target_gene_path'], "operations":[],
                    "stage": stage, "state": 'stop', "search_type": "species_identification_type",
                    "species_identification_dict": species_identification_dict}

        # 用户选择non-target
        else:
            stage += 1
            non_target_search_info_list = []
            non_target_path_list = []
            for non_tar_ in user_select:
                search_info, non_target_path = ncbi_search(non_tar_, "genome",target_type="non_target")  # 检索物种参考基因组
                non_target_search_info_list.append(search_info)
                non_target_path_list.append(non_target_path)
            user_response = f"""Based on the non-target information you selected, I have downloaded the following reference genomes:{non_target_search_info_list} so we will proceed with primer design based solely on the target information."""
            user_response += """Thus far, through our communication, we have confirmed the information regarding both target and non-target. Next, I will inform the Controller to communicate with you to determine the next steps related to primer design."""
            non_target_path_list = [ x[0] for x in non_target_path_list]
            return {"response": user_response, "non_target_path": non_target_path_list, "operations":[],
                    "stage": stage, "state": 'stop', "search_type": "species_identification_type",
                    "species_identification_dict": species_identification_dict}


def ncbi_search(term, db_name, target_type=None):
    # 参考基因组检索
    path = "/reference_data"
    if db_name == "genome":
        if target_type=="target":
            target_path = "target"
            temp_path = "target_zip_file"
        elif target_type=="non_target":
            target_path = "nontarget"
            temp_path = "non_target_zip_file"
        #print(term+"all name")
        handle = Entrez.esearch(db="taxonomy", term=term+" [all names]")
        record = Entrez.read(handle)
        print(record)
        ids = record["IdList"]
        id = ids[0]
        url = f"datasets download genome taxon {id} --reference --assembly-source 'RefSeq' --assembly-level complete --filename /reference_data/{temp_path}/{term.replace(' ', '_')}.zip"
        print(url)
        try:
            os.system(url)
        except:
            print('*****************download error, download again*****************')
            os.system(f"rm /reference_data/{temp_path}/{term.replace(' ', '_')}.zip")
            irl = f"datasets download genome taxon {id} --assembly-source 'RefSeq' --assembly-level complete --filename /reference_data/{temp_path}/{term.replace(' ', '_')}.zip"
            os.system(url)

        file_path = f"/reference_data/{temp_path}/{term.replace(' ', '_')}.zip"
        with zipfile.ZipFile(file_path, 'r') as zip_ref:
            zip_ref.extractall(f"/reference_data/{temp_path}")
        
        #TODO:当存在多个GCF是要先根据level选择complete，然后选择最先的gene
        target_file_path = glob.glob(f"/reference_data/{temp_path}/ncbi_dataset/data/GCF*")[0]
        target_file = glob.glob(target_file_path + '/*')[0]

        file_to_copy = f"/reference_data/{temp_path}/{term.replace(' ', '_')}.fna"
        print(f"Source file: {target_file}")
        print(f"File exists: {os.path.exists(target_file)}")
        print(f"Destination file: {file_to_copy}")
        shutil.copy(target_file, file_to_copy)
        print(f"finish{target_file} --> {file_to_copy}")
        return "Target gene file has been downloaded.", [file_to_copy]  # 返回下载好的数据位置

    # 靶基因检索
    elif db_name == "gene":
        print('+++++++++++++++++++++++',term)
        #return "ok", "test_path"  # 返回下载好的数据位置
        term_str = f"""({term[0]}) AND {term[1]} [Gene Name]"""
        temp_path = "target"
        # 先下载物种参考基因组
        print("term_str",term_str)
        print(term[0])
        print(term[1])
        handle = Entrez.esearch(db="taxonomy", term=term[0])
        record = Entrez.read(handle)
        ids = record["IdList"]
        try:
            id = ids[0]
        except:
            return  "The current search fails. Users are asked to search for relevant target genes on NCBI and other websites according to their needs, and then upload them.", []  # 返回下载好的数据位置
        url = f"datasets download genome taxon {id} --reference --assembly-source 'RefSeq' --assembly-level complete --filename /reference_data/{temp_path}/{term[0].replace(' ', '_')}.zip"
        print(url)
        try:
            os.system(url)
        except:
            print('*****************download error, download again*****************')
            os.system(f"rm /reference_data/{temp_path}/{term[0].replace(' ', '_')}.zip")
            url = f"datasets download genome taxon {id} --assembly-source 'RefSeq' --assembly-level complete --filename /reference_data/{temp_path}/{term[0].replace(' ', '_')}.zip"
            os.system(url)        
 

        file_path = f"/reference_data/{temp_path}/{term[0].replace(' ', '_')}.zip"
        with zipfile.ZipFile(file_path, 'r') as zip_ref:
            zip_ref.extractall(f"/reference_data/{temp_path}")

        target_file_path = glob.glob(f"/reference_data/{temp_path}/ncbi_dataset/data/GCF*")[0]
        target_file = glob.glob(target_file_path+'/*')[0]

        file_to_copy = f"/reference_data/{temp_path}/{term[0].replace(' ', '_')}_{term[1].replace(' ','_')}.fna"
        print(f"Source file: {target_file}")
        print(f"File exists: {os.path.exists(target_file)}")
        print(f"Destination file: {file_to_copy}")
        shutil.copy(target_file, file_to_copy)
        print(f"finish{target_file} --> {file_to_copy}")
        

        #获取靶基因的start和end
        print("target query:", term_str)
        handle = Entrez.esearch(db=db_name, term=term_str, retmax=10)
        record = Entrez.read(handle)
        print("target response:",record)
        ids = record["IdList"]
        print("record id:",len(ids))
        for taxonomy_id in ids:
            esummary_handle = Entrez.esummary(db="gene", id=taxonomy_id)
            esummary_record = Entrez.read(esummary_handle)
            esummary_handle.close()
            # print(esummary_record)
            Name = esummary_record['DocumentSummarySet']['DocumentSummary'][0]['Name']
        
            if term[1].lower().replace(' ','') not in Name.lower().replace(' ',''):  # 确保用户需要的靶基因在其中
                
                continue
            #print(esummary_record)
            Description = esummary_record['DocumentSummarySet']['DocumentSummary'][0]['Description']
            try:
                Assembly_accessions = esummary_record['DocumentSummarySet']['DocumentSummary'][0]['GenomicInfo'][0]['ChrAccVer']
            except:
                Assembly_accessions = "None"

            End_ = esummary_record['DocumentSummarySet']['DocumentSummary'][0]['GenomicInfo'][0]['ChrStart']
            Start_ = esummary_record['DocumentSummarySet']['DocumentSummary'][0]['GenomicInfo'][0]['ChrStop']
            str_=f"Name:{Name}\nDescription{Description}\nAssembly_accessions:{Assembly_accessions}\nStart_:{Start_}\nEnd_:{End_}"
            print(str_)

        #根据参考基因组和靶基因信息，获取目的基因片段
        with open(file_to_copy,'r') as all_file:
            content = all_file.read()
            target_gene = content[int(Start_):int(End_)]
        all_file.close() 
        target_gene_file_path = f"/reference_data/dna_file/{term[0]}_{term[1]}.fna"
        with open(target_gene_file_path, 'w') as f_w:
            target_gene_str_ = f">{Name} | {term[0]} | {Assembly_accessions}"
            target_gene_str_ += target_gene
            f_w.write(target_gene_str_)
        f_w.close()

        return str_, [target_gene_file_path]  # 返回下载好的数据位置
        

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


