import os
import json
import pandas as pd
import openai
import re


# cancer search tool
def cancer_search(instruction, stage):
    # description cancer search tool
    cancer_search_system_prompt = """The current Search Agent has entered the cancer retrieval phase and 
needs to interact with the user to obtain the genomic and corresponding drug resistance locus information for the cancer. 
The main task here is to interact with the user, guiding them to select the Gene name, Primary site, Site subtype 1, 
Primary histology, and Histology subtype 1. Based on the selected content and filtering through the Cosmic database, 
it acquires the gene name and corresponding SNP information, and also retrieves the genomic sequence information from NCBI using the gene name."""
    # cancer name --> gene name list
    if stage == 1:
        stage += 1
        with open("/reference_data/cancer_directory_filter.json", "r") as f:
            Tumor_catalog = json.load(f)
        f.close()
        # cancer_name_catalog = Tumor_catalog[cancer_name]
        response = "I will provide a catalog of tumors for you ,you can select the type of  Primary site, Site subtype 1, Primary histology, Histology subtype 1 about cancer which you insterst in.  If the cancer you want to study is not in our database, you can also choose to upload the target genes corresponding to the cancer yourself."
        return {"response": response, "operations":[{"title":"cancer_list","options":Tumor_catalog,"type":["single","input"],"key":"cancer_select","value":[]},
        {"title":"upload gene file", "options":[], "type":["file"], "key":"gene_path", "value":[]}], "stage": stage, "search_type": "cancer_type",
                "cancer_dict": {}, "state": 'continue'}

    # user select --> gene sequence and snp
    else:
        
        stage += 1
        gene_path = instruction['instruction']['gene_path']
        if len(gene_path) != 0:
            response = 'File upload accepted'
        
            return {"response": response, "operations":[],"target_gene_info": gene_path, "search_type": "cancer_type",
                    "search_system_prompt": cancer_search_system_prompt,"stage": stage, "state": 'stop'}
        # user_select_list = instruction
        name_list = instruction['instruction']['cancer_select']
        print('name_list', name_list)
        data = pd.read_csv("/reference_data/Cosmicfilter.tsv", header=0, sep='\t')
        primary_site = name_list[0].lower().strip().replace(' ', '_')
        site_subtype_1 = name_list[1].lower().strip().replace(' ', '_')
        primary_histology = name_list[2].lower().strip().replace(' ', '_')
        histology_subtype_1 = name_list[3].lower().strip().replace(' ', '_')


        # catalog list
        filtered_data = data[
            (data['Primary site'].str.lower() == primary_site) &
            (data['Site subtype 1'].str.lower() == site_subtype_1) &
            (data['Primary histology'].str.lower() == primary_histology) &
            (data['Histology subtype 1'].str.lower() == histology_subtype_1)
            ]
        name_mutation_csv = filtered_data.loc[:, ['Gene name', 'Chromosome', 'Gene start',
                                                  'Gene end', 'Start position in gene', 'End position in gene']]
        result_data = name_mutation_csv.to_dict()

        # file_path_list.append(fna_path)
        if len(name_mutation_csv) == 0:
            stage = stage - 2
            with open("/reference_data/cancer_directory.json", "r") as f:
                Tumor_catalog = json.load(f)
            f.close()

            response = "We are sorry that the relevant target genes were not found in our database for the type of disease you selected. You can reselect a type of disease.If the cancer you want to study is not in our database, you can also choose to upload the target genes corresponding to the cancer yourself"
            return {"response": response, "operations":[{"title":"cancer_select",
            "options":Tumor_catalog,"type":["single","input"],"key":"cancer_select","value":[]},{"title":"upload file", "options":[], "type":["file"], "key":"gene_path", "value":[]}],
             "stage": stage, "search_type": "cancer_type",
                    "search_system_prompt": cancer_search_system_prompt, "state": 'continue'}

        else:
            response = 'I will provide the snp information'
            # save snp
            # target_gene_dict[gene_seq]=SNP
            return {"response": response, "operations":[],"target_gene_info": result_data, "search_type": "cancer_type",
                    "search_system_prompt": cancer_search_system_prompt, "state": 'stop',"stage": stage}

                    
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

