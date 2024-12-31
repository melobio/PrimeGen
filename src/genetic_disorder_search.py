
import pandas as pd
import openai
import os
import json
import re
from llm_utils.utils import get_llm_chat_completion
import logging
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)
# genetic disorder search tool
def genetic_disorder_search(instruction, stage):
    # description genetic disorder search tool
    genetic_disorder_system_prompt = """
    Current search agents are involved in search projects related to human genetic diseases.First, it requires user interaction to extract the entered human genetic disease name and retrieve the gene names related to the genetic disease on Omim.
And feed back the retrieved information, such as genetic disease phenotypes, related genes and other information.
Further user interaction will select the genetic disease and associated genes required for study. 
    """
    # genetic disorder name --> gene name list
    if stage == 1:
        # Use LLM to extract disease
        experiment_type_prompt = """Based on the user's contextual dialogue, you need to perform information extraction. The information to be extracted includes:  
        1.disease_name (human genetic diseases name,Return null if no diseases was mentioned in the user conversation)   
        The results should be returned to me in JSON format.  
        Example is shown below:
        '''json
        {  
          "disease_name": 'Alzheimer',   
        }
        '''
        """
        conversation = instruction['conversation']
        conversation_str = "\n\n".join([i["role"] + ':\n' + i['content'] for i in conversation][-2:])
        print('**********************conversation_str',conversation_str)
        experiment_info_response = get_llm_chat_completion(
            messages=[
                {"role": "system", "content": experiment_type_prompt},
                {"role": "user", "content": conversation_str}
            ],
            response_format={"type": "json_object"},
        )
        total_tokens = experiment_info_response.usage.total_tokens
        print(f"************Total tokens used: {total_tokens}***************")
        logging.info(f'*************Total tokens used*******：\n{total_tokens}')
        disease_json_data = extract_json_response(experiment_info_response)
        disease_name = disease_json_data['disease_name']
        print('**************disease_name***********',disease_name)
        response = f'I will provide the name of the genetic disease and the corresponding gene name related to {disease_name}'
        stage += 1
        # OMIM dataset
        # Extract gene names based on disorder name
        gene_result = get_disease_gene(disease_name)
        history_data = {}
        # Determine if this is first search or repeat search based on history_data
        if 'data' in instruction['instruction'].keys() and 'history_data' in instruction['instruction']['data'].keys():
            history_data = instruction['instruction']['data']['history_data']

        if isinstance(gene_result, str):
            response =  f"Sorry, the target gene of the disease {disease_name} was not found in our database. You can refer to the phenotype description in OMIM and ask the question again."
            return {"response": response, "stage": 1,
                    "search_type": "genetic_disorder_type", "data": {'history_data': history_data},
                    "search_system_prompt": genetic_disorder_system_prompt,
                    "state": 'continue'}
            # return {"response": gene_result, "operations":[{"title":f"Sorry, the target gene of the disease {disease_name} was not found in our database. You can refer to the phenotype description in OMIM and ask the question again or You can choose to upload the gene file yourself or study other diseases.", "options":[], "type":["file"], "key":"gene_path", "value":[]}], "stage": stage,
            #         "search_type": "genetic_disorder_type", "data":{'history_data':history_data}, "search_system_prompt": genetic_disorder_system_prompt,
            #         "state": 'continue'}
        else:
            # Return gene names, stage, type
            return {"response": response,"operations":[{"title":"gene and disease information","options":gene_result,"type":["multi","input"],"key":"gene_select","value":[]}], "stage": stage,
                    "search_type": "genetic_disorder_type", "data":{'history_data':history_data},"search_system_prompt": genetic_disorder_system_prompt,
                    "state": 'continue'}
    # user select --> gene name list --> gene sequence and snp
    elif stage == 2:
        stage += 1
        response = f'The step of determining the DNA sequence is complete.'
        history_data = instruction['instruction']['data']['history_data']
        if 'gene_path' in str(instruction['instruction']):
            gene_path = instruction['instruction']['gene_path']
            return {"response": response, "operations":[],
                    "data":{"target_gene_info": gene_path,"history_data":history_data}, "search_type": "genetic_disorder_type",
                    "search_system_prompt": genetic_disorder_system_prompt, "stage": stage,"state": 'stop'}
        
        target_gene_names_list = instruction['instruction']['gene_select']
        temp_gene = target_gene_names_list[0]
        if ',' in temp_gene:
            target_gene_names_list = temp_gene.split(',')
        # Open ClinVar dataset
        gene_info_df = pd.read_csv('/reference_data/gene_info_gencodev46.txt', sep='\t')
        print('target_gene_names_list', target_gene_names_list)
        if history_data:
            for ge in target_gene_names_list:
                history_data['gen_name_list'].append(ge)
        else:
            history_data = {'gen_name_list':target_gene_names_list}

        df_list = []
        for gene_name in target_gene_names_list:
            temp_df = gene_info_df[gene_info_df['gene'] == gene_name]
            df_list.append(temp_df)
        # Return target gene position information on chromosome
        res = pd.concat(df_list)
        if len(res) == 0:
            response = f'sorry, the gene {gene_name} seq not found'
            return {"response": response,
                    "data": {"target_gene_list": target_gene_names_list, "history_data": history_data},
                    "search_type": 'genetic_disorder_type',
                    "search_system_prompt": genetic_disorder_system_prompt, "stage": stage, "state": 'continue'}
        res.drop_duplicates()
        gene_names = list(res['gene'])
        end_list = list(res['end'])
        start_list = list(res['start'])
        chr_list = list(res['chrom'])
        # Define file paths
        fasta_file = "/reference_data/GCF_000001405.40_GRCh38.p14_genomic_filter.fna"
        dna_file_path = '/reference_data/dna_file'
        result = []
        len_list =[]
        # Open file and read content
        with open(fasta_file, 'r') as file:
            file_content = file.read().splitlines()
            file.close()
        for ii in range(len(gene_names)):
            na = gene_names[ii]
            temp_seq = get_dna_seq(chr_list[ii][3:], start_list[ii], end_list[ii], file_content)
            len_list.append(len(temp_seq))
            temp_path = os.path.join(dna_file_path, f'{na}.fasta')
            with open(temp_path, 'w') as f:
                f.write(f'>{na} {chr_list[ii]} {start_list[ii]} {end_list[ii]} \n')
                f.write(temp_seq)
                f.close()
            result.append(temp_path)
        if 'dna_seqs' in history_data.keys():
            for seq in result:
                history_data['dna_seqs'].append(seq)
        else:
            history_data['dna_seqs'] = result
        # history_data['dna_seqs'] = result
        response = f'The gene {target_gene_names_list} has been identified in relation to the current genetic disease. If you would like to search for additional genetic information related to another disease, please provide the name of the disease you want to explore. If not, we will proceed to the next stage of the primer design process.'
        return {"response": response,
                "data":{"target_gene_list": target_gene_names_list,"history_data":history_data},
                "search_type": 'genetic_disorder_type',
                "search_system_prompt": genetic_disorder_system_prompt, "stage": stage,"state": 'continue'}
    else:
        # Use LLM to determine whether to continue searching, if continue then call itself
        experiment_type_prompt = """Based on the user's contextual dialogue, you need to perform information extraction. The information to be extracted includes:  
                1.disease_name (human genetic diseases,Return null if no diseases was mentioned in the user conversation)   
                The results should be returned to me in JSON format.  
                Example is shown below:
                '''json
                {  
                  "disease_name": 'Alzheimer',   
                }
                '''
                """
        conversation = instruction['conversation']
        conversation_str = "\n\n".join([i["role"] + ':\n' + i['content'] for i in conversation][-1])
        experiment_info_response = get_llm_chat_completion(
            messages=[
                {"role": "system", "content": experiment_type_prompt},
                {"role": "user", "content": conversation_str}
            ],
            response_format={"type": "json_object"},
        )
        disease_json_data = extract_json_response(experiment_info_response)
        disease_name = disease_json_data['disease_name']
        history_data = instruction['instruction']['data']['history_data']
        print("xxxxxxxxxxdisease_namexxxxxxxxxxxxxx",disease_name)
        response = f'The step of determining the DNA sequence is complete.'
        if disease_name:
            result_response = genetic_disorder_search(instruction, stage=1)
            return result_response
        else:
            return {"response": response, "operations": [],'all_seqs':history_data['dna_seqs'],
                "data": {'history_data': history_data},
                "search_type": 'genetic_disorder_type',
                "search_system_prompt": genetic_disorder_system_prompt, "stage": stage, "state": 'stop'}

def get_dna_seq(chr_num, start, end,file_content):
    chr_name = f'chromosome {chr_num}'
    # Chromosome name in fasta file
    all_chr_name = 'Homo sapiens ' + chr_name + ', GRCh38.p14 Primary Assembly'
    print(all_chr_name)
    # Store extracted content
    extracted_content = ""
    # Whether chr is found
    found_chr = False
    # Search for lines containing "all_chr_name" in memory and extract content
    for i, line in enumerate(file_content):
        if not found_chr and all_chr_name in line:
            found_chr = True
        elif found_chr:
            if line.startswith('>'):
                # Stop extracting if it's a new title line
                break
            extracted_content += line
    temp_seq = extracted_content[int(start):int(end)].upper()
    return temp_seq
def extract_json_response(response):
    # Extract JSON content from response
    experiment_info_dict = response.choices[0].message.content
    try:
        #if json.loads(experiment_info_dict):
        return json.loads(experiment_info_dict)
    except:
        match = re.search(r'```json\n(.*?)\n```', experiment_info_dict, re.DOTALL)
        if match:
            experiment_info_json_str = match.group(1)
            experiment_info_json_data = json.loads(experiment_info_json_str)
        else:
            print("GPT JSON construction error, situation shown below:\n",experiment_info_dict)
            exit()
        return experiment_info_json_data

def get_disease_gene(disease_name):
    path = '/reference_data/pheno2gene.csv'
    gene_df = pd.read_csv(path)
    if disease_name is None:
        return "disease_name cannot be None"
    if '"' in disease_name:
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
