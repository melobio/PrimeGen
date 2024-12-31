
import pandas as pd
import openai
import os
import json
import re
from llm_utils.utils import get_llm_chat_completion
from search_app import get_ner_entities, get_meta_info, check_synonym
from ncbi_function import sort_genomes
import zipfile
from pathlib import Path
import logging
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)
# whole_genome search tool
def whole_genome_search(instruction, stage):
    # description genetic disorder search tool
    whole_genome_test_system_prompt = """. 
    """
    if stage == 1:
        virus_identification_dict = {
                    "accession":[],
                    "taxid":[],
                    "genus": [],
                    "species": [],
                    "organism": [],
                    "infraspecific":[],
                    "isolate":[],
                    "assembly_level":[],
                    "ftp_path":[],
                    "group":[],
                    "genome_path":"",
                    "target_cds_path": "", 
                    "target_number":"" 
                    }
        # Use LLM to extract virus name
        experiment_type_prompt = """According to the user's contextual dialogue, you need to extract and summary the information needed for biological sequence retrieval.
You just need to summarize the extracted information and reply to me with a complete biological information description. You can think of it as rewriting a clear query. The result needs to be returned in json format.

Note:
1) The biological information description must be complete and professional. If it is a variant, strain, isolate, etc., just extract it directly. Remember not to use brackets or other additional descriptions.
2) Attention! If the user does not mention it, please reply null, do not guess the user's intention!
3) Strain must be abbreviated as "str.", substrain must be abbreviated as "substr.", subspecies must be abbreviated as "subsp.", and varieties must be abbreviated as "var.".
4) When extracting user information, do not modify punctuation marks, especially periods, brackets, dashes, underscores, qual sign, space etc.

Example1:
User: I want to conduct a multiplex PCR whole genome sequencing experiment for Influenza A virus (A/pika/Qinghai/BI/2007(H5N1)).
Assistant: {"virus_name":["Influenza A virus (A/pika/Qinghai/BI/2007(H5N1))"]}
----------
Example2:
User: I plan to perform a multiplex PCR-based whole genome sequencing experiment on the Vaccinia virus. 
Assistant: {"virus_name":["Vaccinia virus"]}
----------
Example3:
User: My goal is to conduct whole genome sequencing of the Acanthocystis turfacea Chlorella virus MO0605SPH using multiplex PCR.
Assistant: {"virus_name":["Acanthocystis turfacea Chlorella virus MO0605SPH"]}
----------
Example4:
User: I intend to perform a multiplex PCR experiment for the whole genome sequencing of Circoviridae sp. ctfvP2. 
Assistant: {"virus_name":["Circoviridae sp. ctfvP2"]}"""

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

        disease_json_data = extract_json_response(experiment_info_response)
        virus_name = disease_json_data['virus_name'][0]
        print('**************disease_name***********',virus_name)

        response = f'I will provide the sequenece of the {virus_name}'


        # Currently only consider users who fully say the virus name
        # TODO: Fuzzy matching and vector search
        # SPECIES_SEARCH_HOST = "10.49.60.23"
        # SPECIES_SEARCH_PORT = "8119"
        SPECIES_SEARCH_HOST = os.getenv('SPECIES_SEARCH_HOST')
        SPECIES_SEARCH_PORT = os.getenv('SPECIES_SEARCH_PORT')
        base_url = f"http://{SPECIES_SEARCH_HOST}:{SPECIES_SEARCH_PORT}"
        logging.info(f'**********virus_name**********：\n{virus_name}')
        virus_info_dict = get_seq_main(base_url, virus_name)
        logging.info(f'**********virus_info_dict**********：\n{virus_info_dict}')
        operations_list=[]
        for i,(x,y) in enumerate(virus_info_dict.items()):
            operations_list.append({"widget":"table",
                                    "title":f"The {x} is unclear, which may not be accurate. Please confirm again.(required)", 
                                    "options":y, 
                                    "type":["single","required"],
                                    "key":f"virus_select_{i+1}", 
                                    "value":[]})
        response = """It appears that the virus you're looking to retrieve has multiple variants. Could you please take a look at the table below and select the specific samples you'd like to include in the WGS experiment?"""
        logging.info(f'**********operations**********：\n{operations_list}')
        return {"response": response,
                            "operations": operations_list,
                            "stage": stage+1,
                            "state": 'continue',
                            "search_type": "whole_genome_type",
                            "data":{"species_identification_dict": virus_identification_dict} 
                            }

    # user select --> gene name list --> gene sequence and snp
    elif stage == 2:
        stage = stage + 1
        operation_result = [value[0] for key, value in instruction['instruction'].items() if key.startswith("virus_select_")]
        term = operation_result[0]["organism_name"]
        accession_id = operation_result[0]["accession"]
        print(f'term:{term}\n')
        print(f'accession:{accession_id}\n')
        # Get fna data path
        fna_files_list = get_ncbi_data(accession_id, term)
        # Let user upload SNP site information
        print('fna_files_list',fna_files_list)
        response = f'The sequence files corresponding to the pathogen you selected are as follows. Please confirm the model category for your subsequent primer design. If you need to upload a SNP file for design, please select the SNP type and upload the file.'
        return {"response": response,
                "operations": [{"title": "primer design type","options": ["Tiling","SNP"], "type": ["single","select"], "key": "design_type", "value": []},
                               {"title": "upload SNP file", "options": [], "type": ["file","option"], "key": "snp_file", "value": []}],
                "data":{"target_fna_list": fna_files_list},
                "search_type": 'whole_genome_type',
                "search_system_prompt": whole_genome_test_system_prompt, "stage": stage, "state": 'continue'}
    else:
        fna_files_list = instruction['instruction']['data']['target_fna_list']
        design_type = instruction['instruction']['design_type']
        snp_file = instruction['instruction']['snp_file']
        response ='The sequence search process has been completed and we will proceed to the subsequent primer design process.'
        return {"response": response,
                "data": {"target_fna_list": fna_files_list,"design_type":design_type,"snp_file":snp_file},
                "search_type": 'whole_genome_type',
                "search_system_prompt": whole_genome_test_system_prompt, "stage": stage, "state": 'stop'}

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


def get_seq_main(base_url, name):
    # Find corresponding name in database, could be genus, species, organism etc.
    ner_result = get_ner_entities(base_url, name)
    print(ner_result)
    # logging.info(f'**********ner_result**********：\n{ner_result}')
    # For COVID it's at organism level, so build dict since get_meta_info backend service has fixed format requirements
    extract_meta_info_dict={"genus":'', "species":'',"organism":'',"infraspecific":'',"isolate":''}
    extract_meta_info_dict['organism']=ner_result[0].get('text')
    print(extract_meta_info_dict)

    # Entity detection will find many species, same name but different isolates etc., organize them for user selection theoretically
    meta_result = get_meta_info(base_url, [extract_meta_info_dict])[0]["text"][0]

    title_list = ['accession','representative', 'taxid', 'genus', 'species', 'organism_name', 'infraspecific_name',  'isolate', 'assembly_level', 'ftp_path', 'group', 'name', 'origin', 'contig_n50', 'scaffold_n50', 'total_sequence_length']
    new_title_list = ['accession','representative','organism_name', 'isolate', 'assembly_level',  'contig_n50', 'scaffold_n50']
    related_info_dict = {}
    related_info_list = []
    for single_info in meta_result:
        info_dict={}
        single_info = single_info[:13] + single_info[-3:]
        for x,y in zip(title_list,single_info):
            if x in new_title_list:
                info_dict[x] = y

        related_info_list.append(info_dict)
    related_info_dict[extract_meta_info_dict["organism"]] = related_info_list

    sort_dict = sort_genomes(related_info_dict[extract_meta_info_dict["organism"]])
    # print(sort_dict)
    related_info_dict[extract_meta_info_dict["organism"]] = sort_dict[:50]
    return related_info_dict

def get_ncbi_data(accession_id, name):
    name = '_'.join(name.split(' '))
    file_path = f"/reference_data/dna_file/{name}.zip"
    output_path = f"/reference_data/dna_file/{name}"
    download_url = f"datasets download genome accession {accession_id} --filename {file_path}"
    os.system(download_url)

    # Extract zip file to specified "make file str" path
    with zipfile.ZipFile(file_path, 'r') as zip_ref:
        zip_ref.extractall(output_path)
    
    fna_files = []  # List to store file paths
    root_path = Path(output_path)
    for fna_file in root_path.rglob('*.fna'):
        fna_files.append(str(fna_file.absolute()))
        break
    return fna_files

