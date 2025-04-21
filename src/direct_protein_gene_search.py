import uvicorn
from Bio import Entrez
import os
import pandas as pd
import requests
import openai
import json
import re
import logging
import traceback
import subprocess
import zipfile
import tempfile
import shutil
from bioservices import UniProt
import io
from NCBI_API_use import retrieve_ncbi_ids, get_gene_detailed_info, download_fasta_files
import Levenshtein
from llm_utils.utils import get_llm_chat_completion
from protein_mutation_search import search_uniprot_one, user_select_protein_ncbi_download2, get_nucleotide_info

Entrez.email = "htang26186@gmail.com"
Entrez.api_key = "68d985ea43f2d14fe207179f9d1722720c08"#"db4b31edb00e220c3dd378908403eddbc308"
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)

# direct search tool
def direct_protein_gene_search(instruction, stage):

    
    # stage1: Identify user needs.
    if stage == 1:   
        direct_protein_gene_search_dict = {
            "Experiment_Direction": "direct protein gene search",
            "original_query": "",
            "species": "",
            "protein": "",
            "gene": [],
            "database_query": "",
            "class_type": ""
        }
        

        if 'gene_path' in instruction['instruction'] and instruction['instruction']['gene_path'] != []:
            response = 'Gene file successfully uploaded. Sequence file retrieval completed.'
            new_file_list = instruction['instruction']['gene_path']
            return {"response": response, "operations":[], "stage": 5,
                    "search_type": "direct_protein_gene_search_type",
                    "data":{"direct_protein_gene_search_dict": direct_protein_gene_search_dict, "gene_seq_path": new_file_list}, "state": 'stop'}
        
        experiment_type_prompt = """You are a bioinformatics assistant specializing in standardizing user queries related to proteins and genes. Your task is to provide the corresponding official or scientific species names, protein names, and gene names of a oral expressed target protein/gene in a given query and convert them into a standardized ncbi format for database searches.

# Instructions:
1. replace the original species word with its scientific Latin name (e.g., mouse → Mus musculus, human → Homo sapiens). If the user has not mentioned it, please provide the full species name, as a bioinformatician.
2. replace the original protein word with its official or scientific name (e.g., p53 → TP53 tumor suppressor protein). please extract it directly and do not make excessive adjustments or rewirte.
   * If the user has not mentioned it, returns null.
3. replace the original gene word with its official or scientific gene name (e.g., V segment → IGHV (Immunoglobulin Heavy Chain Variable).
    * Example: kappa light chain gene → IGKV.
    * If the user has not mentioned it, returns null.
4. You need to determine the question type of the user's query. There are only three types of questions: species, protein, and gene.
5. For fields that are missing from user input, you need to make your own judgment to supplement them. For example, if it is a scene species, you can supplement it. For missing gene fields, you need to consider carefully because there are many genes, unless it is a scene with a small number of genes.
The result should be returned to me in JSON format.

# Entity quality verfication: 
**1) For the extracted protein, species or gene name, do not append excessive explanation in the format of bracket**

# Example 1:
User query: "I need information about mouse antibody V segment"
{
  "original_query": "mouse antibody V segment",
  "standardized_query": {
    "species": "Mus musculus",
    "protein": "immunoglobulin variable segment",
    "gene": ["IGHV", "IGKV", "IGLV"],
    "database_query": "(Mus musculus[Organism]) AND (immunoglobulin variable segment[Protein] OR IGHV[Gene] OR IGKV[Gene] OR IGLV[Gene])",
    "class_type": "gene"
  }
}
# Example 2:
User query: "Looking for human p53 protein information"
{
  "original_query": "human p53 protein information",
  "standardized_query": {
    "species": "Homo sapiens",
    "protein": "P53 protein",
    "gene": null,
    "database_query": "(Homo sapiens[Organism]) AND (P53 protein[Protein])",
    "class_type": "protein"
  }
}
# Example 3:
User query: "Need data on E. coli lac operon"
{
  "original_query": "E. coli lac operon",
  "standardized_query": {
    "species": "Escherichia coli",
    "protein": "lac operon",
    "gene": null,
    "database_query": "(Escherichia coli[Organism]) AND (lac operon[Protein])",
    "class_type": "protein"
  }
}
# Example 4:
User query: "How to design primers for E.coli"
{
  "original_query": "How to design primers for E.coli",
  "standardized_query": {
    "species": "Escherichia coli",
    "protein": null,
    "gene": null,
    "database_query": "Escherichia coli[Organism]",
    "class_type": "species"
  }
}

# Example 5:
User query: "Give details about zebrafish actin genes"
{
  "original_query": "Give details about zebrafish actin genes",
  "standardized_query": {
    "species": "Danio rerio",
    "protein": "Actin",
    "gene": ["acta1", "actb1", "actc1"],
    "database_query": "(Danio rerio[Organism]) AND (Actin[Protein] OR acta1[Gene] OR actb1[Gene] OR actc1[Gene])",
    "class_type": "gene"
  }
}
# Example 6:
User query: "I want Arabidopsis flowering proteins"
{
  "original_query": "I want Arabidopsis flowering proteins",
  "standardized_query": {
    "species": "Arabidopsis thaliana",
    "protein": "Flowering proteins",
    "gene": null,
    "database_query": "(Arabidopsis thaliana[Organism]) AND (Flowering proteins[Protein])",
    "class_type": "protein"
  }
}
# Example 7:
User query: "Looking for rice genome sequences"
{
  "original_query": "Looking for rice genome sequences",
  "standardized_query": {
    "species": "Oryza sativa",
    "protein": null,
    "gene": null,
    "database_query": "Oryza sativa[Organism]",
    "class_type": "species"
  }
}
# Example 8:
User query: "Looking for Whey protein"
{
  "original_query": "Looking for Whey protein",
  "standardized_query": {
    "species": "Bos taurus",
    "protein": Whey protein,
    "gene": null,
    "database_query": "(Bos taurus[Organism]) AND (Whey [Protein])",
    "class_type": "protein"
  }
}
# Example 9:
User query: "prime design for Restriction Enzyme EcoRI"
{
  "original_query": "prime design for Restriction Enzyme EcoRI",
  "standardized_query": {
    "species": "Escherichia coli",
    "protein": "EcoRI restriction enzyme",
    "gene": null,
    "database_query": "(Escherichia coli[Organism]) AND (EcoRI restriction enzyme [Protein])",
    "class_type": "protein"
  }
}
# Example 10:
User query: "prime design for T4 DNA Ligase"
{
  "original_query": "prime design for T4 DNA Ligase",
  "standardized_query": {
    "species": "Enterobacteria phage T4",
    "protein": "T4 DNA Ligase",
    "gene": null,
    "database_query": "(Enterobacteria phage T4[Organism]) AND (T4 DNA Ligase [Protein])",
    "class_type": "protein"
  }
}
"""
        # check user send new query
        user_query = None
        if 'query_again' in instruction['instruction'] and instruction['instruction']['query_again'] is not None and instruction['instruction']['query_again'][0] is not None:
            user_query = instruction['instruction']['query_again'][0]
            logging.info(f"Using the new query provided by the user: {user_query}")
        
        # add conversation
        conversation = instruction['conversation']
        conversation_str = "\n\n".join([i["role"] + ':\n' + i['content'] for i in conversation[:-1]])
        
        # add user new query
        if user_query:
            conversation_str += f"Because of the ambiguity of the user's query, the database search yielded no results. User gave a new query:\n{user_query}"
        

        logging.info(f'-------stage1 input system:\n{experiment_type_prompt}')
        logging.info(f'-------stage1 input str:\n{conversation_str}')
        
        # Information extraction by llm prompt
        experiment_info_response = get_llm_chat_completion(
            messages=[
                {"role": "system", "content": experiment_type_prompt},
                {"role": "user", "content": conversation_str}
            ],
            response_format={"type": "json_object"},
        )
        logging.info(f'-------stage1 json extract llm response:\n{experiment_info_response}')
        # extract json data
        experiment_info_json_data = extract_json_response(experiment_info_response)
        logging.info(f'-------stage1 json extract parse:\n{experiment_info_json_data}')

        # update direct_protein_gene_search_dict
        direct_protein_gene_search_dict["original_query"] = experiment_info_json_data["original_query"] or ""
        direct_protein_gene_search_dict["species"] = experiment_info_json_data["standardized_query"]["species"] or ""
        direct_protein_gene_search_dict["protein"] = experiment_info_json_data["standardized_query"]["protein"] or ""
        gene_value = experiment_info_json_data["standardized_query"]["gene"] or []
        if isinstance(gene_value, str):
            direct_protein_gene_search_dict["gene"].append(gene_value)
        elif isinstance(gene_value, list):
            direct_protein_gene_search_dict["gene"] += gene_value
        direct_protein_gene_search_dict["database_query"] = experiment_info_json_data["standardized_query"]["database_query"] or ""
        direct_protein_gene_search_dict["class_type"] = experiment_info_json_data["standardized_query"]["class_type"] or "" 
        
        # Generate professional English description to inform the user about the extracted content
        extracted_info_message = f"Based on your query, we have extracted the following information:\n"
        extracted_info_message += f"- Your requirement: {direct_protein_gene_search_dict['original_query']}\n"
        extracted_info_message += f"- Species: {direct_protein_gene_search_dict['species'] or 'Not specified'}\n"
        extracted_info_message += f"- Protein: {direct_protein_gene_search_dict['protein'] or 'Not specified'}\n"
        extracted_info_message += f"- Gene(s): {', '.join(direct_protein_gene_search_dict['gene']) if direct_protein_gene_search_dict['gene'] else 'Not specified'}\n"
        extracted_info_message += f"\nWe will use this information to perform a targeted database search."

        # Provide options for the user to confirm or modify regardless of extraction results
        option_list = [
            {"title": 'Species', "options": [], "type": ["single", "input", "text"], "key": 'Species',
             "value": [direct_protein_gene_search_dict["species"]], "widget":"query"},
            {"title": 'Protein', "options": [], "type": ["single", "input", "text"], "key": 'Protein',
             "value": [direct_protein_gene_search_dict["protein"]],"widget":"query"},
            {"title": 'Gene Name', "options": [], "type": ["single", "input", "text"], "key": 'Gene_name',
             "value": [direct_protein_gene_search_dict["gene"][0] if direct_protein_gene_search_dict["gene"] else None],"widget":"query"}
        ]
        
        response = f"""{extracted_info_message}
        To enhance search accuracy, please verify or modify the extracted information. 
        If the information is correct, you may proceed by clicking submit. 
        If adjustments are needed, please modify the respective fields."""

        logging.info(f'-------response:\n{response}')
        logging.info(f'-------option_list:\n{option_list}')
        logging.info(f'-------direct_protein_gene_search_dict:\n{direct_protein_gene_search_dict}')
        stage += 1
        return {"response": response, "operations":option_list,
                    "stage": stage, "search_type": "direct_protein_gene_search_type", #new type: direct_protein_gene_search
                   "data":{"direct_protein_gene_search_dict": direct_protein_gene_search_dict},  "state": 'continue'}
    
    # stage2: Based on the requirements, use NCBI API to retrieve information and provide feedback to users.
    elif stage == 2:
        logging.info(f'-------start stage2 [direct]')
        direct_protein_gene_search_dict = instruction['instruction']['data']['direct_protein_gene_search_dict']
        raw_query = direct_protein_gene_search_dict["original_query"]
        database_query = direct_protein_gene_search_dict["database_query"]
        chang_flag = False
        # update direct_protein_gene_search_dict by user modify
        if instruction['instruction']["operations"][0]['value'][0] != ['']:
            if instruction['instruction']["operations"][0]['value'][0] == 'null':
                direct_protein_gene_search_dict["species"] = ''
            else:
                direct_protein_gene_search_dict["species"] = instruction['instruction']["operations"][0]['value'][0]
            chang_flag = True
        else:
            direct_protein_gene_search_dict["species"] = instruction['instruction']["Species"]
        
        if instruction['instruction']["operations"][1]['value'][0] != ['']:
            if instruction['instruction']["operations"][1]['value'][0] == 'null':
                direct_protein_gene_search_dict["protein"] = ''
            else:
                direct_protein_gene_search_dict["protein"] = instruction['instruction']["operations"][1]['value'][0]
            chang_flag = True
        else:
            direct_protein_gene_search_dict["protein"] = instruction['instruction']["Protein"]
        
        if instruction['instruction']["operations"][2]['value'][0] != ['']:
            if instruction['instruction']["operations"][2]['value'][0] == 'null':
                direct_protein_gene_search_dict["gene"] = ''
            else:
                direct_protein_gene_search_dict["gene"] = instruction['instruction']["operations"][2]['value'][0]
            chang_flag = True
        else:
            direct_protein_gene_search_dict["gene"] = instruction['instruction']["gene"]

        type = direct_protein_gene_search_dict["class_type"]
        if chang_flag:
            database_query=""
            database_quer_species=""
            database_query_protein=""
            database_query_gene=""

            if direct_protein_gene_search_dict['species'] != '':
                database_quer_species = f'({direct_protein_gene_search_dict["species"]}[Organism])'
            if direct_protein_gene_search_dict['protein'] != '':
                database_query_protein = f'({direct_protein_gene_search_dict["protein"]}[Protein])'
            if direct_protein_gene_search_dict['gene'] != '':
                database_query_gene= f'({direct_protein_gene_search_dict["gene"]}[Gene])'
            database_query = ' AND '.join([i for i in[database_quer_species, database_query_protein, database_query_gene] if i!=''])
            logging.info(f'-------stage2 [direct]: user new database_query:{database_query}')

        
        logging.info(f'-------stage2 [direct]: \ndirect_protein_gene_search_dict{direct_protein_gene_search_dict}')
        logging.info(f'-------stage2 [direct]: database_query:{database_query}')
        
        # Protein priority search using uniprot
        if type == "protein":
            # Note that search_uniprot_one only takes [0]
            gene_list = direct_protein_gene_search_dict["gene"]
            # if gene_list is not empty, search_gene is the first gene
            if gene_list:
                search_gene = gene_list[0]
            else:
                search_gene = None
            
            # search uniprot
            _, download_infor = search_uniprot_one(direct_protein_gene_search_dict["protein"], direct_protein_gene_search_dict["species"], search_gene)
            # if no uniprot result, use easy query to search
            if len(download_infor) == 0:
                esay_database_query = get_easy_query(database_query)
                download_infor = query_uniprot(esay_database_query)
                # easy query result is not empty
                if len(download_infor) != []:
                    bioserver_list = []
                    for entry_id in download_infor:
                        _, download_infor_entry_id = search_uniprot_one([],[],[],entry_id,type='entry')
                        bioserver_list += download_infor_entry_id

                    download_infor = bioserver_list
                    logging.info(f'-------stage2 [direct]: search_uniprot_one entry result: {download_infor}')
                    # easy query entry id is empty
                    if len(download_infor) == []:
                        download_infor = retrieve_ncbi_ids(database_query, db_type='gene')
                        if len(download_infor) == 0:
                            # Use easy query to ensure retrieval fallback.
                            esay_database_query = get_easy_query(database_query)
                            logging.info(f'--esay_database_query-----{esay_database_query}')
                            download_infor = retrieve_ncbi_ids(esay_database_query, db_type='gene')
                            # If still no gene found, use nucleotide database.
                            if len(download_infor) == 0:
                                logging.info(f'-------stage2 retrieval fallback: no gene found, use nucleotide database')
                                download_infor = retrieve_ncbi_ids(esay_database_query, db_type='nucleotide')
                                direct_protein_gene_search_dict["class_type"] = 'nucleotide'
                                
                                #if len(download_infor) == 0:
                                #    logging.info(f'-------stage2 retrieval fallback: no nucleotide found, use esay_database_query')
                                #    download_infor = retrieve_ncbi_ids(esay_database_query, db_type='nucleotide')
                    
                    else:
                        download_infor_list = []
                        seq_list = []
                        for d_i in download_infor:
                            dict_ = {}
                            dict_['Protein ID'] = d_i.get('Protein ID', 'None')
                            dict_['Protein Name'] = d_i.get('Protein Name', 'None')
                            dict_['Organism'] = d_i.get('Organism', 'None')
                            dict_['Length'] = str(d_i.get('Length', 'None'))
                            function_info = d_i.get('Function information', {})
                            if isinstance(function_info, dict):
                                dict_['FUNCTION'] = function_info.get('FUNCTION', 'None')
                            else:
                                dict_['FUNCTION'] = 'None'
                            dict_['Protein Name fasta'] = d_i.get('Protein Name fasta', 'None')
                            dict_['DNA List'] = str(d_i.get('DNA List', 'None'))
                            seq_list.append(d_i.get('Protein ID', '') + '***' + d_i.get('Protein Sequence', ''))
                            # dict_['Protein Sequence'] = d_i['Protein Sequence']
                            download_infor_list.append(dict_)
                        direct_protein_gene_search_dict['uniprot_download'] = True
                        logging.info(f'-------stage2 [direct]: 搜索uniprot结果: {download_infor_list}')

                        seq_path = '/reference_data/dna_file/temp_seq.fasta'
                        with open(seq_path, 'w') as f:
                            for se in seq_list:
                                f.write(se + '\n')
                            f.close()
                        logging.info(f'-------stage2 seq_path: {seq_path}')
                        response = f'I have found the following protein information in the Uniprot. Please select the most relevant protein target for your research needs:'
                        return {"response": response, "operations": [
                            {"title": "protein information", "options": download_infor_list[:10],
                             "type": ["single", "input"],
                             "key": "gene_select", "value": [], "widget": "table"}],
                                "stage": stage + 1, "search_type": "direct_protein_gene_search_type",
                                "data": {"direct_protein_gene_search_dict": direct_protein_gene_search_dict,
                                         "seq_path": seq_path}, "state": 'continue'}
                
                else:
                    download_infor = retrieve_ncbi_ids(database_query, db_type='gene')
                    if len(download_infor) == 0:
                        # Use easy query to ensure retrieval fallback.
                        esay_database_query = get_easy_query(database_query)
                        logging.info(f'--esay_database_query-----{esay_database_query}')
                        download_infor = retrieve_ncbi_ids(esay_database_query, db_type='gene')
                        # If still no gene found, use nucleotide database.
                        if len(download_infor) == 0:
                            logging.info(f'-------stage2 retrieval fallback: no gene found, use nucleotide database')
                            download_infor = retrieve_ncbi_ids(esay_database_query, db_type='nucleotide')
                            direct_protein_gene_search_dict["class_type"] = 'nucleotide'
                            
                            #if len(download_infor) == 0:
                            #    logging.info(f'-------stage2 retrieval fallback: no nucleotide found, use esay_database_query')
                            #    download_infor = retrieve_ncbi_ids(esay_database_query, db_type='nucleotide')
            
            else:
                download_infor_list=[]
                seq_list = []
                for d_i in download_infor:
                    dict_={}
                    dict_['Protein ID'] = d_i.get('Protein ID', 'None')
                    dict_['Protein Name'] = d_i.get('Protein Name', 'None')
                    dict_['Organism'] = d_i.get('Organism', 'None')
                    dict_['Length'] = str(d_i.get('Length', 'None'))
                    function_info = d_i.get('Function information', {})
                    if isinstance(function_info, dict):
                        dict_['FUNCTION'] = function_info.get('FUNCTION', 'None')
                    else:
                        dict_['FUNCTION'] = 'None'
                    dict_['Protein Name fasta'] = d_i.get('Protein Name fasta', 'None')
                    dict_['DNA List'] = str(d_i.get('DNA List', 'None'))
                    seq_list.append(d_i.get('Protein ID', '') + '***' + d_i.get('Protein Sequence', ''))
                    # dict_['Protein Sequence'] = d_i['Protein Sequence']
                    download_infor_list.append(dict_)
                direct_protein_gene_search_dict['uniprot_download'] = True
                logging.info(f'-------stage2 [direct]: 搜索uniprot结果: {download_infor_list}')
                
                seq_path = '/reference_data/dna_file/temp_seq.fasta'
                with open(seq_path, 'w') as f:
                    for se in seq_list:
                        f.write(se + '\n')
                    f.close()
                logging.info(f'-------stage2 seq_path: {seq_path}')
                response = f'I have found the following protein information in the NCBI. Please select the most relevant protein target for your research needs:'
                return {"response": response, "operations": [
                    {"title": "protein information", "options": download_infor_list[:10], "type": ["single", "input"],
                     "key": "gene_select", "value": [], "widget": "table"}],
                        "stage": stage+1, "search_type": "direct_protein_gene_search_type",
                        "data": {"direct_protein_gene_search_dict": direct_protein_gene_search_dict,
                                 "seq_path": seq_path}, "state": 'continue'}

        # use ncbi api to search gene
        else:
            download_infor = retrieve_ncbi_ids(database_query, db_type='gene')
            if len(download_infor) == 0:
                # Use easy query to ensure retrieval fallback.
                esay_database_query = get_easy_query(database_query)
                download_infor = retrieve_ncbi_ids(esay_database_query, db_type='gene')
                # If still no gene found, use nucleotide database.
                if len(download_infor) == 0:
                    logging.info(f'-------stage2 retrieval fallback: no gene found, use nucleotide database')
                    download_infor = retrieve_ncbi_ids(esay_database_query, db_type='nucleotide')
                    direct_protein_gene_search_dict["class_type"] = 'nucleotide'

                    
                    #if len(download_infor) == 0:
                    #    logging.info(f'-------stage2 retrieval fallback: no nucleotide found, use esay_database_query')
                    #    download_infor = retrieve_ncbi_ids(esay_database_query, db_type='nucleotide')

        logging.info(f'-------stage2 [direct]: search result length: {len(download_infor)}')
        
        if len(download_infor) == 0:
            rewrite_attempts = 0
            max_attempts = 5
            
            while rewrite_attempts < max_attempts:
                rewrite_prompt = """You are an expert in the field of bioinformatics. 
                Under the user query, the extracted content is shown in the following JSON. 
                After searching the NCBI API, no searched content was found. 
                Please re-understand the user's needs and search objects, 
                ensure that the search query of the NCBI API is correct, 
                and return a new search JSON according to the format.
                """
                conversation = instruction['conversation']
                conversation_str = "\n\n".join([i["role"] + ':\n' + i['content'] for i in conversation])
                conversation_str += f"""\n\nThe extracted content is shown in the following JSON:\n{direct_protein_gene_search_dict}\n\n
                This is not retrieved by the NCBI API, and may be missing information or vague description. Please imitate the following JSON format and re-complete the extracted information.
                Json format:\n
                {{
                "original_query": "I want to sequence Tag DNA Polymerase",
                "standardized_query": {{
                    "species": "Thermus aquaticus",
                    "protein": " Tag DNA Polymerase",
                    "gene": 'polA',
                    "database_query": "(Thermus aquaticus[Organism]) AND (Taq DNA polymerase[Protein] OR polA[Gene])",
                    "class_type": "protein"
                }}
                }}
                """
                experiment_info_response = get_llm_chat_completion(
                    messages=[
                        {"role": "system", "content": rewrite_prompt},
                        {"role": "user", "content": conversation_str}
                    ],
                    response_format={"type": "json_object"},
                )
                
                # Extract and update the dictionary with new search parameters
                experiment_info_json_data = extract_json_response(experiment_info_response)
                logging.info(f"----------self try again({rewrite_attempts+1}):", experiment_info_json_data)
                direct_protein_gene_search_dict["species"] = experiment_info_json_data["standardized_query"]["species"] or  ""
                # direct_protein_gene_search_dict["protein"] = experiment_info_json_data["standardized_query"]["protein"] or ""
                direct_protein_gene_search_dict["gene"]=[]
                gene_value = experiment_info_json_data["standardized_query"]["gene"] or []
                if isinstance(gene_value, str):
                    direct_protein_gene_search_dict["gene"].append(gene_value)
                elif isinstance(gene_value, list):
                    direct_protein_gene_search_dict["gene"] += gene_value
                direct_protein_gene_search_dict["database_query"] = experiment_info_json_data["standardized_query"]["database_query"] or ""
                
                # Try searching again with the new parameters
                if rewrite_attempts > 3:
                    download_infor = retrieve_ncbi_ids(direct_protein_gene_search_dict["database_query"], "protein")
                else:
                    download_infor = retrieve_ncbi_ids(direct_protein_gene_search_dict["database_query"], direct_protein_gene_search_dict["class_type"])
                
                if len(download_infor) > 0:
                    # If we found results, break out of the loop
                    logging.info(f'-------stage2 [direct]: 重试次数{rewrite_attempts+1}次，搜索结果数量{len(download_infor)}')
                    break
                
                rewrite_attempts += 1
            
            if rewrite_attempts >= max_attempts:
                logging.info(f'-------stage2 [direct]: 重试次数超过最大次数')
                stage = 1  
                response = f"""After multiple attempts, we were unable to find any results matching your query in the NCBI database. I sincerely apologize for the failure in downloading the sequence. As an alternative, you may upload an existing gene file, 
                and I will gladly proceed with primer design for you directly. Alternatively, you may follow the instructions below to submit a new question or restate your previous question, 
                allowing me to attempt the retrieval again. 
                Additionally, you can simply type **restart** in the dialogue box to reset our conversation and avoid any unintended influences from previous interactions."""
                return {"response": response, "operations":[{"title":"upload gene file", "options":[], "type":["file",'require-one'], "key":"gene_path", "value":[]},
                                                            {"title": 'Possibly due to ambiguity in your question, the database search yielded no results. Please clarify your requirements again.', "options": [], "type": ["single", "input", "text","require-one"], "key": 'query_again',"value": [None], "widget":"query"}], 
                        "stage": stage, "search_type": "direct_protein_gene_search_type",
                        "data":{"direct_protein_gene_search_dict": direct_protein_gene_search_dict,'restart':True}, "state": 'continue'}
            else:
                logging.info(f'-------stage2 [direct]: 重试次数未超过最大次数')
                logging.info(f'-------stage2 retrydownload_information:{download_infor}')
                logging.info(f'-------stage2 [direct]: 开始获取详细信息')
                stage += 1
                download_infor_list = []
                for id in download_infor:
                  
                    logging.info(f'-------stage2 [direct]: 获取详细信息: {id}')
                    try:
                        gene_details = get_gene_detailed_info(id, output_file= f"/reference_data/dna_file/{raw_query.replace(' ','_').replace('/','_')}_info.csv", db_type=direct_protein_gene_search_dict["class_type"])
                    except Exception as e:
                        import traceback
                        logging.info(f"获取详细信息失败: {traceback.format_exc()}")
                       
                        import time
                        logging.info("API可能暂时不可用，等待15秒后重试...")
                        time.sleep(10)
                        try:
                            gene_details = get_gene_detailed_info(id, output_file=f"/reference_data/dna_file/{raw_query.replace(' ','_').replace('/','_')}_info.csv", db_type=direct_protein_gene_search_dict["class_type"])
                            logging.info(f"重试成功获取详细信息: {id}")
                        except Exception as e:
                            import traceback
                            logging.info(f"重试后仍然失败: {traceback.format_exc()}")
                            logging.info("跳过此ID，继续处理其他ID")
                            continue

                    if gene_details:
                        gene_details['ncbi_id'] = id
                        download_infor_list.append(gene_details)
                logging.info(f'-------stage2 [direct]: 获取详细信息: {download_infor_list}')
                seq_path =''
                if download_infor_list==[]:
                    
                    _, download_infor = search_uniprot_one(direct_protein_gene_search_dict["protein"], direct_protein_gene_search_dict["species"], direct_protein_gene_search_dict["gene"])
                    logging.info(f'-------stage2 [direct]: 搜索uniprot结果: {download_infor_list}')
                    download_infor_list=[]
                    seq_list = []
                    for d_i in download_infor:
                        dict_={}
                        dict_['Protein ID'] = d_i.get('Protein ID', 'None')
                        dict_['Protein Name'] = d_i.get('Protein Name', 'None')
                        dict_['Organism'] = d_i.get('Organism', 'None')
                        dict_['Length'] = str(d_i.get('Length', 'None'))
                        function_info = d_i.get('Function information', {})
                        if isinstance(function_info, dict):
                            dict_['FUNCTION'] = function_info.get('FUNCTION', 'None')
                        else:
                            dict_['FUNCTION'] = 'None'
                        dict_['Protein Name fasta'] = d_i.get('Protein Name fasta', 'None')
                        dict_['DNA List'] = str(d_i.get('DNA List', 'None'))
                        seq_list.append(d_i.get('Protein ID', '') + '***' + d_i.get('Protein Sequence', ''))
                        # dict_['Protein Sequence'] = d_i['Protein Sequence']
                        download_infor_list.append(dict_)
                    direct_protein_gene_search_dict['uniprot_download'] = True
                    logging.info(f'-------stage2 [direct]: 搜索uniprot结果: {download_infor_list}')

                    
                    seq_path  = '/reference_data/dna_file/temp_seq.fasta'
                    with open(seq_path,'w') as f:
                        for se in seq_list:
                            f.write(se+'\n')
                        f.close()
                    logging.info(f'-------stage2 seq_path: {seq_path}')
                if direct_protein_gene_search_dict["class_type"] == "protein":
                    response = f'I have found the following protein information in the NCBI. Please select the most relevant protein target for your research needs:'
                elif direct_protein_gene_search_dict["class_type"] == "gene":
                    response = f'I have found the following gene information in the NCBI. Please select the most relevant genetic element for your research needs:'
                else:  # species
                    response = f'I have found the following taxonomic information in the NCBI. Please select the most relevant organism for your research needs:'
                return {"response": response, "operations":[{"title":"protein information","options":download_infor_list[:10],"type":["single","input"],"key":"gene_select","value":[],"widget":"table"}],
                        "stage": stage, "search_type": "direct_protein_gene_search_type",
                        "data":{"direct_protein_gene_search_dict": direct_protein_gene_search_dict,"seq_path":seq_path},"state": 'continue'}

        else:
            logging.info(f'-------stage2 [direct]: 搜索结果数量大于0')
            logging.info(f'-------stage2 [direct]: 开始获取详细信息: {download_infor}')
            stage += 1
            download_infor_list = []
            class_type = direct_protein_gene_search_dict["class_type"]
            if class_type == 'gene':
                db_type = 'gene'
            else:
                db_type = 'nucleotide'
            for id in download_infor:
                try:
                    
                    gene_details = get_gene_detailed_info(id, output_file= f"/reference_data/dna_file/{raw_query.replace(' ','_').replace('/','_')}_info.csv", db_type=db_type)
                    
                except Exception as e:
                    import traceback
                    logging.info(f"获取详细信息失败: {traceback.format_exc()}")
                    
                    import time
                    logging.info("API可能暂时不可用，等待15秒后重试...")
                    time.sleep(10)
                    try:
                        gene_details = get_gene_detailed_info(id, output_file=f"/reference_data/dna_file/{raw_query.replace(' ','_').replace('/','_')}_info.csv", db_type=direct_protein_gene_search_dict["class_type"])
                        logging.info(f"重试成功获取详细信息: {id}")
                    except Exception as e:
                        import traceback
                        logging.info(f"重试后仍然失败: {traceback.format_exc()}")
                        logging.info("跳过此ID，继续处理其他ID")
                        continue

                if gene_details:
                    gene_details['ncbi_id'] = id
                    download_infor_list.append(gene_details)
            logging.info(f'-------stage2 [direct]: 获取详细信息完成')
            logging.info(f'-------stage2 download_infor_list: {download_infor_list}')
            if download_infor_list==[]:
                stage = 1  
                response = f"""After multiple attempts, we were unable to find any results matching your query in the NCBI database. I sincerely apologize for the failure in downloading the sequence. As an alternative, you may upload an existing gene file, 
                and I will gladly proceed with primer design for you directly. Alternatively, you may follow the instructions below to submit a new question or restate your previous question, 
                allowing me to attempt the retrieval again. 
                Additionally, you can simply type **restart** in the dialogue box to reset our conversation and avoid any unintended influences from previous interactions."""
                return {"response": response, "operations":[{"title":"upload gene file", "options":[], "type":["file",'require-one'], "key":"gene_path", "value":[]},
                                                            {"title": 'Possibly due to ambiguity in your question, the database search yielded no results. Please clarify your requirements again.', "options": [], "type": ["single", "input", "text","require-one"], "key": 'query_again',"value": [None], "widget":"query"}], 
                        "stage": stage, "search_type": "direct_protein_gene_search_type",
                        "data":{"direct_protein_gene_search_dict": direct_protein_gene_search_dict,'restart':True}, "state": 'continue'}
            
            if direct_protein_gene_search_dict["class_type"] == "protein":
                response = f'I have found the following protein information in the NCBI. Please select the most relevant protein target for your research needs:'
            elif direct_protein_gene_search_dict["class_type"] == "gene":
                response = f'I have found the following gene information in the NCBI. Please select the most relevant genetic element for your research needs:'
            else:  # species
                response = f'I have found the following taxonomic information in the NCBI. Please select the most relevant organism for your research needs:'
            return {"response": response, "operations":[{"title":"protein information","options":download_infor_list[:10],"type":["single","input"],"key":"gene_select","value":[],"widget":"table"}],
                    "stage": stage, "search_type": "direct_protein_gene_search_type",
                    "data":{"direct_protein_gene_search_dict": direct_protein_gene_search_dict, "seq_path":''}, "state": 'continue'}
        


    elif stage == 3:
        logging.info(f'-------start stage3 [direct]')
        logging.info(f'-------stage3 [direct]: instruction: {instruction["instruction"]}')
        direct_protein_gene_search_dict = instruction['instruction']['data']['direct_protein_gene_search_dict']
        logging.info(f'-------stage3 [direct]: direct_protein_gene_search_dict: {direct_protein_gene_search_dict}')
        if 'gene_select' in instruction['instruction']:
            logging.info(f'-------stage3 [direct]: gene_select: {instruction["instruction"]["gene_select"]}')
            data_list = instruction['instruction']['gene_select']  
            download_path = []
            logging.info(f'--data_list key---{data_list[0].keys()}--')
            try:
                if direct_protein_gene_search_dict["class_type"] == "gene" and 'ncbi_id' in data_list[0].keys():

                    ids = []
                    for data_ in data_list:
                        logging.info(f'-------stage3 [direct]: data_: {data_}')
                        ids.append(data_['ncbi_id'])
                    logging.info(f'-------stage3 [direct]: ids: {ids}')
                    download_path = download_fasta_files(ids=ids, output_dir="/reference_data/dna_file/",
                                                         db_type="gene")

                elif direct_protein_gene_search_dict["class_type"] == "protein" or 'Protein ID' in data_list[0].keys():

                    seq_path = instruction['instruction']['data']['seq_path']
                    seq_list = []
                    with open(seq_path,'r') as f:
                        lines = f.readlines()
                        for ll in lines:
                            temp_dict = {}
                            temp_list = ll.split('***')
                            temp_dict['Protein ID'] = temp_list[0]
                            temp_dict['Protein sequence'] = temp_list[1]
                            seq_list.append(temp_dict)
                        f.close()
                    logging.info(f'-------seq_list: ids: {seq_list}')
                    aa_dna_content, download_path = user_select_protein_ncbi_download2(data_list,seq_list)
                    logging.info(f'-------stage3 download_path: {download_path}')
                else:
                    ids = []
                    for data_ in data_list:
                        logging.info(f'-------stage3 [direct]: data_: {data_}')
                        ids.append(data_['ncbi_id'])
                    download_path = download_fasta_files(ids=ids, output_dir="/reference_data/dna_file/", db_type=direct_protein_gene_search_dict["class_type"], details_list=data_list)

                #download_path = download_fasta_files(ids=ids, output_dir="/reference_data/dna_file/", db_type=direct_protein_gene_search_dict["class_type"])
                if download_path:
                    seqs = ''
                    with open(download_path[0], 'r') as f:
                        lines = f.readlines()
                        for ll in lines:
                            if '>' in ll:
                                continue
                            seqs += ll.replace('\n', '').strip()
                    seq_len = len(seqs)
                    response = f'Gene search has completed and the length of the gene sequence ranges from 1 to {seq_len}. Please further provide the target region for the designed primers, and we will design primers for this region.'
                    stage += 1
                    return {"response": response,  "stage": stage,
                        "search_type": "direct_protein_gene_search_type",
                        "operations": [
                                {"title": "start", "options": [], "type": ["input", "number","required"], "key": "start_pos",
                                 "value": [], "widget": "query"},
                                {"title": "end", "options": [], "type": ["input", "number","required"], "key": "end_pos",
                                 "value": [], "widget": "query"}
                            ],
                        "data":{"direct_protein_gene_search_dict": direct_protein_gene_search_dict,"gene_seq_path": download_path}, "state": 'continue'}
                else:
                    response = f"""I sincerely apologize for the failure in downloading the sequence. As an alternative, you may upload an existing gene file, 
                and I will gladly proceed with primer design for you directly. Alternatively, you may follow the instructions below to submit a new question or restate your previous question, 
                allowing me to attempt the retrieval again. 
                Additionally, you can simply type "restart" in the dialogue box to reset our conversation and avoid any unintended influences from previous interactions."""
                    stage = 1
                    return {"response": response, "operations":[{"title":"upload gene file", "options":[], "type":["file",'require-one'], "key":"gene_path", "value":[]},
                                                                {"title": 'Possibly due to ambiguity in your question, the database search yielded no results. Please clarify your requirements again.', "options": [], "type": ["single", "input", "text","require-one"], "key": 'query_again',"value": [None], "widget":"query"}], 
                            "stage": stage, "search_type": "direct_protein_gene_search_type",
                        "data":{"direct_protein_gene_search_dict": direct_protein_gene_search_dict,'restart':True}, "state": 'continue'}
                
            except Exception as e:
                import traceback
                logging.info(f"下载失败: {traceback.format_exc()}")
                response = f"""I sincerely apologize for the failure in downloading the sequence. As an alternative, you may upload an existing gene file, 
                and I will gladly proceed with primer design for you directly. Alternatively, you may follow the instructions below to submit a new question or restate your previous question, 
                allowing me to attempt the retrieval again. 
                Additionally, you can simply type "restart" in the dialogue box to reset our conversation and avoid any unintended influences from previous interactions."""
                stage = 1
                return {"response": response, "operations":[{"title":"upload gene file", "options":[], "type":["file",'require-one'], "key":"gene_path", "value":[]},
                                                            {"title": 'Possibly due to ambiguity in your question, the database search yielded no results. Please clarify your requirements again.', "options": [], "type": ["single", "input", "text","require-one"], "key": 'query_again',"value": [None], "widget":"query"}], 
                        "stage": stage, "search_type": "direct_protein_gene_search_type",
                        "data":{"direct_protein_gene_search_dict": direct_protein_gene_search_dict,'restart':True}, "state": 'continue'}
                        
        elif 'gene_path' in instruction['instruction']:
            logging.info(f'-------stage3 [direct]: gene_path: {instruction["instruction"]["gene_path"]}')
            new_file_list = instruction['instruction']['gene_path']
            seqs = ''
            with open(new_file_list[0], 'r') as f:
                lines = f.readlines()
                for ll in lines:
                    if '>' in ll:
                        continue
                    seqs += ll.replace('\n', '').strip()
            seq_len = len(seqs)
            response = f'Gene search has completed and the length of the gene sequence ranges from 1 to {seq_len}. Please further provide the target region for the designed primers, and we will design primers for this region.'
            stage += 1
            return {"response": response,
                    "operations": [
                                {"title": "start", "options": [], "type": ["input", "number"], "key": "start_pos",
                                 "value": [], "widget": "query"},
                                {"title": "end", "options": [], "type": ["input", "number"], "key": "end_pos",
                                 "value": [], "widget": "query"}
                            ], "stage": stage,
                   "search_type": "direct_protein_gene_search_type",
                   "data":{"direct_protein_gene_search_dict": direct_protein_gene_search_dict,
                           "gene_seq_path": new_file_list}, "state": 'continue'}
    elif stage == 4:
        response = 'The step of determining the DNA sequence is complete.'
        stage+=1
        new_file_list = instruction['instruction']['data']['gene_seq_path']
        direct_protein_gene_search_dict = instruction['instruction']['data']['direct_protein_gene_search_dict']
        direct_protein_gene_search_dict['start_pos'] = instruction['instruction']['start_pos'][0]
        direct_protein_gene_search_dict['end_pos'] = instruction['instruction']['end_pos'][0]
        return {"response": response,  "operations":[],"stage": stage,
                "search_type": "direct_protein_gene_search_type",
                "data":{"direct_protein_gene_search_dict": direct_protein_gene_search_dict,"gene_seq_path": new_file_list},
                "state": 'stop'}

def get_easy_query_old(s):
    import re
    
    s = re.sub(r'\[.*?\]', '', s)
    
    s = re.sub(r'[()]', '', s)  
    s = re.sub(r'\s+AND\s+', ' ', s)  
  
    return ' '.join(s.split()).strip()

def get_easy_query(s):
    
    s = re.sub(r'\[.*?\]', '', s)
    
   
    organism_match = re.search(r'\(([^()]+)\)\s+AND\s+\((.*)\)', s)
    if organism_match:
        organism = organism_match.group(1).strip()
        query_terms = organism_match.group(2).strip()
        
        return f"{organism} ({query_terms})"
    
    
    s = re.sub(r'\s+AND\s+', ' ', s)  
    s = re.sub(r'\s+OR\s+', ' OR ', s)  
    return ' '.join(s.split()).strip()

def extract_json_response(response):
    
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
            logging.info(f"GPT构造JSON错误，情况如下所示：\n{experiment_info_dict}")
            exit()
        return experiment_info_json_data
    
def query_uniprot(query_term: str, columns: str = "id,length,accession, gene_names", frmt: str = "tsv") -> pd.DataFrame:

    u = UniProt()
    result_text = u.search(query_term, columns=columns, frmt=frmt)#, limit=20)
    
    result_df = pd.read_csv(io.StringIO(result_text), sep="\t")
    if len(result_df)>6:
        result_df = result_df.head(6)
    
    if result_df.empty:
        print("未找到匹配的UniProt条目")
        return []
    
   
    print(f"找到 {len(result_df)} 个UniProt条目")
    
   
    if 'Entry' in result_df.columns:
        return result_df['Entry Name'].tolist()
    else:
        print("警告：结果中没有'Entry'列")
        return []



def is_match(str1, str2, threshold=0.7):

    if not str1 and not str2:  
        return True, 1.0
    if not str1 or not str2:   
        return False, 0.0

    distance = Levenshtein.distance(str1, str2)

    max_len = max(len(str1), len(str2))
    similarity = 1 - (distance / max_len)
    
    return similarity >= threshold, similarity