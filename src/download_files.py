
import pandas as pd 
import requests
import os
from Bio import Entrez
import json
import re

def download_files(instruction,stage):
    download_system_prompt = ''
    if stage == 1:
        stage += 1
        protein_name = instruction
        protein_info, data_list = search_uniprot_one(protein_name)
        if data_list ==[]:
            response = f'No {instruction} relevant sequence was found on the NCBI website. Sorry for the inconvenience, you can change the name to get the information you want.'
            return {"response": response, "operations":[], "data_list": data_list, "stage": stage,
                    "search_type": "download_type", "search_system_prompt": download_system_prompt,
                    "state": 'stop'}
        else:
            response = f'I will provide the protein information related to {instruction}'
            return {"response": response, "operations":[{"title":"information","options":data_list[:50],"type":["single","input"],"key":"gene_select","value":[]}], "data_list": data_list, "stage": stage,
                    "search_type": "download_type", "search_system_prompt": download_system_prompt,
                    "state": 'continue'}
        
        
            
    else:
        stage+=1
        data_list = instruction['instruction']['gene_select']
        
        # user_select_protien_ncbi_download=""
        aa_dna_content, download_path = user_select_protein_ncbi_download(data_list)
        response = f"Protein gene selection has been completed as needed, the details are as follows: {aa_dna_content} "
        return {"response": response, "operations":[], "protein_gene_seq_path": download_path,
                "search_type": "download_type", "search_system_prompt": download_system_prompt, "stage": stage,
                "state": 'stop'}

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
    
def get_nucleotide_info(name, aa_seq):
    # 搜索nucleotide数据库
    query = f"{name}"
    handle = Entrez.esearch(db="nucleotide", term=query, retmax=100000)
    record = Entrez.read(handle)
    handle.close()
    nucleotide_ids = record["IdList"]
    print('nucleotide_ids', nucleotide_ids)

    # 根据uniport提供的序列和名称来判断，取出真正的DNA序列
    for nucleotide_id in nucleotide_ids:
        handle = Entrez.efetch(db="nucleotide", id=nucleotide_id, retmode="xml")
        records = Entrez.read(handle)
        handle.close()
        flag = False
        for record_info in records[0]["GBSeq_feature-table"]:
            if "CDS" == record_info["GBFeature_key"]:
                product_found = False
                for quals in record_info["GBFeature_quals"]:
                    if quals["GBQualifier_name"] == "translation":
                        pro_seq_d = quals["GBQualifier_value"]
                        if pro_seq_d == aa_seq:
                            loc_split = record_info["GBFeature_location"].split('(')[-1].split(')')[0].split(',')
                            if len(loc_split) == 1:
                                print('GBFeature_location', record_info["GBFeature_location"])
                                loc = [int(x.replace('<', '').replace('>', '')) - 1 for x in
                                       record_info["GBFeature_location"].split('(')[-1].split(')')[0].split('..')]
                                print(loc)
                                dna_seq = records[0]["GBSeq_sequence"][loc[0]:loc[1]]
                            else:
                                dna_seq = ''
                                print(loc_split)
                                loc = [x.split('..') for x in loc_split]
                                print(loc)
                                for loc_ in loc:
                                    print(loc_)
                                    dna_seq += records[0]["GBSeq_sequence"][int(loc_[0]):int(loc_[1])]
                            dna_define = records[0]["GBSeq_definition"]
                            dna_name = records[0]["GBSeq_accession-version"]
                            flag = True
                            break
                        else:
                            dna_define = 'None'
                            dna_seq = 'None'
                            dna_name = 'None'
                            flag = False

                if flag == True:
                    break
            else:
                dna_define = 'None'
                dna_seq = 'None'
                dna_name = 'None'
                flag = False

    return dna_seq, dna_define, dna_name