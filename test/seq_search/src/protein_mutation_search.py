import uvicorn
from Bio import Entrez
import os
import pandas as pd
import requests
import openai
import json
import re

Entrez.email = "912837656@qq.com"
Entrez.api_key = "24e5b06635f7b32dbd9b89c178a3caf3f009"#"db4b31edb00e220c3dd378908403eddbc308"

# protein engineering search tool
def protein_mutation_search(instruction, stage):
    # description protein engineering search tool
    protein_engineering_system_prompt = """The current search agent is involved in search projects related to various enzyme mutant libraries. 
First it requires user interaction to extract the entered protein name and retrieve the relevant sequence on the UniProt. 
And feed back the retrieved sequence information, such as function, species, length and other information. 
Further user interaction will select the desired protein. Finally, based on the user's selection, 
the corresponding protein's gene sequence is downloaded from NCBI. 
This sequence is then provided to the user for use in multiplex PCR methods, facilitating tNGS library construction and sequencing."""

    # protein name --> protein information
    if stage == 1:
        stage += 1
          
        protein_mutation_dict = {
            "Experiment_Direction": "protein mutation",
            "protein_name": "",
            "Organism":"",
           	"Gene_Names":"",
            "Length":"",
        }
        
        experiment_type_prompt = """Based on the user's contextual dialogue, you need to perform information extraction. The information to be extracted includes:  
1.protein_name (If the protein name mentioned by the user in the conversation is an abbreviation, please return its full name.Return null if the protein name was not mentioned in the user conversation)   

The results should be returned to me in JSON format.  
Example is shown below:
'''json
{  
  "protein_name": 'luciferase',   
}
'''
"""
        conversation = instruction['conversation']
        conversation_str = "\n\n".join([i["role"] + ':\n' + i['content'] for i in conversation])
        experiment_info_response = openai.ChatCompletion.create(
            messages=[{"role": "system", "content": experiment_type_prompt},
                      {"role": "user", "content": conversation_str}
                      ],
            deployment_id="XMGI-Chat1-GPT4"
        )
        experiment_info_json_data = extract_json_response(experiment_info_response)
        protein_mutation_dict["protein_name"] = experiment_info_json_data["protein_name"]
        
        protein_name = protein_mutation_dict["protein_name"]
        print('protein_name',protein_mutation_dict["protein_name"])
        response = f'In order to search more accurately, if you know the Organism, Gene Names, Length and other information related to the protein {protein_name}, please provide it to us.'
        return {"response": response, "operations":[],
                    "stage": stage, "search_type": "protein_mutation_type",
                    "protein_mutation_dict": protein_mutation_dict, "state": 'continue'} 

        
        
    elif stage == 2:
    
        stage+=1
        protein_mutation_dict = instruction['instruction']['protein_mutation_dict']
        protein_name = protein_mutation_dict["protein_name"]
        
        experiment_type_prompt = """Based on the user's contextual dialogue, you need to perform information extraction. If the user provides an abbreviation, please restore it to the full name. The information to be extracted includes:  
1.Organism (the Species information or Organism information mentioned by users,Return null if it was not mentioned in the user conversation)
2.Gene_Names (the gene name. Return null if it was not mentioned in the user conversation)
3.Length (the protein length. Return null if it was not mentioned in the user conversation)

The results should be returned to me in JSON format.  
Example is shown below:
'''json
{  
  "Organism": 'Escherichia coli',
  "Gene_Names": 'ssuD', 
  "Length": 384,   
}
'''
"""
        conversation = instruction['conversation']
        conversation_str = "\n\n".join([i["role"] + ':\n' + i['content'] for i in conversation])
        experiment_info_response = openai.ChatCompletion.create(
            messages=[{"role": "system", "content": experiment_type_prompt},
                      {"role": "user", "content": conversation_str}
                      ],
            deployment_id="XMGI-Chat1-GPT4"
        )
        
        
        experiment_info_json_data = extract_json_response(experiment_info_response)
        
        protein_mutation_dict["Organism"] = experiment_info_json_data["Organism"]
        protein_mutation_dict["Gene_Names"] = experiment_info_json_data["Gene_Names"]
        protein_mutation_dict["Length"] = experiment_info_json_data["Length"]
        
        print('experiment_info_json_data',experiment_info_json_data)

        temp_O = protein_mutation_dict["Organism"]
        temp_n = protein_mutation_dict["Gene_Names"]
        temp_l = protein_mutation_dict["Length"]
        
        #protein_name = protein_name+' '+temp_O+' '+temp_n+' '+temp_l
        #protein_name = protein_name.replace('null','')
        print('protein_name',protein_name)
        
        protein_info, data_list = search_uniprot_one(protein_name,temp_O,temp_n)
        print('搜索结果数量', len(data_list))
        
        if len(data_list) == 0:
            response = f'Sorry, the current search results for the protein {protein_name} on uniprot are empty. You can manually upload the plasmid-modified gene file related to the protein {protein_name} or adjust the name and search again.'
            return {"response": response, "operations":[{"title":"upload gene file", "options":[], "type":["file"], "key":"gene_path", "value":[]}], "stage": stage,
                     "search_type": "protein_mutation_type",
                    "protein_mutation_dict": protein_mutation_dict, "state": 'continue'}

        else:
            # print('xxxxxprotein_info',len(protein_info))
            response = f'I will provide the protein information related to {protein_name}'
            return {"response": response, "operations":[{"title":"protein information","options":data_list[:50],"type":["single","input"],"key":"protein_select","value":[]}],
                    "stage": stage, "search_type": "protein_mutation_type",
                    "protein_mutation_dict": protein_mutation_dict, "state": 'continue'}     
    # user select proetin number --> gene sequence --> continue primer design
    elif stage == 3:
        stage+=1
        protein_mutation_dict = instruction['instruction']['protein_mutation_dict']
        
        if 'protein_select' in instruction['instruction']:
        
            data_list = instruction['instruction']['protein_select'] # user select
            #gene_path = instruction['instruction']['gene_path']
            # user_select_protien_ncbi_download=""
            try:
                aa_dna_content, download_path = user_select_protein_ncbi_download(data_list)
            except Exception as e:
                raise Exception('文件下载失败')
            response = f"Protein gene selection has been completed as needed, the details are as follows: {aa_dna_content} "
            
            if download_path == []:
                return {"response": 'Sorry, the gene sequence corresponding to the protein you selected was not found.You can manually upload the plasmid-modified gene file or select a different protein.',
                "operations":[{"title":"upload file", "options":[], "type":["file"], "key":"gene_path", "value":[]}],
                        "protein_gene_seq_path": '',
                        "search_type": "protein_mutation_type", 
                        "protein_mutation_dict": protein_mutation_dict,"stage": stage,
                        "state": 'stop'}
            else:
                response+=" In order to consider the full length of the gene sequence when designing primers, please click the upload button below to upload your plasmid-modified files, and keep the number of file uploads and file downloads equal."
                return {"response": response, "protein_gene_seq_path": download_path,
                "operations":[{"title":"upload file", 
                "options":[], "type":["file"], "key":"gene_path", "value":[]}],
                        "search_type": "protein_mutation_type", "stage": stage,
                        "protein_mutation_dict": protein_mutation_dict, "state": 'continue'}
                        
        elif 'gene_path' in instruction['instruction']:
        
            response = 'Sequence file retrieval completed'
            new_file_list = instruction['instruction']['gene_path']
            return {"response": response, "protein_gene_seq_path": new_file_list, "operations":[],"stage": stage,
                "search_type": "protein_mutation_type",
                "protein_mutation_dict": protein_mutation_dict, "state": 'stop'}
        
        
    else:
        response = ''
        stage+=1
        
        protein_mutation_dict = instruction['instruction']['protein_mutation_dict']
        new_file_list = instruction['instruction']['gene_path']
        return {"response": response, "protein_gene_seq_path": new_file_list, "operations":[],"stage": stage,
                "search_type": "protein_mutation_type",
                "protein_mutation_dict": protein_mutation_dict, "state": 'stop'}

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

def search_uniprot_one(protein_name,Organism,Gene_Names):
    if Organism ==None and Gene_Names == None:
        get_info_url = f"http://rest.uniprot.org/uniprotkb/search?query=%28protein_name%3A{protein_name}%29+AND+%28reviewed%3Atrue%29&format=json"
    
    elif Organism == None and Gene_Names != None:
    
        #get_info_url = f"http://rest.uniprot.org/uniprotkb/search?query={term}%28reviewed%3Atrue%29&format=json"
        get_info_url = f"http://rest.uniprot.org/uniprotkb/search?query=%28protein_name%3A{protein_name}%29+AND+%28gene%3A{Gene_Names}%29&format=json"
        #http://rest.uniprot.org/uniprotkb/search?query=%28protein_name%3Aluciferase%29+AND+%28gene%3AssuD%29&format=json
    elif Organism!= None and Gene_Names == None: 
        
        get_info_url = f"http://rest.uniprot.org/uniprotkb/search?query=%28protein_name%3A{protein_name}%29+AND+%28organism_name%3A%22{Organism}%22%29&format=json"
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
