import openai
import re
import json
import os
from search_app import *
from ncbi_function import *

from llm_utils.utils import get_llm_chat_completion


def extract_json_response(response):
    # 从response中提取JSON内容
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
            print("GPT构造JSON错误，情况如下所示：\n",experiment_info_dict)
            experiment_info_json_data = None  # 默认值，避免未绑定变量错误
        return experiment_info_json_data


def ner2metainfo(extract_information,base_url,get_ner_entities):
    # [NER] Based on the species information, genus, species, organism, infraspecific, and isolate information are extracted.
    extract_meta_info_list=[]
    for extract_info in extract_information:
        extract_meta_info_dict={"genus":'', "species":'',"organism":'',"infraspecific":'',"isolate":''}
        ner_result = get_ner_entities(base_url, extract_info)
        print("[NER information]: ", ner_result) 
        

        for ner_ in ner_result:
            label = ner_.get('label').split(',')
            if len(label) != 1:
                if 'infraspecific' in label:
                    extract_meta_info_dict['infraspecific']=ner_.get('text')
                elif 'isolate' in label:
                    extract_meta_info_dict['isolate']=ner_.get('text')
                elif 'genues' in label:
                    extract_meta_info_dict['genus']=ner_.get('text')
                elif 'species' in label:
                    extract_meta_info_dict['species']=ner_.get('text')
                else:
                    extract_meta_info_dict['organism']=ner_.get('text')
            else:
                extract_meta_info_dict[ner_.get('label')]=ner_.get('text')

        extract_meta_info_list.append(extract_meta_info_dict)
    print("[current Meta information]: ",extract_meta_info_list)
    return extract_meta_info_list

def extract_metainfor(extract_meta_info_list,base_url):
    # [Meta data] Extract meta information based on NER data
    title_list = ['accession','representative', 'taxid', 'genus', 'species', 'organism_name', 'infraspecific_name',  'isolate', 'assembly_level', 'ftp_path', 'group', 'origin','contig_n50','scaffold_n50','total_sequence_length']
    related_info_dict = {}
    for meta_info_dict in extract_meta_info_list:
        related_info_list = []
        related_info = get_meta_info(base_url, [meta_info_dict])[0]["text"][0]
        for single_info in related_info:
            info_dict={}
            single_info = single_info[:12] + single_info[-3:]
            for x,y in zip(title_list,single_info):
                info_dict[x]=y
            related_info_list.append(info_dict)
        related_info_dict[meta_info_dict["organism"]] = related_info_list

    
    species_num = len(related_info_dict)
    target_num_list = [len(y) for _,y in related_info_dict.items()]
    return species_num,target_num_list,related_info_dict

def accession2seqs(term, accession_id, taxid_list):
    #s_n, 'genome', t
    species_file_path, species_info = ncbi_search(term, 'genome', accession_id=accession_id, target_type="target", search_type="species_cds")
    species_list = [] 
    species_list.append((species_file_path,species_info))
    species_path_list = [x[0] for x in species_list]
    species_info_list = [x[1] for x in species_list]
    #taxid_list = [related_info_dict[term][0]['taxid']]      
    core_cds_path_list=[]
    for species_path,taxid  in zip(species_path_list,taxid_list):
        print(f'开始{species_path}')
        core_cds_path_list.append(get_core_gene(file_path=[species_path],taxid=[taxid], type_cds="have_cds"))
    return species_info_list, core_cds_path_list


def get_coregene_response(species_info_list, related_info_dict, core_cds_path_list, species_identification_dict, term):
    species_info_str = "Downloading sequence information is as follow:\n"
    for index, species_info in enumerate(species_info_list):
        if len(species_info_list)>1:
            species_info_str += f'Species {index+1}'
        species_info_str = "**Reference Genome Information**\n"+f'**Accession**:\t{related_info_dict[term][0]["accession"]}\n'+f'**Genus Name**:\t{related_info_dict[term][0].get("genus","")}\n'+f'**Species Nmae**:\t{related_info_dict[term][0].get("species","")}\n'+f'**Organism Name**:\t{related_info_dict[term][0].get("organism_name","")}\n'+f'**Infraspecific Name**:\t{related_info_dict[term][0].get("infraspecific_name","").replace("strain=","")}\n'+f'**Isolate Name**:\t{related_info_dict[term][0].get("isolate","")}\n'+f'**Assembly Level**:\t{related_info_dict[term][0].get("assembly_level","")}\n'+f'**Link**:\thttps://www.ncbi.nlm.nih.gov/search/all/?term={species_info[0]}\n'+f'**File Path**:\t{related_info_dict[term][0].get("ftp_path","")}\t'

    strain_response = "The current sequence search requirements have been fulfilled. Below is the Specific Core Sequence file:"
    species_identification_dict["target_number"] = 'multi'
    species_identification_dict["target_cds_path"] = core_cds_path_list

    response_str = """Next, I will inform the Controller to communicate with you regarding the next steps for primer design."""
    user_response = species_info_str+'\n\n' + response_str + '\n\n' + strain_response
    print('species_identification_dict:\n', species_identification_dict)
    return user_response, core_cds_path_list, species_identification_dict
    

# species identification search tool
def species_identification_search(instruction, stage):
    SPECIES_SEARCH_HOST = os.getenv('SPECIES_SEARCH_HOST')
    SPECIES_SEARCH_PORT = os.getenv('SPECIES_SEARCH_PORT')
    base_url = f"http://{SPECIES_SEARCH_HOST}:{SPECIES_SEARCH_PORT}"
    
    if stage == 1:
        print('----------first stage :', instruction)
        print(instruction)
         
        if 'search_type' not in instruction['instruction']:
            conversation = instruction['conversation']
            # Initialized information extraction dictionary for species identification
            species_identification_dict = {
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

            # Try to extract information in the user's language: experiment type, species name, strain name
            experiment_type_prompt = """According to the user's contextual dialogue, you need to extract and summary the information needed for biological sequence retrieval.
You just need to summarize the extracted information and reply to me with a complete biological information description. You can think of it as rewriting a clear query. The result needs to be returned in json format.

Note:
1) The biological information description must be complete and professional. If it is a variant, strain, isolate, etc., just extract it directly. Remember not to use brackets or other additional descriptions.
2) Attention! If the user does not mention it, please reply null, do not guess the user's intention!
3) Strain must be abbreviated as "str.", substrain must be abbreviated as "substr.", subspecies must be abbreviated as "subsp.", and varieties must be abbreviated as "var.".
4) When extracting user information, do not modify punctuation marks, especially periods, brackets, dashes, underscores, qual sign, space etc.

Example1:
User:Species identification of Salmonella enterica subspecies enterica serovar Uganda varieties 15+
Assistant: {"species_information":["Salmonella enterica subspecies enterica serovar Uganda varieties 15+"]}
----------
Example2:
User:Species identification of Salmonella enterica subspecies enterica serovar Uganda varieties 15+ and Francisella tularensis subsp. novicida.
Assistant: {"species_information":["Salmonella enterica subspecies enterica serovar Uganda varieties 15+","Francisella tularensis subsp. novicida"]}
----------
Example3:
User: I want to conduct a multiplex PCR species identitifcation experiment for Influenza A virus (A/pika/Qinghai/BI/2007(H5N1)).
Assistant: {"species_information":["Influenza A virus (A/pika/Qinghai/BI/2007(H5N1))"]}
----------
User: I want to conduct a multiplex PCR species identitifcation experiment for Helicobacter pylori NCTC 11637 = CCUG 17874 = ATCC 43504 = JCM 12093.
Assistant: {"species_information":["Helicobacter pylori NCTC 11637 = CCUG 17874 = ATCC 43504 = JCM 12093"]}
"""

            conversation_str = "\n\n".join([i["role"] + ':\n' + i['content'] for i in conversation])
        
            #JSON MODE extract bio-information
            experiment_info_response = get_llm_chat_completion(
                messages=[
                    {"role": "system", "content": experiment_type_prompt},
                    {"role": "user", "content": "##user's contextual dialogue:\n"+"'''\n"+conversation_str+"\n'''"}
                ],
                response_format={"type": "json_object"},
            )
            print('[GPT4 extract bio-information]: ', experiment_info_response) 
            extract_information = json.loads(experiment_info_response.choices[0].message.content or '{}')["species_information"]
            print("[Json extract bio-information]: ", extract_information)

            # [NER] Based on the species information, genus, species, organism, infraspecific, and isolate information are extracted.
            extract_meta_info_list = ner2metainfo(extract_information,base_url,get_ner_entities)
            print("[NER result]：", extract_meta_info_list)
            species_num,target_num_list,related_info_dict = extract_metainfor(extract_meta_info_list,base_url)
            #print("commend information: ",related_info_dict)
            #If NER succeeds, it returns no value info, which means there may be synonyms.
            is_all_empty = all(not value for key, value in related_info_dict.items())
            if is_all_empty:
                print('empty information!!!!start synonyms ')
                synonym_name = [k for k,_ in related_info_dict.items()]
                synonym_result = check_synonym(base_url, synonym_name)
                print('[synonym results]:', synonym_result)
                extract_synonym_information = synonym_result[0]['text']
                extract_synonym_meta_info_list = ner2metainfo(extract_synonym_information,base_url,get_ner_entities)
                species_num,target_num_list,related_info_dict = extract_metainfor(extract_synonym_meta_info_list,base_url)
                #print("commend information: ",related_info_dict)
                print('target_num_list: ',target_num_list)
                print('species_num: ', species_num)
            else:
                # [Meta data] Extract meta information based on NER data
                species_num,target_num_list,related_info_dict = extract_metainfor(extract_meta_info_list,base_url)
                #print("commend information: ",related_info_dict)
                print('target_num_list: ',target_num_list)
                print('species_num: ', species_num)
                
            # A large amount of information below the species level is obtained and the user confirms it.
            operations_list=[]
            if sum(target_num_list) > species_num:
                for i,(n,(x,y)) in enumerate(zip(target_num_list, related_info_dict.items())):
                    if n>1:
                        commend_y = commend_seqs_filte(y)
                        if commend_y==[]:
                            continue
                        operations_list.append({"widget":"table",
                                                "title":f"The {x} is unclear, which may not be accurate. Please confirm again.(required)", 
                                                "options":commend_y, 
                                                "type":["single","required"],
                                                "key":f"species_select_{i+1}", 
                                                "value":[]})
                if operations_list==[]:
                    response = f"Sorry, the current {extract_information[0]} sequence does not have **complete genome** data and cannot be used to design primers. Please provide new biological sequence information."
                    return {"response": response,
                        "operations": [],
                        "stage": stage,
                        "state": 'continue',
                        #"search_type": "controller_type", 
                        "data":{"species_identification_dict": species_identification_dict} 
                        }
                else:
                    response = """We would like to inform you that based on the database search results, the species name remains ambiguous. Please complete the confirmation below, and we will proceed with further experimental analysis once the information has been clarified."""
                    return {"response": response,
                            "operations": operations_list,
                            "stage": stage,
                            "state": 'continue',
                            "search_type": "species_identification_type", 
                            "data":{"species_identification_dict": species_identification_dict} 
                            }

            #获得唯一数据，可直接accession下载序列
            else:
                print(related_info_dict)
                term=list(related_info_dict.keys())[0]
                accession_id = related_info_dict[term][0]['accession']
                
                species_info_list,core_cds_path_list = accession2seqs(term, accession_id, [related_info_dict[term][0]['taxid']])
                print('完成 core cds：',core_cds_path_list)
                user_response, core_cds_path_list, species_identification_dict = get_coregene_response(species_info_list, related_info_dict, core_cds_path_list, species_identification_dict,term)
                
                species_identification_dict["accession"] = [accession_id]
                species_identification_dict["taxid"] = [related_info_dict[term][0]['taxid']]
                species_identification_dict["genus"] = [related_info_dict[term][0]['genus']]
                species_identification_dict["species"] = [related_info_dict[term][0]['species']]
                species_identification_dict["organism"] = [related_info_dict[term][0]['organism_name']]
                species_identification_dict["infraspecific"] = [related_info_dict[term][0]['infraspecific_name']]
                species_identification_dict["isolate"] = [related_info_dict[term][0]['isolate']]
                species_identification_dict["assembly_level"] = [related_info_dict[term][0]['assembly_level']]
                species_identification_dict["ftp_path"] = [related_info_dict[term][0]['ftp_path']]
                species_identification_dict["group"] = [related_info_dict[term][0]['group']]
                species_identification_dict["target_cds_path"] = core_cds_path_list#[related_info_dict[term][0]['core_cds_path_list']]
                species_identification_dict["ftp_path"] = [related_info_dict[term][0]['ftp_path']]

                return {"response": user_response, 
                        "stage": stage+2,
                        "operations":[],
                        "core_cds_path": core_cds_path_list, 
                        "state": 'stop',
                        "search_type": "species_identification_type",
                        "data":{"species_identification_dict": species_identification_dict}
                    }     
        
        #user choice data
        else:
            print(instruction)
            print('user select data')
            operation_result = [value[0] for key, value in instruction['instruction'].items() if key.startswith("species_select_")]
            term = operation_result[0]["organism_name"]
            accession_id = operation_result[0]["accession"]
            tax_id= operation_result[0]["taxid"]
            species_identification_dict = instruction['instruction']["species_identification_dict"]
            species_info_list,core_cds_path_list = accession2seqs(term, accession_id, [tax_id])
            print('完成 core cds：',core_cds_path_list)
            related_info_dict = {term:operation_result}
            print('related_info_dict: ',related_info_dict)
            user_response, core_cds_path_list, species_identification_dict = get_coregene_response(species_info_list, related_info_dict, core_cds_path_list, species_identification_dict,term)
            species_identification_dict["accession"] = [accession_id]
            species_identification_dict["taxid"] = [related_info_dict[term][0]['taxid']]
            species_identification_dict["genus"] = [related_info_dict[term][0]['genus']]
            species_identification_dict["species"] = [related_info_dict[term][0]['species']]
            species_identification_dict["organism"] = [related_info_dict[term][0]['organism_name']]
            species_identification_dict["infraspecific"] = [related_info_dict[term][0]['infraspecific_name']]
            species_identification_dict["isolate"] = [related_info_dict[term][0]['isolate']]
            species_identification_dict["assembly_level"] = [related_info_dict[term][0]['assembly_level']]
            species_identification_dict["ftp_path"] = [related_info_dict[term][0]['ftp_path']]
            species_identification_dict["group"] = [related_info_dict[term][0]['group']]
            species_identification_dict["target_cds_path"] = core_cds_path_list
            species_identification_dict["ftp_path"] = [related_info_dict[term][0]['ftp_path']]

            return {"response": user_response, 
                    "stage": stage+2,
                    "operations":[],
                    "core_cds_path": core_cds_path_list, 
                    "state": 'stop',
                    "search_type": "species_identification_type",
                    "data":{"species_identification_dict": species_identification_dict}
                }