import openai
import re
import json
import os
from search_app import *
from ncbi_function import *


def extract_json_response(response):
    # Extract the JSON content from the response
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
            print("GPT constructs JSON incorrectly, as shown below:\n",experiment_info_dict)
            experiment_info_json_data = None  # Default value to avoid unbound variable errors
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
    species_file_path, species_info = ncbi_search(term, 'genome', accession_id=accession_id, target_type="target", search_type="species_cds") # type: ignore
    species_list = [] 
    species_list.append((species_file_path,species_info))
    species_path_list = [x[0] for x in species_list]
    species_info_list = [x[1] for x in species_list]
    return species_info_list, species_path_list


def get_coregene_response(species_info_list, related_info_dict, species_path_list, wgs_dict, term):
    species_info_str = "Downloading sequence information is as follow:\n"
    for index, species_info in enumerate(species_info_list):
        if len(species_info_list)>1:
            species_info_str += f'Species {index+1}'
        species_info_str = "**Reference Genome Information**\n"+f'**Accession**:\t{related_info_dict[term][0]["accession"]}\n'+f'**Genus Name**:\t{related_info_dict[term][0].get("genus","")}\n'+f'**Species Nmae**:\t{related_info_dict[term][0].get("species","")}\n'+f'**Organism Name**:\t{related_info_dict[term][0].get("organism_name","")}\n'+f'**Infraspecific Name**:\t{related_info_dict[term][0].get("infraspecific_name","").replace("strain=","")}\n'+f'**Isolate Name**:\t{related_info_dict[term][0].get("isolate","")}\n'+f'**Assembly Level**:\t{related_info_dict[term][0].get("assembly_level","")}\n'+f'**Link**:\thttps://www.ncbi.nlm.nih.gov/search/all/?term={species_info[0]}\n'+f'**File Path**:\t{related_info_dict[term][0].get("ftp_path","")}\t'

    strain_response = "The current sequence search requirements have been fulfilled. Below is the Specific Core Sequence file:"
    wgs_dict["genome_path"] = species_path_list

    response_str = """Next, I will inform the Controller to communicate with you regarding the next steps for primer design."""
    user_response = species_info_str+'\n\n' + response_str + '\n\n' + strain_response
    print('wgs_dict:\n', wgs_dict)
    return user_response, species_path_list, wgs_dict
    

# whole genome sequencing search tool
def whole_genome_sequencing_search(instruction, stage):
    SPECIES_SEARCH_HOST = os.getenv('SPECIES_SEARCH_HOST')
    SPECIES_SEARCH_PORT = os.getenv('SPECIES_SEARCH_PORT')
    base_url = f"http://{SPECIES_SEARCH_HOST}:{SPECIES_SEARCH_PORT}"
    
    if stage == 1:
        print('----------first stage :', instruction)
        print(instruction)
         
        if 'search_type' not in instruction['instruction']:
            conversation = instruction['conversation']
            # Initialized information extraction dictionary for species identification
            wgs_dict = {
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
User: I want to conduct a multiplex PCR species identitifcation experiment for Influenza A virus (A/pika/Qinghai/BI/2007(H5N1)).
Assistant: {"species_information":["Influenza A virus (A/pika/Qinghai/BI/2007(H5N1))"]}
----------
Example2:
User: I want to conduct a multiplex PCR species identitifcation experiment for Acholeplasma laidlawii phage L2. Assistant: {"species_information":["Acholeplasma laidlawii phage L2"]}
----------
Example3:
User: I want to conduct a multiplex PCR species identitifcation experiment for Pacmanvirus str. A23. Assistant: {"species_information":["Pacmanvirus str. A23"]}
----------
Example4:
User: I want to conduct a multiplex PCR species identitifcation experiment for Pseudomonas phage vB_Paer_Ps25. Assistant: {"species_information":["Pseudomonas phage vB_Paer_Ps25"]}
"""

            conversation_str = "\n\n".join([i["role"] + ':\n' + i['content'] for i in conversation])
        
            #JSON MODE extract bio-information
            experiment_info_response = openai.ChatCompletion.create(
                messages=[{"role": "system", "content": experiment_type_prompt},
                          {"role": "user", "content": "##user's contextual dialogue:\n"+"'''\n"+conversation_str+"\n'''"}],
                response_format={"type": "json_object"},
                deployment_id="gpt-4o",
            )
            print('[GPT4 extract bio-information]: ', experiment_info_response) 
            extract_information = json.loads(experiment_info_response["choices"][0]["message"]["content"])["species_information"]
            print("[Json extract bio-information]: ", extract_information)
            
            # [NER] Based on the species information, genus, species, organism, infraspecific, and isolate information are extracted.
            extract_meta_info_list = ner2metainfo(extract_information,base_url,get_ner_entities)
            print("[NER result]：", extract_meta_info_list)
            species_num,target_num_list,related_info_dict = extract_metainfor(extract_meta_info_list,base_url)
            #print("commend information: ",related_info_dict)
            #If NER succeeds, it returns no value info, which means there may be synonyms.
            is_all_empty = all(not value for _, value in related_info_dict.items())
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
                        "data":{"wgs_dict": wgs_dict} 
                        }
                else:
                    response = """We would like to inform you that based on the database search results, the species name remains ambiguous. Please complete the confirmation below, and we will proceed with further experimental analysis once the information has been clarified."""
                    return {"response": response,
                            "operations": operations_list,
                            "stage": stage,
                            "state": 'continue',
                            "search_type": "wgs_type", 
                            "data":{"wgs_dict": wgs_dict} 
                            }

            #Obtaining whole genome data
            else:
                print(related_info_dict)
                term=list(related_info_dict.keys())[0]
                accession_id = related_info_dict[term][0]['accession']
                
                species_info_list, species_path_list = accession2seqs(term, accession_id, [related_info_dict[term][0]['taxid']])
                print('get reference gene：',species_path_list)
                
                user_response, species_path_list, wgs_dict = get_coregene_response(species_info_list, related_info_dict, species_path_list, wgs_dict,term)
                
                wgs_dict["accession"] = [accession_id]
                wgs_dict["taxid"] = [related_info_dict[term][0]['taxid']]
                wgs_dict["genus"] = [related_info_dict[term][0]['genus']]
                wgs_dict["species"] = [related_info_dict[term][0]['species']]
                wgs_dict["organism"] = [related_info_dict[term][0]['organism_name']]
                wgs_dict["infraspecific"] = [related_info_dict[term][0]['infraspecific_name']]
                wgs_dict["isolate"] = [related_info_dict[term][0]['isolate']]
                wgs_dict["assembly_level"] = [related_info_dict[term][0]['assembly_level']]
                wgs_dict["ftp_path"] = [related_info_dict[term][0]['ftp_path']]
                wgs_dict["group"] = [related_info_dict[term][0]['group']]
                wgs_dict["genome_path"] = species_path_list#[related_info_dict[term][0]['core_cds_path_list']]
                wgs_dict["ftp_path"] = [related_info_dict[term][0]['ftp_path']]
                
                user_response+="""\nTo design primers that comprehensively cover all strains and variants of pathogenic viruses, please provide the following SNP information:
- **SNP list**: Include the chromosome number, genomic position (e.g., chr3:21563), reference allele, and variant allele.
- **Additional details** (optional): Any annotations (e.g., functional impact, associated diseases) or frequency data, if available.
- **Preferred format**: CSV or Excel with columns for:
  - Chromosome
  - Genomic position
  - Reference allele
  - Variant allele
  - Frequency (if available)
  - Annotation (if available)\n
You can also choose not to upload the SNP file, and we will use the bioinformatics algorithm to calculate the SNP content.
For algorithmic SNP identification, please provide a mutation frequency threshold to filter significant variants.
- **Frequency threshold**: A minimum mutation frequency (e.g., 0.05) to include SNPs occurring in at least 5% of the population."""
                return {"response": user_response, 
                        "operations":[{"title":"upload SNP file", "options":[], "type":["file"], "key":"snp_path", "value":[]}],
                        "stage": stage+1,
                        "operations":[],
                        "core_cds_path": species_path_list, 
                        "state": 'continue',
                        "search_type": "wgs_type",
                        "data":{"wgs_dict": wgs_dict}
                    }     
        
        #user choice data
        else:
            print(instruction)
            print('user select data')
            operation_result = [value[0] for key, value in instruction['instruction'].items() if key.startswith("species_select_")]
            term = operation_result[0]["organism_name"]
            accession_id = operation_result[0]["accession"]
            tax_id= operation_result[0]["taxid"]
            wgs_dict = instruction['instruction']["wgs_dict"]
            species_info_list,species_path_list = accession2seqs(term, accession_id, [tax_id])
            print('get reference gene：',species_path_list)
            related_info_dict = {term:operation_result}
            print('related_info_dict: ',related_info_dict)
            user_response, species_path_list, wgs_dict = get_coregene_response(species_info_list, related_info_dict, species_path_list, wgs_dict,term)
            wgs_dict["accession"] = [accession_id]
            wgs_dict["taxid"] = [related_info_dict[term][0]['taxid']]
            wgs_dict["genus"] = [related_info_dict[term][0]['genus']]
            wgs_dict["species"] = [related_info_dict[term][0]['species']]
            wgs_dict["organism"] = [related_info_dict[term][0]['organism_name']]
            wgs_dict["infraspecific"] = [related_info_dict[term][0]['infraspecific_name']]
            wgs_dict["isolate"] = [related_info_dict[term][0]['isolate']]
            wgs_dict["assembly_level"] = [related_info_dict[term][0]['assembly_level']]
            wgs_dict["ftp_path"] = [related_info_dict[term][0]['ftp_path']]
            wgs_dict["group"] = [related_info_dict[term][0]['group']]
            wgs_dict["genome_path"] = species_path_list
            wgs_dict["ftp_path"] = [related_info_dict[term][0]['ftp_path']]
            
            user_response+="""\nTo design primers that comprehensively cover all strains and variants of pathogenic viruses, please provide the following SNP information:
- **SNP list**: Include the chromosome number, genomic position (e.g., chr3:21563), reference allele, and variant allele.
- **Additional details** (optional): Any annotations (e.g., functional impact, associated diseases) or frequency data, if available.
- **Preferred format**: CSV or Excel with columns for:
  - Chromosome
  - Genomic position
  - Reference allele
  - Variant allele
  - Frequency (if available)
  - Annotation (if available)\n
You can also choose not to upload the SNP file, and we will use the bioinformatics algorithm to calculate the SNP content.
For algorithmic SNP identification, please provide a mutation frequency threshold to filter significant variants.
- **Frequency threshold**: A minimum mutation frequency (e.g., 0.05) to include SNPs occurring in at least 5% of the population.
"""
            return {"response": user_response, 
                    "operations":[{"title":"upload SNP file", "options":[], "type":["file"], "key":"snp_path", "value":[]}], 
                    "stage": stage+1,
                    "operations":[],
                    "core_cds_path": species_path_list, 
                    "state": 'continue',
                    "search_type": "wgs_type",
                    "data":{"wgs_dict": wgs_dict}
                }
        
    #User upload snp file 
    if stage == 2:
        stage+=1
        wgs_dict = instruction['instruction']['data']['wgs_dict']
        new_file_list = instruction['instruction']['snp_path']
        
        #User chooses to provide SNP file
        if new_file_list != []:
            response = f'Next, we will conduct a whole genome sequencing primer design experiment based on the complete genome and the SNP file provided by the user.'
            return {"response": response,
                "operations": [],
                "stage": stage+1,
                "search_type": "wgs_type",
                "data":{"protein_mutation_dict": wgs_dict,"wgs_seq_path": new_file_list},
                "state": 'continue'}
        
        #The user chooses to provide automated analysis to obtain SNPs
        else:
            print("Under development")

    
    #Save snp file and 
    #Ask the user whether to provide new experimental objects.
    if stage == 3:
        response = f'In addition to the data we discussed above, do you have more experimental goals you want to study?'
            return {"response": response,
                "operations": [],
                "stage": stage+1,
                "search_type": "wgs_type",
                "data":{"protein_mutation_dict": wgs_dict,"wgs_seq_path": new_file_list},
                "state": 'continue'}