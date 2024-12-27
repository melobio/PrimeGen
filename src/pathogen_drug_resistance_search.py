import pandas as pd
import os
import openai
import re
import json

from llm_utils.utils import get_llm_chat_completion
import logging
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

#root_path = '/home/x/x-multiplex-pcr/.data/xsearch/reference_data'
#增加药品名称，增加可选择SNP
def pathogen_drug_resistance_search(instruction,stage):

    pathogen_drug_resistance_prompt = """ 
    Currently, Search Agent has entered pathogen resistance search. 
    It requires interaction with the user to obtain the name of the pathogen species, 
    and by querying the pathogen library, 
    it can find the SNP and sequence information of the corresponding species for 
    subsequent design of primers for drug resistance analysis.
    """
    # 第一轮对话，提取实体，提取出用户提到的物种信息，并与现有耐药库中的物种进行对比。
    if stage == 1:
        stage += 1
        pathogen_drug_resistance_dict = {
            "Experiment_Direction": "pathogen drug resistance",
            "species_name": "",
            "exist": "",
        }

        conversation = instruction['conversation']
        
        #response = f'I will provide resistance information and site information related to species {instruction}'
        #species_name = instruction.lower()
        dfpatho = pd.read_csv('/reference_data/pathodrugdb.csv', header=0, sep=',')
        patholist = list(set(list(dfpatho['Species'])))
        
        experiment_type_prompt = """Based on the user's contextual dialogue, you need to perform information extraction. The information to be extracted includes:  
1.species_name (If the species name mentioned by the user in the conversation is an abbreviation, please return its full name. If the species name mentioned by the user belongs to one of our species lists"""+f'{patholist}'+""", please return it with the species name in the species list.
 Return null if the species name was not mentioned in the user conversation)  

The results should be returned to me in JSON format.  
Example is shown below:
'''json
{  
  "species_name": 'null',  
}
'''
"""
        
        conversation_str = "\n\n".join([i["role"] + ':\n' + i['content'] for i in conversation])
        experiment_info_response = get_llm_chat_completion(
            messages=[{"role": "system", "content": experiment_type_prompt},
                      {"role": "user", "content": conversation_str}
                      ],
            response_format={"type": "json_object"},
        )
        experiment_info_json_data = extract_json_response(experiment_info_response)
        temp_species_name = experiment_info_json_data["species_name"]
        exist_flag = False
        for ss in patholist:
            if temp_species_name.lower() == ss.lower():
                exist_flag = True
        pathogen_drug_resistance_dict["exist"] = exist_flag
        pathogen_drug_resistance_dict["species_name"] = temp_species_name
        logging.info(f'**********species_name**********：\n{pathogen_drug_resistance_dict["species_name"]}  {pathogen_drug_resistance_dict["exist"]}')
        if pathogen_drug_resistance_dict["exist"] and pathogen_drug_resistance_dict["species_name"] != 'null':
        
          
            species_name = pathogen_drug_resistance_dict["species_name"]
            dfresult = dfpatho[dfpatho['Species'].str.lower() == species_name.lower()]
            print('species_name','dfresult',dfresult)
            response = f'Based on the species {species_name} you provided, we found the following information in our resistance database.'
            df_unique = dfresult.drop_duplicates(subset=['Drug'])
  
            gen_list = list(df_unique['Gene'])
            pos_list = list(df_unique['Position'])
            Drug_list = list(df_unique['Drug'])

            result = []
            for ii in range(len(gen_list)):
                temp_dict ={}
                temp_dict['gene'] = gen_list[ii]
                temp_dict['position'] = pos_list[ii]
                temp_dict['Drug'] = Drug_list[ii]
                result.append(temp_dict)
            return {"response": response, "operations":[{"title":"Resistance information","options":result,"type":["single","input"],"key":"drug_select","value":[]}],"stage": stage,
                            "state": 'continue', "search_type": "pathogen_drug_resistance_type",
                            "data":{"pathogen_drug_resistance_dict": pathogen_drug_resistance_dict}}
                            
        elif not pathogen_drug_resistance_dict["exist"] and pathogen_drug_resistance_dict["species_name"] != 'null':
         
            response = f"Sorry, the species name you provided is not in our database. We currently only support resistance analysis of the following species {patholist}. You can choose one of the species for research or you can choose to upload the resistance information of the species you are interested in ( Must contain gene sequence, SNP site information)."
       
            return {"response": response, "operations":[],"stage": stage-1,
                            "state": 'continue', "search_type": "",
                            "data":{"pathogen_drug_resistance_dict": pathogen_drug_resistance_dict}}
        else:
            response = f"Please provide the name of the species you want to study so that we can find its corresponding resistance information for you." 
            return {"response": response, "operations":[],"stage": stage-1,
                            "state": 'continue', "search_type": '',
                            "data":{"pathogen_drug_resistance_dict": pathogen_drug_resistance_dict}}
             
    elif stage == 2:
        stage+=1
        pathogen_drug_resistance_dict = instruction['instruction']['data']['pathogen_drug_resistance_dict']
        if pathogen_drug_resistance_dict["exist"] and pathogen_drug_resistance_dict["species_name"] != 'null':
            drug_select = instruction['instruction']['drug_select']
            species_name = pathogen_drug_resistance_dict['species_name']
            
            #drug_name = instruction[:-1]
            response = 'Based on your final choice, I give the following DNA sequence and snp information'
            dfpatho = pd.read_csv('/reference_data/pathodrugdb.csv', header=0, sep=',')
            dfresult = dfpatho[dfpatho['Species'].str.lower() == species_name.lower()]
            print('drug_select',drug_select)
            dfresult = dfresult[dfresult['Drug'].isin(drug_select)]
            print('xxxx',dfresult)
            save_path = f'/reference_data/{species_name}_drug_filter_result.csv'
            dfresult.to_csv(save_path, index=False)
            if species_name.lower() == 'mycobacterium tuberculosis':
                
                MTBseq = '/reference_data/MTB.fa'
                
                return {"response": response,"operations":[],
                        "data":{"filter_result_path":save_path,'sequences':MTBseq},
                    "search_type": "pathogen_drug_resistance_type", "stage":stage,"search_system_prompt": pathogen_drug_resistance_prompt,
                    "state": 'stop'}
                    
            else:
                return {"response": response,"operations":[],
                        "data":{"filter_result_path":save_path,'sequences':''},"search_type": "pathogen_drug_resistance_type",
                    "stage":stage,"search_system_prompt": pathogen_drug_resistance_prompt,"state": 'stop'}
        else:
        
            pass


def extract_json_response(response):
    # 从response中提取JSON内容
    experiment_info_dict = response.choices[0].message.content
    try:
        # if json.loads(experiment_info_dict):
        return json.loads(experiment_info_dict)
    except:
        match = re.search(r'```json\n(.*?)\n```', experiment_info_dict, re.DOTALL)
        if match:
            experiment_info_json_str = match.group(1)
            experiment_info_json_data = json.loads(experiment_info_json_str)
        else:
            print("GPT构造JSON错误，情况如下所示：\n", experiment_info_dict)
            exit()
        return experiment_info_json_data