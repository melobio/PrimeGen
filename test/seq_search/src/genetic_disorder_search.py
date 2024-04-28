
import pandas as pd
import openai
import os



# genetic disorder search tool
def genetic_disorder_search(instruction, stage):
    # description genetic disorder search tool
    genetic_disorder_system_prompt = """
    Current search agents are involved in search projects related to human genetic diseases.First, it requires user interaction to extract the entered human genetic disease name and retrieve the gene names related to the genetic disease on Omim.
And feed back the retrieved information, such as genetic disease phenotypes, related genes and other information.
Further user interaction will select the genetic disease and associated genes required for study. Finally, according to the user’s choice,
The gene sequence of the corresponding gene and the position information on the genome were downloaded from NCBI.
    """

    # genetic disorder name --> gene name list
    if stage == 1:
        response = f'I will provide the name of the genetic disease and the corresponding gene name related to {instruction}'
        stage += 1
        # cancer_name = instruction
        # OMIM dataset
        # 根据 disorder name来抽取gene names
        gene_result = get_disease_gene(instruction)
        if isinstance(gene_result, str):
            return {"response": gene_result, "operations":[{"title":"Sorry, the target gene of the disease was not found in our database. You can choose to upload the gene file yourself or study other diseases.", "options":[], "type":["file"], "key":"gene_path", "value":[]}], "stage": stage,
                    "search_type": "genetic_disorder_type", "search_system_prompt": genetic_disorder_system_prompt,
                    "state": 'continue'}
        else:
            # 返回gene names 、 stage、type
            return {"response": response,"operations":[{"title":"gene and disease information","options":gene_result,"type":["single","input"],"key":"gene_select","value":[]}], "stage": stage,
                    "search_type": "genetic_disorder_type", "search_system_prompt": genetic_disorder_system_prompt,
                    "state": 'continue'}
    # user select --> gene name list --> gene sequence and snp
    elif stage == 2 :
    
        stage += 1
        if 'gene_path' in str(instruction['instruction']):
            gene_path = instruction['instruction']['gene_path']
            return {"response": response, "operations":[],"target_gene_info": gene_path, "search_type": "genetic_disorder_type",
                    "search_system_prompt": cancer_search_system_prompt, "stage": stage,"state": 'stop'}
        
        target_gene_names_list = instruction['instruction']['gene_select']
        # 打开 ClinVar dataset
        gene_info_df = pd.read_csv('/reference_data/gene_info.tsv', sep='\t')
        print('target_gene_names_list', target_gene_names_list)
        # target_gene_dict = {}
        df_list = []
        for gene_name in target_gene_names_list:
            temp_df = gene_info_df[gene_info_df['gene'] == gene_name]
            df_list.append(temp_df)
        # 返回目标基因在染色体上的位置信息
        target_gene_dict = pd.concat(df_list).to_dict()
        response = f'I will provide the location information of the target gene on the whole genome'
        return {"response": response, "operations":[],"target_gene_info": target_gene_dict, "search_type": genetic_disorder_type,
                "search_system_prompt": genetic_disorder_system_prompt, "stage": stage,"state": 'stop'}

def get_disease_gene(disease_name):
    path = '/reference_data/pheno2gene.csv'
    gene_df = pd.read_csv(path)
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
