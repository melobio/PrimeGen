
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
    elif stage == 2:

        stage += 1
        response = f'The step of determining the DNA sequence is complete.'
        if 'gene_path' in str(instruction['instruction']):
            gene_path = instruction['instruction']['gene_path']
            return {"response": response, "operations":[],"target_gene_info": gene_path, "search_type": "genetic_disorder_type",
                    "search_system_prompt": genetic_disorder_system_prompt, "stage": stage,"state": 'stop'}
        
        target_gene_names_list = instruction['instruction']['gene_select']
        # 打开 ClinVar dataset
        gene_info_df = pd.read_csv('/reference_data/gene_info.tsv', sep='\t')

        df_list = []
        for gene_name in target_gene_names_list:
            temp_df = gene_info_df[gene_info_df['gene'] == gene_name]
            df_list.append(temp_df)

        res = pd.concat(df_list)
        res.drop_duplicates()
        gene_names = list(res['gene'])
        end_list = list(res['end'])
        start_list = list(res['start'])
        chr_list = list(res['chrom'])

        fasta_file = "/reference_data/GCF_000001405.40_GRCh38.p14_genomic_filter.fna"
        dna_file_path = '/reference_data/dna_file'
        result = []
        len_list =[]
        # 打开文件并读入内容
        with open(fasta_file, 'r') as file:
            file_content = file.read().splitlines()
            file.close()
        for ii in range(len(gene_names)):
            na = gene_names[ii]
            temp_seq = get_dna_seq(chr_list[ii][3:], start_list[ii], end_list[ii], file_content)
            len_list.append(len(temp_seq))
            temp_path = os.path.join(dna_file_path, f'{na}.fasta')
            with open(temp_path, 'w') as f:
                f.write(f'>{na}\n')
                f.write(temp_seq)
                f.close()
            result.append(temp_path)

        response = f'I will provide the gene files,The length of the gene sequence ranges from 1 to {len_list[0]}. Please further provide the target region for the designed primers, and we will design primers for this region.'
        return {"response": response, "operations":[{"title":"start","options":[],"type":["input"],"key":"start_pos","value":[]},
                                                {"title": "end","options": [],"type": ["input"],"key": "end_pos","value": []}],
                "target_gene_path": result,"max_length": len_list[0], "search_type": 'genetic_disorder_type',
                "search_system_prompt": genetic_disorder_system_prompt, "stage": stage,"state": 'continue'}
    else:
        result = instruction['instruction']['target_gene_path']
        start_pos = instruction['instruction']['start_pos']
        end_pos = instruction['instruction']['end_pos']
        response = f'The step of determining the DNA sequence is complete.'
        return {"response": response, "operations": [], "start_end_pos": [start_pos,end_pos],"gene_file_path":result,
                "search_type": 'genetic_disorder_type',
                "search_system_prompt": genetic_disorder_system_prompt, "stage": stage, "state": 'stop'}
def get_dna_seq(chr_num, start, end,file_content):
    chr_name = f'chromosome {chr_num}'
    # fasta文件中的染色体名称
    all_chr_name = 'Homo sapiens ' + chr_name + ', GRCh38.p14 Primary Assembly'
    # 用于存储提取的内容
    extracted_content = ""
    found_chr = False
    # 在内存中搜索包含"all_chr_name"的行并提取内容
    for i, line in enumerate(file_content):
        if not found_chr and all_chr_name in line:
            found_chr = True
        elif found_chr:
            if line.startswith('>'):
                # 如果是新的标题行，则停止提取
                break
            extracted_content += line
    temp_seq = extracted_content[int(start):int(end)].upper()
    return temp_seq
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
    if len(disease_gene_list) == 0:
        return f'sorry, No causative gene related to {disease_name} has been found,You can select other genetic diseases for analysis.'
    else:
        return disease_gene_list
