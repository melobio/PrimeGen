
import pandas as pd
import os

import requests
import openai
import json
import re




def snp_primer_design(instruction,stage):
    if stage ==1:
        
        stage+=1
        response = '''The next critical step is to "determine primer design parameters". This step involves defining optimal conditions for primer design, which includes specifying maximum and minimum amplicon lengths, desired annealing temperatures, and GC content ranges. Additionally, it is important to determine the maximum length of the primers. These parameters will form the basis for designing specific primers. This step is critical to ensure the specificity and efficiency of subsequent amplification and detection processes, helping you move toward your experimental goals. You can refer to the following method to give parameter values: {"max_amp_len": 300,
   "min_amp_len":250,
  "temperature":65,
   "min_GC":40,
   "max_GC":60,
   "max_len":23, }, these are the default values of each parameter.'''
        return {'response':response,"operations":[], 'state': 'continue', "primer_type": "snp_primer_design_type",
                     "stage":stage}
     
    elif stage ==2:
    
        stage+=1
        primer_design_dict = {
                "Experiment_Direction": "primer design",
                "search_type":"",
                "max_amp_len": "",
                "min_amp_len":"",
               	"temperature":"",
                "min_GC":"",
                "max_GC":"",
                "max_len":"",
            }
            
        experiment_type_prompt = 'Based on the user`s contextual dialogue, you need to perform information extraction. The information to be extracted includes:  1.max_amp_len (the Maximum amplicon length, The default value is 300)   2.min_amp_len (the Minimum amplicon length, The default value is 250)   3.temperature (the Annealing temperature, The default value is 65)   4.min_GC (the Minimum GC content of a primer, The default value is 40)   5.max_GC (the Maximum GC content of a primer, The default value is 60)   6.max_len (the Maximum length of a primer , The default value is 23)   Notice: If the user wants to take default parameters, each field returns the default value, otherwise it needs to be extracted from the user reply.Please output a string of a JSON object. Example is shown below:{"max_amp_len": 300,"min_amp_len":250,"temperature":65,"min_GC":40,"max_GC":60,"max_len":23}'
        conversation = instruction['conversation'][-4:]
        #print('instruction',instruction)
        search_type = instruction['instruction']['search_responses']['search_type']
        
        conversation_str = "\n\n".join([i["role"] + ':\n' + i['content'] for i in conversation])
        experiment_info_response = openai.ChatCompletion.create(
            messages=[{"role": "system", "content": experiment_type_prompt},
                      {"role": "user", "content": conversation_str}
                      ],
            deployment_id="XMGI-Chat1-GPT4"
        )
        #experiment_info_json_data = extract_json_response(experiment_info_response)
        experiment_info_dict = experiment_info_response["choices"][0]["message"]['content']
        print('experiment_info_dict',experiment_info_dict)
        experiment_info_json_data = json.loads(experiment_info_dict)
        primer_design_dict["max_amp_len"] = experiment_info_json_data["max_amp_len"]
        primer_design_dict["min_amp_len"] = experiment_info_json_data["min_amp_len"]
        primer_design_dict["temperature"] = experiment_info_json_data["temperature"]
        primer_design_dict["min_GC"] = experiment_info_json_data["min_GC"]
        primer_design_dict["max_GC"] = experiment_info_json_data["max_GC"]
        primer_design_dict["max_len"] = experiment_info_json_data["max_len"]
        primer_design_dict["search_type"] = search_type
        print('primer_design_dict',primer_design_dict)
        response = f'Please confirm whether the following parameter values and experiment categories are correct. If correct, please reply yes, otherwise reply no and provide correct parameter values. primer_design_dict {primer_design_dict}'
        return {'response':response,"operations":[], 'state': 'continue', "primer_type": "snp_primer_design_type",
                     'primer_design_dict': primer_design_dict,"stage":stage}
        
    elif stage ==3:
        print('进入stage 3')
        stage+=1
        primer_design_dict = instruction['instruction']['primer_design_dict']
        print('primer_design_dict',primer_design_dict)
        conversation = instruction['conversation'][-6:]
        experiment_type_prompt = """Based on the user's contextual dialogue,  you need to perform information extraction. The information to be extracted includes:
1.flag (If the user replies yes, it returns TRUE, if it is no, it returns FALSE.)
The results should be returned to me in JSON format.  
Example is shown below:
'''json
{  
  "flag": "TRUE",
  
}
''' 
"""
        conversation_str = "\n\n".join([i["role"] + ':\n' + i['content'] for i in conversation])
        experiment_info_response = openai.ChatCompletion.create(
            messages=[{"role": "system", "content": experiment_type_prompt},
                      {"role": "user", "content": conversation_str}
                      ],
            deployment_id="XMGI-Chat1-GPT4"
        )
        experiment_info_json_data = extract_json_response(experiment_info_response)
        temp_flag = experiment_info_json_data['flag']
        print('temp_flag',temp_flag)
        if temp_flag =='FALSE':
        
            response = ''
            return {'response':response,"operations":[], 'state': 'continue', "primer_type": "snp_primer_design_type",
                     'primer_design_dict': primer_design_dict,"stage":stage-2}
        
        
        parameters = [primer_design_dict["max_amp_len"],primer_design_dict["min_amp_len"],primer_design_dict["temperature"],
        primer_design_dict["min_GC"],primer_design_dict["max_GC"],primer_design_dict["max_len"]]
        
        
        snp_primer_design_prompt = """
here is the designed primer file.
                    """
        out_path = '/reference_data/primer_result'
        dna_file_path = '/reference_data/dna_file'
        snp_file_path = '/reference_data/snp_file'
        # if search_type == "genetic_disorder_type":
        search_type = primer_design_dict["search_type"]
    
        print('search_type',search_type)

        # #获取引物设计所必需的参数
        # parameters = instruction['parameters']
        # parameters = []
        if search_type == 'genetic_disorder_type':
            print(f'当前结果来自：{search_type}')
            temp_dict = instruction['instruction']['search_responses']['target_gene_info']
    
            res = pd.DataFrame.from_dict(temp_dict)
    
            print('res',res)
            gene_names = list(res['gene'])
            end_list = list(res['end'])
            start_list = list(res['start'])
            chr_list = list(res['chrom'])
            # 定义文件路径
            fasta_file = "/reference_data/GCF_000001405.40_GRCh38.p14_genomic_filter.fna"
            result = []
            # 打开文件并读入内容
            with open(fasta_file, 'r') as file:
                file_content = file.read().splitlines()
                file.close()
            for ii in range(len(gene_names)):
                na = gene_names[ii]
                temp_seq = get_dna_seq(chr_list[ii][3:], start_list[ii], end_list[ii],file_content)
                temp_path = os.path.join(dna_file_path,f'{na}.fasta')
                with open(temp_path,'w') as f:
                    f.write(f'>{na}\n')
                    f.write(temp_seq)
                    f.close()
    
                fasta_path = temp_path
                snp_path = ''
                if parameters ==[]:
                    cmd = 'python /reference_data/plasmid.py --fasta_path {}  --out_path {} --gene_name {}'.format(fasta_path, out_path,na)
                else:
                    cmd = 'python /reference_data/plasmid.py --fasta_path {}  --out_path {} --gene_name {} --max_amp_len {} --min_amp_len {} --temperature {} --min_GC {} --max_GC {} --max_len {}'.format(
                        fasta_path, out_path, na,parameters[0],parameters[1],parameters[2],parameters[3],parameters[4],parameters[5])
                print('-' * 30, 'cmd: \n', cmd)
                process = os.popen(cmd)  # return file
                output = process.read()
                print(output)
                process.close()
                for ff in os.listdir(out_path):
                    if ff.endswith('.csv') and na in ff:
                        result.append(os.path.join(out_path,ff))
    
    
        elif search_type == 'protein_mutation_type':
            print(f'当前结果来自：{search_type} ')
            path_list = instruction['instruction']['search_responses']['protein_gene_seq_path']
            #path_list = instruction['instruction']['upload_file']
            #if len(path_list_select) != len(path_list):
                #print('上传文件数量与原文件数量不相等。')
                #print('path_list_select',path_list_select)
                #print('path_list', path_list)
            name_list = []
            for pa in path_list:
                temp_name = os.path.basename(pa)
                temp_name = temp_name.split('.')[0]
                temp_name = temp_name.replace('(','').replace(')','')
                name_list.append(temp_name)
                # cmd = 'python /reference_data/plasmid.py --fasta_path {}  --out_path {} --gene_name {}'.format(pa,
                #                                                                                               out_path, temp_name)
                if parameters ==[]:
                    cmd = 'python /reference_data/plasmid.py --fasta_path {}  --out_path {} --gene_name {}'.format(pa, out_path,temp_name)
                else:
                    cmd = 'python /reference_data/plasmid.py --fasta_path {}  --out_path {} --gene_name {} --max_amp_len {} --min_amp_len {} --temperature {} --min_GC {} --max_GC {} --max_len {}'.format(
                        pa, out_path, temp_name,parameters[0],parameters[1],parameters[2],parameters[3],parameters[4],parameters[5])
                print('-' * 30, 'cmd: \n', cmd)
                process = os.popen(cmd)  # return file
                output = process.read()
                print(output)
                process.close()
            result =[]
            for temp_name in name_list:
                for file in os.listdir(out_path):
                    if temp_name in file and '.csv' in file:
                        result.append(os.path.join(out_path,file))
            
        elif search_type == 'pathogen_drug_resistance_type':
            print(f'当前结果来自：{search_type} ')
            res_path = instruction['instruction']['search_responses']['filter_result_path']
            res_df = pd.read_csv(res_path)
    
            gene_list = list(set(list(res_df['Gene'])))
            print('gene_list',gene_list)
            for gen in gene_list:
                temp = res_df[res_df['Gene']==gen]
                snp_list = list(temp['Position'])
                seq_temp = list(temp['Sequence'])[0]
                snp_path = os.path.join(snp_file_path,gen+'_var.csv')
                fa_path = os.path.join(dna_file_path,gen+'_seq.fasta')
                freq = [1 for x in range(len(snp_list))]
                snp_df = pd.DataFrame({'START':snp_list,'STOP':snp_list,'FREQ':freq})
                snp_df = snp_df.drop_duplicates()
                snp_df.to_csv(snp_path,index=False)
                with open(fa_path,'w') as f:
                    f.write('>'+gen+'\n')
                    f.write(seq_temp)
                    f.close()
                # cmd = 'python /reference_data/plasmid.py --fasta_path {}  --snp_path {} --out_path {} --gene_name {}'.format(fa_path,snp_path,out_path,gen)
                if parameters ==[]:
                    cmd = 'python /reference_data/plasmid.py --fasta_path {} --snp_path {} --out_path {} --gene_name {}'.format(fa_path,snp_path, out_path,gen)
                else:
                    cmd = 'python /reference_data/plasmid.py --fasta_path {} --snp_path {} --out_path {} --gene_name {} --max_amp_len {} --min_amp_len {} --temperature {} --min_GC {} --max_GC {} --max_len {}'.format(
                        fa_path,snp_path, out_path,gen,parameters[0],parameters[1],parameters[2],parameters[3],parameters[4],parameters[5])
                print('-' * 30, 'cmd: \n', cmd)
                process = os.popen(cmd)  # return file
                output = process.read()
                print(output)
                process.close()
                time.sleep(3)
            result = []
            for temp_name in gene_list:
                for file in os.listdir(out_path):
                    if temp_name in file and '.csv' in file and 'primer_design' in file:
                        result.append(os.path.join(out_path, file))
    
        elif search_type == 'cancer_type':
            print(f'当前结果来自：{search_type}')
            res_dict = instruction['instruction']['search_responses']['target_gene_info']
            res_df = pd.DataFrame.from_dict(res_dict)
            gene_list = list(set(list(res_df['Gene name'])))
            print('gene_list', gene_list)
            # 定义文件路径
            fasta_file = "/reference_data/GCF_000001405.40_GRCh38.p14_genomic_filter.fna"
    
            # 打开文件并读入内容
            with open(fasta_file, 'r') as file:
                file_content = file.read().splitlines()
                file.close()
            for gen in gene_list:
                temp_df = res_df[res_df['Gene name']==gen]
                snp_path = os.path.join(snp_file_path, gen + '_var.csv')
                snp_list = list(temp_df['Start position in gene'])
                freq = [1 for x in range(len(snp_list))]
                snp_df = pd.DataFrame({'START': snp_list, 'STOP': snp_list, 'FREQ': freq})
                snp_df = snp_df.drop_duplicates()
                snp_df.to_csv(snp_path, index=False)
    
                temp_seq = get_dna_seq(list(temp_df['Chromosome'])[0], list(temp_df['Gene start'])[0], list(temp_df['Gene end'])[0], file_content)
                temp_path = os.path.join(dna_file_path,gen+'_seq.fa')
                with open(temp_path, 'w') as f:
                    f.write(f'>{gen}\n')
                    f.write(temp_seq)
                    f.close()
    
                fasta_path = temp_path
                print(f'{gen}基因序列提取完成。开始设计引物')
                # snp_path = ''
                # cmd = 'python /reference_data/plasmid.py --fasta_path {} --snp_path {} --out_path {} --gene_name {}'.format(fasta_path,snp_path,
                #                                                                                                out_path, gen)
                if parameters ==[]:
                    cmd = 'python /reference_data/plasmid.py --fasta_path {} --snp_path {} --out_path {} --gene_name {}'.format(fasta_path,snp_path,out_path,gen)
                else:
                    cmd = 'python /reference_data/plasmid.py --fasta_path {} --snp_path {} --out_path {} --gene_name {} --max_amp_len {} --min_amp_len {} --temperature {} --min_GC {} --max_GC {} --max_len {}'.format(
                        fasta_path,snp_path, out_path,gen,parameters[0],parameters[1],parameters[2],parameters[3],parameters[4],parameters[5])
                print('-' * 30, 'cmd: \n', cmd)
                process = os.popen(cmd)  # return file
                output = process.read()
                print(output)
                process.close()
            result = []
            for temp_name in gene_list:
                for file in os.listdir(out_path):
                    if temp_name in file and '.csv' in file and 'primer_design' in file:
                        result.append(os.path.join(out_path, file))
        elif search_type == 'species_identification_type':
    
            print(f'当前结果来自：{search_type}')
            if type == '':
                target_path = instruction['instruction']['search_responses']['target_gene_info'][0]
                nontarget_path  = instruction['search_responses']['target_gene_info'][1]
                cmd = f'python /reference_data/run_primer_design.py --target_path {target_path} --nontarget_path {nontarget_path}'
                print('-' * 30, 'cmd: \n', cmd)
                process = os.popen(cmd)  # return file
                output = process.read()
                print(output)
                process.close()
            #直接上传单个基因序列文件进行引物设计。
            elif type =='  ':
                files_list= instruction['search_responses']['upload']
                output_path = out_path
                result = primer3_design(files_list,output_path)
    
    
        if len(result)==0:
            responses = {'response': 'Sorry, primer design failed.',"operations":[], 'primer': result, 'state': 'stop',
                         "primer_type": "snp_primer_design_type",
                         'primer_design_prompt': '',"stage":stage}
        else:
            responses = {'response':snp_primer_design_prompt,"operations":[],'primer': result, 'state': 'stop', "primer_type": "snp_primer_design_type",
                         'primer_design_prompt': snp_primer_design_prompt,"stage":stage}
    
        return responses
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
#def extract_json_response(response):
    # 从response中提取JSON内容
#    experiment_info_dict = response["choices"][0]["message"]['content']
#    print('experiment_info_dict',experiment_info_dict)
#    experiment_info_json_data = json.loads(experiment_info_dict)
    #match = re.search(r'```json\n(.*?)\n```', experiment_info_dict, re.DOTALL)
    #if match:
    #    experiment_info_json_str = match.group(1)
    #    experiment_info_json_data = json.loads(experiment_info_json_str)
    #else:
    #    print("GPT构造JSON错误，情况如下所示：\n",experiment_info_dict)
    #    exit()
 #  return experiment_info_json_data
    
#针对单个基因序列文件设计引物
def primer3_design(fa_path_list,output_path):
    result = []
    for fa_path in fa_path_list:
        seqs = ''
        with open(fa_path,'r') as f:
            lines = f.readlines()
            for ll in lines:
                if '>' in ll:
                    continue
                else:
                    seqs =seqs + ll
        i_length = len(seqs)
        filename = os.path.basename(fa_path).split('.')[0]
        primers = primer3.bindings.design_primers(
                seq_args={
                    "SEQUENCE_ID": filename,
                    "SEQUENCE_TEMPLATE": seqs,
                    "SEQUENCE_INCLUDED_REGION": [0, int(i_length)],
                },
                global_args={
                    "PRIMER_OPT_SIZE": 20,
                    "PRIMER_MIN_SIZE": 18,
                    "PRIMER_MAX_SIZE": 22,
                    "PRIMER_OPT_TM": 60.0,
                    "PRIMER_MIN_TM": 58.0,
                    "PRIMER_MAX_TM": 63.0,
                    "PRIMER_PAIR_MAX_DIFF_TM": 2.0,
                    "PRIMER_MIN_GC": 40.0,
                    "PRIMER_MAX_GC": 60.0,
                    "PRIMER_MAX_POLY_X": 3,
                    "PRIMER_PICK_RIGHT_PRIMER": 1,
                    "PRIMER_PICK_LEFT_PRIMER": 1,
                    "PRIMER_PICK_INTERNAL_OLIGO": 0,
                    "PRIMER_PRODUCT_SIZE_RANGE": "75-150",
                    "PRIMER_GC_CLAMP": 1,
                    "PRIMER_NUM_RETURN": 4,
                },
            )
                    # Create a name specific to the bacterium and CDS name
        filenameTSV = filename + ".tsv"
        # Save the dictionary as a .tsv file
        res = os.path.join(output_path, filenameTSV)
        with open(res, "w") as f:
            for key in primers.keys():
                f.write("%s\t %s\n" % (key, primers[key]))
        result.append(res)
    return result

def get_dna_seq(chr_num, start, end,file_content):
    chr_name = f'chromosome {chr_num}'
    # fasta文件中的染色体名称
    all_chr_name = 'Homo sapiens ' + chr_name + ', GRCh38.p14 Primary Assembly'
    print(all_chr_name)
    # 用于存储提取的内容
    extracted_content = ""

    # 是否找到了chr
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