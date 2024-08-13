from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

import pandas as pd
import os

import requests
import openai
import json
import re
import time
import primer3



def snp_primer_design(instruction,stage):
    if stage == 1:

        stage += 1
        response = '''The next critical step is to "determine primer design parameters". This step involves defining optimal conditions for primer design, which includes specifying maximum and minimum amplicon lengths, desired annealing temperatures, and GC content ranges. Additionally, it is important to determine the maximum length of the primers. These parameters will form the basis for designing specific primers. This step is critical to ensure the specificity and efficiency of subsequent amplification and detection processes, helping you move toward your experimental goals.'''
        return {'response':response,"operations":[{"key": "min_amp_len","type": ["input","required"],"title": "Min Amp length","value": [250],"options": []},
        {"key": "max_amp_len","type": ["input","required"],"title": "Max Amp length","value": [300],"options": []},
        {"key": "temperature","type": ["input","required"],"title": "temperature","value": [65],"options": []},
        {"key": "min_GC","type": ["input","required"],"title": "Min GC%","value": [40],"options": []},
        {"key": "max_GC","type": ["input","required"],"title": "Max GC%","value": [60],"options": []},
        {"key": "max_primer_len","type": ["input","required"],"title": "Max primer length","value": [23],"options": []},
        ], 'state': 'continue', "primer_type": "snp_primer_design_type","stage":stage}


    elif stage ==2:
        print('进入stage 2')
        primer_design_dict = {
        "Experiment_Direction": "primer design",
        "search_type":"",
        "max_amp_len": "",
        "min_amp_len":"",
        "temperature":"",
        "min_GC":"",
        "max_GC":"",
        "max_primer_len":"",}
        primer_design_dict["max_amp_len"] = instruction['instruction']["max_amp_len"]
        primer_design_dict["min_amp_len"] = instruction['instruction']["min_amp_len"]
        primer_design_dict["temperature"] = instruction['instruction']["temperature"]
        primer_design_dict["min_GC"] = instruction['instruction']["min_GC"]
        primer_design_dict["max_GC"] = instruction['instruction']["max_GC"]
        primer_design_dict["max_primer_len"] = instruction['instruction']["max_primer_len"]
        search_type = instruction['instruction']['search_responses']['search_type']
        primer_design_dict["search_type"] = search_type

        #exit()
        stage += 1
        parameters = [primer_design_dict["max_amp_len"][0],primer_design_dict["min_amp_len"][0],primer_design_dict["temperature"][0],
        primer_design_dict["min_GC"][0],primer_design_dict["max_GC"][0],primer_design_dict["max_primer_len"][0]]

        snp_primer_design_prompt = """
here is the designed primer file.
                    """
        out_path = '/reference_data/primer_result'
        dna_file_path = '/reference_data/dna_file'
        snp_file_path = '/reference_data/snp_file'
        search_type = primer_design_dict["search_type"]
        # #获取引物设计所必需的参数
        if search_type == 'genetic_disorder_type':
            gene_file_path = instruction['instruction']['search_responses']['gene_file_path']
            start_end_pos = instruction['instruction']['search_responses']['start_end_pos']

            result = []
            for temp_path in gene_file_path:
                fasta_path = temp_path
                na = os.path.basename(fasta_path).split('.')[0]
                start_pos = start_end_pos[0]
                end_pos = start_end_pos[1]
                script_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'plasmid.py')
                snp_path = ''
                if parameters ==[]:

                    cmd = 'python {} --fasta_path {}  --out_path {} --gene_name {} --start_pos {} --end_pos {}'.format(script_path,fasta_path, out_path,na,start_pos,end_pos)
                else:
                    cmd = 'python {} --fasta_path {}  --out_path {} --gene_name {} --max_amp_len {} --min_amp_len {} --temperature {} --min_GC {} --max_GC {} --max_len {} --start_pos {} --end_pos {}'.format(
                        script_path,fasta_path, out_path, na,parameters[0],parameters[1],parameters[2],parameters[3],parameters[4],parameters[5],start_pos,end_pos)
                print('-' * 30, 'cmd: \n', cmd)
                process = os.popen(cmd)  # return file
                output = process.read()
                print(output)
                process.close()
                for ff in os.listdir(out_path):
                    if ff.endswith('.csv') and na in ff and 'pool' in ff:
                        temp_path = os.path.join(out_path, ff)
                        primer_file = modify_primer_file(gene_file_path[0], temp_path)
                        result.append(primer_file)


        elif search_type == 'protein_mutation_type':
            print(f'当前结果来自：{search_type} ')
            path_list = instruction['instruction']['search_responses']['protein_gene_seq_path']
            protein_mutation_dict= instruction['instruction']['search_responses']['protein_mutation_dict']
            start_pos = protein_mutation_dict['start_pos']
            end_pos = protein_mutation_dict['end_pos']
            name_list = []
            for pa in path_list:
                temp_name = os.path.basename(pa)
                temp_name = temp_name.split('.')[0]
                temp_name = temp_name.replace('(','').replace(')','')
                name_list.append(temp_name)
                script_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'plasmid.py')
                # cmd = 'python /reference_data/plasmid.py --fasta_path {}  --out_path {} --gene_name {}'.format(pa,
                #                                                                                               out_path, temp_name)
                if parameters ==[]:
                    cmd = 'python {} --fasta_path {}  --out_path {} --gene_name {} --start_pos {} --end_pos {}'.format(script_path ,pa, out_path,
                                                                                                                                               temp_name,start_pos,end_pos)
                else:
                    cmd = 'python {} --fasta_path {}  --out_path {} --gene_name {} --max_amp_len {} --min_amp_len {} --temperature {} --min_GC {} --max_GC {} --max_len {} --start_pos {} --end_pos {}'.format(
                        script_path, pa, out_path, temp_name,parameters[0],parameters[1],parameters[2],parameters[3],parameters[4],parameters[5],start_pos,end_pos)
                print('-' * 30, 'cmd: \n', cmd)
                run_cmd(cmd)

            result =[]
            for temp_name in name_list:
                for file in os.listdir(out_path):
                    if temp_name in file and '.csv' in file:
                        temp_path = os.path.join(out_path,file)
                        gene_file = path_list[0]
                        primer_file = modify_primer_file(gene_file,temp_path)
                        result.append(primer_file)

        elif search_type == 'pathogen_drug_resistance_type':

            res_path = instruction['instruction']['search_responses']['filter_result_path']
            res_df = pd.read_csv(res_path)

            gene_list = list(set(list(res_df['Gene'])))

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
                script_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'plasmid.py')
                if parameters ==[]:
                    cmd = 'python {} --fasta_path {} --snp_path {} --out_path {} --gene_name {}'.format(script_path,fa_path,snp_path, out_path,gen)
                else:
                    cmd = 'python {} --fasta_path {} --snp_path {} --out_path {} --gene_name {} --max_amp_len {} --min_amp_len {} --temperature {} --min_GC {} --max_GC {} --max_len {}'.format(
                        script_path,fa_path, snp_path, out_path,gen,parameters[0],parameters[1],parameters[2],parameters[3],parameters[4],parameters[5])
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

            res_dict = instruction['instruction']['search_responses']['target_gene_info']
            res_df = pd.DataFrame.from_dict(res_dict)
            gene_list = list(set(list(res_df['Gene name'])))
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
                if len(snp_df) ==1:
                    start_pos = int(list(snp_df['START'])[0])
                    temp_len = int(parameters[0] - parameters[5]*2)
                    target_gene_part = temp_seq[start_pos-temp_len:start_pos+temp_len]
                    i_length = len(target_gene_part)
                    primers = primer3.bindings.design_primers(
                        seq_args={
                            "SEQUENCE_ID": gen,
                            "SEQUENCE_TEMPLATE": target_gene_part,
                            "SEQUENCE_INCLUDED_REGION": [0, int(i_length)],
                        },
                        global_args={
                            "PRIMER_OPT_SIZE": 20,
                            "PRIMER_MIN_SIZE": 18,
                            "PRIMER_MAX_SIZE": 23,
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
                    filenameTSV = gen + "_primer_design.tsv"
                    # Save the dictionary as a .tsv file
                    res = os.path.join(out_path, filenameTSV)
                    with open(res, "w") as f:
                        for key in primers.keys():
                            f.write("%s\t %s\n" % (key, primers[key]))
                    # result.append(res)
                else:
                    script_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'plasmid.py')
                    if parameters ==[]:
                        cmd = 'python {} --fasta_path {} --snp_path {} --out_path {} --gene_name {}'.format(script_path,fasta_path,snp_path,out_path,gen)
                    else:
                        cmd = 'python {} --fasta_path {} --snp_path {} --out_path {} --gene_name {} --max_amp_len {} --min_amp_len {} --temperature {} --min_GC {} --max_GC {} --max_len {}'.format(
                            script_path,fasta_path,snp_path, out_path,gen,parameters[0],parameters[1],parameters[2],parameters[3],parameters[4],parameters[5])

                    process = os.popen(cmd)  # return file
                    output = process.read()
                    process.close()
            result = []
            for temp_name in gene_list:
                for file in os.listdir(out_path):
                    if len(snp_df) == 1:
                        if temp_name in file  and 'design.tsv' in file:
                            result.append(os.path.join(out_path, file))
                    else:
                        if temp_name in file and 'pool' in file:
                            result.append(os.path.join(out_path, file))
        elif search_type == 'species_identification_type':
            search_result_dict = instruction['instruction']['search_responses']['species_identification_dict']
            name_list = search_result_dict['organism_names']

            result = []
            #如果是基因文件，则直接设计引物
            gene_path = search_result_dict['target_cds_path']
            output_path = '/reference_data/primer_result'
            if type(gene_path) == str:
                gene_path = [gene_path]
            result = primer3_design(gene_path, parameters, output_path, name_list)


        if len(result)==0:
            responses = {'response': 'Sorry, primer design failed.',"operations":[], 'primer': result, 'state': 'stop',
                         "primer_type": "snp_primer_design_type",
                         'primer_design_prompt': '',"stage":stage}
        else:
            responses = {'response':snp_primer_design_prompt,"operations":[],'primer': result, 'state': 'stop', "primer_type": "snp_primer_design_type",
                         'primer_design_prompt': snp_primer_design_prompt,"stage":stage}
        return responses

def run_cmd(cmd):
    import subprocess

    # 使用 subprocess.Popen 来运行命令
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    # 实时输出标准输出和标准错误
    while True:
        output = process.stdout.readline()
        if output == '' and process.poll() is not None:
            break
        if output:
            print(output.strip())

    # 检查是否有错误输出
    err = process.stderr.read()
    if err:
        print("Error:\n", err)

def modify_primer_file(fasta_path,primer_file):
    seq_raw = ''
    with open(fasta_path,'r')as f:
        lines = f.readlines()
        for ll in lines:
            if '>'in ll:
                continue
            seq_raw+= ll.replace('\n','')
    reference_sequence = seq_raw.lower()
    df = pd.read_csv(primer_file)
    fp_list = list(df['fP'])
    rp_list = list(df['rP'])
    start =[]
    insert_start = []
    insert_end = []
    end = []
    for ii in range(len(fp_list)):
        pos_all = find_primer_positions(reference_sequence, fp_list[ii], rp_list[ii])
        start.append(pos_all[0])
        insert_start.append(pos_all[1])
        insert_end.append(pos_all[2])
        end.append(pos_all[3])
    new_df = df.loc[:,['amp_id','fP','rP','amp','insert','fP_full','rP_full']]
    new_df.insert(3,'start',start)
    new_df.insert(4,'insert_start',insert_start)
    new_df.insert(5,'insert_end',insert_end)
    new_df.insert(6,'end',end)
    new_df.to_csv(primer_file)
    return primer_file

def reverse_complement(seq):
    complement = {'a': 't', 't': 'a', 'c': 'g', 'g': 'c'}
    return "".join(complement[base] for base in reversed(seq))


def find_primer_positions(reference_seq, forward_primer, reverse_primer):
    # 找到左引物的起始位置
    forward_start = reference_seq.find(forward_primer)
    # 获取右引物的互补反向序列
    reverse_complement_primer = reverse_complement(reverse_primer)
    # 找到右引物的起始位置
    reverse_start = reference_seq.find(reverse_complement_primer)
    # 如果找不到引物，返回None
    if forward_start == -1 or reverse_start == -1:
        return None

    # 计算各个位置
    start = forward_start + 1
    insert_start = start + len(forward_primer)
    insert_end = reverse_start
    end = insert_end + len(reverse_primer)

    # 返回四个位置
    return start, insert_start, insert_end, end
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

#针对单个基因序列文件设计引物
def primer3_design(fa_path_list, parameters, output_path, name_list=[[1]]):
    result = []
    #exit()
    for name_info, fa_path in zip(name_list, fa_path_list):
        #print(name_info)
        #continue
        name = name_info[0]
        seqs = ''
        with open(fa_path,'r') as f:
            lines = f.readlines()
            for ll in lines:
                if '>' in ll:
                    continue
                else:
                    seqs =seqs + ll.replace('\n','')
        #print(seqs)
        i_length = len(seqs)
        primer_num = 100
        filename = os.path.basename(fa_path).split('.')[0]
        primers = primer3.bindings.design_primers(
                seq_args={
                    "SEQUENCE_ID": filename,
                    "SEQUENCE_TEMPLATE": seqs,
                    "SEQUENCE_INCLUDED_REGION": [0, int(i_length)],
                },
                global_args={
                    "PRIMER_OPT_SIZE": parameters[-1],
                    "PRIMER_MIN_SIZE": parameters[-1]-2,
                    "PRIMER_MAX_SIZE": parameters[-1]+2,
                    "PRIMER_OPT_TM": parameters[2],
                    "PRIMER_MIN_TM": parameters[2]-3,
                    "PRIMER_MAX_TM": parameters[2]+3,
                    "PRIMER_PAIR_MAX_DIFF_TM": 2.0,
                    "PRIMER_MIN_GC": parameters[3],
                    "PRIMER_MAX_GC": parameters[4],
                    "PRIMER_MAX_POLY_X": 3,
                    "PRIMER_PICK_RIGHT_PRIMER": 1,
                    "PRIMER_PICK_LEFT_PRIMER": 1,
                    "PRIMER_PICK_INTERNAL_OLIGO": 0,
                    "PRIMER_PRODUCT_SIZE_RANGE": f"{parameters[1]}-{parameters[0]}",
                    "PRIMER_GC_CLAMP": 1,
                    "PRIMER_NUM_RETURN": primer_num,
                },
            )
        # Create a name specific to the bacterium and CDS name
        # Create a name specific to the bacterium and CDS name
        PRIMER_PAIR_0_PENALTY = []
        PRIMER_PAIR_0_COMPL_ANY_TH = []
        PRIMER_PAIR_0_COMPL_END_TH= []
        PRIMER_PAIR_0_PRODUCT_SIZE= []
        PRIMER_LEFT_0_SEQUENCE= []
        PRIMER_RIGHT_0_SEQUENCE= []
        PRIMER_LEFT_0= []
        PRIMER_RIGHT_0= []
        PRIMER_LEFT_0_TM= []
        PRIMER_RIGHT_0_TM= []
        PRIMER_LEFT_0_GC_PERCENT= []
        PRIMER_RIGHT_0_GC_PERCENT= []
        PRIMER_LEFT_0_SELF_ANY_TH= []
        PRIMER_RIGHT_0_SELF_ANY_TH= []
        PRIMER_LEFT_0_SELF_END_TH= []
        PRIMER_RIGHT_0_SELF_END_TH= []
        PRIMER_LEFT_0_HAIRPIN_TH= []
        PRIMER_RIGHT_0_HAIRPIN_TH= []
        PRIMER_LEFT_0_END_STABILITY= []
        PRIMER_RIGHT_0_END_STABILITY= []

        for ii in range(primer_num):
            PRIMER_PAIR_0_PENALTY.append(primers[f'PRIMER_PAIR_{ii}_PENALTY'])
            PRIMER_PAIR_0_COMPL_ANY_TH.append(primers[f'PRIMER_PAIR_{ii}_COMPL_ANY_TH'])
            PRIMER_PAIR_0_COMPL_END_TH.append(primers[f'PRIMER_PAIR_{ii}_COMPL_END_TH'])
            PRIMER_PAIR_0_PRODUCT_SIZE.append(primers[f'PRIMER_PAIR_{ii}_PRODUCT_SIZE'])
            PRIMER_LEFT_0_SEQUENCE.append(primers[f'PRIMER_LEFT_{ii}_SEQUENCE'])
            PRIMER_RIGHT_0_SEQUENCE.append(primers[f'PRIMER_RIGHT_{ii}_SEQUENCE'])
            PRIMER_LEFT_0.append(primers[f'PRIMER_LEFT_{ii}'])
            PRIMER_RIGHT_0.append(primers[f'PRIMER_RIGHT_{ii}'])
            PRIMER_LEFT_0_TM.append(primers[f'PRIMER_LEFT_{ii}_TM'])
            PRIMER_RIGHT_0_TM.append(primers[f'PRIMER_RIGHT_{ii}_TM'])
            PRIMER_LEFT_0_GC_PERCENT.append(primers[f'PRIMER_LEFT_{ii}_GC_PERCENT'])
            PRIMER_RIGHT_0_GC_PERCENT.append(primers[f'PRIMER_RIGHT_{ii}_GC_PERCENT'])
            PRIMER_LEFT_0_SELF_ANY_TH.append(primers[f'PRIMER_LEFT_{ii}_SELF_ANY_TH'])
            PRIMER_RIGHT_0_SELF_ANY_TH.append(primers[f'PRIMER_RIGHT_{ii}_SELF_ANY_TH'])
            PRIMER_LEFT_0_SELF_END_TH.append(primers[f'PRIMER_LEFT_{ii}_SELF_END_TH'])
            PRIMER_RIGHT_0_SELF_END_TH.append(primers[f'PRIMER_RIGHT_{ii}_SELF_END_TH'])
            PRIMER_LEFT_0_HAIRPIN_TH.append(primers[f'PRIMER_LEFT_{ii}_HAIRPIN_TH'])
            PRIMER_RIGHT_0_HAIRPIN_TH.append(primers[f'PRIMER_RIGHT_{ii}_HAIRPIN_TH'])
            PRIMER_LEFT_0_END_STABILITY.append(primers[f'PRIMER_LEFT_{ii}_END_STABILITY'])
            PRIMER_RIGHT_0_END_STABILITY.append(primers[f'PRIMER_RIGHT_{ii}_END_STABILITY'])
        raw_all_param_df = pd.DataFrame({"PRIMER_PAIR_PENALTY":PRIMER_PAIR_0_PENALTY,"PRIMER_PAIR_COMPL_ANY_TH":PRIMER_PAIR_0_COMPL_ANY_TH,
                                     "PRIMER_PAIR_COMPL_END_TH":PRIMER_PAIR_0_COMPL_END_TH,"PRIMER_PAIR_PRODUCT_SIZE":PRIMER_PAIR_0_PRODUCT_SIZE,
                                     "PRIMER_LEFT_SEQUENCE":PRIMER_LEFT_0_SEQUENCE,"PRIMER_RIGHT_SEQUENCE":PRIMER_RIGHT_0_SEQUENCE,
                                     "PRIMER_LEFT":PRIMER_LEFT_0,"PRIMER_RIGHT":PRIMER_RIGHT_0,
                                     "PRIMER_LEFT_TM":PRIMER_LEFT_0_TM,"PRIMER_RIGHT_TM":PRIMER_RIGHT_0_TM,
                                     "PRIMER_LEFT_GC_PERCENT":PRIMER_LEFT_0_GC_PERCENT,"PRIMER_RIGHT_GC_PERCENT":PRIMER_RIGHT_0_GC_PERCENT,
                                     "PRIMER_LEFT_SELF_ANY_TH":PRIMER_LEFT_0_SELF_ANY_TH,"PRIMER_RIGHT_SELF_ANY_TH":PRIMER_RIGHT_0_SELF_ANY_TH,
                                     "PRIMER_LEFT_SELF_END_TH":PRIMER_LEFT_0_SELF_END_TH,"PRIMER_RIGHT_SELF_END_TH":PRIMER_RIGHT_0_SELF_END_TH,
                                     "PRIMER_LEFT_END_STABILITY":PRIMER_LEFT_0_END_STABILITY,"PRIMER_RIGHT_END_STABILITY":PRIMER_RIGHT_0_END_STABILITY,
                                     "PRIMER_LEFT_HAIRPIN_TH":PRIMER_LEFT_0_HAIRPIN_TH,"PRIMER_RIGHT_HAIRPIN_TH":PRIMER_RIGHT_0_HAIRPIN_TH})


        all_param_df = raw_all_param_df[(raw_all_param_df['PRIMER_LEFT_HAIRPIN_TH']<=0.0) & (raw_all_param_df["PRIMER_RIGHT_HAIRPIN_TH"]<=0.0)]

        if len(all_param_df)==0:
            left_sorted_df = raw_all_param_df.sort_values(by='PRIMER_LEFT_HAIRPIN_TH').head(50)
            right_sorted_df = raw_all_param_df.sort_values(by='PRIMER_RIGHT_HAIRPIN_TH').head(50)
            intersection_indices = set(left_sorted_df.index) & set(right_sorted_df.index)
            all_param_df = raw_all_param_df.loc[intersection_indices]

        filenameTSV = filename + ".csv"
        # Save the dictionary as a .tsv file
        res = os.path.join(output_path, filenameTSV)
        top_100_df = all_param_df.head(100)
        order_top100_df = top_100_df.sort_values('PRIMER_PAIR_PENALTY')#.reset_index(drop=True, inplace=True)#,inplace=True)
        order_top100_df.reset_index(drop=True, inplace=True)
        order_top100_df.to_csv(res,index=False)

        #排查二聚体情况
        top_50_index = detect_dimer(order_top100_df)
        top50_df = order_top100_df.loc[top_50_index]
        #排查mismatch情况
        mismatch_list = detect_mismatch(top50_df, seqs)
        filtermismatch_top50_df = top50_df[[bool(value) for value in mismatch_list]]
        top50_df = filtermismatch_top50_df
        result_primer = top50_df.head(1)
        result_primer.insert(0,'Organism Name',name)
        #print(result_primer)
        filenameTSV = name.replace(' ','_')+'_primer_' + ".csv"
        res = os.path.join(output_path, filenameTSV)
        result_primer.to_csv(res,index=False)
        test_pd = pd.read_csv(res)

        result.append(res)
        #TODO:通过数据库全比对，提高引物设计的特异性和成功率，计算引物与非target的匹配数量，选择最小的推荐
        #get_best_primer(top50_df,filename)
    return result

#3'端的50%序列内，不能错配
def detect_mismatch(primer3_df, target_cds):
    print('+++',primer3_df)
    left_primer = primer3_df['PRIMER_LEFT_SEQUENCE'].tolist()
    left_primer_start = primer3_df['PRIMER_LEFT'].apply(lambda x: x[0]).tolist()
    left_primer_length = primer3_df['PRIMER_LEFT'].apply(lambda x: x[1]).tolist()
    right_primer = primer3_df['PRIMER_RIGHT_SEQUENCE'].tolist()
    right_primer_start = primer3_df['PRIMER_RIGHT'].apply(lambda x: x[0]).tolist()
    right_primer_length = primer3_df['PRIMER_RIGHT'].apply(lambda x: x[1]).tolist()

    left_primer_list = [seq[-3:] for seq in left_primer]
    right_primer_list = [seq[-1:-4:-1] for seq in right_primer]


    left_target_list = [target_cds[i:i+l] for i,l in zip(left_primer_start,left_primer_length)]
    print('left target', left_target_list)
    target_left_match_list = [i[-3:] for i in left_target_list]
    right_target_list = [target_cds[i-l+1:i+1] for i,l in zip(right_primer_start,right_primer_length)]
    print('right target', right_target_list)
    target_right_match_list = [i[:3] for i in right_target_list]


    #target_right_match_list = [target_cds[i[0]:i[1]] for i in target_cds[i-l+1:i+1] for i,l in zip(right_primer_start,right_primer_length)]]
    #exit()
    right_index = []
    for x,y,z,u in zip(left_primer_list,right_primer_list,target_left_match_list,target_right_match_list):

        if x==z and Seq(u) == Seq(y).complement():
            print('Match')
            right_index.append(1)
            print('+++++++++++++++')
            continue
        print('No Match')
        right_index.append(0)
        print('+++++++++++++++')

    print(right_index)
    #exit()
    return right_index

def detect_dimer(df, method="single"):
    f_list = list(df['PRIMER_LEFT_SEQUENCE'])
    r_list = list(df['PRIMER_RIGHT_SEQUENCE'])

    if method=='multi':
        all_list = f_list+r_list
        for ii in range(len(all_list)-1):
            s1 = all_list[ii]
            for jj in range(ii+1,len(all_list)):
                s2 = all_list[jj]
                match_list.append((s1,s2))
                dimer_1F_1R = primer3.calcHeterodimer(s1, s2)
                free_energy_list.append(dimer_1F_1R.dg)
                structure_list.append(dimer_1F_1R.structure_found)

    match_list =[]
    free_energy_list = []
    structure_list = []
    index_list =[]
    for index,(s1,s2) in enumerate(zip(f_list,r_list)):
        match_list.append((s1,s2))
        dimer_1F_1R = primer3.calcHeterodimer(s1, s2)
        free_energy_list.append(dimer_1F_1R.dg)
        structure_list.append(dimer_1F_1R.structure_found)
        index_list.append(index)
    print('---------',index_list)
    result_df = pd.DataFrame({'original_index':index_list, 'F_match_R':match_list,'free_energy':free_energy_list,'structure_found':structure_list})
    result_df.sort_values('free_energy',inplace=True)
    top50_df = result_df.head(50)
    result_df.to_csv('/reference_data/二聚体分析_两两匹配.csv',index=False)

    return top50_df['original_index']



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

def getglCMD(kmer_file, ofile):
    cmd = f'/home/YJhou/workspace/pcr_paper/search_base/fasta-36.3.8i/bin/glsearch36 {kmer_file} /home/YJhou/workspace/pcr_paper/search_base/blast_db/refseq.fa -E 30000 -m 8 -n -T 30 --V > {ofile}'
    #cmd = 'time blastn -task blastn-short'
    #cmd += f' -query {kmer_file} -db {db} -num_threads 100 -perc_identity {ident} -qcov_hsp_perc {qcov_hsp_perc} -evalue 10000 -outfmt "6 qseqid sseqid qstart qend sstart send length pident mismatch gapopen evalue bitscore staxids sscinames qlen"  > {ofile}'
    return cmd

def calculate_optimal_primer(blast_result_path):
    #with open('blast_result_path','r') as f:
    blast_df = pd.read_csv(blast_result_path, sep='\t', header=None,
                                   names=['qseqid', 'sseqid', 'qstart', 'qend', 'sstart', 'send', 'length', 'pident', 'mismatch', 'gapopen',
                                          'evalue', 'bitscore', 'staxids', 'sscinames', 'qlen'])
    print('blast特异性比对结果：\n', blast_df)
    exit()


def get_best_primer(filter_primer_df, name):
    filter_primer_path = '/reference_data/primer_blast_result'
    #with open(f'{filter_primer_path}/{name}_primer.fna', 'w') as f:
    records = []
    for index, row in filter_primer_df.iterrows():
        f_seq = row['PRIMER_LEFT_SEQUENCE']
        r_seq = row['PRIMER_RIGHT_SEQUENCE']
        f_l = int(len(f_seq)*1)#向下取整
        r_l = int(len(r_seq)*1)
        f_record = SeqRecord(Seq(f_seq[-f_l:]), id=f"lcl|{index+1}_1", description="Forward primer chunk")
        r_record = SeqRecord(Seq(r_seq[:r_l]), id=f"lcl|{index+1}_2", description="Reverse primer chunk")
        records.append(f_record)
        records.append(r_record)

    output_file = f'{filter_primer_path}/{name}_primer.fna'
    SeqIO.write(records, output_file, 'fasta')
    #TODO:补充序列全比对过滤，确保引物与模板没有mismatch
    #TODO:
    ofile = f'{filter_primer_path}/{name}_primer_glsearch36_output'
    primer_file_path = output_file
    cmd = getglCMD(primer_file_path, ofile)
    os.system(cmd)
    calculate_optimal_primer(ofile)