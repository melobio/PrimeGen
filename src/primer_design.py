from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

import pandas as pd
import os, sys, gzip, re
import requests
import openai
import json
import time
import primer3
import logging
import math
#logger = logging.getLogger(__name__)
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)
def snp_primer_design(instruction,stage):
    if stage == 1: 
        stage += 1
        response = '''The next crucial step is to "determine the parameters for primer design". This includes defining optimal conditions such as the maximum and minimum amplicon lengths, desired annealing temperatures, and GC content ranges. These parameters will serve as the foundation for designing specific primers, ensuring the specificity and efficiency of subsequent amplification and detection processes, ultimately guiding you toward your experimental goals.'''
        return {'response':response,"operations":[{"key": "min_amp_len","type": ["input","required"],"title": "Min Amp length","value": [210],"options": []},
        {"key": "max_amp_len","type": ["input","required"],"title": "Max Amp length", "value": [270],"options": []},
        {"key": "temperature","type": ["input","required"],"title": "temperature","value": [60], "options": []},
        {"key": "min_GC","type": ["input","required"],"title": "Min GC","value": [0.4],"options": []},
        {"key": "max_GC","type": ["input","required"],"title": "Max GC","value": [0.65],"options": []},
        {"key": "max_primer_len","type": ["input","required"],"title": "Max primer length","value": [23],"options": []},
        ], 'state': 'continue', "primer_type": "snp_primer_design_type","stage":stage}

    elif stage ==2:
        print('Entering stage 2')
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
        print('Received search:', instruction['instruction']['search_responses'])
        print('Received parameters:', instruction['instruction'])
        #exit()
        stage += 1
        print('primer_design_dict',primer_design_dict)
        parameters = [primer_design_dict["min_amp_len"][0],primer_design_dict["max_amp_len"][0],primer_design_dict["temperature"][0],
        primer_design_dict["min_GC"][0],primer_design_dict["max_GC"][0],primer_design_dict["max_primer_len"][0]]
        
        snp_primer_design_prompt = """
here is the designed primer file.
                    """
        out_path = '/reference_data/primer_result'
        dna_file_path = '/reference_data/dna_file'
        snp_file_path = '/reference_data/snp_file'
        
        search_type = primer_design_dict["search_type"]
    
        print('search_type',search_type)

        # Get required parameters for primer design
        if search_type == 'genetic_disorder_type':
            print(f'Current results from: {search_type}')
            gen_name_list = instruction['instruction']['search_responses']['data']['history_data']['gen_name_list']
            temp_file_name = '_'.join(gen_name_list)
            records = []
            with gzip.open('/reference_data/gencode.v46.primary_assembly.annotation.gtf.gz', 'rt') as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    tmp = line.strip().split("\t")
                    if tmp[2] == 'transcript':
                        chrom, ori, start, end, strand, info = tmp[0], tmp[1], tmp[3], tmp[4], tmp[6], tmp[8]
                        infod = dict()
                        gene = re.search(r'gene_name\s"(.*?)"', info).group(1)
                        transcript = re.search(r'transcript_id\s"(.*?)"', info).group(1)
                        level = re.search(r'transcript_support_level\s"(.*?)"', info).group(
                            1) if 'transcript_support_level' in info else ""
                        tags = " | ".join(re.findall(r'tag\s"(.*?)"', info)) if 'tag' in info else ""
                        records.append([chrom, ori, start, end, strand, gene, transcript, level, tags])

            trans_tmp = pd.DataFrame(records,
                                     columns=['chrom', 'ori', 'start', 'end', 'strand', 'gene', 'transcript', 'level',
                                              'tags'])

            trans = trans_tmp[trans_tmp['gene'].isin(gen_name_list)]
                # pd.merge(disease, trans_tmp.drop_duplicates(), left_on='Gene', right_on='gene', how='left')
            trans = trans[trans.tags.map(lambda x: "MANE_Select" in x)]
            mane = set(trans.transcript)

            records = []
            with gzip.open(f'/reference_data/gencode.v46.primary_assembly.annotation.gtf.gz', 'rt') as f:  # Modify path based on location
                for line in f:
                    if line.startswith('#'):
                        continue
                    tmp = line.strip().split("\t")
                    if tmp[2] == 'exon':  # exon
                        chrom, ori, start, end, strand, info = tmp[0], tmp[1], tmp[3], tmp[4], tmp[6], tmp[8]
                        infod = dict()
                        gene = re.search(r'gene_name\s"(.*?)"', info).group(1)
                        exon = re.search(r'exon_number\s(\d+);', info).group(1)
                        transcript = re.search(r'transcript_id\s"(.*?)"', info).group(1)
                        if transcript in mane:
                            records.append([chrom, start, end, strand, gene, transcript, exon])

            cds = pd.DataFrame(records,
                               columns=['chrom', 'start', 'end', 'strand', 'gene', 'transcript', 'exon_number'])
            cds = cds.drop_duplicates()
            cds_disease = cds[cds['gene'].isin(gen_name_list)]
            cds_disease['start'] = cds_disease['start'].astype('int')
            cds_disease['end'] = cds_disease['end'].astype('int')
            wd = '/reference_data'
            merge_num = 580
            # Merge CDS with interval xx
            cds_disease[['chrom', 'start', 'end', 'strand', 'gene', 'transcript', 'exon_number']].to_csv(
                f'{wd}/cds.bed', sep="\t", index=False, header=False)
            os.system(f'bedtools sort -i {wd}/cds.bed > {wd}/cds.sort.bed')
            os.system(
                f'bedtools merge -d {merge_num} -o collapse -c 4,5,7 -i {wd}/cds.sort.bed -delim : > {wd}/cds.merge.bed')

            # Add header
            cds = pd.read_csv(f'{wd}/cds.merge.bed', sep="\t", header=None)
            cds.columns = ['chrom', 'start', 'end', 'strand', 'gene', 'exon']
            cds.to_csv(f'{wd}/hg38_disease.cds_pos.tsv', index=False, sep="\t")
            os.system(f'rm {wd}/cds.bed {wd}/cds.sort.bed {wd}/cds.merge.bed')

            # Extend
            max_amp_len = 270
            chr2seq = dict()
            for r in SeqIO.parse(f'{wd}/GCF_000001405.40_GRCh38.p14_genomic.fna', 'fasta'):
                chr2seq[r.id] = str(r.seq)
            cds = pd.read_csv(f'{wd}/hg38_disease.cds_pos.tsv', sep="\t")
            expand = 100  # Number of bases to extend CDS boundaries (adjust as needed) #72
            final_fasta_file = f'{wd}/{temp_file_name}_cds.fasta'
            with open(final_fasta_file, 'w') as out:
                for i, r in cds.iterrows():
                    header = "_".join([r['gene'], r['strand'], str(r['exon'])])
                    if (r['end'] - r['start'] + 1) <= max_amp_len:
                        start = int(r['start']) - 1 - expand - math.ceil(max_amp_len / 2)
                        end = int(r['end']) + expand + math.floor(max_amp_len / 2)
                    else:
                        start, end = int(r['start']) - 1 - expand, int(r['end']) + expand

                    seq = chr2seq[r['chrom']][start:end]
                    seq = seq.lower()
                    out.write(f'>{header}\n{seq}\n')

            script_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'primer_design_agent.py')
            cmd = f'python {script_path} --fasta_path {final_fasta_file} --file_name {temp_file_name}'
            print("cmd: ",cmd)
            process = os.popen(cmd)  # return file
            output = process.read()
            print(output)
            process.close()
            result =[]
            out_path = '/reference_data/primer_result/output/'
            for ff in os.listdir(out_path):
                if ff.endswith('.csv') and temp_file_name+'_pool' in ff:
                    temp_path = os.path.join(out_path, ff)
                    result.append(temp_path)

        elif search_type == 'protein_mutation_type':
            print(f'Current results from: {search_type} ')
            logging.info(f'Current results from: {search_type} ')
            path_list = instruction['instruction']['search_responses']['data']['protein_gene_seq_path']
            protein_mutation_dict = instruction['instruction']['search_responses']['data']['protein_mutation_dict']
            start_pos = protein_mutation_dict['start_pos']
            end_pos = protein_mutation_dict['end_pos']
            name_list = []
            print('path_list',path_list)
            for pa in path_list:
                temp_name = os.path.basename(pa)
                temp_name = temp_name.split('.')[0]
                temp_name = temp_name.replace('(', '').replace(')', '')
                name_list.append(temp_name)
                script_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'plasmid.py')
                if parameters == []:
                    cmd = 'python {} --fasta_path {}  --out_path {} --gene_name {} --start_pos {} --end_pos {}'.format(
                        script_path, pa, out_path,
                        temp_name, start_pos, end_pos)
                    print('mmmmmmmmm',cmd)
                else:
                    cmd = 'python {} --fasta_path {}  --out_path {} --gene_name {} --min_amp_len {} --max_amp_len {} --temperature {} --min_GC {} --max_GC {} --max_len {} --start_pos {} --end_pos {}'.format(
                        script_path, pa, out_path, temp_name, parameters[0], parameters[1], parameters[2],
                        parameters[3], parameters[4], parameters[5], start_pos, end_pos)
                    print('xxxxxxxxxxxxxx',cmd)
                logging.info('-' * 30+'cmd: \n'+cmd)
                process = os.popen(cmd)  # return file
                output = process.read()
                print(output)
                process.close()


            result = []
            for temp_name in name_list:
                for file in os.listdir(out_path):
                    if temp_name in file and '.csv' in file:
                        temp_path = os.path.join(out_path, file)
                        # gene_file = path_list[0]
                        # primer_file = modify_primer_file(gene_file, temp_path)
                        result.append(temp_path)

        elif search_type == 'pathogen_drug_resistance_type':
            print(f'Current results from: {search_type} ')
            res_path = instruction['instruction']['search_responses']['data']['filter_result_path']
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
                script_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'plasmid.py')
                if parameters ==[]:
                    cmd = 'python {} --fasta_path {} --snp_path {} --out_path {} --gene_name {}'.format(script_path,fa_path,snp_path, out_path,gen)
                else:
                    cmd = 'python {} --fasta_path {} --snp_path {} --out_path {} --gene_name {} --min_amp_len {} --max_amp_len {} --temperature {} --min_GC {} --max_GC {} --max_len {}'.format(
                        script_path,fa_path,snp_path, out_path,gen,parameters[0],parameters[1],parameters[2],parameters[3],parameters[4],parameters[5])
                logging.info('-' * 30+'cmd: \n'+cmd)
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

        elif search_type == 'whole_genome_type':
            print(f'Current results from: {search_type}')
            seqs_file = instruction['instruction']['search_responses']['data']['target_fna_list'][0]
            design_tpye = instruction['instruction']['search_responses']['data']['design_type'][0]

            snp_file = instruction['instruction']['search_responses']['data']['snp_file'][0]
            file_name = os.path.basename(seqs_file).split('.')[0]
            out_path = '/reference_data/primer_result/output/'
            script_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'primer_design_agent_snp_whole_panel.py')
            cmd = f'python {script_path} --fasta_path {seqs_file} --snp_file_path {snp_file} --file_name {file_name} --output_path {out_path}'
            print('cmd',cmd)
            process = os.popen(cmd)  # return file
            output = process.read()
            print(output)
            process.close()
            result = []
            for file in os.listdir(out_path):
                if file_name in file and 'pool' in file and '.csv' in file:
                    result.append(os.path.join(out_path, file))

        elif search_type == 'cancer_type':
            print(f'Current results from: {search_type}')
            res_dict = instruction['instruction']['search_responses']['data']['target_gene_info']
            res_df = pd.DataFrame.from_dict(res_dict)
            gene_list = list(set(list(res_df['Gene name'])))
            print('gene_list', gene_list)
            #
            fasta_file = "/reference_data/GCF_000001405.40_GRCh38.p14_genomic_filter.fna"

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
                print('Total length of this gene:', len(temp_seq))
                temp_path = os.path.join(dna_file_path,gen+'_seq.fa')
                with open(temp_path, 'w') as f:
                    f.write(f'>{gen}\n')
                    f.write(temp_seq)
                    f.close()
    
                fasta_path = temp_path
                print(f'Gene sequence extraction complete. Starting primer design')
                if len(snp_df) ==1:
                    start_pos = int(list(snp_df['START'])[0])
                    temp_len = int(int(parameters[0]) - int(parameters[5])*2)
                    target_gene_part = temp_seq[start_pos-temp_len:start_pos+temp_len]
                    i_length = len(target_gene_part)
                    print('Starting primer3 design')
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
                            "PRIMER_MIN_GC": int(parameters[3]),
                            "PRIMER_MAX_GC": int(parameters[4]),
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
                        cmd = 'python {} --fasta_path {} --snp_path {} --out_path {} --gene_name {} --min_amp_len {} --max_amp_len {} --temperature {} --min_GC {} --max_GC {} --max_len {}'.format(
                            script_path,fasta_path,snp_path, out_path,gen,parameters[0],parameters[1],parameters[2],parameters[3],parameters[4],parameters[5])
                    logging.info('-' * 30+'cmd: \n'+cmd)
                    process = os.popen(cmd)  # return file
                    output = process.read()
                    print(output)
                    process.close()
            result = []
            for temp_name in gene_list:
                for file in os.listdir(out_path):
                    if temp_name in file  and 'primer_design' in file and('.csv'in file or '.tsv'in file):
                        result.append(os.path.join(out_path, file))
                        
        elif search_type == 'species_identification_type':
            search_result_dict = instruction['instruction']['search_responses']["data"]['species_identification_dict']
            print(f'Current results from: {search_type}')
            logging.info(f'Current results from: \n{search_type}')
            #print('Search results: \n',search_result_dict)
            logging.info(f'Search results: \n{search_result_dict}')
            name_list = search_result_dict['organism']
            #print(search_result_dict['Target_cds_path'])
            result = []
            # If it's a gene file, design primers directly
            gene_path = search_result_dict['target_cds_path']
            output_path = '/reference_data/primer_result'
            if type(gene_path) == str:
                gene_path = [gene_path]
            result = primer3_design(gene_path, parameters, output_path, name_list)
    
    
        if len(result)==0:
            responses = {'response': 'Sorry, primer design failed.',"operations":[],
                         "data":{'primer': result}, 'state': 'stop',
                         "primer_type": "snp_primer_design_type",
                         'primer_design_prompt': '',"stage":stage}
        else:
            responses = {'response':snp_primer_design_prompt,"operations":[],"data":{'primer': result}, 'state': 'stop', "primer_type": "snp_primer_design_type",
                         'primer_design_prompt': snp_primer_design_prompt,"stage":stage}    
        return responses

    #elif stage=='3':

def redesign_primer(instruction,stage):
    #
    if stage == 1:
        stage += 1
        response = '''The current stage is the primer redesign stage. You need to provide the following three files. We will analyze the contents in the tables and redesign the primers that do not meet the requirements.'''
        return {'response': response, "operations": [
            {"key": "summary_depth", "type": ["file"], "title": "upload summary_depth file", "value": [],
             "options": []},
            {"key": "efficiency", "type": ["file"], "title": "upload efficiency file", "value": [],
             "options": []},
            {"key": "first_design_file", "type": ["file"], "title": "upload first_design pkl file", "value": [],
             "options": []},
          ], 'state': 'continue', "primer_type": "redesign_primer_type", "stage": stage}

    elif stage == 2:
        stage += 1
        # Analyze primers that do not meet the requirements and find the corresponding reference sequences
        summary_depth_file = instruction['instruction']['summary_depth']
        efficiency_file = instruction['instruction']['efficiency']
        first_design_file = instruction['instruction']['first_design_file']
        uniformity_bad_primer_list,eiffency_bad_primer_list = get_bad_primer(summary_depth_file+efficiency_file)
        temp_list =[]
        pri_list = []
        for temp_amp in uniformity_bad_primer_list:
            if temp_amp in eiffency_bad_primer_list:
                temp_list.append(temp_amp)
                pri_list.append('fr')

        redesign_df = pd.DataFrame({"amp_id":temp_list,"primer_type":pri_list})
        redesign_amp_info_path = '/reference_data/primer_result/redesign_amp_info.csv'
        redesign_df.to_csv(redesign_amp_info_path,index=False)
        response = 'The following are the primers that we screened out and need to be redesigned. You can download and view them. We will redesign these primers later.'

        responses = {'response': response, "operations": [], "data": {'redesign_primer': redesign_amp_info_path,'first_design_file':first_design_file},
                     'state': 'continue', "primer_type": "redesign_primer_type",
                    "stage": stage}
        return responses

    elif stage == 3:
        redesign_df_path = instruction['instruction']['data']['redesign_primer']
        first_design_file = instruction['instruction']['data']['first_design_file'][0]
        redesign_df = pd.read_csv(redesign_df_path)
        redesign_num = len(redesign_df)
        # redesign_amp_dict = {}
        # for k , v in zip(list(redesign_df['amp_id']),list(redesign_df['primer_type'])):
        #     redesign_amp_dict[k] = v

        result_name = f'redesign_panel_{redesign_num}'
        script_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'primer_agent_redesign_getPrimer.py')
        cmd = f'python {script_path} --first_design_file {first_design_file} --bad_primer_path {redesign_df_path} --file_name {result_name}'
        print("cmd: ", cmd)
        process = os.popen(cmd)  # return file
        output = process.read()
        print(output)
        process.close()
        result = []
        output_panel = '/reference_data/primer_result/output/'
        for ff in os.listdir(output_panel):
            if result_name in ff:
                result.append(os.path.join(output_panel, ff))

        response = f'According to the experimental result file you provided, we filtered out {redesign_num} pairs of primers that did not meet the requirements. We redesigned them and the design results are as follows:'
        responses = {'response': response, "operations": [], "data": {'redesign_result': result},
                     'state': 'stop', "primer_type": "redesign_primer_type",
                    "stage": stage}
        return responses

def run_cmd(cmd):

    import subprocess
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    while True:
        output = process.stdout.readline()
        if output == '' and process.poll() is not None:
            break
        if output:
            print(output.strip())
    err = process.stderr.read()
    if err:
        print("Error:\n", err)

def extract_json_response(response):
    # Extract JSON content from the response
    experiment_info_dict = response.choices[0].message.content
    match = re.search(r'```json\n(.*?)\n```', experiment_info_dict, re.DOTALL)
    if match:
        experiment_info_json_str = match.group(1)
        experiment_info_json_data = json.loads(experiment_info_json_str)
    else:
        print("GPT constructs JSON error, the situation is as follows：\n",experiment_info_dict)
        exit()
    return experiment_info_json_data

def primer3_design(fa_path_list, parameters, output_path, name_list=[[1]]):
    result = []

    print(fa_path_list)
    logging.info(f"Current fna path: {fa_path_list}")
    print(parameters)
    logging.info(f"Current primer3 parameters: {parameters}")
    #print([type(y) for y in parameters])
    
    for name_info, fa_path in zip(name_list, fa_path_list):
        #print(name_info)
        #continue
        name = name_info
        
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
        logging.info(f"seqs_length:{i_length}")


        primer_num = 100
        filename = os.path.basename(fa_path).split('.')[0]
        logging.info(f"filename：{filename}")
        if i_length<300:
            if int(parameters[-1])+2>=int(i_length*0.1):
                PRIMER_PRODUCT_SIZE_RANGE = f"30-{int(i_length*0.5)}"
            else:
                PRIMER_PRODUCT_SIZE_RANGE = f"{int(i_length*0.1)}-{int(i_length*0.2)}"

            logging.info(f"core gene small, dynamic PRIMER_PRODUCT_SIZE_RANGE, {int(i_length*0.3)}-{int(i_length*0.6)}")
            primers = primer3.bindings.design_primers(
                        seq_args={
                            "SEQUENCE_ID": filename,
                            "SEQUENCE_TEMPLATE": seqs,
                            "SEQUENCE_INCLUDED_REGION": [0, int(i_length)],
                            },
                        global_args={
                            "PRIMER_OPT_SIZE": 18,#parameters[-1],
                            "PRIMER_MIN_SIZE": 15,#int(parameters[-1])-2,
                            "PRIMER_MAX_SIZE": 30,#int(parameters[-1])+2,
                            "PRIMER_OPT_TM": parameters[2],
                            "PRIMER_MIN_TM": int(parameters[2])-3,
                            "PRIMER_MAX_TM": int(parameters[2])+3,
                            "PRIMER_PAIR_MAX_DIFF_TM": 2.0,
                            "PRIMER_MIN_GC": float(parameters[3])*100,
                            "PRIMER_MAX_GC": float(parameters[4])*100,
                            "PRIMER_MAX_POLY_X": 3,
                            "PRIMER_PICK_RIGHT_PRIMER": 1,
                            "PRIMER_PICK_LEFT_PRIMER": 1,
                            "PRIMER_PICK_INTERNAL_OLIGO": 0,
                            "PRIMER_PRODUCT_SIZE_RANGE": PRIMER_PRODUCT_SIZE_RANGE,#f"{int(i_length*0.1)}-{int(i_length*0.2)}",
                            "PRIMER_GC_CLAMP": 1,
                            "PRIMER_NUM_RETURN": primer_num,
                            },
                        )
        else:
            logging.info(f"传入的参数：{parameters}")
            logging.info(f"参数类型：{[type(y) for y in parameters]}")
            primers = primer3.bindings.design_primers(
                    seq_args={
                        "SEQUENCE_ID": filename,
                        "SEQUENCE_TEMPLATE": seqs,
                        "SEQUENCE_INCLUDED_REGION": [0, int(i_length)],
                    },
                    global_args={
                        "PRIMER_OPT_SIZE": parameters[-1],
                        "PRIMER_MIN_SIZE": int(parameters[-1])-2,
                        "PRIMER_MAX_SIZE": int(parameters[-1])+2,
                        "PRIMER_OPT_TM": parameters[2],
                        "PRIMER_MIN_TM": int(parameters[2])-3,
                        "PRIMER_MAX_TM": int(parameters[2])+3,
                        "PRIMER_PAIR_MAX_DIFF_TM": 2.0,
                        "PRIMER_MIN_GC": float(parameters[3])*100,
                        "PRIMER_MAX_GC": float(parameters[4])*100,
                        "PRIMER_MAX_POLY_X": 3,
                        "PRIMER_PICK_RIGHT_PRIMER": 1,
                        "PRIMER_PICK_LEFT_PRIMER": 1,
                        "PRIMER_PICK_INTERNAL_OLIGO": 0,
                        "PRIMER_PRODUCT_SIZE_RANGE": f"{int(parameters[0])}-{int(parameters[1])}",
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
        
        logging.info(f"--------------------------------primer3 result:\n{primers}")
        for ii in range(int(primers["PRIMER_PAIR_NUM_RETURNED"])):
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
        
        logging.info(f'Primer3 output:\n{raw_all_param_df}') 
        all_param_df = raw_all_param_df[(raw_all_param_df['PRIMER_LEFT_HAIRPIN_TH']<=0.0) & (raw_all_param_df["PRIMER_RIGHT_HAIRPIN_TH"]<=0.0)]
        logging.info(len(all_param_df))

        if len(all_param_df)==0:
            left_sorted_df = raw_all_param_df.sort_values(by='PRIMER_LEFT_HAIRPIN_TH').head(50)
            logging.info(f'left----------\n{left_sorted_df}')
            right_sorted_df = raw_all_param_df.sort_values(by='PRIMER_RIGHT_HAIRPIN_TH').head(50)
            logging.info(f'right----------\n{right_sorted_df}')
            intersection_indices = list(set(left_sorted_df.index) & set(right_sorted_df.index))
            all_param_df = raw_all_param_df.loc[intersection_indices]
        logging.info(f'HAIRPIN filter results:\n{all_param_df}')

        filenameTSV = filename + ".csv"
        # Save the dictionary as a .tsv file
        res = os.path.join(output_path, filenameTSV)
        top_100_df = all_param_df.head(100)
        order_top100_df = top_100_df.sort_values('PRIMER_PAIR_PENALTY')#.reset_index(drop=True, inplace=True)#,inplace=True)
        order_top100_df.reset_index(drop=True, inplace=True)
        order_top100_df.to_csv(res,index=False)
        logging.info("Order PENALTY&HAIRPIN filter df:\n{order_top100_df}")
        logging.info("Order PENALTY&HAIRPIN filter df index:\n{order_top100_df.index}")
        
        #排查二聚体情况
        top_50_index = detect_dimer(order_top100_df)
        logging.info('+'*25)
        logging.info(f"Dimer result index:\n{top_50_index}")
        top50_df = order_top100_df.loc[top_50_index]
        logging.info(top50_df)
        #排查mismatch情况
        mismatch_list = detect_mismatch(top50_df, seqs)
        filtermismatch_top50_df = top50_df[[bool(value) for value in mismatch_list]]
        top50_df = filtermismatch_top50_df
        logging.info('-'*25)
        logging.info(f'free energy filter df：\n{top50_df}')
        logging.info(f'best primer:\n{top50_df.head(1)}')
        logging.info(f'left primer:\n{top50_df.head(1)["PRIMER_LEFT_SEQUENCE"]}')
        logging.info(f'right primer:\n{top50_df.head(1)["PRIMER_RIGHT_SEQUENCE"]}')
        logging.info('-'*25)
        result_primer = top50_df.head(5)
        result_primer.insert(0,'Organism Name',name)
        filenameTSV = name.replace(' ','_')+'_primer_' + ".csv"
        res = os.path.join(output_path, filenameTSV)
        result_primer.to_csv(res,index=False)

        result.append(res)
    return result


def analysis_target_effiency(summary_df, efficiency_df):
    summary_df_summary = summary_df.parse(sheet_name='summary')
    sample_name_list = summary_df_summary[summary_df_summary["map_target_rate(%)"] < 90]["sample_name"].to_list()
    all_bad_primer_dict = {}
    for name in sample_name_list:
        # print(name)
        name_efficiency_df = efficiency_df.parse(sheet_name=name)
        filter_name_efficiency_df = name_efficiency_df[name_efficiency_df["primer_efficiency(%)"] < 50]  # .to_list()
        for index, info_ in filter_name_efficiency_df.iterrows():
            if info_["amp_name"] not in all_bad_primer_dict:
                all_bad_primer_dict[info_["amp_name"]] = []
            all_bad_primer_dict[info_["amp_name"]].append(info_["primer_efficiency(%)"])
    # print('all_bad_primer_dict',all_bad_primer_dict)
    # 如果3次重复只有有两次都符合，则也算符合条件
    new_primer_dict = {}
    for k in all_bad_primer_dict.keys():
        if len(all_bad_primer_dict[k]) == 1:
            continue
        else:
            new_primer_dict[k] = all_bad_primer_dict[k]
    # print('new_primer_dict', new_primer_dict)
    effiency_bad_primer_result = {key: min(value) for key, value in new_primer_dict.items()}
    # print('effiency_bad_primer_result', effiency_bad_primer_result)

    return effiency_bad_primer_result


def analysis_uniformity(summary_df):
    summary_df_amp = summary_df.parse(sheet_name='summary')
    sample_name_list = summary_df_amp[summary_df_amp["uniformity(>0.1x)"] < 90]["sample_name"].to_list()
    name_efficiency_df = summary_df.parse(sheet_name="amp_unique_region_depth")
    all_bad_primer_list = []
    for name in sample_name_list:
        # print(name)
        uni = name_efficiency_df[name].to_list()
        mean = sum(uni) / len(uni)
        filter_name_efficiency_df = name_efficiency_df[name_efficiency_df[name] < mean * 0.1]
        # print('filter_name_efficiency_df',filter_name_efficiency_df)
        bad_primer = filter_name_efficiency_df.loc[:, ['name']+sample_name_list]
        all_bad_primer_list.append(bad_primer)

    all_bad_primer_df = pd.concat(all_bad_primer_list)
    all_bad_primer_df.drop_duplicates()
    # exit()

    #     # continue
    #     filter_name_efficiency_df = name_efficiency_df[name_efficiency_df["primer_efficiency(%)"]<50]
    #     for index, info_ in filter_name_efficiency_df.iterrows():
    #         if info_["amp_name"] not in all_bad_primer_dict:
    #             all_bad_primer_dict[info_["amp_name"]] = []
    #         all_bad_primer_dict[info_["amp_name"]].append(info_["primer_efficiency(%)"])
    # # exit()
    # effiency_bad_primer_result = {key: min(value) for key, value in all_bad_primer_dict.items()}
    # print('analysis_uniformity', all_bad_primer_df)
    return all_bad_primer_df


def get_bad_primer(report_file_list):
    for file_path in report_file_list:
        if 'summary' in file_path:
            if 'xlsx' in file_path:
                summary_df = pd.ExcelFile(file_path)
            elif 'csv' in file_path:
                summary_df = pd.read_csv(file_path)
            else:
                raise ValueError("[File error]: Only support .csv or .xlsx file.")

        elif 'efficiency' in file_path:
            if 'xlsx' in file_path:
                efficiency_df = pd.ExcelFile(file_path)
            elif 'csv' in file_path:
                efficiency_df = pd.read_csv(file_path)
            else:
                raise ValueError("[File error]: Only support .csv or .xlsx file.")

        elif 'dimer' in file_path:
            if 'xlsx' in file_path:
                dimer_df = pd.read_excel(file_path)
            elif 'csv' in file_path:
                dimer_df = pd.read_csv(file_path)
            else:
                raise ValueError("[File error]: Only support .csv or .xlsx file.")
        else:
            raise ValueError(
                "[File error]: Only summary, efficiency, and dimer files are accepted, and their names must contain the above fields.")
    uniformity_bad_primer_df = analysis_uniformity(summary_df)
    uniformity_bad_primer_list = list(uniformity_bad_primer_df['name'])
    eiffency_bad_primer_result = analysis_target_effiency(summary_df, efficiency_df)
    eiffency_bad_primer_list = eiffency_bad_primer_result.keys()

    #dimer analysis
    return uniformity_bad_primer_list,eiffency_bad_primer_list


#3'端的50%序列内，不能错配
def detect_mismatch(primer3_df, target_cds):
    print('+++',primer3_df)
    left_primer = primer3_df['PRIMER_LEFT_SEQUENCE'].tolist()
    left_primer_start = primer3_df['PRIMER_LEFT'].apply(lambda x: x[0]).tolist()
    left_primer_length = primer3_df['PRIMER_LEFT'].apply(lambda x: x[1]).tolist()
    right_primer = primer3_df['PRIMER_RIGHT_SEQUENCE'].tolist()
    right_primer_start = primer3_df['PRIMER_RIGHT'].apply(lambda x: x[0]).tolist()
    right_primer_length = primer3_df['PRIMER_RIGHT'].apply(lambda x: x[1]).tolist()
    print('left primer: ',left_primer)
    print('left start: ',left_primer_start)
    print('left length: ',left_primer_length)
    print('right primer: ',right_primer)
    print('right start: ',right_primer_start)
    print('right length: ',right_primer_length)
    left_primer_list = [seq[-3:] for seq in left_primer]
    right_primer_list = [seq[-1:-4:-1] for seq in right_primer]
    print('3‘right: ', right_primer_list)
    print('3‘left: ', left_primer_list)

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
        print('左引物：',x)
        print('左模板：',z)
        print('右引物：',y)
        print('右模板：',u)
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
def calculate_optimal_primer(blast_result_path):
    #with open('blast_result_path','r') as f:
    blast_df = pd.read_csv(blast_result_path, sep='\t', header=None, 
                                   names=['qseqid', 'sseqid', 'qstart', 'qend', 'sstart', 'send', 'length', 'pident', 'mismatch', 'gapopen', 
                                          'evalue', 'bitscore', 'staxids', 'sscinames', 'qlen'])
    print('blast特异性比对结果：\n', blast_df)
    exit()
    

def get_best_primer(filter_primer_df, name):
    filter_primer_path = '/reference_data/primer_blast_result'
    print(filter_primer_df)
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
    db = '/home/YJhou/workspace/pcr_paper/search_base/blast_db/refseq'
    cmd = getglCMD(primer_file_path, ofile)
    print(cmd)
    os.system(cmd)
    print(f'完成引物blast，存于{ofile}')
    calculate_optimal_primer(ofile)

