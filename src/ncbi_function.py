import re
import os
import json
import sys
import traceback
import time
import subprocess
import zipfile
import shutil
import glob
import pandas as pd
from Bio import SeqIO, SeqRecord
from Bio.Seq import Seq
import logging

logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

def commend_seqs_filte(information_list):
    commend_seqs_dict_list=[]
    for infor_dict in information_list:
        if infor_dict["assembly_level"] != "Complete Genome":
            continue
        if infor_dict["representative"] == "reference genome":
            print('get reference genome')
            commend_seqs_dict_list.append(infor_dict)
    
    # When there is no reference genome, get Complete genome
    if len(commend_seqs_dict_list) == 0:
        print('No reference genome')
        for infor_dict in information_list:
            if infor_dict["assembly_level"] != "Complete Genome":
                continue
            else:
                commend_seqs_dict_list.append(infor_dict)
        print(commend_seqs_dict_list)                               
        sorted_commend_seqs_dict_list = sorted(commend_seqs_dict_list, key=lambda x: x['total_sequence_length'], reverse=True) # Complete genome needs to be sorted by total_sequence_length
        commend_seqs_dict_list = sorted_commend_seqs_dict_list[:20] # Recommend TOP20

    return commend_seqs_dict_list # Return list if complete exists, empty list if no complete

def get_gca_gene(data,data_info_list,genome_dict):
    for data_info in data:
        if "contigN50" not in data_info["assemblyStats"]:
            continue
        if "assemblyLevel" not in data_info["assemblyInfo"]:
            continue

        accession = data_info["accession"]
        assemblyLevel = data_info["assemblyInfo"]["assemblyLevel"]
        contigN50 = data_info["assemblyStats"]["contigN50"]
        data_info_list.append((accession, assemblyLevel, contigN50))

        if assemblyLevel == "Complete Genome":
            genome_dict['complete'][accession] = (assemblyLevel, contigN50)
            complete_flag = True
        else:
            genome_dict['noncomplete'][accession] = (assemblyLevel, contigN50)
            complete_flag = False

    return data_info_list, complete_flag, genome_dict

def get_assembly_data_report_info(data_info):
    accession = data_info["accession"]
    assemblyLevel = data_info["assemblyInfo"]["assemblyLevel"]
    if "scaffoldN50" not in data_info["assemblyStats"]:
        scaffoldN50 = 'None'
    else:
        scaffoldN50 = data_info["assemblyStats"]["scaffoldN50"]
    if "annotationInfo" not in data_info:
        releaseDate = ""
    else:
        releaseDate = data_info["annotationInfo"]["releaseDate"]
    assemblyName = data_info["assemblyInfo"]["assemblyName"]
    organismName = data_info["organism"]["organismName"]
    if "infraspecificNames" not in data_info["organism"]:
        strainName = 'None'
    else:
        if "strain" not in data_info["organism"]["infraspecificNames"]:
            strainName = 'None'
        else:
            strainName = data_info["organism"]["infraspecificNames"]["strain"]
    if "checkmInfo" not in data_info:
        completeness=""
        contamination=""
    else:
        completeness = data_info["checkmInfo"]["completeness"]
        if "contamination" not in data_info["checkmInfo"]:
            contamination = 'None'
        else:
            contamination = data_info["checkmInfo"]["contamination"]
    return (accession, assemblyLevel, scaffoldN50, releaseDate, assemblyName, organismName, strainName, completeness, contamination)

def get_best_genome(path):
    # For GCF, consider scaffoldN50 and assemblyLevel - maximize scaffoldN50 and require Complete Genome
    # For GCA, consider contigN50 and assemblyLevel - maximize contigN50 with or without complete
    with open(path, 'r') as f:
        lines = f.readlines()
        data = [json.loads(line) for line in lines]

    data_info_list = []
    genome_dict = {'complete':{}, 'noncomplete':{}}
    complete_flag = False
    print("start to find GCF")
    for data_info in data:
        if "refseqCategory" in data_info["assemblyInfo"]:
            if 'reference' in data_info["assemblyInfo"]["refseqCategory"] or 'representative' in data_info["assemblyInfo"]["refseqCategory"]:
                # Return directly
                information_ = get_assembly_data_report_info(data_info)
                return (information_[0],(information_[1:]))
            
        else:
            if 'GCF' not in data_info["accession"]:
                continue
            information_ = get_assembly_data_report_info(data_info)
            data_info_list.append(information_)
            if information_[1] == "Complete Genome":
                genome_dict['complete'][information_[0]] = information_[1:]#(assemblyLevel, scaffoldN50,releaseDate, assemblyName, organismName, strainName, completeness, contamination)
                complete_flag = True
                print(genome_dict)
            else:
                genome_dict['noncomplete'][information_[0]] = information_[1:]#(assemblyLevel, scaffoldN50,releaseDate, assemblyName, organismName, strainName, completeness, contamination)
    
    print('----------extract information:\n', data_info_list)
    if data_info_list == []:
        print("NO GCF, refind GCA")
        data_info_list,complete_flag,genome_dict = get_gca_gene(data,data_info_list,genome_dict)
        
    if complete_flag:
        best_complete_genome = sorted(genome_dict['complete'].items(), key=lambda x: x[1][1], reverse=True)[0] if genome_dict['complete'] else None
        return best_complete_genome

    else:
        best_noncomplete_genome = sorted(genome_dict['noncomplete'].items(), key=lambda x: x[1][1], reverse=True)[0] if genome_dict['noncomplete'] else None
        return best_noncomplete_genome
    

def extract_gene(fna_content, gene_name):
    pattern = r">(.*?)\n([ATGC\n]+)"
    matches = re.findall(pattern, fna_content, re.MULTILINE)
    
    for header, sequence in matches:
        if f"[gene={gene_name}]" in header:
            seq = sequence.replace('\n','')
            return f">{header}\n{seq}"
    return ""

def getCMD(kmer_file, klen, db, ofile, max_target_seqs=100000000,cov_hsp_perc=50, ident=15, subtree_path=""):
    if klen <= 50:
        cmd = 'time blastn -task blastn-short'
    elif klen <= 500:
        cmd = 'time blastn -task blastn'
    else:
        cmd = 'time blastn -task megablast'
    cmd += f' -query {kmer_file} -db {db} -max_target_seqs {max_target_seqs} -num_threads 100 -perc_identity {ident} -qcov_hsp_perc {cov_hsp_perc} -evalue 10000  -outfmt "6 qseqid sseqid qstart qend sstart send length pident mismatch gapopen evalue bitscore staxids sscinames qlen"  > {ofile}'
    print(cmd)
    return cmd

# Get taxon subtree taxid using taxonkit
def get_taxon_subtree(taxid):
    # Call taxonkit list command
    result = subprocess.run(["/app/bin/taxonkit", "list", "--ids", str(taxid)],
            capture_output=True,
            text=True,
            check=True)
    # Process output
    output = result.stdout.split('\n')
    subtree = [int(i.replace(' ','')) for i in output if i!='']
    return subtree

# Different species and strains
def filter_cds_blast_output(cds_file_path, file_path, subtree, identity_threshold=30, type_='general_speices_blast'):
    cds_count = 0
    sequences = {}
    for record in SeqIO.parse(cds_file_path, "fasta"):
        record_id = str(record.id)
        record_seq = str(record.seq)
        cds_count += 1
        sequences[record_id] = len(record_seq)

    # Parse Blast results
    column_names = ["qseqid", "sseqid", "qstart", "qend", "sstart", "send", "length", "pident", "mismatch", "gapopen", "evalue", "bitscore", "staxids", "sscinames", "qlen"]
    blast_output = pd.read_csv(file_path, sep="\t", header=None, names=column_names)
    
    # Make distinction
    if type_ == "general_speices_blast":
        filter_blast_list = []
        # Filter same species
        filtered_subtree_blast = blast_output[~blast_output["staxids"].isin(subtree)]
        
        qseqid_species_counts = filtered_subtree_blast.groupby(['qseqid', 'staxids']).size().reset_index(name='count')
        print("qseqid_species_counts:",qseqid_species_counts)
        
        qseqid_total_counts = qseqid_species_counts.groupby('qseqid')['count'].sum().reset_index(name='total_count')
        print("qseqid_total_counts:",qseqid_total_counts)
        
        sorted_qseqid_total_counts = qseqid_total_counts.sort_values(by='total_count', ascending=True)
        print("sorted_qseqid_total_counts", sorted_qseqid_total_counts)
        
        count_dict = dict(zip(sorted_qseqid_total_counts['qseqid'], sorted_qseqid_total_counts['total_count']))
        
        with open(file_path+'_blast_count.json','w') as f:
            json.dump(count_dict, f, indent=4)
        print("count_dict_result path:",file_path+'_blast_count.json')
        
        if cds_count>len(count_dict):
            print(cds_count)
            print(len(count_dict))
            print('Best specific CDS exists')
            filter_blast_list = filtered_subtree_blast["qseqid"].unique().tolist()
            print('Number to filter:',len(filter_blast_list))
            return filter_blast_list, ''
            
        else:
            filter_dict={}
            min_total_count = sorted_qseqid_total_counts['total_count'].min()
            filter_blast_list = sorted_qseqid_total_counts[sorted_qseqid_total_counts['total_count'] == min_total_count]['qseqid'].tolist()
            for id_ in filter_blast_list:
                filter_dict[id_]=sequences[id_]
            max_seq_id = max(filter_dict, key=filter_dict.get)
            max_value = filter_dict[max_seq_id]
            print(max_seq_id, max_value)
            return '', max_seq_id

def get_cds_blast_result(cds_name_list, inputfile, outputfile):
    with open(inputfile, 'r') as input_handle, open(outputfile, 'w') as output_handle:
        for record in SeqIO.parse(input_handle, 'fasta'):
            name = record.id
            if name not in cds_name_list:
                SeqIO.write(record, output_handle, 'fasta')

def get_best_core_gene_length(outputfile, best_core_gene_path):
    sequences = list(SeqIO.parse(outputfile, "fasta")) 
    longest_seq = max(sequences, key=lambda record: len(record.seq))
    with open(best_core_gene_path, "w") as output_file:
        SeqIO.write(longest_seq, output_file, "fasta")

def get_best_core_gene_match_count(cds_file_path, match_count, best_core_gene_path):
    for record in SeqIO.parse(cds_file_path, "fasta"):
        if record.id == match_count:
            with open(best_core_gene_path, "w") as output_file:
                SeqIO.write(record, output_file, "fasta")

def genKmer(seq,  kmer_len):
    kmer = set()
    for i in range(0, len(seq)):
        if i + kmer_len > len(seq):     # Skip if last window cannot reach kmer length
            break
        kmer.add(seq[i:(i + kmer_len)])
    return kmer

def getkmer(infile, outfile):
    try:
        seq_list = []
        with open(outfile, 'w') as out_fasta:
            average_len = 150
            for r in SeqIO.parse(infile, 'fasta'):     # Read genome file, need to iterate since some genomes are incomplete with many sequences
                seq = str(r.seq)
                seq_list.append((seq, r.id))

            # Sort by seq length and select top 3 longest elements
            top_3_longest = sorted(seq_list, key=lambda x: len(x[0]),reverse=True)[:3]

            for seq in top_3_longest:
                print(len(seq[0]))
                kmers = genKmer(seq[0], average_len)
                name = seq[1]
                for idx, kmer in enumerate(kmers):
                    if 'lcl' in name:
                        record = SeqRecord.SeqRecord(Seq(kmer), id=f"{name}_cds_kmer_{idx}", description="")
                    else:
                        record = SeqRecord.SeqRecord(Seq(kmer), id=f"lcl|{name}_cds_kmer_{idx}", description="")

                    SeqIO.write(record, out_fasta, 'fasta')
    except:
        traceback.print_exc()
        sys.exit(1)

def get_core_gene(file_path, taxid=[],type_num="single", type_cds="",type_experiment=""):
    print('***************Start getting core gene***************')
    print(file_path)
    print(type_experiment)
    
    if type_cds == "have_cds":
        kmer_file = file_path[0]
        
    else:
        # When no cds or gtf file, need to do kmer processing
        infile = file_path[0]
        kmer_file = ""
        outfile = infile+'_chunk_kmer.fasta'
        getkmer(infile, outfile)
        
    # kmer length
    start_time = time.time()
    len_list=[]
    for r in SeqIO.parse(kmer_file, 'fasta'):
        seqlen = len(str(r.seq))
        len_list.append(seqlen)
    klen = int(sum(len_list)/len(len_list))
    print('Average fragment length:', klen)
        
    # blast database
    db = "/reference_data/blastdb_refs/refseq"
    
    # blast output
    path_lst = kmer_file.split('/')
    output_path = '/'.join(path_lst[:-1])
    name = path_lst[-1]
    blast_name = name.split('.fna')[0]+'_blast_output'
    blast_ofile = output_path+'/'+blast_name

    print('kmer_file:',kmer_file)
    cmd = getCMD(kmer_file, int(klen), db, blast_ofile, max_target_seqs=100000000,cov_hsp_perc=30, ident=30 )
    print('blast command:', cmd)
    os.system(cmd)
    print('blast time used:', time.time()-start_time)
    print("[cds blast output] ", blast_ofile)    
    
    start_time = time.time()
    subtree = get_taxon_subtree(taxid[0])
    filter_blast_list, match_count = filter_cds_blast_output(kmer_file, blast_ofile, subtree, 30)
    print("Time used for filtering by species branch and identity:",time.time()-start_time)
    print("Number to filter by species branch:\n", len(filter_blast_list))
    
    if match_count == '':
        # Filter cds
        start_time = time.time()
        outputfile = output_path+'/'+blast_name+'_core_cds.fna'
        get_cds_blast_result(filter_blast_list, kmer_file, outputfile)
        print("[output csd path] ",outputfile)

        best_core_gene_path = "/reference_data/target/"+name.split('.fna')[0] + '_core_cds.fna'
        
        # from length
        get_best_core_gene_length(outputfile,best_core_gene_path)
        
    else:
        best_core_gene_path = "/reference_data/target/"+name.split('.fna')[0] + '_core_cds.fna'
        get_best_core_gene_match_count(kmer_file, match_count, best_core_gene_path)
    print("[best core gene path] ",best_core_gene_path)
    print("Time used to get core gene:",time.time()-start_time) 

    return best_core_gene_path

def ncbi_search(term, db_name, accession_id=None, target_type=None, search_type=None):
    # Task is to download whole genome and CDS files
    if db_name == "genome":
        # target gene scenario
        if target_type=="target":
            temp_path = "target_zip_file"
        
        # non-target gene scenario
        elif target_type=="non_target":
            temp_path = "non_target_zip_file"
        
        # taxid
        id = accession_id
        
        datasets_limitation = ["--assembly-source RefSeq --reference", "--assembly-source RefSeq", "--reference'",""]
        failed_list = []
        
        # cds download
        if search_type == "species_cds":
            term = term.replace('/','_').replace("(","_").replace(")","_").replace(' ','_')
            make_file_str = f"/reference_data/{temp_path}/{term}"
            print('-------------------------------:',term)
            for idx in range(len(datasets_limitation)):
                # Record error sequences
                if len(failed_list) != 0:
                    print(f'*****************{failed_list} download error, try change limitation *****************')
                    # All searches failed
                if idx == len(datasets_limitation):
                    return f'{term} gene search failed.',''
                # Select search mode
                datasets_limi_str_ = datasets_limitation[idx]
                # Current sequence already downloaded
                if os.path.exists(f"{make_file_str}/{term}.zip"):
                    file_path = f"{make_file_str}/{term}.zip"
                # Current sequence not downloaded, start downloading
                else:
                    print("Start creating folder")
                    os.system(f"mkdir {make_file_str}")
                    url = f"/app/bin/datasets download genome accession {id} {datasets_limi_str_} --include cds --filename {make_file_str}/{term}.zip"
                    print(url)
                    os.system(url)
                    if not os.path.exists(f"{make_file_str}/{term}.zip"):
                        failed_list += datasets_limitation[idx]
                        continue
                    else:
                        print('--------Finished genome download')
                        file_path = f"{make_file_str}/{term}.zip"
                        break

        elif search_type == "dna_species_gene" or search_type == "rna_species_gene":
            term = term.replace('/','_').replace("(","_").replace(")","_").replace(' ','_')
            make_file_str = f"/reference_data/{temp_path}/{term}"
            if search_type == "dna_species_gene":
                url = f"/app/bin/datasets download gene taxon {id} --include gene --filename {make_file_str}/{term}.zip"
            else:
                url = f"/app/bin/datasets download genome taxon {id} --filename {make_file_str}/{term}.zip"
            print(url)
            
            os.system(url)
            if not os.path.exists(f"{make_file_str}/{term}.zip"):
                print('----------Current file has no gene file, try downloading genome file')
                for idx in range(len(datasets_limitation)):
                    # Record error sequences
                    if len(failed_list) != 0:
                        print(f'*****************{failed_list} download error, try change limitation *****************')
                    # All searches failed
                    if idx == len(datasets_limitation):
                        return f'{term} gene search failed.',''
                    # Select search mode
                    datasets_limi_str_ = datasets_limitation[idx]
                    # Current sequence already downloaded
                    if os.path.exists(f"{make_file_str}/{term}.zip"):
                        file_path = f"{make_file_str}/{term}.zip"
                    else:
                        os.system(f"mkdir {make_file_str}")
                        url = f"/app/bin/datasets download genome taxon {id} {datasets_limi_str_} --filename {make_file_str}/{term}.zip"
                        os.system(url)
                        if not os.path.exists(f"{make_file_str}/{term}.zip"):
                            failed_list += datasets_limitation[idx]
                            continue
                        # Download successful, break download mode loop
                        else:
                            print('----------genome download complete')
                            file_path = f"{make_file_str}/{term}.zip"
                            break
            else:
                print('----------gene download complete')

        else:
            print(f"Search type problem: {search_type}")
            pass

        # Extract zip to specified "make file str" path
        with zipfile.ZipFile(file_path, 'r') as zip_ref:
            zip_ref.extractall(make_file_str)
        
        # Locate best genome by reading report
        best_gene_info = get_best_genome(f"{make_file_str}/ncbi_dataset/data/assembly_data_report.jsonl")
        # Get best genome
        logging.info(f"best_gene_info: {best_gene_info}")
        path = f"{make_file_str}/ncbi_dataset/data/{best_gene_info[0]}/*"
        logging.info(f"genomic_source_path: {path}")
        genomic_source = glob.glob(f"{make_file_str}/ncbi_dataset/data/{best_gene_info[0]}/*")
        logging.info(f"genomic_source_path: {genomic_source}")
        genomic_source_path = genomic_source[0]
        logging.info(f"genomic_source: {genomic_source_path}")

        # Genome copy destination
        genomic_copy_path = f"/reference_data/{temp_path}/{term.replace(' ', '_')}.fna"
        print(f"Source file: {genomic_source_path}")
        print(f"File exists: {os.path.exists(genomic_source_path)}")
        print(f"Destination file: {genomic_copy_path}")
        
        # Complete genome file copy task for later use
        shutil.copy(genomic_source_path, genomic_copy_path)
        print(f"finish{genomic_source_path} --> {genomic_copy_path}")

        # If user only needs species level content, return whole genome file now
        # TODO: Need to further refine task types
        if search_type in ["species_genome", "dna_species_gene", "rna_species_gene"]:
            return [genomic_copy_path]

        # Non-species level, still need to mine cds gene information
        make_cds_file_str = f"/reference_data/{temp_path}/{term}/cds"
        # Check if cds file downloaded, download cds file using accession number from previous step
        if not os.path.exists(make_cds_file_str):
            os.system(f"mkdir {make_cds_file_str}")

        target_gene_list = []
        # File transfer
        if search_type=="species_cds":
            # cds file path
            cds_path = genomic_source_path
            
            # New cds file path for transfer
            if not os.path.exists(f"/reference_data/{temp_path}/{term}/cds/ncbi_dataset/data/{best_gene_info[0]}/"):
                os.makedirs(f"/reference_data/{temp_path}/{term}/cds/ncbi_dataset/data/{best_gene_info[0]}/")
            new_cds_path = f"/reference_data/{temp_path}/{term}/cds/ncbi_dataset/data/{best_gene_info[0]}/{term}_cds_from_genomic.fna"
        
        else:
            if not os.path.exists(make_cds_file_str):
                os.system(f"mkdir {make_cds_file_str}")
                get_cds_url = f"/app/bin/datasets download genome accession {best_gene_info[0]} --include cds --filename {make_cds_file_str}/{term}_cds.zip"
                os.system(get_cds_url)
            # Extract cds file
            with zipfile.ZipFile(f"{make_cds_file_str}/{term}_cds.zip", 'r') as zip_ref:
                zip_ref.extractall(make_cds_file_str)
            
            # cds file path
            cds_path = f"/reference_data/{temp_path}/{term}/cds/ncbi_dataset/data/{best_gene_info[0]}/cds_from_genomic.fna"
            # New cds file path for transfer
            new_cds_path = f"/reference_data/{temp_path}/{term}/cds/ncbi_dataset/data/{best_gene_info[0]}/{term}_cds_from_genomic.fna"
        
        # File rename using name as distinction
        print(f"Path transfer: cp {cds_path} {new_cds_path}")
        os.system(f"cp {cds_path} {new_cds_path}")
        
        return [new_cds_path, best_gene_info]

# Sort function
def sort_genomes(genomes):
    # Define sorting keys (priority order)
    def sorting_key(genome):
        # Priority 1: representative is reference genome
        representative_priority = 0 if genome['representative'] == 'reference genome' else 1

        # Priority 2: assembly_level is Complete Genome
        assembly_level_priority = 0 if genome['assembly_level'] == 'Complete Genome' else 1

        # Priority 3: accession is GCF
        accession_priority = 0 if genome['accession'].startswith('GCF') else 1

        # Priority 4: scaffold_n50 and contig_n50
        # For GCF, sort by scaffold_n50 descending
        if genome['accession'].startswith('GCF'):
            n50_priority = -genome['scaffold_n50']
        else:  # For GCA, sort by contig_n50 descending
            n50_priority = -genome['contig_n50']

        # Return tuple, first element is main sorting criteria, second is secondary, etc.
        return (representative_priority, assembly_level_priority, accession_priority, n50_priority)

    # Sort dictionary list by defined sorting rules
    return sorted(genomes, key=sorting_key)