from Bio import Entrez
import sys
from typing import List, Optional
import os
from Bio import SeqIO
import xml.etree.ElementTree as ET
import pandas as pd
from io import StringIO
import subprocess
import zipfile
from Bio.SeqFeature import SeqFeature
import logging
from collections import defaultdict
import re
from llm_utils.utils import get_llm_chat_completion
import json

Entrez.email = ""
Entrez.api_key = ""
def retrieve_ncbi_ids(query_name, db_type="gene"):
    """
    Retrieve gene or protein IDs from NCBI based on name
    :param query_name: Search term (such as gene name or protein name)
    :param db_type: Specify the database to search ("gene", "protein", or "species")
    :param organism: Optional species filter (e.g. "Homo sapiens")
    :return: List of IDs
    """
    logging.info('start retrieve')
    
    if db_type == "gene":
       
        search_db = "gene"
        fetch_db = "gene"
        
        validation_keyword = "Annotation:"
        query_name = f"({query_name}) AND alive[prop]"
    else:
        search_db = "nucleotide"
        fetch_db = "nucleotide"
        
    logging.info(f"Using database: {search_db}, query: {query_name}")
    handle = Entrez.esearch(db=search_db, term=query_name, sort="relevance")
    record = Entrez.read(handle)
    logging.info(record)
    ids = record["IdList"]
    handle.close()
    logging.info(f"Found IDs: {ids}")
    
  
    new_ids = []
    for id in ids:
        try:
            handle = Entrez.efetch(db=fetch_db, id=id, rettype="gb", retmode="text")
        except Exception as e:
            logging.info(f"ID {id} efetch failed, skipping")
            continue
        record_text = handle.read()
        handle.close()
        
        
        if db_type == "gene":
           
            if validation_keyword in record_text:
                logging.info(f"Valid ID: {id}")
                new_ids.append(id)
        else:
            
            try:
               
                record = SeqIO.read(StringIO(record_text), "genbank")
                
                
                logging.info(record.annotations)
                molecule_type = record.annotations.get("molecule_type", "")
                logging.info(f"ID {id} sequence type: {molecule_type}")

                if not molecule_type:
                    locus_line = re.search(r'^LOCUS\s+\S+\s+\d+\s+bp\s+(\S+)', record_text, re.MULTILINE)
                    if locus_line:
                        molecule_type = locus_line.group(1).upper()
                
                logging.info(f"ID {id} sequence type: {molecule_type}")
                
                if db_type == "protein" or db_type == "species":
                    new_ids.append(id)
                else:
                    new_ids.append(id)
                    
            except Exception as e:
                logging.info(f"Error parsing ID {id}: {str(e)}")

    
    logging.info(f"Number of valid IDs after filtering: {len(new_ids)}")
    logging.info(f"Valid IDs: {new_ids}")
    return new_ids


def parse_gene_info2(genbank_str, db_type="gene"):
    """
    Parse GenBank string data and extract specified fields
    
    :param genbank_str: GenBank format data string
    :param db_type: Database type, 'gene' or 'nucleotide'/'protein'/'species'
    :return: Dictionary containing extracted fields or None (if parsing fails)
    """
    if db_type == "gene":
        gene_info = defaultdict(lambda: "N/A")
        qualifiers = defaultdict(dict)
        
        try:
            record = SeqIO.read(StringIO(genbank_str), "genbank")
        except Exception as e:
            import traceback
            logging.error(f"Parse failed: {e}\nContent: {genbank_str[:500]}")
            logging.error(f"Parse failed: {traceback.format_exc()}")
            return None

        gene_info["Accession"] = record.id
        gene_info["Description"] = record.description
        gene_info["Length(bp)"] = len(record.seq)
        gene_info["Org-ref_taxname"] = record.annotations.get("organism", "N/A")

        for feature in record.features:
            ft_type = feature.type.upper()
            qual = feature.qualifiers
            
            if ft_type == "SOURCE":
                gene_info["OrgName_lineage"] = qual.get("db_xref", ["N/A"])[0].split(":")[1] if "db_xref" in qual and qual["db_xref"] else "N/A"
                gene_info["Strain"] = qual.get("strain", ["N/A"])[0]
                
            elif ft_type == "GENE":
                gene_info["Gene-ref_locus"] = qual.get("gene", ["N/A"])[0]
                gene_info["Gene-track_geneid"] = qual.get("gene_id", ["N/A"])[0]
                
                synonyms = qual.get("gene_synonym", [])
                gene_info["Gene-ref_syn"] = "|".join(synonyms) if synonyms else "N/A"
                
                gene_info["Gene-ref_desc"] = qual.get("note", ["N/A"])[0]
                
            elif ft_type == "CDS":
                gene_info["ProteinID"] = qual.get("protein_id", ["N/A"])[0]
                gene_info["ProteinName"] = qual.get("product", ["N/A"])[0]
                
                gene_info["Seq-interval_from"] = feature.location.start
                gene_info["Seq-interval_to"] = feature.location.end
                
                gene_info["Gene-commentary_accession"] = qual.get("transcript_id", ["N/A"])[0]
                
            if "note" in qual and qual["note"]:
                gene_info["Gene-commentary_text"] = qual["note"][0]
        
        gene_info["Seq-interval_from"] = str(gene_info["Seq-interval_from"])
        gene_info["Seq-interval_to"] = str(gene_info["Seq-interval_to"])
        
        required_fields = [
            "Gene-track_geneid", "Gene-ref_locus", "Gene-ref_syn", 
            "Org-ref_taxname", "Gene-ref_desc", "Gene-commentary_text", 
            "OrgName_lineage", "Gene-commentary_accession", 
            "Seq-interval_from", "Seq-interval_to", "Strain"
        ]
        for field in required_fields:
            if field not in gene_info:
                gene_info[field] = "N/A"
        
        return dict(gene_info)
    
    else:
        nucleotide_info = {
            "Accession": "na",
            "Description": "na",
            "Length(bp)": "na",
            "Organism": "na",
            "DNA_sequence": "na",
            "Gene_name": "na",
            "Product": "na",
            "Protein_id": "na",
            "db_xref": "na",
            "References": "na"
        }

        try:
            record = SeqIO.read(StringIO(genbank_str), "genbank")

            nucleotide_info["Accession"] = record.id
            nucleotide_info["Description"] = record.description
            nucleotide_info["Organism"] = record.annotations.get("organism", "na")

            if hasattr(record, 'seq') and isinstance(record.seq, Seq):
                nucleotide_info["Length(bp)"] = str(len(record.seq))
            else:
                nucleotide_info["Length(bp)"] = "na"

            dna_sequence = "na"
            try:
                seq_str = str(record.seq)
                if len(seq_str) > 100:
                    dna_sequence = seq_str[:100] + "..."
                else:
                    dna_sequence = seq_str
            except:
                dna_sequence = "na"
            nucleotide_info["DNA_sequence"] = dna_sequence
            gene_names = []
            product_names = []
            protein_ids = []
            db_xrefs = []
            features = getattr(record, 'features', [])
            for feature in features:
                qualifiers = feature.qualifiers if hasattr(feature, 'qualifiers') else {}
                if "gene" in qualifiers:
                    gene_names.extend(qualifiers.get("gene", []))
                if "product" in qualifiers:
                    product_names.extend(qualifiers.get("product", []))
                if "protein_id" in qualifiers:
                    protein_ids.extend(qualifiers.get("protein_id", []))
                if "db_xref" in qualifiers:
                    db_xrefs.extend(qualifiers.get("db_xref", []))
            for key, values, field in zip(
                    ["Gene_name", "Product", "Protein_id", "db_xref"],
                    [gene_names, product_names, protein_ids, db_xrefs],
                    ["gene", "product", "protein_id", "db_xref"]
            ):
                if values:
                    unique_values = list(dict.fromkeys(values))
                    first_two = unique_values[:2]
                    nucleotide_info[key] = "|".join(first_two)
                else:
                    nucleotide_info[key] = "na"
            references = []
            for ref in record.annotations.get("references", []):
                this_ref = {
                    "title": getattr(ref, "title", "na"),
                    "authors": getattr(ref, "authors", "na"),
                    "journal": getattr(ref, "journal", "na")
                }
                references.append(this_ref)
            nucleotide_info["References"] = references if references else "na"
            return nucleotide_info

        except Exception as e:
            return nucleotide_info

def extract_json_response(response):
    experiment_info_dict = response.choices[0].message.content
    try:
        return json.loads(experiment_info_dict)
    except:
        match = re.search(r'```json\n(.*?)\n```', experiment_info_dict, re.DOTALL)
        if match:
            experiment_info_json_str = match.group(1)
            experiment_info_json_data = json.loads(experiment_info_json_str)
        else:
            logging.info(f"GPT JSON error: {experiment_info_dict}")
            exit()
        return experiment_info_json_data
def parse_gene_info(genbank_str, db_type="gene"):
    """
    Parse GenBank string data and extract specified fields

    :param genbank_str: GenBank format data string
    :param db_type: Database type, 'gene' or 'nucleotide'/'protein'/'species'
    :return: Dictionary containing extracted fields or None (if parsing fails)
    """
    if db_type == "gene":
        user_input = genbank_str
        logging.info(f'----user_input---: {user_input}')
        experiment_type_prompt = '''Given the genetic information stored in the string gb_str privided by user, extract the specified data fields and return the results as a Python dictionary.
        
        Fields to Extract (required_fields):
        Accession: Extract the value after "ID:" (e.g., 2064).
        Length(bp): Calculate the genomic length in base pairs using the chromosome coordinates in the Annotation field (e.g., 39728658 - 39688094 + 1).
        Description: Extract the full gene description from the "Name" field in the "Official Symbol" line (e.g., erb-b2 receptor tyrosine kinase 2 [Homo sapiens (human)]).
        Gene-track_geneid: Extract the MIM number (e.g., 164870).
        Gene-ref_locus: Extract the "Official Symbol" 
        Gene-ref_syn: Extract all entries from the "Other Aliases" field, separated by commas (e.g., CD340, HER-2, [...]).
        Org-ref_taxname: Extract the species name from the brackets in the "Name" field (e.g., Homo sapiens).
        Gene-ref_desc: Extract the gene description without the species bracket (e.g., erb-b2 receptor tyrosine kinase 2).
        Gene-commentary_text: Combine all entries from the "Other Designations" field, separated by semicolons (e.g., receptor tyrosine-protein kinase erbB-2; c-erb B2/neu protein; [...).
        OrgName_lineage: If present, extract the taxonomic lineage (e.g., "Homo sapiens, Primates, [...]"); mark as N/A if absent.
        Gene-commentary_accession: Use the same value as Accession.
        Seq-interval_from and Seq-interval_to: Extract the raw start and end coordinates from the Annotation field (e.g., 39688094 and 39728658).
        Strain: Mark as N/A (no strain information is provided in the input).
        Output Format:
        Return a JSON with all keys from required_fields, formatted as follows:
        {  
        "Accession": "...",  
        "Length(bp)": ...,  
        "Description": "...",  
        "Gene-track_geneid":"",
        "Gene-ref_locus":"",
        "Gene-ref_syn":"",
        "Org-ref_taxname":"..",
        "Gene-ref_desc":"..",
        "Gene-commentary_text":"..",
        "OrgName_lineage":"..",
        "Gene-commentary_accession":"..",
        "Seq-interval_from and Seq-interval_to":"..",
        "Strain":"...",
        }
        Notes:
        If a field is not explicitly present in the input string, use N/A.
        Ensure numerical values (e.g., Length(bp) or coordinates) are stored as integers or strings as appropriate.
        must return a JSON
        '''
        experiment_info_response = get_llm_chat_completion(
            messages=[
                {"role": "system", "content": experiment_type_prompt},
                {"role": "user", "content": user_input}
            ],
            response_format={"type": "json_object"},
        )
        gene_info = extract_json_response(experiment_info_response)
        logging.info(f'----gene_info---: {gene_info}')

        return gene_info

    else:
        nucleotide_info = {
            "Accession": "na",
            "Description": "na",
            "Length(bp)": "na",
            "Organism": "na",
            "DNA_sequence": "na",
            "Gene_name": "na",
            "Product": "na",
            "Protein_id": "na",
            "db_xref": "na",
            "References": "na"
        }

        try:
            record = SeqIO.read(StringIO(genbank_str), "genbank")

            nucleotide_info["Accession"] = record.id
            nucleotide_info["Description"] = record.description
            nucleotide_info["Organism"] = record.annotations.get("organism", "na")

            if hasattr(record, 'seq') and isinstance(record.seq, Seq):
                nucleotide_info["Length(bp)"] = str(len(record.seq))
            else:
                nucleotide_info["Length(bp)"] = "na"

            dna_sequence = "na"
            try:
                seq_str = str(record.seq)
                if len(seq_str) > 100:
                    dna_sequence = seq_str[:100] + "..."
                else:
                    dna_sequence = seq_str
            except:
                dna_sequence = "na"
            nucleotide_info["DNA_sequence"] = dna_sequence
            gene_names = []
            product_names = []
            protein_ids = []
            db_xrefs = []
            features = getattr(record, 'features', [])
            for feature in features:
                qualifiers = feature.qualifiers if hasattr(feature, 'qualifiers') else {}
                if "gene" in qualifiers:
                    gene_names.extend(qualifiers.get("gene", []))
                if "product" in qualifiers:
                    product_names.extend(qualifiers.get("product", []))
                if "protein_id" in qualifiers:
                    protein_ids.extend(qualifiers.get("protein_id", []))
                if "db_xref" in qualifiers:
                    db_xrefs.extend(qualifiers.get("db_xref", []))
            for key, values, field in zip(
                    ["Gene_name", "Product", "Protein_id", "db_xref"],
                    [gene_names, product_names, protein_ids, db_xrefs],
                    ["gene", "product", "protein_id", "db_xref"]
            ):
                if values:
                    unique_values = list(dict.fromkeys(values))
                    first_two = unique_values[:2]
                    nucleotide_info[key] = "|".join(first_two)
                else:
                    nucleotide_info[key] = "na"
            references = []
            for ref in record.annotations.get("references", []):
                this_ref = {
                    "title": getattr(ref, "title", "na"),
                    "authors": getattr(ref, "authors", "na"),
                    "journal": getattr(ref, "journal", "na")
                }
                references.append(this_ref)
            nucleotide_info["References"] = references if references else "na"
            return nucleotide_info

        except Exception as e:
            return nucleotide_info
    
def retrieve_gene_info(gene_id, rettype="gb", retmode="text", db_type="gene"):
    """
    Retrieve detailed gene records (GenBank format) using eFetch
    :param gene_id: NCBI gene ID (e.g.: 1017)
    :param rettype: Return data type (e.g.: "gb" for GenBank format, "xml" for XML format)
    :param retmode: Return format (default "text")
    :param db_type: Database type ("gene", "protein", "nucleotide", or "species")
    :return: Sequence record or raw data string
    """
    if db_type == "gene":
        fetch_db = "gene"
    else:
        fetch_db = "nucleotide"
    logging.info(f'----genbank_str---: {fetch_db}  {gene_id}  {rettype}  {retmode}')
    handle = Entrez.efetch(db=fetch_db, id=gene_id, rettype=rettype, retmode=retmode)
    return handle.read()

def get_gene_detailed_info(gene_id, output_file=None, db_type="gene"):
    """
    Get detailed gene information and save (optional)
    :param gene_id: NCBI gene ID
    :param output_file: Optional file path to save results as JSON or CSV
    :param db_type: Database type ("gene", "protein", "nucleotide", or "species")
    :return: Gene information dictionary
    """

    if db_type == "gene":
        genbank_str = retrieve_gene_info(gene_id, rettype="gb", retmode="text", db_type=db_type)
        logging.info(f'----genbank_str---: {genbank_str}')
        try:
            gene_data = parse_gene_info(genbank_str, db_type=db_type)
        except Exception as e:
            logging.info(f"Parse failed: trying protein database")
            genbank_str = retrieve_gene_info(gene_id, rettype="gb", retmode="text", db_type="nucleotide")
            gene_data = parse_gene_info(genbank_str, db_type="nucleotide")

        if output_file:
            _, ext = os.path.splitext(output_file)
            if ext == ".csv":
                pd.DataFrame([gene_data]).to_csv(output_file, index=False, encoding='utf-8-sig')
            elif ext == ".json":
                import json
                with open(output_file, 'w') as f:
                    json.dump(gene_data, f, indent=2)
            else:
                raise ValueError("Only CSV or JSON formats supported")
        
        return gene_data
        

    else:

        genbank_str = retrieve_gene_info(gene_id, rettype="gb", retmode="text", db_type=db_type)
        gene_data = parse_gene_info(genbank_str, db_type=db_type)
        logging.info(f'-------gene_data: {gene_data}')
        if gene_data is None:
            genbank_str = retrieve_gene_info(gene_id, rettype="gb", retmode="text", db_type="gene")
            gene_data = parse_gene_info(genbank_str, db_type="gene")
        return gene_data
    
def download_fasta_files(ids, output_dir, db_type="gene", details_list=[]):
    """
    Download FASTA files for genes or proteins using NCBI dataset tool
    :param ids: List of NCBI IDs
    :param output_dir: Output directory
    :param db_type: Database type ("gene" or "protein")
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    fna_path_list = []
    logging.info(f'-------Downloading FASTA files: {ids}')
    if (db_type == "species" or db_type == "protein"):
        logging.info(f"Using parsed sequence data to create fna files...")
        
        for detail in details_list:
            if 'ncbi_id' in detail and detail['ncbi_id'] in ids:
                ncbi_id = detail['ncbi_id']
                
                logging.info(f'-------Downloading FASTA file: {ncbi_id}')
                matchake_file_str = os.path.join(output_dir,ncbi_id)
                logging.info(f'-------Downloading FASTA file: {matchake_file_str}')
                if not os.path.exists(matchake_file_str):
                    logging.info(f'-------Directory does not exist, creating: {matchake_file_str}')
                    os.system(f"mkdir {matchake_file_str}")

                id_dir = f"{matchake_file_str}/{ncbi_id}"
                if not os.path.exists(id_dir):
                    os.makedirs(id_dir)
                    
                data_dir = f"{id_dir}/ncbi_dataset/data"
                if not os.path.exists(data_dir):
                    os.makedirs(data_dir)
                
                accession = detail.get('Accession', 'unknown')
                organism = detail.get('Organism', 'unknown')
                gene_name = detail.get('Gene_name', 'unknown')
                description = detail.get('Description', 'unknown')
                sequence = detail.get('DNA_sequence', '')
                
                if not sequence:
                    logging.info(f"Warning: ID {ncbi_id} has no DNA sequence data, skipping")
                    continue
                
                fasta_path = f"{data_dir}/{accession}.fna"
                logging.info(f'-------Downloading FASTA file: {fasta_path}')
                if gene_name != 'N/A' and gene_name != 'unknown':
                    header = f">{accession} {gene_name} [organism={organism}] {description}"
                else:
                    header = f">{accession} [organism={organism}] {description}"
                
                with open(fasta_path, "w") as fasta_file:
                    fasta_file.write(f"{header}\n")
                    
                    for i in range(0, len(sequence), 60):
                        fasta_file.write(f"{sequence[i:i+60]}\n")
                
                organism_safe = organism.replace(" ", "_")
                gene_name_safe = gene_name.replace(" ", "_")
                new_filename_path = f"{data_dir}/{organism_safe}_{gene_name_safe}_{accession}.fna"
                
                os.rename(fasta_path, new_filename_path)
                fna_path_list.append(new_filename_path)
                logging.info(f"Successfully created file: {new_filename_path}")
        
        return fna_path_list
    elif db_type == "nucleotide":
        logging.info(f'-------nucleotide download start')
        for ncbi_id in ids:
            logging.info(f'-------Downloading nucleotide FASTA file: {ncbi_id}')
            matchake_file_str = f"{output_dir}/{ncbi_id}"
            gene_file_path = f"{matchake_file_str}/ncbi_dataset/data/{ncbi_id}.fna"
            logging.info(f'-------Downloading nucleotide FASTA file: {gene_file_path}')
            
            if not os.path.exists(matchake_file_str):
                os.makedirs(matchake_file_str)
            
            data_dir = f"{matchake_file_str}/ncbi_dataset/data"
            if not os.path.exists(data_dir):
                os.makedirs(data_dir)
                
            handle = Entrez.efetch(db="nucleotide", id=ncbi_id, rettype="fasta", retmode="text")
            gene_str = handle.read()
            logging.info(f'-------Writing nucleotide FASTA file: {gene_file_path}')
            with open(gene_file_path, "w") as output_handle:
                output_handle.write(gene_str)
            fna_path_list.append(gene_file_path)
        logging.info(f'-------nucleotide download end')
        return fna_path_list
    else:
        for ncbi_id in ids:
            try:
                logging.info(f'-------Downloading FASTA file: {ncbi_id}')
                matchake_file_str = f"{output_dir}/{ncbi_id}"
                logging.info(f'-------Downloading FASTA file: {matchake_file_str}')
                if not os.path.exists(matchake_file_str):
                    logging.info(f'-------Directory does not exist, creating: {matchake_file_str}')
                    os.system(f"mkdir {matchake_file_str}")

                output_path = f'{matchake_file_str}/{ncbi_id}.zip'
                logging.info(f'-------Downloading FASTA file: {output_path}')
                command = [
                    'datasets',
                    'download',
                    'gene',
                    '--id' if db_type == 'protein' else ('gene-id'),
                    ncbi_id,
                    '--include',
                    'gene',
                    '--filename',
                    output_path
                ]
                logging.info(f'-------Download command: {command}')
                try:
                    result = subprocess.run(command, check=True, capture_output=True, text=True)
                except subprocess.CalledProcessError as e:
                    logging.info(f"Error: {e}")
                    logging.info(f"Command Output: { e.stdout}")
                    logging.info(f"Command Error: {e.stderr}")
                logging.info(f"Successfully downloaded {db_type} ID: {ncbi_id}")
                
                logging.info(f'-------Extracting FASTA file: {output_path}')
                with zipfile.ZipFile(output_path, 'r') as zip_ref:
                    zip_ref.extractall(matchake_file_str)

                gene_file_path = f"{matchake_file_str}/ncbi_dataset/data/gene.fna"
                logging.info(f'-------Downloading FASTA file: {gene_file_path}')

                records = list(SeqIO.parse(gene_file_path, "fasta"))
                first_record = records[0]
                logging.info(f'-------Writing FASTA file: {gene_file_path}')
                with open(gene_file_path, "w") as output_handle:
                    SeqIO.write(first_record, output_handle, "fasta")

                record = next(SeqIO.parse(gene_file_path, "fasta"))
                header = record.description

                gene_name_match = re.search(r"\s(\S+)\s\[", header)
                gene_id_match = re.search(r"\[GeneID=(\d+)\]", header)
                organism_match = re.search(r"\[organism=([^\]]+)\]", header)

                if gene_name_match and gene_id_match and organism_match:
                    gene_name = gene_name_match.group(1)
                    gene_id = gene_id_match.group(1)
                    organism = organism_match.group(1).replace(" ", "_")

                    new_filename_path = '/'.join(gene_file_path.split('/')[:-1])+'/'+f"{organism}_{gene_name}_GeneID_{gene_id}.fna"
                    logging.info(f'-------Downloading FASTA file: {new_filename_path}')
                    os.rename(gene_file_path, new_filename_path)
                    fna_path_list.append(new_filename_path)
                    logging.info(f'-------File renamed to: {new_filename_path}')
                else:
                    logging.info(f'-------Could not extract complete information, please check format')


            except subprocess.CalledProcessError as e:
                import traceback
                logging.info(f"Error downloading {db_type} ID: {ncbi_id}\n{e.stderr}")
                logging.info(f"Error downloading {db_type} ID: {ncbi_id}\n{traceback.format_exc()}")
                return None, None

        return fna_path_list

if __name__ == "__main__":
    QUERY_NAME = "Escherichia coli[Organism]"
    OUTPUT_DIR = "./"
    SELECT_DB = "gene"

    try:
        ids = retrieve_ncbi_ids(
            query_name=QUERY_NAME,
            db_type=SELECT_DB,
        )
        
        if not ids:
            raise ValueError(f"No {SELECT_DB} IDs found matching criteria!")

    except Exception as e:
        print(f"检索失败：\n{str(e)}")
        exit(1)
    print('检索到的gene id 为：',ids)
    details_list = []
    for id in ids:
         # 获取详细信息
        gene_details = get_gene_detailed_info(id, output_file= f"./tp53_{id}_info.csv")
        # 添加id字段到gene_details字典中
        
        print(f"{id}:", gene_details)
        if gene_details:
            gene_details['ncbi_id'] = id
            details_list.append(gene_details)
    print(details_list)

        
    # 2. 下载FASTA文件
    try:
        download_fasta_files(
            ids=[ids[2]],
            output_dir=OUTPUT_DIR,
            db_type=SELECT_DB
        )
    except Exception as e:
        print(f"下载失败：\n{str(e)}")
        
    
    
    