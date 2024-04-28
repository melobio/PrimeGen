from Bio import Entrez
from Bio import SeqIO
import json
import time
# API请求有频率限制，TODO：多账号混合、时延机制
Entrez.email = "912837656@qq.com"
Entrez.api_key = "24e5b06635f7b32dbd9b89c178a3caf3f009"#"db4b31edb00e220c3dd378908403eddbc308"

def ncbi_search(term, db_name):
    #参考基因组检索
    if db_name=="genome":
        handle = Entrez.esearch(db=db_name, term=term, retmax=10)
        record = Entrez.read(handle)
        print(record)
        ids = record["IdList"]
        for id in ids:
            esummary_handle = Entrez.esummary(db="genome", id=id)
            esummary_record = Entrez.read(esummary_handle)
            esummary_handle.close()

            # 打印每个记录的摘要信息
            assembly_name = esummary_record[0]["Assembly_Name"]
            handle = Entrez.esearch(db="assembly", term=assembly_name, retmax=100000)
            record = Entrez.read(handle)
            handle.close()
            assembly_ids = record["IdList"]
            handle = Entrez.efetch(db="assembly", id=assembly_ids, rettype="docsum", retmode="xml")
            records = Entrez.read(handle)
            handle.close()
            assembly_accessions = [doc['AssemblyAccession'] for doc in records['DocumentSummarySet']['DocumentSummary']]
            print(assembly_accessions)
            #exit()
            url = f"datasets download genome accession {assembly_accessions[0]} --include gff3,rna,cds,protein,genome,seq-report --filename "#{path}/{term.replace(' ','_')}.zip"
            print(url)
            exit()
            try:
                os.system(url)
            except:
                print('*****************download error, download again*****************')
                os.system(f'rm {path}/{term.replace(" ","_")}.zip')
                os.system(url)

            file_path = f"{path}/{term.replace(' ','_')}.zip"
            with zipfile.ZipFile(file_path, 'r') as zip_ref:
                zip_ref.extractall(f"{path}/")

            target_file_path = glob.glob(f"{path}/ncbi_dataset/data/G*")[0]
            target_file = target_file_path+'/cds_from_genomic.fna'

            #file_to_copy = f"./non_target_cds/{term.replace(' ','_')}_cds_from_genomic.fna"
            file_to_copy = f"/reference_data/target/{term.replace(' ','_')}_cds_from_genomic.fna"
            print(f"Source file: {target_file}")
            print(f"File exists: {os.path.exists(target_file)}")
            print(f"Destination file: {file_to_copy}")
            shutil.copy(target_file, file_to_copy)
            print(f"finish{target_file} --> {file_to_copy}")

        return file_to_copy #返回下载好的数据位置
    
    #靶基因检索
    elif db_name=="gene":
        term_str = f"""({term[0]}) AND {term[1]} [Text Word]"""
        
        #先下载物种参考基因组
        handle = Entrez.esearch(db="genome", term=term[0], retmax=10)
        record = Entrez.read(handle)
        ids = record["IdList"]
        for id in ids:
            esummary_handle = Entrez.esummary(db="genome", id=id)
            esummary_record = Entrez.read(esummary_handle)
            esummary_handle.close()

            # 打印每个记录的摘要信息
            assembly_name = esummary_record[0]["Assembly_Name"]
            handle = Entrez.esearch(db="assembly", term=assembly_name, retmax=100000)
            record = Entrez.read(handle)
            handle.close()
            assembly_ids = record["IdList"]
            handle = Entrez.efetch(db="assembly", id=assembly_ids, rettype="docsum", retmode="xml")
            records = Entrez.read(handle)
            handle.close()
            assembly_accessions = [doc['AssemblyAccession'] for doc in records['DocumentSummarySet']['DocumentSummary']]
            print(assembly_accessions)
            exit()
            url = f"datasets download genome accession {assembly_accessions[0]} --include gff3,rna,cds,protein,genome,seq-report --filename {path}/{term.replace(' ','_')}.zip"
            try:
                os.system(url)
            except:
                print('*****************download error, download again*****************')
                os.system(f'rm {path}/{term.replace(" ","_")}.zip')
                os.system(url)

            file_path = f"{path}/{term.replace(' ','_')}.zip"
            with zipfile.ZipFile(file_path, 'r') as zip_ref:
                zip_ref.extractall(f"{path}/")

            target_file_path = glob.glob(f"{path}/ncbi_dataset/data/G*")[0]
            target_file = target_file_path+'/cds_from_genomic.fna'

            #file_to_copy = f"./non_target_cds/{term.replace(' ','_')}_cds_from_genomic.fna"
            file_to_copy = f"/reference_data/target/{term.replace(' ','_')}_cds_from_genomic.fna"
            print(f"Source file: {target_file}")
            print(f"File exists: {os.path.exists(target_file)}")
            print(f"Destination file: {file_to_copy}")
            shutil.copy(target_file, file_to_copy)
            print(f"finish{target_file} --> {file_to_copy}")
        
        #file_to_copy
        



        #kaishi1
        handle = Entrez.esearch(db=db_name, term=term_str, retmax=10)
        record = Entrez.read(handle)
        print(record)
        ids = record["IdList"]
        print(len(ids))
        for taxonomy_id in ids:
            esummary_handle = Entrez.esummary(db="gene", id=taxonomy_id)
            esummary_record = Entrez.read(esummary_handle)
            esummary_handle.close()
            #print(esummary_record)
            Name = esummary_record['DocumentSummarySet']['DocumentSummary'][0]['Name']
            if term[1] not in Name: #确保用户需要的靶基因在其中
                continue
            print(esummary_record)
            Description =  esummary_record['DocumentSummarySet']['DocumentSummary'][0]['Description']
            Assembly_accessions = esummary_record['DocumentSummarySet']['DocumentSummary'][0]['GenomicInfo'][0]['ChrAccVer']
            Start_ = esummary_record['DocumentSummarySet']['DocumentSummary'][0]['GenomicInfo'][0]['ChrStart']
            End_ = esummary_record['DocumentSummarySet']['DocumentSummary'][0]['GenomicInfo'][0]['ChrStop']
            print(Name)
            print(Description)
            print(Assembly_accessions)
            print(Start_)
            print(End_)
            continue
        exit()
        pass


if __name__ == "__main__":
    ncbi_search("Escherichia coli O26:H11 str. 11368[Organism]","genome")
    #ncbi_search("Severe acute respiratory syndrome coronavirus 2 (SARS-CoV-2)","genome")
    #ncbi_search(["Mycobacterium tuberculosis H37Rv","katG"],"gene")
