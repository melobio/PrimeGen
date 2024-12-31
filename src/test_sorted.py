from operator import itemgetter

# Sorting function
def sort_genomes(genomes):
    # Define sorting keys (priority order)
    def sorting_key(genome):
        # Priority 1: representative is reference genome
        representative_priority = 0 if genome['representative'] == 'reference genome' else 1
        
        # Priority 2: assembly_level is Complete Genome
        assembly_level_priority = 0 if genome['assembly_level'] == 'Complete Genome' else 1
        
        # Priority 3: whether accession starts with GCF
        accession_priority = 0 if genome['accession'].startswith('GCF') else 1
        
        # Priority 4: scaffold_n50 and contig_n50
        # For GCF, sort by scaffold_n50 in descending order
        if genome['accession'].startswith('GCF'):
            n50_priority = -genome['scaffold_n50']
        else:  # For GCA, sort by contig_n50 in descending order
            n50_priority = -genome['contig_n50']
        
        # Return a tuple, first element is main sorting key, second is secondary key, and so on
        return (representative_priority, assembly_level_priority, accession_priority, n50_priority)

    # Sort dictionary list according to defined sorting rules
    return sorted(genomes, key=sorting_key)

# Example data
genomes = [
    #{'accession': 'GCF_009858895.2', 'representative': 'reference genome', 'taxid': 2697049, 'genus': 'Betacoronavirus', 'species': 'Severe acute respiratory syndrome-related coronavirus', 'organism_name': 'Severe acute respiratory syndrome coronavirus 2', 'infraspecific_name': 'na', 'isolate': 'Wuhan-Hu-1', 'assembly_level': 'Complete Genome', 'ftp_path': 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3', 'group': 'viral', 'name': 'Severe acute respiratory syndrome coronavirus 2 Wuhan-Hu-1', 'origin': 'refseq', 'contig_n50': 29903.0, 'scaffold_n50': 29903.0, 'total_sequence_length': 29903.0},
    #{'accession': 'GCA_009948465.1', 'representative': 'na', 'taxid': 2697049, 'genus': 'Betacoronavirus', 'species': 'Severe acute respiratory syndrome-related coronavirus', 'organism_name': 'Severe acute respiratory syndrome coronavirus 2', 'infraspecific_name': 'na', 'isolate': 'WIV04', 'assembly_level': 'Complete Genome', 'ftp_path': 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/948/465/GCA_009948465.1_ASM994846v1', 'group': 'viral', 'name': 'Severe acute respiratory syndrome coronavirus 2 WIV04', 'origin': 'genbank', 'contig_n50': 29891.0, 'scaffold_n50': 29891.0, 'total_sequence_length': 29891.0},
    {'accession': 'GCA_009948525.1', 'representative': 'na', 'taxid': 2697049, 'genus': 'Betacoronavirus', 'species': 'Severe acute respiratory syndrome-related coronavirus', 'organism_name': 'Severe acute respiratory syndrome coronavirus 2', 'infraspecific_name': 'na', 'isolate': 'WIV06', 'assembly_level': 'Complete Genome', 'ftp_path': 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/948/525/GCA_009948525.1_ASM994852v1', 'group': 'viral', 'name': 'Severe acute respiratory syndrome coronavirus 2 WIV06', 'origin': 'genbank', 'contig_n50': 29854.0, 'scaffold_n50': 29854.0, 'total_sequence_length': 29854.0},
    {'accession': 'GCA_009937905.1', 'representative': 'na', 'taxid': 2697049, 'genus': 'Betacoronavirus', 'species': 'Severe acute respiratory syndrome-related coronavirus', 'organism_name': 'Severe acute respiratory syndrome coronavirus 2', 'infraspecific_name': 'na', 'isolate': 'SARS-CoV-2/human/USA/WA-CDC-02982586-001/2020', 'assembly_level': 'Complete Genome', 'ftp_path': 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/937/905/GCA_009937905.1_ASM993790v1', 'group': 'viral', 'name': 'Severe acute respiratory syndrome coronavirus 2 SARS-CoV-2/human/USA/WA-CDC-02982586-001/2020', 'origin': 'genbank', 'contig_n50': 29882.0, 'scaffold_n50': 29882.0, 'total_sequence_length': 29882.0},
    {'accession': 'GCA_009937915.1', 'representative': 'na', 'taxid': 2697049, 'genus': 'Betacoronavirus', 'species': 'Severe acute respiratory syndrome-related coronavirus', 'organism_name': 'Severe acute respiratory syndrome coronavirus 2', 'infraspecific_name': 'na', 'isolate': 'SARS-CoV-2/human/USA/IL-CDC-02983522-001/2020', 'assembly_level': 'Complete Genome', 'ftp_path': 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/937/915/GCA_009937915.1_ASM993791v1', 'group': 'viral', 'name': 'Severe acute respiratory syndrome coronavirus 2 SARS-CoV-2/human/USA/IL-CDC-02983522-001/2020', 'origin': 'genbank', 'contig_n50': 29882.0, 'scaffold_n50': 29882.0, 'total_sequence_length': 29882.0},
    #{'accession': 'GCA_009858895.3', 'representative': 'reference genome', 'taxid': 2697049, 'genus': 'Betacoronavirus', 'species': 'Severe acute respiratory syndrome-related coronavirus', 'organism_name': 'Severe acute respiratory syndrome coronavirus 2', 'infraspecific_name': 'na', 'isolate': 'Wuhan-Hu-1', 'assembly_level': 'Complete Genome', 'ftp_path': 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/858/895/GCA_009858895.3_ASM985889v3', 'group': 'viral', 'name': 'Severe acute respiratory syndrome coronavirus 2 Wuhan-Hu-1', 'origin': 'genbank', 'contig_n50': 29903.0, 'scaffold_n50': 29903.0, 'total_sequence_length': 29903.0}
    # More dictionary data can be added here
]

# Sorted results
sorted_genomes = sort_genomes(genomes)

# Print sorted dictionary data
for genome in sorted_genomes:
    print(genome['accession'], genome['name'])
