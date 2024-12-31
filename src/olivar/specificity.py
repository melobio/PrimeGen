'''
1. redesign -- Specificity scenarios
    Scenarios: Species identification, proteins, genetic diseases/tumors
    Definition of non-target: Based on scenarios to determine what non-target is. For proteins/genetic diseases/tumors, non-target could be other regions of the genome besides primer regions; for species identification, could be other species
    Sources of non-target: (1) User provided; (2) Genome excluding current design sequence (proteins/genetic diseases/tumors); (3) Choose core-nt database as background (or DB composed of representative sequences); core-nt and representative species DB are already downloaded

    Simplest approach: Direct count - 1 (Assuming provided sequence contains primer design region, it will match current primer position, so count-1 represents matches to other positions)


2. Usage, can directly import all functions,
    Call: BLAST_seq(s_raw, non_target_fasta, spec_blast_conf, minus=True)
    s_raw: Primer sequence to analyze
    non_target_fasta: User input non_target in fasta format (workflow: make blastdb > blast > parse result > return count)
    spec_blast_conf: Configuration including output path, blast mode, threads
        Example:
        tmp_id = generate_random_string()
        spec_blast_conf = {
            'title': 'human',
            'db': f'/home/QShe/workspace/pcr/PrimeGen/{tmp_id}',
            'mode': 'precise',
            'thread': '20',
            'blast_out': f'{tmp_id}_blast.out'
        }
    minus: Default simplified mode (True), meaning non_target contains primer design region, so blast result will include that position, therefore -1 represents non-specific matches
'''

import os, sys
from Bio.Blast import NCBIXML
import random
import string

def generate_random_string(length=6):
    chars = string.ascii_letters + string.digits
    first_char = random.choice(string.ascii_letters)
    remaining_chars = ''.join(random.choices(chars, k=length - 1))
    return first_char + remaining_chars

def buildNCBIdb(file, title, db):  ## db: includes output directory, outputs to current directory if no directory specified
    ## Ensure makeblastdb is in PATH or replace with absolute path
    os.system(f'makeblastdb -in {file} -dbtype nucl -title {title} -parse_seqids -hash_index -out {db}')

def runNCBI(s_raw, db, mode, thread, blast_fa, blast_out):
    '''
    Reference Olivar's non-specificity evaluation:
        precise: Similar to Primer-BLAST approach, considers 3' end matching
        rough: Counts any match
    '''
    config = None
    if mode == 'precise':
        # referred to Primer-BLAST paper
        config = {
            'outfmt': 5, # xml
            'evalue': 5000, 
            'reward': 1, 
            'penalty': -1, 
            'gapopen': 2, 
            'gapextend': 1
        }
    elif mode == 'rough':
        # default of blastn-short
        config = {
            'outfmt': 6, # tabular
            'evalue': 10, 
            'reward': 1, 
            'penalty': -3, 
            'gapopen': 5, 
            'gapextend': 2
        }

    config = {k: str(v) for k, v in config.items()}

    with open(blast_fa, 'w') as out:
        out.write('>qseq\n')
        out.write(f'{s_raw}\n')
    
    # cmd = (f"echo -e '>qseq\\n' | blastn -db {db} -out {blast_out} " + 
    cmd = (f"blastn -query {blast_fa} -db {db} -out {blast_out} " + 
        f"-task blastn-short -num_threads {thread} -outfmt {config['outfmt']} -evalue {config['evalue']} " + 
        f"-reward {config['reward']} -penalty {config['penalty']} -gapopen {config['gapopen']} " + 
        f"-gapextend {config['gapextend']}")
    # print(cmd)
    os.system(cmd)

def parseBlast(BLAST_output, mode):
    '''
    Reference Olivar's non-specificity analysis parsing:
        Parse based on blast output mode, return count
    '''
    all_count = []
    all_loc = []
    with open(BLAST_output, 'r') as handle:
        if mode == 'precise':
            all_records = NCBIXML.parse(handle)
            for record in all_records:      ## For each sequence
                # single BLAST query
                count = 0
                loc_dict = {}
                qlen = record.query_letters
                for alignment in record.alignments:        ## Multiple matches
                    loc_list = []
                    for hsp in alignment.hsps:
                        if (hsp.query_end == qlen) and \
                        (hsp.match[-5:].count(' ')<2) and \
                        (hsp.identities>qlen-5):
                            # no 3' overhang
                            # no more than 1 mismatch within 5nt at 3'end
                            # no more than 4 total mismatches
                            loc_list.append((hsp.sbjct_start, hsp.sbjct_end, record.query))
                            count += 1
                    if loc_list:
                        loc_dict[alignment.accession] = loc_list
                all_count.append(count)   ## Record number of matches meeting conditions for each query
                all_loc.append(loc_dict)  ## Record detailed info of matches meeting conditions for each query
        elif mode == 'rough':             ## Direct count, above needs to meet conditions
            '''
            1. qseqid      query or source (e.g., gene) sequence id
            2. sseqid      subject  or target (e.g., reference genome) sequence id
            3. pident      percentage of identical matches
            4. length      alignment length (sequence overlap)
            5. mismatch    number of mismatches
            6. gapopen     number of gap openings
            7. qstart      start of alignment in query
            8. qend        end of alignment in query
            9. sstart      start of alignment in subject
            10. send        end of alignment in subject
            11. evalue      expect value
            12. bitscore    bit score
            '''
            # get first hsp
            hsp = handle.readline()[:-1].split('\t')
            query_name = hsp[0]
            count = 1
            for line in handle:
                hsp = line[:-1].split('\t')
                if hsp[0] == query_name:
                    count += 1
                else:
                    # next query sequence
                    all_count.append(count)
                    query_name = hsp[0]
                    count = 1
            all_count.append(count) # save count of last record

    return all_count, all_loc

def BLAST_seq(s_raw, db, mode='rough', thread=20, minus=True):
    '''
    s_raw: Primer to evaluate
    db
    mode: chosen('precise', 'rough'), Reference Olivar's design, precise similar to primer-blast considering 3' end; rough counts any match
    thread: Number of threads
    minus: Default simplified mode, meaning non_target contains primer design region, so blast result will include that position, therefore -1 represents non-specific matches
    '''
    tmp_id = generate_random_string()
    blast_fa = f'{tmp_id}.fasta'
    blast_out = f'{tmp_id}_blast.out'
    runNCBI(s_raw, db, mode, thread, blast_fa, blast_out)
    all_count, all_loc = parseBlast(blast_out, mode)
    os.system(f'rm {blast_out} {blast_fa}')
    return all_count[0] - 1 if minus else all_count[0]