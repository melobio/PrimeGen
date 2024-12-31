import pandas as pd
import numpy as np
import os
import pdfplumber

def get_protocol_process(template_protocol_name):
    protocol_dict = {}
    protocol_path = '/reference_data/protocol_code'
    file_list = os.listdir(protocol_path)
    ATOPlexDNA_list = []
    ATOPlexHIV_list = []
    ATOPlexMPVX_list = []
    ATOPlexRNA_list =[]
    MGIEasy_fast_list =[]
    MGIEasy_PRCfree_list = []
    Chromosome_copy_variation_list = []
    Exome_digestion_list =[]
    Free_DNA_list = []
    Enzyme_digestion_list = []
    General_DNA_list =[]
    RNA_directionality_list = []
    fastRNA_list = []
    Single_cell_list = []
    for ff in file_list:
        if 'ATOPlexDNA' in ff:
            ATOPlexDNA_list.append(ff)
        elif 'ATOPlexHIV' in ff:
            ATOPlexHIV_list.append(ff)
        elif 'ATOPlexMPVX' in ff:
            ATOPlexMPVX_list.append(ff)
        elif 'ATOPlexRNA' in ff:
            ATOPlexRNA_list.append(ff)
        elif 'MGIEasy_fast_' in ff:
            MGIEasy_fast_list.append(ff)
        elif 'MGIEasy_PRCfree' in ff:
            MGIEasy_PRCfree_list.append(ff)
        elif 'Single Cell' in ff:
            Single_cell_list.append(ff)
            
        elif 'Chromosome Copy Number Variation' in ff:
            Chromosome_copy_variation_list.append(ff)
        elif 'Exome Digestion' in ff:
            Exome_digestion_list.append(ff)
        elif 'Free DNA' in ff:
            Free_DNA_list.append(ff)
        elif 'General DNA' in ff:
            General_DNA_list.append(ff)
        elif 'Enzyme Digestion DNA' in ff:
            Enzyme_digestion_list.append(ff)
        elif 'Directionality' in ff:
            RNA_directionality_list.append(ff)
        elif 'fastRNA' in ff:
            fastRNA_list.append(ff)

            
    protocol_dict['ATOPlex DNA Library Prep Kit'] = ATOPlexDNA_list
    protocol_dict['ATOPlex HIV Library Prep Kit'] = ATOPlexHIV_list
    protocol_dict['ATOPlex MPXV Library Prep Kit'] = ATOPlexMPVX_list
    protocol_dict['ATOPlex RNA Library Prep Kit'] = ATOPlexRNA_list
    protocol_dict['MGIEasy Fast Enzyme Digestion DNA Library Prep Kit'] = MGIEasy_fast_list
    protocol_dict['MGIEasy Fast PCR-FREE Enzyme Digestion Library Prep Kit'] = MGIEasy_PRCfree_list
    protocol_dict['MGICare Chromosome Copy Number Variation Detection Kit'] = Chromosome_copy_variation_list
    protocol_dict['MGIEasy Exome Enzyme Digestion Library Prep Kit'] = Exome_digestion_list
    protocol_dict['MGIEasy Free DNA Library Prep Kit'] = Free_DNA_list
    protocol_dict['MGIEasy Enzyme Digestion DNA Library Prep Kit'] = Enzyme_digestion_list
    protocol_dict['MGIEasy General DNA Library Prep Kit'] = General_DNA_list    
    protocol_dict['MGIEasy RNA Directional Library Prep Kit'] = RNA_directionality_list
    protocol_dict['MGIEasy Fast RNA Library Prep Kit'] = fastRNA_list
    protocol_dict['MGICare Single Cell Chromosome Copy Number Variation Detection Kit'] = Single_cell_list
    
    process_list = []
    temp_list = protocol_dict[template_protocol_name]
    for ii in range(len(temp_list)):
        for xx in temp_list:
            if str(ii+1) in xx:
                process_list.append(xx)
                
    return process_list

def get_beads_volume(Barcode_type,panel_amp_length):
    beads_volume_1st = 0
    beads_volume_2nd = 0
    if Barcode_type == 'single_pcr_barcode':
        if panel_amp_length >=100 and panel_amp_length <200:
            beads_volume_1st = 32.5
            beads_volume_2nd = 27.5
        elif panel_amp_length >=200 and panel_amp_length <300:
            beads_volume_1st = 25
            beads_volume_2nd = 22.5
        else:
            beads_volume_1st = 22.5
            beads_volume_2nd = 20
    elif Barcode_type == 'double_pcr_barcode':
        if panel_amp_length >=60 and panel_amp_length <100:
            beads_volume_1st = 37.5
            beads_volume_2nd = 37.5
        elif panel_amp_length >=100 and panel_amp_length <200:
            beads_volume_1st = 32.5
            beads_volume_2nd = 25
        elif panel_amp_length >=200 and panel_amp_length <300:
            beads_volume_1st = 25
            beads_volume_2nd = 20
        else:
            beads_volume_1st = 22.5
            beads_volume_2nd = 20

    return beads_volume_1st,beads_volume_2nd

def calculating_solution_volume(param_dict,protocol_name):

    
    if protocol_name =='ATOPlex RNA Library Prep Kit':
        AtoplexRNA_param_dict = get_protocol_param_dict(protocol_name)
        sample_nums = int(param_dict['sample_nums'])
        AtoplexRNA_param_dict['input_sample_nums_1'] = sample_nums        
        AtoplexRNA_param_dict['input_Spike_in_Control_PCR_block_2'] = param_dict['spike_in_control_pcr_block']
        AtoplexRNA_param_dict['input_amplicon_multiplicity_2'] = param_dict['amplicon_multiplicity']
        AtoplexRNA_param_dict['input_panel_amp_length_3'] = param_dict['panel_amp_length']
        AtoplexRNA_param_dict['input_Barcode_type_3'] = param_dict['barcode_type']
        beads_volume_1st,beads_volume_2nd = get_beads_volume(param_dict['barcode_type'],param_dict['panel_amp_length'])
        AtoplexRNA_param_dict['input_RT_Master_MIX_volume_1'] = 10*(sample_nums+8)
        AtoplexRNA_param_dict['input_DNA_Beads_volume_3'] = beads_volume_1st*(sample_nums+8)
        AtoplexRNA_param_dict['input_TE_Buffer_volume_3'] = 6.5*(sample_nums+8)
        AtoplexRNA_param_dict['input_Ethanol_80_volume_3'] = 160*(sample_nums+8)*2
        AtoplexRNA_param_dict['input_DNA_Beads_volume_5'] = beads_volume_2nd*(sample_nums+8)
        AtoplexRNA_param_dict['input_TE_Buffer_volume_5'] = 25*(sample_nums+8)
        AtoplexRNA_param_dict['input_Ethanol_80_volume_5'] = 160*(sample_nums+8)*2
        if param_dict['spike_in_control_pcr_block']:
            AtoplexRNA_param_dict['input_PCR_MIX_volume_2'] = 16*(sample_nums+8)
            AtoplexRNA_param_dict['input_PCR_MIX_volume_4'] = 15.5*(sample_nums+8)
        else:
            AtoplexRNA_param_dict['input_PCR_MIX_volume_2'] = 15*(sample_nums+8)
            AtoplexRNA_param_dict['input_PCR_MIX_volume_4'] = 14.5*(sample_nums+8)
        return AtoplexRNA_param_dict
    
    elif protocol_name =='ATOPlex DNA Library Prep Kit':
        AtoplexDNA_param_dict = get_protocol_param_dict(protocol_name)
        sample_nums = int(param_dict['sample_nums'])
        AtoplexDNA_param_dict['input_sample_nums_1'] = sample_nums
        AtoplexDNA_param_dict['input_Spike_in_Control_PCR_block_1'] = param_dict['spike_in_control_pcr_block']
        AtoplexDNA_param_dict['input_amplicon_multiplicity_1'] = param_dict['amplicon_multiplicity']
        AtoplexDNA_param_dict['input_panel_amp_length_2'] = param_dict['panel_amp_length']
        AtoplexDNA_param_dict['input_Barcode_type_2'] = param_dict['barcode_type']
        beads_volume_1st,beads_volume_2nd = get_beads_volume(param_dict['barcode_type'],param_dict['panel_amp_length'])

        AtoplexDNA_param_dict['input_DNA_Beads_volume_2'] = beads_volume_1st*(sample_nums+8)
        AtoplexDNA_param_dict['input_TE_Buffer_volume_2'] = 6.5*(sample_nums+8)
        AtoplexDNA_param_dict['input_Ethanol_80_volume_2'] = 160*(sample_nums+8)*2
        
        AtoplexDNA_param_dict['input_DNA_Beads_volume_4'] = beads_volume_2nd*(sample_nums+8)
        AtoplexDNA_param_dict['input_TE_Buffer_volume_4'] = 25*(sample_nums+8)
        AtoplexDNA_param_dict['input_Ethanol_80_volume_4'] = 160*(sample_nums+8)*2
        if param_dict['spike_in_control_pcr_block']:
            AtoplexDNA_param_dict['input_PCR_MIX_1st_volume_1'] = 16*(sample_nums+8)
            AtoplexDNA_param_dict['input_PCR_MIX_2nd_volume_3'] = 15.5*(sample_nums+8)
        else:
            AtoplexDNA_param_dict['input_PCR_MIX_1st_volume_1'] = 15*(sample_nums+8)
            AtoplexDNA_param_dict['input_PCR_MIX_2nd_volume_3'] = 14.5*(sample_nums+8)
        return AtoplexDNA_param_dict
    
    elif protocol_name =='ATOPlex HIV Library Prep Kit':
        AtoplexHIV_param_dict = get_protocol_param_dict(protocol_name)
        sample_nums = int(param_dict['sample_nums'])
        AtoplexHIV_param_dict['input_sample_nums_1'] = sample_nums        
        AtoplexHIV_param_dict['input_PCR_MIX_volume_1'] = 30*(sample_nums+8)       
        AtoplexHIV_param_dict['input_DNA_Beads_volume_2'] =35*(sample_nums+8)
        AtoplexHIV_param_dict['input_TE_Buffer_volume_2'] = 32*(sample_nums+8)
        AtoplexHIV_param_dict['input_Ethanol_80_volume_2'] = 160*(sample_nums+8)*2
        AtoplexHIV_param_dict['input_Enzyme_MIX_volume_3'] = 15*(sample_nums+8)      
        AtoplexHIV_param_dict['input_DNA_Beads_volume_4'] = 60*(sample_nums+8)
        AtoplexHIV_param_dict['input_TE_Buffer_volume_4'] = 45*(sample_nums+8)
        AtoplexHIV_param_dict['input_Ethanol_80_volume_4'] = 160*(sample_nums+8)*2        
        AtoplexHIV_param_dict['input_Adapter_ligation_volume_5'] = 30*(sample_nums+8)
        AtoplexHIV_param_dict['input_DNA_Beads_volume_6'] = 20*(sample_nums+8)
        AtoplexHIV_param_dict['input_TE_Buffer_volume_6'] = 20*(sample_nums+8)
        AtoplexHIV_param_dict['input_Ethanol_80_volume_6'] = 160*(sample_nums+8)*2
        AtoplexHIV_param_dict['input_DNA_Beads_volume_7'] = 20*(sample_nums+8)
        AtoplexHIV_param_dict['input_TE_Buffer_volume_7'] = 47*(sample_nums+8)
        AtoplexHIV_param_dict['input_Ethanol_80_volume_7'] = 160*(sample_nums+8)*2
        return AtoplexHIV_param_dict
    
    elif protocol_name =='ATOPlex MPXV Library Prep Kit':
        AtoplexMPVX_param_dict = get_protocol_param_dict(protocol_name)
        sample_nums = int(param_dict['sample_nums'])
        AtoplexMPVX_param_dict['input_sample_nums_1'] = sample_nums        
        AtoplexMPVX_param_dict['input_PCR_MIX_volume_1'] = 15*(sample_nums+8)       
        AtoplexMPVX_param_dict['input_PCR_MIX_volume_2'] =15*(sample_nums+8)
        AtoplexMPVX_param_dict['input_DNA_Beads_volume_3'] = 75*(sample_nums+8)
        AtoplexMPVX_param_dict['input_TE_Buffer_volume_3'] = 32*(sample_nums+8)
        AtoplexMPVX_param_dict['input_Ethanol_80_volume_3'] = 160*(sample_nums+8)*2      
        AtoplexMPVX_param_dict['input_Enzyme_MIX_volume_4'] = 15*(sample_nums+8)
        AtoplexMPVX_param_dict['input_Adapter_ligation_volume_5'] = 30*(sample_nums+8)
        AtoplexMPVX_param_dict['input_DNA_Beads_volume_6'] = 20*(sample_nums+8)
        AtoplexMPVX_param_dict['input_TE_Buffer_volume_6'] = 20*(sample_nums+8)
        AtoplexMPVX_param_dict['input_Ethanol_80_volume_6'] = 160*(sample_nums+8)*2
        AtoplexMPVX_param_dict['input_DNA_Beads_volume_7'] = 20*(sample_nums+8)
        AtoplexMPVX_param_dict['input_TE_Buffer_volume_7'] = 47*(sample_nums+8)
        AtoplexMPVX_param_dict['input_Ethanol_80_volume_7'] = 160*(sample_nums+8)*2
        return AtoplexMPVX_param_dict
    elif protocol_name =='MGIEasy Free DNA Library Prep Kit':         
        FreeDNA_param_dict = get_protocol_param_dict(protocol_name)
        sample_nums = int(param_dict['sample_nums'])
        FreeDNA_param_dict['input_sample_nums_1'] = sample_nums        
        FreeDNA_param_dict['input_A_MIX_volume_1'] = 10*(sample_nums+8)       
        FreeDNA_param_dict['input_Adapter_ligation_volume_2'] =25*(sample_nums+8)
        FreeDNA_param_dict['input_DNA_Beads_volume_3'] = 40*(sample_nums+8)
        FreeDNA_param_dict['input_TE_Buffer_volume_3'] = 23*(sample_nums+8)
        FreeDNA_param_dict['input_Ethanol_80_volume_3'] =200*(sample_nums+8)*2          
        FreeDNA_param_dict['input_PCR_MIX_volume_4'] = 29*(sample_nums+8)
        FreeDNA_param_dict['input_DNA_Beads_volume_5'] = 50*(sample_nums+8)
        FreeDNA_param_dict['input_TE_Buffer_volume_5'] = 32*(sample_nums+8)
        FreeDNA_param_dict['input_Ethanol_80_volume_5'] = 200*(sample_nums+8)*2
        return FreeDNA_param_dict
    elif protocol_name =='MGIEasy外显子组酶切文库制备试剂盒':         
        Exome_digestion_param_dict = get_protocol_param_dict(protocol_name)
        sample_nums = int(param_dict['sample_nums'])
        Exome_digestion_param_dict['input_sample_nums_1'] = sample_nums        
        Exome_digestion_param_dict['input_Frag_Master_Mix_volume_1'] = 15*(sample_nums+8)       
        Exome_digestion_param_dict['input_DNA_Beads_volume_2'] = 48*(sample_nums+8)
        Exome_digestion_param_dict['input_TE_Buffer_volume_2'] = 43*(sample_nums+8)
        Exome_digestion_param_dict['input_Ethanol_80_volume_2'] =200*(sample_nums+8)*2          
        Exome_digestion_param_dict['input_A_MIX_volume_3'] = 10*(sample_nums+8)
        Exome_digestion_param_dict['input_Adapter_ligation_volume_4'] = 25*(sample_nums+8)
        Exome_digestion_param_dict['input_TE_Buffer_volume_4'] = 20*(sample_nums+8)        
        Exome_digestion_param_dict['input_DNA_Beads_volume_5'] = 50*(sample_nums+8)
        Exome_digestion_param_dict['input_TE_Buffer_volume_5'] = 21*(sample_nums+8)
        Exome_digestion_param_dict['input_Ethanol_80_volume_5'] = 200*(sample_nums+8)*2
        Exome_digestion_param_dict['input_PCR_MIX_volume_6'] = 31*(sample_nums+8)
        Exome_digestion_param_dict['input_DNA_Beads_volume_7'] = 50*(sample_nums+8)
        Exome_digestion_param_dict['input_TE_Buffer_volume_7'] = 32*(sample_nums+8)
        Exome_digestion_param_dict['input_Ethanol_80_volume_7'] = 200*(sample_nums+8)*2
        return Exome_digestion_param_dict
    
    elif protocol_name =='MGIEasy酶切DNA文库制备试剂盒':         
        Enzyme_digestion_DNA_param_dict = get_protocol_param_dict(protocol_name)
        sample_nums = int(param_dict['sample_nums'])
        Enzyme_digestion_DNA_param_dict['input_sample_nums_1'] = sample_nums        
        Enzyme_digestion_DNA_param_dict['input_Enzyme_Mix_volume_1'] = 15*(sample_nums+8)       
        Enzyme_digestion_DNA_param_dict['input_DNA_Beads_volume_2'] = 48*(sample_nums+8)
        Enzyme_digestion_DNA_param_dict['input_TE_Buffer_volume_2'] = 43*(sample_nums+8)
        Enzyme_digestion_DNA_param_dict['input_Ethanol_80_volume_2'] =200*(sample_nums+8)*2          
        Enzyme_digestion_DNA_param_dict['input_A_MIX_volume_3'] = 10*(sample_nums+8)
        Enzyme_digestion_DNA_param_dict['input_Adapter_ligation_volume_4'] = 25*(sample_nums+8)
        Enzyme_digestion_DNA_param_dict['input_TE_Buffer_volume_4'] = 20*(sample_nums+8)        
        Enzyme_digestion_DNA_param_dict['input_DNA_Beads_volume_5'] = 50*(sample_nums+8)
        Enzyme_digestion_DNA_param_dict['input_TE_Buffer_volume_5'] = 21*(sample_nums+8)
        Enzyme_digestion_DNA_param_dict['input_Ethanol_80_volume_5'] = 200*(sample_nums+8)*2
        Enzyme_digestion_DNA_param_dict['input_PCR_MIX_volume_6'] = 31*(sample_nums+8)
        Enzyme_digestion_DNA_param_dict['input_DNA_Beads_volume_7'] = 50*(sample_nums+8)
        Enzyme_digestion_DNA_param_dict['input_TE_Buffer_volume_7'] = 32*(sample_nums+8)
        Enzyme_digestion_DNA_param_dict['input_Ethanol_80_volume_7'] = 200*(sample_nums+8)*2
        return Enzyme_digestion_DNA_param_dict  
    
    elif protocol_name =='MGIEasy RNA Directional Library Prep Kit':         
        MGIEasy_RNA_directionality_param_dict = get_protocol_param_dict(protocol_name)
        sample_nums = int(param_dict['sample_nums'])
        MGIEasy_RNA_directionality_param_dict['input_sample_nums_1'] = sample_nums        
        MGIEasy_RNA_directionality_param_dict['input_RT_MIX_volume_1'] = 6*(sample_nums+8)         
        MGIEasy_RNA_directionality_param_dict['input_second_strand_MIX_volume_1'] = 30*(sample_nums+8)
        MGIEasy_RNA_directionality_param_dict['input_Fragmentation_Buffer_volume_1'] = 4*(sample_nums+8)
        MGIEasy_RNA_directionality_param_dict['input_DNA_Beads_volume_2'] = 75*(sample_nums+8)
        MGIEasy_RNA_directionality_param_dict['input_TE_Buffer_volume_2'] = 42*(sample_nums+8)
        MGIEasy_RNA_directionality_param_dict['input_Ethanol_80_volume_2'] = 200*(sample_nums+8)*2
        MGIEasy_RNA_directionality_param_dict['input_A_MIX_volume_3'] = 25*(sample_nums+8)
        MGIEasy_RNA_directionality_param_dict['input_ERAT_MIX_volume_3'] = 10*(sample_nums+8)          
        MGIEasy_RNA_directionality_param_dict['input_DNA_Beads_volume_4'] = 50*(sample_nums+8)
        MGIEasy_RNA_directionality_param_dict['input_TE_Buffer_volume_4'] = 22*(sample_nums+8)
        MGIEasy_RNA_directionality_param_dict['input_Ethanol_80_volume_4'] = 200*(sample_nums+8)*2
        MGIEasy_RNA_directionality_param_dict['input_PCR_MIX_volume_5'] = 30*(sample_nums+8)
        MGIEasy_RNA_directionality_param_dict['input_DNA_Beads_volume_6'] = 60*(sample_nums+8)
        MGIEasy_RNA_directionality_param_dict['input_TE_Buffer_volume_6'] = 32*(sample_nums+8)
        MGIEasy_RNA_directionality_param_dict['input_Ethanol_80_volume_6'] = 200*(sample_nums+8)*2
        return MGIEasy_RNA_directionality_param_dict  
    
    elif protocol_name =='MGIEasy Fast RNA文库制备试剂盒':         
        MGIEasy_fastRNA_param_dict = get_protocol_param_dict(protocol_name)
        sample_nums = int(param_dict['sample_nums'])
        MGIEasy_fastRNA_param_dict['input_sample_nums_1'] = sample_nums
        MGIEasy_fastRNA_param_dict['input_Washing_Buffer_volume_1'] = 50*(sample_nums+8)
        MGIEasy_fastRNA_param_dict['input_Beads_volume_1'] = 25*(sample_nums+8)
        MGIEasy_fastRNA_param_dict['input_RT_MIX_volume_2'] = 5*(sample_nums+8)         
        MGIEasy_fastRNA_param_dict['input_second_strand_MIX_volume_2'] = 30*(sample_nums+8)
        MGIEasy_fastRNA_param_dict['input_Fragmentation_Buffer_volume_2'] = 5*(sample_nums+8)
        MGIEasy_fastRNA_param_dict['input_ligation_MIX_volume_3'] = 25*(sample_nums+8)                                         
        MGIEasy_fastRNA_param_dict['input_DNA_Beads_volume_4'] = 30*(sample_nums+8)
        MGIEasy_fastRNA_param_dict['input_TE_Buffer_volume_4'] = 23*(sample_nums+8)
        MGIEasy_fastRNA_param_dict['input_Ethanol_80_volume_4'] =200*(sample_nums+8)*2          
        MGIEasy_fastRNA_param_dict['input_PCR_MIX_volume_5'] = 25*(sample_nums+8)
        MGIEasy_fastRNA_param_dict['input_DNA_Beads_volume_6'] = 60*(sample_nums+8)
        MGIEasy_fastRNA_param_dict['input_TE_Buffer_volume_6'] = 32*(sample_nums+8)
        MGIEasy_fastRNA_param_dict['input_Ethanol_80_volume_6'] = 200*(sample_nums+8)*2
        return MGIEasy_fastRNA_param_dict  
    
    elif protocol_name =='MGIEasy Fast PCR-FREE酶切文库制备试剂盒':         
        MGIEasy_PRCfree_param_dict = get_protocol_param_dict(protocol_name)
        sample_nums = int(param_dict['sample_nums'])
        MGIEasy_PRCfree_param_dict['input_sample_nums_1'] = sample_nums
        MGIEasy_PRCfree_param_dict['input_FS_Enzyme_MIX_volume_1'] = 15*(sample_nums+8)                                       
        MGIEasy_PRCfree_param_dict['input_DNA_Beads_volume_2'] = 48*(sample_nums+8)
        MGIEasy_PRCfree_param_dict['input_TE_Buffer_volume_2'] = 47*(sample_nums+8)
        MGIEasy_PRCfree_param_dict['input_Adapter_ligation_volume_3'] = 30*(sample_nums+8) 
        MGIEasy_PRCfree_param_dict['input_DNA_Beads_volume_4'] = 40*(sample_nums+8)
        MGIEasy_PRCfree_param_dict['input_TE_Buffer_volume_4'] = 50*(sample_nums+8)
        MGIEasy_PRCfree_param_dict['input_TE_Buffer_volume_5'] = 51*(sample_nums+8)
        MGIEasy_PRCfree_param_dict['input_Ethanol_80_volume_5'] = 160*(sample_nums+8)*2
        return MGIEasy_PRCfree_param_dict   
    
    elif protocol_name =='MGICare染色体拷贝数变异检测试剂盒':         
        Chromosome_copy_param_dict = get_protocol_param_dict(protocol_name)
        sample_nums = int(param_dict['sample_nums'])
        Chromosome_copy_param_dict['input_sample_nums_1'] = sample_nums
        Chromosome_copy_param_dict['input_frament_end_repair_volume_1'] = 10*(sample_nums+8)
        Chromosome_copy_param_dict['input_Ligation_Master_Mix_volume_2'] = 35*(sample_nums+8)                                
        Chromosome_copy_param_dict['input_DNA_Beads_volume_3'] = 40*(sample_nums+8)
        Chromosome_copy_param_dict['input_TE_Buffer_volume_3'] = 25*(sample_nums+8)
        Chromosome_copy_param_dict['input_Ethanol_80_volume_3'] =200*(sample_nums+8)*2          
        Chromosome_copy_param_dict['input_PCR_MIX_volume_4'] = 27*(sample_nums+8)
        Chromosome_copy_param_dict['input_DNA_Beads_volume_5'] = 50*(sample_nums+8)
        Chromosome_copy_param_dict['input_TE_Buffer_volume_5'] = 32*(sample_nums+8)
        Chromosome_copy_param_dict['input_Ethanol_80_volume_5'] = 200*(sample_nums+8)*2
        return Chromosome_copy_param_dict 
    
    elif protocol_name =='MGICare单细胞染色体拷贝数变异检测试剂盒':         
        Single_cell_chromosome_param_dict = get_protocol_param_dict(protocol_name)
        sample_nums = int(param_dict['sample_nums'])
        Single_cell_chromosome_param_dict['input_sample_nums_1'] = sample_nums
        Single_cell_chromosome_param_dict['input_Cell_lysis_master_mix_volume_1'] = 5*(sample_nums+8)
        Single_cell_chromosome_param_dict['input_Pre_PCR_master_mix_volume_1'] = 5*(sample_nums+8)
        Single_cell_chromosome_param_dict['input_Post_PCR_master_mix_volume_1'] = 25.8*(sample_nums+8)
        Single_cell_chromosome_param_dict['input_Nuclease_free_water_volume_1'] = 34.2*(sample_nums+8)
        Single_cell_chromosome_param_dict['input_Fragment_MIX_volume_2'] = 6*(sample_nums+8)
        Single_cell_chromosome_param_dict['input_TE_Buffer_volume_2'] = 20*(sample_nums+8)
        Single_cell_chromosome_param_dict['input_DNA_Beads_volume_3'] = 48*(sample_nums+8)
        Single_cell_chromosome_param_dict['input_TE_Buffer_volume_3'] = 43*(sample_nums+8)
        Single_cell_chromosome_param_dict['input_Ethanol_80_volume_3'] =200*(sample_nums+8)*2          
        Single_cell_chromosome_param_dict['input_A_MIX_volume_4'] = 10*(sample_nums+8)
        Single_cell_chromosome_param_dict['input_ligation_mix_volume_5'] = 25*(sample_nums+8) 
        Single_cell_chromosome_param_dict['input_DNA_Beads_volume_6'] = 40*(sample_nums+8)
        Single_cell_chromosome_param_dict['input_TE_Buffer_volume_6'] = 23*(sample_nums+8)
        Single_cell_chromosome_param_dict['input_Ethanol_80_volume_6'] = 200*(sample_nums+8)*2
        Single_cell_chromosome_param_dict['input_PCR_MIX_volume_7'] = 29*(sample_nums+8)
        Single_cell_chromosome_param_dict['input_DNA_Beads_volume_8'] = 50*(sample_nums+8)
        Single_cell_chromosome_param_dict['input_TE_Buffer_volume_8'] = 32*(sample_nums+8)
        Single_cell_chromosome_param_dict['input_Ethanol_80_volume_8'] = 200*(sample_nums+8)*2
        return Single_cell_chromosome_param_dict 
    
    elif protocol_name =='MGIEasy Universal DNA Library Prep Kit':         
        General_DNA_param_dict = get_protocol_param_dict(protocol_name)
        sample_nums = int(param_dict['sample_nums'])
        General_DNA_param_dict['input_sample_nums_1'] = sample_nums        
        General_DNA_param_dict['input_A_MIX_volume_1'] = 10*(sample_nums+8)  
        General_DNA_param_dict['input_Adapter_ligation_volume_2'] = 25*(sample_nums+8) 
        General_DNA_param_dict['input_TE_Buffer_volume_2'] = 20*(sample_nums+8)     
        General_DNA_param_dict['input_DNA_Beads_volume_3'] = 50*(sample_nums+8)
        General_DNA_param_dict['input_TE_Buffer_volume_3'] = 40*(sample_nums+8)
        General_DNA_param_dict['input_Ethanol_80_volume_3'] = 200*(sample_nums+8)*2
        General_DNA_param_dict['input_PCR_MIX_volume_4'] = 31*(sample_nums+8)
        General_DNA_param_dict['input_DNA_Beads_volume_5'] = 50*(sample_nums+8)
        General_DNA_param_dict['input_TE_Buffer_volume_5'] = 32*(sample_nums+8)
        General_DNA_param_dict['input_Ethanol_80_volume_5'] = 200*(sample_nums+8)*2
        return General_DNA_param_dict
    
    elif protocol_name =='MGIEasy Cell-Free DNA Library Prep Kit':         
        FreeDNA_param_dict = get_protocol_param_dict(protocol_name)
        sample_nums = int(param_dict['sample_nums'])
        FreeDNA_param_dict['input_sample_nums_1'] = sample_nums        
        FreeDNA_param_dict['input_A_MIX_volume_1'] = 10*(sample_nums+8)  
        FreeDNA_param_dict['input_Adapter_ligation_volume_2'] = 25*(sample_nums+8)    
        FreeDNA_param_dict['input_DNA_Beads_volume_3'] = 40*(sample_nums+8)
        FreeDNA_param_dict['input_TE_Buffer_volume_3'] = 23*(sample_nums+8)
        FreeDNA_param_dict['input_Ethanol_80_volume_3'] = 200*(sample_nums+8)*2
        FreeDNA_param_dict['input_PCR_MIX_volume_4'] = 29*(sample_nums+8)
        FreeDNA_param_dict['input_DNA_Beads_volume_5'] = 50*(sample_nums+8)
        FreeDNA_param_dict['input_TE_Buffer_volume_5'] = 32*(sample_nums+8)
        FreeDNA_param_dict['input_Ethanol_80_volume_5'] = 200*(sample_nums+8)*2
        return FreeDNA_param_dict    

def calculate_liquid(sample_num,liquid_volume):
    liquid_list =[]
    if sample_num <8:
        per_sample_liquid = liquid_volume/sample_num
        temp_num = sample_num
        for ii in range(8):
            if temp_num>0:
                liquid_list.append(per_sample_liquid)
                temp_num=temp_num-1
            else:
                liquid_list.append(0.0)
        return liquid_list
    else:
        count1 = int(sample_num/8)
        count2 = sample_num%8
        per_sample_liquid = liquid_volume/(8*(count1+1)+count2)

        for ii in range(8):
            if count2>0:
                liquid_list.append((count1+2)*per_sample_liquid)
                count2=count2-1
            else:
                liquid_list.append((count1+1)*per_sample_liquid)
        return liquid_list
    
    
if __name__ == "__main__":
    
    pdf_path = './data/MGIEasy_Fast_酶切⽂库制备试剂套装使用说明书3.0-套装版本号v2.0.pdf'
    pdf_to_txt(pdf_path)
    
    
    