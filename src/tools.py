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
        elif '单细胞' in ff:
            Single_cell_list.append(ff)
            
        elif '染色体拷贝数变异' in ff:
            Chromosome_copy_variation_list.append(ff)
        elif '外显子组酶切' in ff:
            Exome_digestion_list.append(ff)
        elif '游离DNA' in ff:
            Free_DNA_list.append(ff)
        elif '通用DNA' in ff:
            General_DNA_list.append(ff)
        elif '酶切DNA' in ff:
            Enzyme_digestion_list.append(ff)
        elif '方向性' in ff:
            RNA_directionality_list.append(ff)
        elif 'fastRNA' in ff:
            fastRNA_list.append(ff)

            
    protocol_dict['ATOPlex DNA建库试剂盒'] = ATOPlexDNA_list
    protocol_dict['ATOPlex HIV建库试剂盒'] = ATOPlexHIV_list
    protocol_dict['ATOPlex MPXV建库试剂盒'] = ATOPlexMPVX_list
    protocol_dict['ATOPlex RNA建库试剂盒'] = ATOPlexRNA_list
    protocol_dict['MGIEasy Fast酶切DNA文库制备试剂盒'] = MGIEasy_fast_list
    protocol_dict['MGIEasy Fast PCR-FREE酶切文库制备试剂盒'] = MGIEasy_PRCfree_list
    protocol_dict['MGICare染色体拷贝数变异检测试剂盒'] = Chromosome_copy_variation_list
    protocol_dict['MGIEasy外显子组酶切文库制备试剂盒'] = Exome_digestion_list
    protocol_dict['MGIEasy 游离DNA文库制备试剂盒'] = Free_DNA_list
    protocol_dict['MGIEasy酶切DNA文库制备试剂盒'] = Enzyme_digestion_list
    protocol_dict['MGIEasy通用DNA文库制备试剂盒'] = General_DNA_list    
    protocol_dict['MGIEasy RNA方向性文库制备试剂盒'] = RNA_directionality_list
    protocol_dict['MGIEasy Fast RNA文库制备试剂盒'] = fastRNA_list
    protocol_dict['MGICare单细胞染色体拷贝数变异检测试剂盒'] = Single_cell_list
    
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
    
    if protocol_name =='ATOPlex RNA建库试剂盒':
        AtoplexRNA_param_dict = get_protocol_param_dict(protocol_name)
        sample_nums = param_dict['sample_nums']
        AtoplexRNA_param_dict['input_sample_nums_1'] = sample_nums        
        AtoplexRNA_param_dict['input_Spike_in_Control_PCR_block_2'] = param_dict['Spike_in_Control_PCR_block']
        AtoplexRNA_param_dict['input_amplicon_multiplicity_2'] = param_dict['amplicon_multiplicity']
        AtoplexRNA_param_dict['input_panel_amp_length_3'] = param_dict['panel_amp_length']
        AtoplexRNA_param_dict['input_Barcode_type_3'] = param_dict['Barcode_type']
        beads_volume_1st,beads_volume_2nd = get_beads_volume(param_dict['Barcode_type'],param_dict['panel_amp_length'])
        AtoplexRNA_param_dict['input_RT_Master_MIX_volume_1'] = 10*(sample_nums+8)
        AtoplexRNA_param_dict['input_DNA_Beads_volume_3'] = beads_volume_1st*(sample_nums+8)
        AtoplexRNA_param_dict['input_TE_Buffer_volume_3'] = 6.5*(sample_nums+8)
        AtoplexRNA_param_dict['input_Ethanol_80_volume_3'] = 160*(sample_nums+8)*2
        AtoplexRNA_param_dict['input_DNA_Beads_volume_5'] = beads_volume_2nd*(sample_nums+8)
        AtoplexRNA_param_dict['input_TE_Buffer_volume_5'] = 25*(sample_nums+8)
        AtoplexRNA_param_dict['input_Ethanol_80_volume_5'] = 160*(sample_nums+8)*2
        if param_dict['Spike_in_Control_PCR_block']:
            AtoplexRNA_param_dict['input_PCR_MIX_volume_2'] = 16*(sample_nums+8)
            AtoplexRNA_param_dict['input_PCR_MIX_volume_4'] = 15.5*(sample_nums+8)
        else:
            AtoplexRNA_param_dict['input_PCR_MIX_volume_2'] = 15*(sample_nums+8)
            AtoplexRNA_param_dict['input_PCR_MIX_volume_4'] = 14.5*(sample_nums+8)
        return AtoplexRNA_param_dict
    
    elif protocol_name =='ATOPlex DNA建库试剂盒':
        AtoplexDNA_param_dict = get_protocol_param_dict(protocol_name)
        sample_nums = param_dict['sample_nums']
        AtoplexDNA_param_dict['input_sample_nums_1'] = sample_nums        
        AtoplexDNA_param_dict['input_Spike_in_Control_PCR_block_1'] = param_dict['Spike_in_Control_PCR_block']
        AtoplexDNA_param_dict['input_amplicon_multiplicity_1'] = param_dict['amplicon_multiplicity']
        AtoplexDNA_param_dict['input_panel_amp_length_2'] = param_dict['panel_amp_length']
        AtoplexDNA_param_dict['input_Barcode_type_2'] = param_dict['Barcode_type']
        beads_volume_1st,beads_volume_2nd = get_beads_volume(param_dict['Barcode_type'],param_dict['panel_amp_length'])

        AtoplexDNA_param_dict['input_DNA_Beads_volume_2'] = beads_volume_1st*(sample_nums+8)
        AtoplexDNA_param_dict['input_TE_Buffer_volume_2'] = 6.5*(sample_nums+8)
        AtoplexDNA_param_dict['input_Ethanol_80_volume_2'] = 160*(sample_nums+8)*2
        
        AtoplexDNA_param_dict['input_DNA_Beads_volume_4'] = beads_volume_2nd*(sample_nums+8)
        AtoplexDNA_param_dict['input_TE_Buffer_volume_4'] = 25*(sample_nums+8)
        AtoplexDNA_param_dict['input_Ethanol_80_volume_4'] = 160*(sample_nums+8)*2
        if param_dict['Spike_in_Control_PCR_block']:
            AtoplexDNA_param_dict['input_PCR_MIX_1st_volume_1'] = 16*(sample_nums+8)
            AtoplexDNA_param_dict['input_PCR_MIX_2nd_volume_3'] = 15.5*(sample_nums+8)
        else:
            AtoplexDNA_param_dict['input_PCR_MIX_1st_volume_1'] = 15*(sample_nums+8)
            AtoplexDNA_param_dict['input_PCR_MIX_2nd_volume_3'] = 14.5*(sample_nums+8)
        return AtoplexDNA_param_dict
    
    elif protocol_name =='ATOPlex HIV建库试剂盒':
        AtoplexHIV_param_dict = get_protocol_param_dict(protocol_name)
        sample_nums = param_dict['sample_nums']
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
    
    elif protocol_name =='ATOPlex MPXV建库试剂盒':
        AtoplexMPVX_param_dict = get_protocol_param_dict(protocol_name)
        sample_nums = param_dict['sample_nums']
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
    elif protocol_name =='MGIEasy 游离DNA文库制备试剂盒':         
        FreeDNA_param_dict = get_protocol_param_dict(protocol_name)
        sample_nums = param_dict['sample_nums']
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
        sample_nums = param_dict['sample_nums']
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
        sample_nums = param_dict['sample_nums']
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
    
    elif protocol_name =='MGIEasy RNA方向性文库制备试剂盒':         
        MGIEasy_RNA_directionality_param_dict = get_protocol_param_dict(protocol_name)
        sample_nums = param_dict['sample_nums']
        MGIEasy_RNA_directionality_param_dict['input_sample_nums_1'] = sample_nums        
        MGIEasy_RNA_directionality_param_dict['input_RT_MIX_volume_1'] = 6*(sample_nums+8)         
        MGIEasy_RNA_directionality_param_dict['input_second_strand_MIX_volume_1'] = 30*(sample_nums+8)
        MGIEasy_RNA_directionality_param_dict['input_Fragmentation_Buffer_volume_1'] = 4*(sample_nums+8)
        MGIEasy_RNA_directionality_param_dict['input_DNA_Beads_volume_2'] = 75*(sample_nums+8)
        MGIEasy_RNA_directionality_param_dict['input_TE_Buffer_volume_2'] = 42*(sample_nums+8)
        MGIEasy_RNA_directionality_param_dict['input_Ethanol_80_volume_2'] =200*(sample_nums+8)*2          
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
        sample_nums = param_dict['sample_nums']
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
        sample_nums = param_dict['sample_nums']
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
        sample_nums = param_dict['sample_nums']
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
        sample_nums = param_dict['sample_nums']
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
    
    elif protocol_name =='MGIEasy通用DNA文库制备试剂盒':         
        General_DNA_param_dict = get_protocol_param_dict(protocol_name)
        sample_nums = param_dict['sample_nums']
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
    
    elif protocol_name =='MGIEasy 游离DNA文库制备试剂盒':         
        FreeDNA_param_dict = get_protocol_param_dict(protocol_name)
        sample_nums = param_dict['sample_nums']
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
def get_protocol_param_dict(protocol_name):
    
    AtoplexDNA_param_dict = {"input_sample_nums_1":96,
              "input_amplicon_multiplicity_1":50,
              "input_Spike_in_Control_PCR_block_1":False,
              "input_PCR_MIX_1st_volume_1":1560, 
              "input_DNA_Beads_volume_2":2600,
              "input_Barcode_type_2":'single_pcr_barcode',
              "input_TE_Buffer_volume_2":676,
              "input_Ethanol_80_volume_2":33280,
              "input_panel_amp_length_2":300,           
              "input_PCR_MIX_2nd_volume_3":1508,     
              "input_DNA_Beads_volume_4":2340,
              "input_TE_Buffer_volume_4":2600,
              "input_Ethanol_80_volume_4":33280,}
    
    AtoplexHIV_param_dict = {"input_sample_nums_1":96,
          "input_PCR_MIX_volume_1":3120,                    
          "input_DNA_Beads_volume_2":3640,
          "input_TE_Buffer_volume_2":3328,
          "input_Ethanol_80_volume_2":33280,                   
          "input_Enzyme_MIX_volume_3":1560,                   
          "input_DNA_Beads_volume_4":6240,
          "input_TE_Buffer_volume_4":4680,
          "input_Ethanol_80_volume_4":33280,                   
          "input_Adapter_ligation_volume_5":3120,                  
          "input_DNA_Beads_volume_6":2080,
          "input_TE_Buffer_volume_6":4160,
          "input_Ethanol_80_volume_6":33280,                      
          "input_DNA_Beads_volume_7":2080,
          "input_TE_Buffer_volume_7":5200,
          "input_Ethanol_80_volume_7":33280,}
        
    AtoplexRNA_param_dict = {"input_sample_nums_1":96,
          "input_RT_Master_MIX_volume_1":1040,           
          "input_PCR_MIX_volume_2":1560, 
          "input_amplicon_multiplicity_2":50,
          "input_Spike_in_Control_PCR_block_2":False,                             
          "input_DNA_Beads_volume_3":2600,
          "input_TE_Buffer_volume_3":676,
          "input_Ethanol_80_volume_3":33280,   
          "input_Barcode_type_3":'single_pcr_barcode',
           "input_panel_amp_length_3":300,                                  
          "input_PCR_MIX_volume_4":1508,                             
          "input_DNA_Beads_volume_5":2600,                  
          "input_TE_Buffer_volume_5":2600,
          "input_Ethanol_80_volume_5":33280,}
    
    AtoplexMPVX_param_dict = {"input_sample_nums_1":96,
          "input_PCR_MIX_volume_1":1560, 
           "input_PCR_MIX_volume_2":1560,                                                   
          "input_DNA_Beads_volume_3":7800,
          "input_TE_Buffer_volume_3":3328,
          "input_Ethanol_80_volume_3":33280,                               
          "input_Enzyme_MIX_volume_4":1560,                   
          "input_Adapter_ligation_volume_5":5600,                                
          "input_DNA_Beads_volume_6":2080,
          "input_TE_Buffer_volume_6":4160,
          "input_Ethanol_80_volume_6":33280,                               
          "input_DNA_Beads_volume_7":2080,
          "input_TE_Buffer_volume_7":5200,
          "input_Ethanol_80_volume_7":33280,}
    
    FreeDNA_param_dict = {"input_sample_nums_1":96,
          "input_A_MIX_volume_1":1040,               
          "input_Adapter_ligation_volume_2":2600,                
          "input_DNA_Beads_volume_3":5200,
          "input_TE_Buffer_volume_3":2392,
          "input_Ethanol_80_volume_3":41600,              
          "input_PCR_MIX_volume_4":3016,                    
          "input_DNA_Beads_volume_5":5200,
          "input_TE_Buffer_volume_5":3328,
          "input_Ethanol_80_volume_5":41600,}
    
    General_DNA_param_dict = {"input_sample_nums_1":96,
          "input_A_MIX_volume_1":1040,               
          "input_Adapter_ligation_volume_2":2600, 
          "input_TE_Buffer_volume_2":2080,                   
          "input_DNA_Beads_volume_3":5200,
          "input_TE_Buffer_volume_3":4160,
          "input_Ethanol_80_volume_3":41600,                    
          "input_PCR_MIX_volume_4":3224,                     
          "input_DNA_Beads_volume_5":5200,
          "input_TE_Buffer_volume_5":3328,
          "input_Ethanol_80_volume_5":41600,}
    
    Enzyme_digestion_DNA_param_dict = {"input_sample_nums_1":96,
          "input_Enzyme_MIX_volume_1":1560,         
          "input_DNA_Beads_volume_2":4992,
          "input_TE_Buffer_volume_2":4992,
          "input_Ethanol_80_volume_2":33280,  
           "input_A_MIX_volume_3":1040,  
          "input_Adapter_ligation_volume_4":2600, 
          "input_TE_Buffer_volume_4":2080,                                                  
          "input_DNA_Beads_volume_5":5200,
          "input_TE_Buffer_volume_5":2184,
          "input_Ethanol_80_volume_5":33280,
           "input_PCR_MIX_volume_6":3200, 
           "input_DNA_quality_6": 50,                          
           "input_output_DNA_quality_6": 300,
           "input_DNA_Beads_volume_7":5200,
          "input_TE_Buffer_volume_7":3328,
          "input_Ethanol_80_volume_7":33280, }
    
    Single_cell_chromosome_param_dict = {"input_sample_nums_1":96,
          "input_Cell_lysis_master_mix_volume_1":520, 
          "input_Pre_PCR_master_mix_volume_1":520,
          "input_Post_PCR_master_mix_volume_1":2684,
           "input_Nuclease_free_water_volume_1":3557,
           "input_Fragment_MIX_volume_2":624,
          "input_TE_Buffer_volume_2":2080,                                  
          "input_DNA_Beads_volume_3":5200,
          "input_TE_Buffer_volume_3":4992,
          "input_Ethanol_80_volume_3":41600,  
           "input_A_MIX_volume_4":1040,  
          "input_ligation_mix_volume_5":2600,      
          "input_DNA_Beads_volume_6":5200,
          "input_TE_Buffer_volume_6":2600,
          "input_Ethanol_80_volume_6":41600,                                         
           "input_PCR_MIX_volume_7":3016, 
           "input_DNA_Beads_volume_8":5200,
          "input_TE_Buffer_volume_8":4992,
          "input_Ethanol_80_volume_8":41600, }
    
    Chromosome_copy_param_dict = {"input_sample_nums_1":96,
          "input_frament_end_repair_volume_1":1040, 
           "input_Ligation_Master_Mix_volume_2":3640,                       
          "input_DNA_Beads_volume_3":4160,
          "input_TE_Buffer_volume_3":2600,
          "input_Ethanol_80_volume_3":41600,                                         
           "input_PCR_MIX_volume_4":2808,
           "input_DNA_Beads_volume_5":5200,                       
          "input_TE_Buffer_volume_5":3424,
          "input_Ethanol_80_volume_5":41600, } 
    
    MGIEasy_fast_param_dict = {"input_sample_nums_1":96,
          "input_FS_Enzyme_MIX_volume_1":1560,  
          "input_DNA_quality_1":50,                     
          "input_DNA_Beads_volume_2":4992,
          "input_TE_Buffer_volume_2":4680,
          "input_Adapter_ligation_volume_3":3120,                                                         
          "input_DNA_Beads_volume_4":2080,
          "input_TE_Buffer_volume_4":2080,
          "input_Ethanol_80_volume_4":33280,
           "input_PCR_MIX_volume_5":2600,                     
           "input_UDB_MIX_volume_5":624, 
           "input_DNA_quality_5": 50, 
           "input_DNA_NG_5": 300, 
           "input_DNA_Beads_volume_6":3952,
          "input_TE_Buffer_volume_6":3328,
          "input_Ethanol_80_volume_6":33280, } 
    
    MGIEasy_PRCfree_param_dict = {"input_sample_nums_1":96,
          "input_FS_Enzyme_MIX_volume_1":1560,             
          "input_DNA_Beads_volume_2":4992,
          "input_TE_Buffer_volume_2":4680,                       
          "input_Adapter_ligation_volume_3":3120,                                                         
          "input_DNA_Beads_volume_4":5200,
          "input_TE_Buffer_volume_4":4160,
           "input_TE_Buffer_volume_5":5304,
           "input_Ethanol_80_volume_5":33280,   } 
    
    MGIEasy_fastRNA_param_dict = {"input_sample_nums_1":96,
          "input_Washing_Buffer_volume_1":10400,  
           "input_Beads_volume_1":2600,                     
          "input_RT_MIX_volume_2":520,
          "input_second_strand_MIX_volume_2":3120,
          "input_Fragmentation_Buffer_volume_2":520, 
           "input_ligation_MIX_volume_3":2600,                       
          "input_DNA_Beads_volume_4":3120,
          "input_TE_Buffer_volume_4":2392,
          "input_Ethanol_80_volume_4":41600,
           "input_PCR_MIX_volume_5":3120,                     
           "input_DNA_Beads_volume_6":4160,
          "input_TE_Buffer_volume_6":3328,
          "input_Ethanol_80_volume_6":41600, }
    
    MGIEasy_RNA_directionality_param_dict = {"input_sample_nums_1":96,                    
          "input_RT_MIX_volume_1":624,
          "input_second_strand_MIX_volume_1":3120,
          "input_Fragmentation_Buffer_volume_1":416,                                                                 
          "input_DNA_Beads_volume_2":7800,
          "input_TE_Buffer_volume_2":4368,
          "input_Ethanol_80_volume_2":41600,
           "input_A_MIX_volume_3":2608,
          "input_ERAT_MIX_volume_3":1040,                                             
           "input_DNA_Beads_volume_4":5200,
          "input_TE_Buffer_volume_4":2392,
          "input_Ethanol_80_volume_4":41600,                                   
           "input_PCR_MIX_volume_5":3120,                                             
           "input_DNA_Beads_volume_6":6240,
          "input_TE_Buffer_volume_6":3328,
          "input_Ethanol_80_volume_6":41600, }    
    
    Exome_digestion_param_dict = {"input_sample_nums_1":96,
          "input_Frag_Master_Mix_volume_1":1560,                                              
          "input_DNA_Beads_volume_2":6240,
          "input_TE_Buffer_volume_2":4472,
          "input_Ethanol_80_volume_2":41600,                        
          "input_A_MIX_volume_3":1040,                        
          "input_Adapter_ligation_volume_4":2600, 
          "input_TE_Buffer_volume_4":2080,                        
          "input_DNA_Beads_volume_5":5200,
          "input_TE_Buffer_volume_5":2184,
          "input_Ethanol_80_volume_5":33280,
          "input_PCR_MIX_volume_6":3200,                       
          "input_DNA_Beads_volume_7":5200,
          "input_TE_Buffer_volume_7":3328,
          "input_Ethanol_80_volume_7":33280, 
          "input_PCR_MIX_volume_8":5824, 
          "input_DNA_Beads_volume_9":10400,
          "input_TE_Buffer_volume_9":3328,
          "input_Ethanol_80_volume_9":41600, }
    
    protocol_dict = {"ATOPlex DNA建库试剂盒":AtoplexDNA_param_dict,
                "ATOPlex HIV建库试剂盒":AtoplexHIV_param_dict,
                "ATOPlex MPXV建库试剂盒":AtoplexMPVX_param_dict,
                "MGIEasy 游离DNA文库制备试剂盒":FreeDNA_param_dict,
                "ATOPlex RNA建库试剂盒":AtoplexRNA_param_dict,
                "MGIEasy外显子组酶切文库制备试剂盒":Exome_digestion_param_dict,
                "MGIEasy Fast RNA文库制备试剂盒":MGIEasy_fastRNA_param_dict,
                "MGIEasy RNA方向性文库制备试剂盒":MGIEasy_RNA_directionality_param_dict,
                "MGIEasy Fast PCR-FREE酶切文库制备试剂盒":MGIEasy_PRCfree_param_dict,
                "MGIEasy Fast酶切DNA文库制备试剂盒":MGIEasy_fast_param_dict,
                "MGICare染色体拷贝数变异检测试剂盒":Chromosome_copy_param_dict,
                "MGIEasy酶切DNA文库制备试剂盒":Enzyme_digestion_DNA_param_dict,
                "MGIEasy通用DNA文库制备试剂盒":General_DNA_param_dict,
                "MGICare单细胞染色体拷贝数变异检测试剂盒":Single_cell_chromosome_param_dict}
    
    return protocol_dict[protocol_name]


def pdf_to_txt(pdf_path):
    
    target_folder = os.path.dirname(pdf_path)
    filename = os.path.basename(pdf_path)
    # 使用pdfplumber打开PDF文件
    with pdfplumber.open(pdf_path) as pdf:
        # 初始化一个空字符串来保存文本内容
        text = ""
        # 遍历PDF中的每一页
        for page in pdf.pages:
            # 提取页面的文本并添加到text变量中
            text += page.extract_text()
            text += "\n\n"  # 添加换行符以分隔不同页面的内容
    # # 构建目标TXT文件的路径，文件名保持不变，只是扩展名改为.txt
    txt_file_path = os.path.join(target_folder, filename.replace(".pdf", ".txt"))
    # 将文本内容写入TXT文件
    with open(txt_file_path, "w", encoding="utf-8") as txt_file:
        txt_file.write(text)

    print(f"已转换文件: {filename} -> {txt_file_path}")
    print(text)
    return text

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
    
    
    