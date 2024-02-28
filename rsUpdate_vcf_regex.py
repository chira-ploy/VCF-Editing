#!usr/bin/env python

#Author = Chiranan Khantham

from fuc import pyvcf
import pandas as pd
import re

from fuc import pyvcf
import pandas as pd 
import re

def split_location(location):
    pattern = r'(\w+):(\d+)-(\d+)'
    match = re.match(pattern, location)
    if match:
        chromosome = match.group(1)
        position = match.group(2)
    else:
        return " "
            
    return chromosome, position


def extract_rsID(variation):
    pattern = r'(\brs\d+)(?:,(\w+))?'
    match = re.match(pattern, variation)
    if match:
        rsID = match.group(1)       
        return rsID
    else:
        return " "
    
def rsIDupdate_vcf(SPnumber):
    
    #read vcf file and construct in vcf frame format
    vcf_file = pyvcf.VcfFrame.from_file(f'{SPnumber}_vep.vcf')
    
    #read txt file and select only the columns we need
    txt_file = pd.read_csv(f'{SPnumber}_vep.txt', sep='\t', header=0)
    txt_working = txt_file[['Location','Allele', 'Existing_variation']]
    txt_working = txt_working.copy()
    
    #implement functions to get Chrom, POS, rsID
    txt_working.loc[:, 'CHROM'] = txt_working['Location'].apply(lambda x: split_location(x)[0])
    txt_working.loc[:, 'POS'] = txt_working['Location'].apply(lambda x: split_location(x)[1])
    txt_working.loc[:, 'rsID'] = txt_working['Existing_variation'].apply(lambda x: extract_rsID(x))

    
    #remove duplicate before merge step
    txt_working_noDup = txt_working.drop_duplicates(keep='first')
    
    #change the data type before merge step
    txt_working_noDup.loc[:, 'CHROM'] = txt_working_noDup['CHROM'].astype(str)
    txt_working_noDup.loc[:, 'POS'] = txt_working_noDup['POS'].astype(int)
    txt_working_noDup.loc[:, 'rsID'] = txt_working_noDup['rsID'].astype(str)
    
    #merge
    merge_df = vcf_file.df.merge(txt_working_noDup, how='left', left_on=['CHROM', 'POS', 'ALT'], 
                              right_on=['CHROM', 'POS', 'Allele'])
    
    #select the vcf columns, change column name, and fill NaN with .
    update_df = merge_df[['CHROM', 'POS', 'rsID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']]
    update_df = update_df.rename(columns={'rsID':'ID'})
    update_df['ID'] = update_df['ID'].fillna('.')
    
    return update_df

def get_header(vcf_filename):
    header_lines = []
    with open(vcf_filename, 'r') as file:
        for line in file:
            if line.startswith('#CHROM'):
                break
            header_lines.append(line.strip())
            
    header_line_str = '\n'.join(header_lines)
    header_line_str += '\n#'
    
    return header_line_str
    
    
def write_vcf(vcf_filename, update_df, new_file_path):
    header_line_str = get_header(vcf_filename)
    
    with open(new_file_path,'w') as file:
        file.write(header_line_str)
        
    update_df.to_csv(new_file_path, mode='a', sep='\t', index=False)
    

def main(SPnumber):
    update_df = rsIDupdate_vcf(SPnumber)
    write_vcf(f'{SPnumber}_vep.vcf', update_df, f'{SPnumber}_vep_rsUpdate.vcf')

