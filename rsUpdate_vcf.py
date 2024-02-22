#!usr/bin/env python

#Author = Chiranan Khantham

from fuc import pyvcf
import pandas as pd
import os

def split_location(loc_str):
    """Splits a location string into CHR and position values.

    Args:
        loc_str: The location string to split.

    Returns:
        A tuple containing the extracted CHR and position values,
        or None if the string is invalid.
    """

    parts = loc_str.split(":")
    if len(parts) == 2:
        chr_val = str(parts[0])
        position_parts = parts[1].split("-")
        if len(position_parts) == 2:
            position_val = str(position_parts[0])
            return chr_val, position_val

    return None

def extract_id(variant):
    """Extracts the first ID starting with "rs" from a variation string.

    Args:
        variant: The variation string (e.g., "rs10904045,COSV59118495", "rs2472749").

    Returns:
        The first ID starting with "rs" or None if no such ID is found.
    """

    for part in variant.split(","):
        if part.startswith("rs"):
            return part
    return None

def match_and_update(sp_num):

    #Construct VcfFrame from a VCF file
    vcf = pyvcf.VcfFrame.from_file(f'{sp_num}_vep.vcf')
    txt_df = pd.read_csv(f'{sp_num}_vep.txt', sep='\t', header=0)
    txt_df = txt_df[['Location', 'Allele', 'REF_ALLELE', 'Existing_variation']]
    txt_df['Existing_variation'] = txt_df['Existing_variation'].astype(str)

    # create new columns, de-duplicate txt file and cast pos to integer
    # Split 'Location' using : into 2 new columns 'CHROM', 'POS'
    txt_df = txt_df.assign(
        CHROM=txt_df['Location'].apply(lambda x: split_location(x)[0]),
        POS=txt_df['Location'].apply(lambda x: split_location(x)[1]),
    )
    txt_df['ID'] = txt_df['Existing_variation'].apply(extract_id)
    txt_df_deduped = txt_df.drop_duplicates().reset_index(drop=True)
    txt_df_deduped['POS'] = txt_df_deduped['POS'].astype(int)

    # merge the two dataframes on the CHROM, POS and ALT/Allele columns
    df = vcf.df.merge(txt_df_deduped,how='left',left_on=['CHROM','POS','ALT'], right_on=['CHROM','POS','Allele'])
    df = df[['CHROM','POS', 'ID_y','REF','ALT','QUAL','FILTER','INFO']]
    df = df.rename({'ID_y':'ID'},axis=1)
    df.fillna('.',inplace=True)

    return df

def get_header(vcf_file):
    """
    Get the header of the original VCF file.
    """
    header_info = []
    with open(vcf_file, 'r') as vcf_original:
        for line in vcf_original:
            if line.startswith('#CHROM'):
                break
            header_info.append(line)
    header_info_str = ''.join(header_info)
    header_info_str += '#'

    return header_info_str

def write_to_vcf_files(vcf_original_path, vcf_updated_path, vcf_df):
    """
    Create a vcf file with from a Pandas dataframe.
    """
    vcf_headers = get_header(vcf_original_path)

    with open(vcf_updated_path,'w') as vcf_updated:
        vcf_updated.write(vcf_headers)

    vcf_df.to_csv(vcf_updated_path, sep='\t', mode='a', index=False)

def main(sp_num):

    final_df = match_and_update(sp_num)
    write_to_vcf_files(f'{sp_num}_vep.vcf',f'{sp_num}_updated.vcf', final_df)
