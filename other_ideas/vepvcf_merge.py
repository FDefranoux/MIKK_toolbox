import vcf
import pandas as pd
import os
import time
import argparse
import sys
<<<<<<< HEAD

vcf_file = 'vcf/chr1.vcf.gz'
vep_file = 'vep_tox/chr1_vep.txt'
sample_file = 'vcf/cram_file_to_line_ids.txt'
=======
from utils_bigdata import read_bigdata, sample_ID_set

vcf_file = 'datafiles/chr1.vcf.gz'
vep_file = 'datafiles/vep_chr1.txt'
sample_file = 'datafiles/cram_file_to_line_ids.txt'
>>>>>>> 187a5d9d0b0f8c472df8a1a9d8e8b3688630f96c


def arg_parsing():
    parser = argparse.ArgumentParser(description='VCF-VEP_merging')

    parser.add_argument('-f', '--vcf_file', help='VCF file to be used for '
                        + 'analysis', type=str,
                        default='/nfs/research/birney/projects/indigene/datafreeze/14-03-2019/vcf/medaka_inbred_panel_ensembl_new_reference_release_94.vcf.gz')

    parser.add_argument('-p', '--vep_file', help='VEP file to be used for '
                        + 'analysis', type=str,
                        default='/nfs/research/birney/projects/indigene/datafreeze/14-03-2019/vep/vep_medaka_new_reference_release_94.txt.gz')

    parser.add_argument('-t', '--sample_file', help='Text file to help '
                        + 'translating the common line names into the vcf one.',
                        type=str, default='vcf/cram_file_to_line_ids.txt')

    args = parser.parse_args()
    print(str(time.asctime()))
    print('The arguments are: ', vars(args), '\n')
    sys.stderr.write('\n' + str(time.asctime()) + ' -- Arguments: ' + str(vars(args)) + '\n')
    return args

<<<<<<< HEAD

def read_bigdata(file_name, cols=None, sep='\t'):
    if cols:
        target_cols = cols.keys()
        target_colnames = cols.values()
        file = pd.read_table(file_name, sep=sep, header=None, usecols=target_cols, names=target_colnames, comment='#')
    else:
        file = pd.read_table(file_name, sep=sep, header=None, comment='#')
    return file


def vcf_genome_id(vcf_file, samples_dict, chunk_size, vcf_info):
    cols = samples_dict.keys()
    names = samples_dict.values()
    if 'Genome_vepID.csv' in os.listdir():
        print('Removing previous file')
        os.remove('Genome_vepID.csv')
    pd.DataFrame(list(names) + ['ID']).T.to_csv('Genome_vepID.csv', header=None, index=False)
    for chunk in pd.read_table(vcf_file, sep='\t', usecols=cols,
                               names=names, comment='#', header=None,
                               chunksize=chunk_size):
        chunk.replace(':.*', '', regex=True, inplace=True)
        genome_vepid = pd.merge(chunk, vcf_info.astype(object), on=['#CHROM', 'POS'])
        genome_vepid.to_csv('Genome_vepID.csv', mode='a', header=False, index=False)
    print('CSV saving ended')



=======
################################################### FUNCTIONS FOR ID CALCULATION
>>>>>>> 187a5d9d0b0f8c472df8a1a9d8e8b3688630f96c
def variant_type(vcf):
    vcf_b = vcf[['REF', 'ALT']].copy()
    vcf_b['TYPE'] = ''
    vcf_b.loc[vcf['REF'].str.fullmatch('[ATGC]')
              & vcf['ALT'].str.fullmatch('[ATGC]'), 'TYPE'] = 'SNV'
    vcf_b.loc[vcf['REF'].str.fullmatch('[ATGC]') == False, 'TYPE'] = 'DEL'
    vcf_b.loc[vcf['ALT'].str.fullmatch('[ATGC]') == False, 'TYPE'] = 'INS'
    # Needs to be at the end not to be undone by DEL and INS
    vcf_b.loc[vcf['ALT'].str.contains(','), 'TYPE'] = 'MULT'
    return vcf_b['TYPE']


def vcf_vepID(vcf, method=None):
    bis = vcf[['ALT', 'POS', 'REF', '#CHROM']].copy()
    bis['ALT'] = vcf['ALT'].str.replace(pat=',', repl='/', regex=True)
    # bis['ALT'] = bis['ALT'].str.replace(pat='\*', repl='-', regex=True)
    if method == 'indel':
        # We need to remove the ref nucleotide to match VEP syntax
        bis['ALT'] = [s[1:] if len(s) > 1 else '-' for s in bis['ALT']]
        bis['REF'] = [s[1:] if len(s) > 1 else '-' for s in bis['REF']]
        bis['POS'] = bis['POS'].astype(int) + 1
    elif method == 'multiple':
        # We need to remove the ref nucleotide to match VEP syntax
        #   in each variant reported in ALT (separated by / in vcf file)
        alt = bis['ALT'].str.split('/', expand=True)\
                        .replace(regex=r'^[ATGC]', value='')\
                        .replace(regex='', value='-').fillna('#')
        bis['ALT'] = alt.apply(lambda x: '/'.join(x), axis=1)\
                        .replace(regex='/#', value='')
        bis['REF'] = bis['REF'].str[1:].replace(regex='', value='-')
        bis['POS'] = bis['POS'].astype(int) + 1

    bis['ID'] = bis['#CHROM'].astype(str) + '_' + bis['POS'].astype(str)\
        + '_' + bis['REF'] + '/' + bis['ALT']
    return bis['ID']


<<<<<<< HEAD
def ID_calculation(vcf):
    # Setting up types
    vcf.loc[:, 'TYPE'] = variant_type(vcf)
    # Normal way of calculating the ID for SNVs and MULT
    vcf.loc[vcf['TYPE'].isin(['MULT', 'SNV']),
            'ID'] = vcf_vepID(vcf[vcf['TYPE'].isin(['MULT', 'SNV'])])
    # Adapted calculation for insertion and deletions
    vcf.loc[vcf['TYPE'].isin(['DEL', 'INS']),
            'ID'] = vcf_vepID(vcf[vcf['TYPE'].isin(['DEL', 'INS'])],
                                  method='indel')
    return vcf


def variant_check(vep_id, vcf):
=======
def variant_check(vcf, vep_id):
>>>>>>> 187a5d9d0b0f8c472df8a1a9d8e8b3688630f96c
    recognized_variants = vcf[(vcf['ID'].isin(vep_id))]
    if len(vcf) == len(recognized_variants):
        print('Horray, all variant recognized!')
    else:
        print('More calculation needed...Please wait')
        vcf.loc[(vcf['ID'].isin(vep_id)) == False,
                    'ID'] = vcf_vepID(vcf[(vcf['ID'].isin(vep_id)) == False],
                                      method='multiple')
        recognized_variants = vcf[(vcf['ID'].isin(vep_id))]
<<<<<<< HEAD
        if len(vcf) == len(recognized_variants):
=======
        vcf_nonrecog = vcf[(vcf_info['ID'].isin(vep_id)) == False]
        if len(vcf) == len(recognized_variants) & len(vcf_nonrecog) == 0:
>>>>>>> 187a5d9d0b0f8c472df8a1a9d8e8b3688630f96c
            print('Horray, all variant recognized!')
        else:
            print('Problems with some variants not solved...' )
            y = len(recognized_variants)*100 / len(vcf)
            print(f'{y}% of the variant recognized')
<<<<<<< HEAD


def alternative_variant_table(vcf, start_pos = '0', ref_seq='', sequence=True):
    # Recuperate each ALT possible per position in differents columns
    alt = vcf.pop('ALT').astype(str)
    alt = alt.str.split(',', expand=True)
    # Merging the corresponding position/reference per index
    merged = pd.merge(vcf, alt, right_index=True, left_index=True)
    # Renaming columns to have a number for each alternative position
    merged['REF_copy'] = merged['REF'].copy()
    rename_dict = {n: f'{n+1}' for n in range(alt.shape[1])}
    if not sequence:
        del(alt)
    rename_dict['REF_copy'] = '0'  # Reference columns is equivalent to sequence 0
    merged.rename(columns=rename_dict, inplace=True)
    # Melting dataframe to have long format dataframe
    melted = pd.melt(merged, id_vars=vcf.columns,
                     var_name='ALT_num', value_name='single_alt')
    melted.dropna(axis=0, how='any', inplace=True)

    # Creation of column containing sequence from POS to POS+1, completing with
    # reference sequence when necessary.
    if sequence:
        overlapping=False
        seq_ls_ref = []
        for variant in melted.index:
            alt = melted.loc[variant, 'single_alt']  # Alternative sequence
            ref = melted.loc[variant, 'REF']  # Reference sequence
            end = int(melted.loc[variant, 'POS+1'])  # Next variant position
            pos = int(melted.loc[variant, 'POS'])
            # Start of reference sequence to complete variant sequence
            if len(ref) < len(alt):  # Case for insertions
                start = pos + 1
            else:  # Case for deletions
                start = pos + len(ref)
            # Appending the sequences
            if '*' in alt:
                seq_ls_ref.append('')
            elif start > end:  # If the alternative seq is including next variant
                overlapping = True
                seq_ls_ref.append(alt[:end-pos])  # Truncation of the alt seq
            else:  # If it is not, we complete the seq with ref sequence
                seq_ls_ref.append(alt + ref_seq[start-start_pos-1:end-start_pos-1])
        melted['SEQ'] = seq_ls_ref  # Adding corresponding seq to the table
        if overlapping:
            print('Beware, Overlapping variants in the sequences!')
    return melted


def main(vcf_file, vep_file, sample_file):

    # Reading the vcf file
    vcf_info = read_bigdata(vcf_file, cols={0: '#CHROM', 1: 'POS', 3: 'REF', 4: 'ALT'})
    # print(vcf_info.head(5))
    # vcf_info = read(vcf_file)
    # vcf_info = vcf_reader(vcf_file, samples=None)
    print(vcf_info.head(5), str(time.asctime()), '\n\n', flush=True)

    # Adjusting the columns in VEP and selecting just chr = 1 + counts
    vep_id = read_bigdata(vep_file, cols={0: '#Uploaded_variation'})
    vep_id = vep_id['#Uploaded_variation']  #.str.replace(pat='*', repl='-', regex=True)

    # Calculations
    vcf_new = ID_calculation(vcf_info)
    variant_check(vep_id, vcf_new) # Verification of the IDs between the files
    #print(vcf_new.head(), str(time.asctime()), '\n\n', flush=True)
    # Checks:
    vcf_nonrecog = vcf_info[(vcf_info['ID'].isin(vep_id)) == False]
    if len(vcf_nonrecog) != 0:
        vcf_nonrecog.to_csv('Vcf_nonrecog.csv')

    del(vep_id)
    vcf_info.drop(['REF', 'ALT', 'TYPE'], axis=1, inplace=True)

    # Drop REF ? ALT?



    # # Genome load
    sample_file = pd.read_table(sample_file)
    samples = sample_file['cram_file'].tolist()
    # VEP_VCF:0 add something to verify the samples are well attributed
    samples_dict = {0: '#CHROM', 1: 'POS'}
    samples_dict.update({n+9: sample for n, sample in enumerate(samples)})
    vcf_genome_id(vcf_file, samples_dict, 100000, vcf_info)
=======
            vcf_nonrecog.to_csv('Vcf_nonrecog.csv')


def ID_calculation(vcf, vep_id):
    # Setting up types
    vcf.loc[:, 'TYPE'] = variant_type(vcf)
    # Normal way of calculating the ID for SNVs and MULT
    vcf.loc[vcf['TYPE'].isin(['MULT', 'SNV']),
            'ID'] = vcf_vepID(vcf[vcf['TYPE'].isin(['MULT', 'SNV'])])
    # Adapted calculation for insertion and deletions
    vcf.loc[vcf['TYPE'].isin(['DEL', 'INS']),
            'ID'] = vcf_vepID(vcf[vcf['TYPE'].isin(['DEL', 'INS'])],
                                  method='indel')
    variant_check(vcf, vep_id) # Verification of the IDs between the files
    return vcf


######################################################### CREATING NEW DATAFILES
def vcf_genome_id(vcf_file, samples_dict, chunk_size, vcf_info):
    cols = samples_dict.keys()
    names = samples_dict.values()
    if 'Genome_vepID.csv' in os.listdir():
        print('Removing previous file')
        os.remove('Genome_vepID.csv')
    # Saving the headers
    pd.DataFrame(list(names) + ['ID']).T.to_csv('Genome_vepID.csv', header=None, index=False)
    for chunk in pd.read_table(vcf_file, sep='\t', usecols=cols,
                               names=names, comment='#', header=None,
                               chunksize=chunk_size):
        chunk.replace(':.*', '', regex=True, inplace=True)
        genome_vepid = pd.merge(chunk, vcf_info.astype(object), on=['#CHROM', 'POS'])
        genome_vepid.to_csv('Genome_vepID.csv', mode='a', header=False, index=False)
    print('CSV saving ended')






def main(vcf_file, vep_file, save_genome=True):
    # Reading the files
    vcf_info = read_bigdata(vcf_file, cols={0: '#CHROM', 1: 'POS',
                                            3: 'REF', 4: 'ALT'})
    vep_id = read_bigdata(vep_file,
                          cols={0: '#Uploaded_variation'})['#Uploaded_variation']

    # Calculations for new VCF ID
    vcf_info = ID_calculation(vcf_info, vep_id)
    vcf_info.drop(['REF', 'ALT', 'TYPE'], axis=1, inplace=True)

### HERE
    # Genome load
    samples_df = sample_ID_set()
    cols_dict = {0: '#CHROM', 1: 'POS'}
    cols_dict.update({n+9: sample for n, sample in enumerate(samples_df['cram_file'].tolist())})

    #Reads in chunk and save chunk by chunk not keeping in memory
    vcf_genome_id(vcf_file, cols_dict, 100000, vcf_info)






>>>>>>> 187a5d9d0b0f8c472df8a1a9d8e8b3688630f96c
    # vcf_genome = read_bigdata(vcf_file, cols=samples_dict)
    # vcf_genome = pd.merge(vcf_info, vcf_genome, on=['#CHROM', 'POS'])
    # del(vcf_info)
    #
    # # Arrange alt to have the right number
    # for sample in samples:
    #     vcf_genome[sample] = vcf_genome[sample].str[:3]
    # print(vcf_genome.head(5), str(time.asctime()), '\n\n', flush=True)
    # melt = alternative_variant_table(vcf_genome, sequence=False)
    # del(vcf_genome)
    # print(melt.head(5), '\n\n', flush=True)
    # melt[melt == './.'] = '0/0'

    # # Final merging and saving dataframe
    # final = pd.merge(vcf_genome, vep_id, on='ID')
    # print('Final assessment of correspondance: ', len(vep_id) - len(final))
    # print(final.head(), '\n\n', flush=True)
    # final.to_csv('VepVcfMerged_melted.csv')

    # # Arrange the SNPs allele information as usable values
    # melt.to_csv('VepVcfMerged_melted.csv')

if __name__ == '__main__':
    main(**vars(arg_parsing()))
