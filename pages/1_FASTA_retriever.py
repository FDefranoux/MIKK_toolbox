import os
import numpy as np
import pandas as pd
import streamlit as st
from streamlit_dynamic_filters import DynamicFilters
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

from pages.utils.alignment_utils import render_alignment, view_alignment, save_bokeh_figure
st.set_page_config(
    page_title="MIKK Toolbox",
    page_icon="🐟",
    layout="wide",
)


# ##### F2 files

# VCF_FILE = '/home/fanny/Documents/Work/HeartRate_Tox/files_for_interactiveplots/region12.vcf.gz'
# FASTA_DIR = '/home/fanny/Documents/Work/fasta_medaka_panel/chr_fasta_ref/'
# SAMPLE_FILE = '/home/fanny/Documents/Work/fasta_medaka_panel/cram_file_f2s_philip_filters.txt'  # Need at least line and cram_file columns

##### FO files

VCF_FILE = '/home/fanny/Documents/Work/fasta_medaka_panel/region12_F0.vcf.gz'
FASTA_DIR = '/home/fanny/Documents/Work/fasta_medaka_panel/chr_fasta_ref/'
SAMPLE_FILE = '/home/fanny/Documents/Work/Database_MIKK/igv_tests/igv_flask/igv_flexlmm/variant_to_sequence/vcf/cram_file_to_line_ids.txt'  # Need at least line and cram_file columns



def samples_definition(sample_df, samples_ls=[], new_name='sample'):
    """
    Description: Permits to translate name of the samples used in the VCF file
    to the one commonly used.
    Required: At least two columns one `cram_file` and one `line`.

    Args:
        - sample_df (table): dataframe.
        - samples_ls, optionnal (list): List of samples names to be used and
        translated.

    Returns:
        - dictionnary of equivalent names.
    """
    if samples_ls:
        samples = sample_df.loc[(sample_df['sample'].isin(samples_ls))
                                | (sample_df.index.isin(samples_ls))].to_dict()
    else:
        samples = sample_df.to_dict()
    return samples[new_name]


@st.cache_data()
def get_vcf(vcf_file, chr, start, end):
    """
    Description: Function that will retrieve a subset of the VCF file according to parameters.
    Required: At least one of the package `pysam` or `pyvcf`.

    Args:
        - chr : region name corresponding to the chromosome.
        - start, end: position of the start and end of targetted subregion.

    Returns:
        - pd.DataFrame of the VCF information for this subset.
    """

    # Get VCF 
    vcf_df = pd.DataFrame()
    try:
        from pysam import VariantFile
        reader = VariantFile(vcf_file)
        # st.write(reader.header)
        for rec in reader.fetch(str(chr),max(0, int(start)-1), int(end)):
            info_pos = {'POS':rec.pos, 'REF':rec.ref, 'ALT':'/'.join(rec.alts)}
            info_pos.update({rec.samples[sample].name:str(rec.samples[sample]['GT'])[1:-1].replace('None, None', '0, 0') for _, sample in enumerate(rec.samples)})
            vcf_df= pd.concat([vcf_df, pd.DataFrame(info_pos, index=[rec.pos])])
    except ModuleNotFoundError:
        import vcf
        reader = vcf.Reader(filename=vcf_file)
        for rec in reader.fetch(str(chr),  max(0, int(start)-1), int(end)):
            info_pos = {'POS':rec.POS, 'REF':rec.REF, 'ALT':'/'.join(rec.ALT)}
            info_pos.update({rec.samples[sample].sample:str(rec.samples[sample]['GT'])[1:-1] for sample, _ in enumerate(rec.samples)})
            vcf_df= pd.concat([vcf_df, pd.DataFrame(info_pos, index=[rec.POS])])
    # st.write(vcf_df)
    return vcf_df


def descriptive_info(vcf_df):
    """
    Description: Assign variant type and return descriptive table for the VCF.

    Args:
        - vcf_df (table): dataframe.

    Returns:
        - dictionnary of some descriptive measures.
    """
    vcf = vcf_df.drop('POS+1', axis=1).melt(id_vars=['REF', 'POS', 'ALT'], var_name='Sample', value_name='Genotype')
    vcf['type'] = 'None'
    vcf.loc[(vcf['REF'].isin(['A', 'T', 'G', 'C'])) & (vcf['ALT'].isin(['A', 'T', 'G', 'C'])), 'type'] = 'SNP'
    vcf.loc[(vcf['REF'].str.len() > 1), 'type'] = 'DEL'
    vcf.loc[(vcf['ALT'].str.len() > 1), 'type'] = 'INS'
    if len(st.session_state['sample_dict'].keys()) > 0:
        vcf['id'] = vcf['Sample'].map(st.session_state['sample_dict'])
    else:
        vcf['id'] = vcf['Sample']
    vcf.drop(['Sample'], axis=1, inplace=True)
    return {'Genotypes per type of variants': vcf.groupby(['POS', 'type', 'Genotype'])['id'].nunique().astype(int).unstack().round(0), 
            'Genotypes per samples': vcf.groupby(['id', 'Genotype'])['POS'].nunique().astype(int).unstack().round(0)}


@st.cache_data()
def recup_sequence(chr, start_pos, end_pos, vcf_file, fasta_dir, samples_ls=[]):
    """
    Description: Function that will retrieve sequences according to the VCF file.

    Args:
        - chr : region name corresponding to the chromosome.
        - start_pos, end_pos: position of the start and end of targetted subregion.
        - vcf_file: path where to find VCF file containing required information.
        - fasta_dir: path of the directory where to find fasta files of the reference.
        - samples_ls, optionnal (list): List of samples names to filter on.
        
    Returns:
        - dictionnary of sequence for each targetted samples.
        - tuple signalling the new range of the sequence (can be modified in case of INDELs and/or previous SNPs)
        - dictionnary of descriptive measures. 
    """
    # Adjusting position number to python indexing
    if start_pos > 0:
        start_pos = start_pos-1
        end_pos = end_pos-1
    if start_pos > end_pos: st.error('Start position greater than end position')

    # Recuperation of the VCF file and position of SNP before
    vcf_before = pd.DataFrame()
    before_pos = start_pos
    while vcf_before.empty:
        before_pos = before_pos - 1000
        if before_pos < 1:
            vcf_before = get_vcf(vcf_file, chr, 0, start_pos)
            break
        vcf_before = get_vcf(vcf_file, chr, before_pos, start_pos)

    vcf_df = get_vcf(vcf_file, chr, start_pos, end_pos)    
    # st.write(vcf_df.head())
    if vcf_df.shape[0] < 0 : st.error('Region not recognized in VCF')
    # assert vcf_df.shape[0] > 0, 'No SNPs have been recognized in the VCF for this region!'

    try:
        if not vcf_before.empty:
            pos_before = int(vcf_before[vcf_before['POS'].astype(int) < start_pos].iloc[-1]['POS']) -1 # Last variant before the sequence of interest 
            len(vcf_before.iloc[-1]['ALT'])
            if start_pos <= pos_before + max(len(vcf_before.iloc[-1]['REF'])+1, len(vcf_before.iloc[-1]['ALT'])+1):
                vcf_df = pd.concat([pd.DataFrame(vcf_before.iloc[-1]).T, vcf_df], ignore_index=True).drop_duplicates()
            else:
                pos_before = start_pos
        else:
            pos_before = start_pos
    except Exception as err:
        st.exception(err)
        pos_before = start_pos

    # Filtering for the samples
    if not samples_ls: 
        samples_ls = vcf_df.drop(['POS', 'REF', 'ALT'], axis=1).columns.tolist()
    if set(vcf_df.columns) & set(samples_ls) != set(samples_ls):
        st.warning('The following samples have not been recognized in the VCF\n'+ '\n '.join(list(set(samples_ls) -  set(vcf_df.columns))) )
    samples_ls = [samp for samp in samples_ls if samp in vcf_df.columns]
    assert len(samples_ls) > 0, 'No samples has been recognized in the VCF!'
    vcf_df = vcf_df[['REF', 'ALT', 'POS'] + samples_ls]  # Filtering to asked samples
    st.write(vcf_df.rename(columns=st.session_state['sample_dict']))

    # Formatting the df
    if vcf_df[vcf_df != '0, 0'].drop(['POS', 'REF', 'ALT'], axis=1).dropna(how='all').dropna(how='all', axis=1).empty:  # If no SNPs found
        # Reference sequence from initial
        ref_seq = ref_sequence_at_pos(fasta_dir, chr, start_pos, end_pos)
        dict_seq = {st.session_state['ref_line'] :ref_seq}
        return dict_seq, (start_pos, end_pos), None

    else:
        # st.write(vcf_df)
        vcf_df = vcf_df[(vcf_df['POS'].astype(int) >= pos_before) & (vcf_df['POS'].astype(int) < end_pos)]
        # vcf_df[['REF', 'ALT']] = vcf_df['ID'].str.split('_', expand=True)[2].str.split('/', expand=True)[[0,1]]
        # if len(vcf_df[['REF', 'POS', 'ALT']].dropna()) != 0:
        #     end_pos = int(vcf_df.iloc[0]['POS'])
        vcf_df.replace('./.', '0/0', inplace=True, regex=False)
        vcf_df.dropna(how='any', inplace=True)
        first_variant = vcf_df['POS'].sort_values().iloc[0]

        # Add the reference sequence
        vcf_df['HdR'] = '0, 0'  
        samples_ref = samples_ls + ['HdR']
        st.session_state['sample_dict']['HdR'] = 'HdR'
        vcf_df['POS+1'] = vcf_df['POS'].tolist()[1:] + [end_pos+1]  # give the position of following mut

        # Sequence table
        ref_seq = ref_sequence_at_pos(fasta_dir, chr, pos_before, end_pos)
        alt_table = alternative_sequence_table(vcf_df[['ALT', 'REF', 'POS', 'POS+1']], pos_before, end_pos, ref_seq) 
        alt_table.sort_values(by=['POS', 'ALT_num'], inplace=True)

        # Reconstitution for each sample sequences per allele
        allele1 = vcf_df[samples_ref].replace(', .', '', regex=True).astype(int)
        allele2 = vcf_df[samples_ref].replace('., ', '', regex=True).astype(int)

        # Formating the sequences
        dict_seq = {}
        if len(st.session_state['sample_dict'].keys()) < 2:
            st.session_state['sample_dict'].update({key:key for key in samples_ls})
        for sample in samples_ref:
            sequence1 = ref_seq[:first_variant - pos_before -1] + variant_to_seq(vcf_df['POS'], allele1[sample],
                                        alt_table)
            sequence2 = ref_seq[:first_variant - pos_before -1] + variant_to_seq(vcf_df['POS'], allele2[sample],
                                        alt_table)
            # st.write(sequence1, sequence2)
            if not sample == st.session_state['ref_line']:
                dict_seq.update({st.session_state['sample_dict'][sample] + '_all1':sequence1})
                dict_seq.update({st.session_state['sample_dict'][sample] + '_all2':sequence2})
            else: 
                dict_seq.update({st.session_state['sample_dict'][sample]: sequence1})
        
        # desc_dict = descriptive_info(vcf_df)

        return dict_seq, (pos_before, end_pos), descriptive_info(vcf_df)



@st.cache_data()
def ref_sequence_at_pos(fasta_dir, chr, start_pos, end_pos):
    """
    Recuperation of the sequence from a fasta file at given chromosome and
    positions.

    Args:
        - fasta_dir (path): path of the directory containing the fasta file for
        each chromosomes.
        - chr (int): chromosome to target.
        - start_pos (int): positions from which we want vcf informations.
        - end_pos (int): positions to which we want vcf informations.

    Returns:
        ref_seq (string): Sequence at given chromosome and positions.

    Exception:
    - If the first line from the fasta file do not correspond to the chromosome
    required, will just send a warning and continue working.
    - If the position is out of bonds, returns the sequence until the
    end of the chromosome.
    """

    # Select the chromosome specific fasta file
    fasta_file = os.path.join(fasta_dir, f'{str(chr)}.fa')
    # Open the fasta file
    with open(fasta_file, 'r', encoding="utf-8") as fasta_content:
        # Verification that the first line correspond to the chr demanded
        chr_read = fasta_content.readlines(1)[0][1:-1]
        if chr_read != str(chr):
            st.warning('Wrong Fasta file!, reading fasta file from chromosome'
                  + f'{chr_read}')
        # Select the sequence as list of lines
        fasta = fasta_content.readlines()
        # Join the list in one string (removing the line jumps)
        fasta = ''.join([f[:-1] if '\n' in f else f for f in fasta])
        # Check that the end_pos is included in the fasta length
        if start_pos > len(fasta):
            st.error('Positon not existing for this chromosome !')
            ref_seq=''
        if end_pos > len(fasta):
            st.warning('End-Position Out of Bonds !')
            ref_seq = fasta[start_pos:]
        else:
            ref_seq = fasta[start_pos:end_pos]
    return ref_seq


@st.cache_data()
def alternative_sequence_table(vcf, start_pos, end_pos, ref_seq, sequence=True):
    """
    From vcf files, recuperates the alternatives sequences for each position.
    Starts from the alternative sequence and complete the sequence with the
    referenced one until next variant position.

    Args:
        - vcf (pandas.DataFrame): vcf table containing the variants in target
        chromosome at required positions. Need at least ALT, REF, POS columns.
        - start_pos (int): position from which sequence starts.
        - end_pos (int): position were sequence needs to end.
        - ref_seq (str): reference sequence to be used to complete the
        alterntive sequences.
        - sequence (bool): To get the theoretical sequence until next position

    Returns:
        melted (Pandas.DataFrame): Long-format table containing sequences for
        each alternative variant from its position to the next one.
    """
    # Recuperate each ALT possible per position in differents columns
    alt = vcf['ALT'].astype(str)
    # alt = alt.str[1:-1]
    alt.replace(' ', '', regex=True, inplace=True)
    alt = alt.str.split(',', expand=True)
    # Merging the corresponding position/reference per index
    merged = pd.merge(vcf[['POS', 'REF', 'POS+1']], alt, right_index=True,
                      left_index=True)
    # Renaming columns to have a number for each alternative position
    merged['REF_copy'] = merged['REF'].copy()
    rename_dict = {n: f'{n+1}' for n in range(alt.shape[1])}
    rename_dict['REF_copy'] = '0'  # Reference columns is equivalent to sequence 0
    merged.rename(columns=rename_dict, inplace=True)
    # Melting dataframe to have long format dataframe
    melted = pd.melt(merged, id_vars=['POS', 'POS+1', 'REF'],
                     var_name='ALT_num', value_name='single_alt')
    melted.dropna(axis=0, how='any', inplace=True)
    # Creating columns that render position + alternative sequence number
    melted['ALT_POS'] = melted['POS'].astype(str) + '_'\
        + melted['ALT_num'].astype(str)

    # Creation of column containing sequence from POS to POS+1, completing with
    # reference sequence when necessary.
    overlapping=False
    if sequence:
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
            st.warning('Beware, Overlapping variants in the sequences! Sequence could end up erroneous')

    return melted


@st.cache_data()
def variant_to_seq(pos, genot, alt_table):
    """
    For a table containing all variants alleles, with the corresponding
    alternative sequence table, returns the sequence for required samples.

    Args:
        - pos (pandas.Serie): POS column of the vcf table containing the
        variants alleles in target chromosome at required positions.
        - genot (pandas.Series): Serie containing the allele number for each
        variant
        - alt_table (pandas.DataFrame, long format): output of the function
        alternative_sequence_table().

    Returns:
        Sequence (str)
    """
    #  Combining allele number and position as ALT_POS in alt_table

    alt_pos1 = pos.astype(str) + '_' + genot.astype(str)
    #  Recuperating all related sequece in alt_table
    all1_seq_ls = alt_table[alt_table['ALT_POS'].isin(alt_pos1)]['SEQ']
    #print(alt_table[alt_table['ALT_POS'].isin(alt_pos1)].to_markdown())
    seq_allele = ''.join(all1_seq_ls)

    return seq_allele


def alignement_res(dict_seqs, ref='HdR', start=0, title=''):
    """
    For a dictionnary of sequences. Compute the alignments pairwise then display multiple alignements.

    Args:
        - dict_seqs (dictionnary): Containing the id of sequence and the sequences. 
        - ref: Id of the reference sequence to be use for general alignments
        - start: Permits to align with the sequence coordinates.
        - title: Title for plot and saving.
    """
    ref_seq = dict_seqs.pop(ref) 

    dict_seqs = {key:val for key, val in sorted(dict_seqs.items(), key=lambda x: len(x[1]), reverse=True)}    # Ordering the sequences from len
    # align_scores = dict()

    #First run of alignments to get the longest refseq aligned sequence and saving alignment score.
    for key, seq in dict_seqs.items():
        refalign_seq = render_alignment(ref_seq, seq)
        # align_scores[key], refalign_seq = alignements(ref_seq, seq)
        # st.write(refalign_seq, align_scores[key])
        refalign_seq = refalign_seq.replace(' ', '').split('\n')
        if len(refalign_seq[0]) > len(ref_seq):
            ref_seq = refalign_seq[0].replace(' ', '') 

    #Second run of alignments to get the alignments from the longest refseq aligned sequence.
    results = '>' + ref + '\n' + ref_seq + '\n' 
    list_seq = []
    dict_seqs = {key:val for key, val in sorted(dict_seqs.items(), key=lambda x: x[0], reverse=True)} # Ordering the sequences from line name
    for key, seq in dict_seqs.items():
        align_seq = render_alignment(ref_seq, seq)
        # _, align_seq = alignements(ref_seq, seq)
        align_seq = align_seq.replace(' ', '').split('\n')#.replace(' ', '')
        results += '>' + key + '_align-score:' + '\n' + align_seq[1] + '\n' #str(round(align_scores[key], 2)) 
            
        list_seq.append(SeqRecord(Seq(align_seq[1]), id=key))
    # Saving
    save_cols = st.columns(2)
    save_cols[1].download_button('Download alignments (text format)', results, file_name=f'FastaMedakaPanel{title}_alignments.txt', use_container_width=True)

    # Preparing alignment figure 
    try:
        align = MultipleSeqAlignment([SeqRecord(Seq(align_seq[0]), id=ref + '_REF')] + list_seq)
        p = view_alignment(align, plot_width=1100, start=start, ref=ref + '_REF')
        param_saving = dict(figure=p, filename = f'FastaMedakaPanel{title}_alignment_figure.html', title=title)
        st.bokeh_chart(p)
        save_cols[0].button('Download alignment figure (html format)', on_click=save_bokeh_figure, kwargs=param_saving, use_container_width=True)

    except Exception as err:
        st.exception(err)
    


def save_seqs_as_fasta(dict_seqs):
    """
    Description: Function to save sequence in fasta format.
    """
    seqs = [f'>{key}\n{item}'  for key, item in dict_seqs.items()]
    return '\n'.join(seqs)
    

def streamlit_params():
    """
    Description: Permits to initialize the parameters and session_state to run the streamlit app.
    """
    if 'submit_fasta' not in st.session_state:
        st.info('Submit the data you want to analyse on the left-side sidebar.')
        st.session_state['submit_fasta'] = False

    if len(set(['output_dir', 'vcf_file', 'fasta_dir', 'sample_file']) - set(st.session_state.keys())) > 0: 
        with st.expander('Files parameters').form('General params', border=False):
            st.header('Enter paths manually')
            params = {}
            for var in ['output_dir', 'vcf_file', 'fasta_dir', 'sample_file']:
                if (var not in st.session_state) or (st.session_state[var] == ''):
                    params[var] = st.text_input(var.replace('_', ' ').replace('dir', 'directory').title(), value='')
            submit_params = st.form_submit_button('Submit', use_container_width=True)
            if submit_params:
                st.session_state.update(params)
    if not 'sample_dict' in st.session_state:
        st.session_state['sample_dict'] = {}
    if 'sample_file' in st.session_state:
        if st.session_state['sample_file'] != '':
            with st.sidebar.expander('Filter parameters').form('Filter Form', border=False):
                sample_df = pd.read_table(st.session_state['sample_file'], index_col='cram_file')
                dynamic_filters = DynamicFilters(sample_df, filters=sample_df.columns.tolist())
                dynamic_filters.display_filters(location='columns', num_columns=2)
                submit_filter_samples = st.form_submit_button('Submit filters')
            if submit_filter_samples:
                st.session_state['sample_dict'] = samples_definition(dynamic_filters.filter_df())
                st.write(st.session_state['sample_dict'])
                st.session_state['sample_ls'] = st.session_state['sample_dict'].keys()


    with st.sidebar.form('Genome browser for MIKK panel'):
        sample_dict = st.session_state['sample_dict']
        cols = st.columns(3)
        chr = cols[0].number_input('Chromosome', 1, 25, value=1)
        start = cols[1].number_input('Start position', 0, step=1000)
        end = cols[2].number_input('End position', 1000, step=1000)
        st.session_state['ref_line'] = 'HdR'
        submit = st.form_submit_button('Submit parameters')
        if submit:
            st.session_state.update({'chr': chr, 'start_pos': start, 'end_pos':end, 
                                        'samples_ls':list(sample_dict.keys()), 'submit_fasta':True})


def main():
    st.header('FASTA retriever tool')
    streamlit_params()
    if st.session_state['submit_fasta']: 
        st.write(st.session_state['sample_dict'])
        seqs, new_range, desc_dict = recup_sequence(chr=st.session_state['chr'], start_pos=st.session_state['start_pos'], end_pos=st.session_state['end_pos'], samples_ls=st.session_state['samples_ls'], vcf_file=st.session_state['vcf_file'], fasta_dir=st.session_state['fasta_dir'])
        if (len(seqs.keys()) == 1) & ('HdR' in seqs.keys()):
            st.warning(f'No SNPs found in this region for those filters. Displaying the reference sequence (HdR) for the region {new_range[0]+1} and {new_range[1]+1}bp: \n')
            # st.code(seqs['HdR'])
        else:
            tab_seq, tab_plot, tab_desc = st.tabs(["Sequences", "Alignment plot", "Descriptive table"])
            title=f'_{st.session_state["chr"]}_{st.session_state["start_pos"]}_{st.session_state["end_pos"]}'

            with tab_seq:
                st.success(f'The sequence represents the region between the variants at {new_range[0]+1} and {new_range[1]+1}bp\n\n')
                for line, seq in seqs.items():
                    st.markdown(f'### {line}')
                    st.code(seq)
                st.download_button('Download sequence (fasta format)', save_seqs_as_fasta(seqs), file_name=f'FastaMedakaPanel{title}_alignments.txt', use_container_width=True)
            
            with tab_plot:
                st.success(f'The plot represents the alignments in the region between the variants at {new_range[0]+1} and {new_range[1]+1}bp\n\n')
                alignement_res(seqs, ref=st.session_state['ref_line'], start=new_range[0], title=title)
            
            with tab_desc:
                # st.write(desc_dict)
                desc_cols = st.columns(len(desc_dict.keys()))
                for n, (key, table) in enumerate(desc_dict.items()):
                    desc_cols[n].markdown(f'##### {key}')
                    desc_cols[n].table(table)
            

if __name__ == '__main__':
    main()

