import sys
import os
import numpy as np
import pandas as pd
import socket
import subprocess
from io import BytesIO
# import seaborn as sns
# import matplotlib.pyplot as plt
# from matplotlib.collections import PathCollection
# import matplotlib.patches as mpatches
import streamlit as st
import plotly.express as px
from streamlit_plotly_events import plotly_events
from st_btn_select import st_btn_select
import glob
import argparse
from plotly.graph_objs._figure import Figure
from _io import TextIOWrapper, BufferedReader
import pandas as pd
import sys
import plotly.io as pio
from streamlit_dynamic_filters import DynamicFilters
pio.templates.default = "simple_white"
st.set_page_config(layout="wide")

#IDEA: Save specific loci in wanted file + Load the selection you want

# @st.cache_data()
def initialization(folder, phenotype_file, covariates_file):

    # Initiation of the state session
    st.session_state['submitted'] = True
    st.session_state['folder'] = folder

    # Reading tables: 
    gwas_file = os.path.join(folder, st.session_state['gwas_phenotype'] , f'input_{st.session_state["chr"]}_{st.session_state["gwas_phenotype"] }.gwas.tsv.gz')
    print('\n\n', gwas_file)
    df_gwas = pd.read_table(gwas_file)
    df_gwas['log_p'] = -np.log10(df_gwas['lrt_p'])
    df_gwas['selected'] = 'False'
    st.session_state['gwas_df'] = df_gwas

    pheno_df = pd.read_table(phenotype_file)
    cov_df = pd.read_table(covariates_file)
    df = pd.merge(pheno_df, cov_df, on='#IID')
    df['vcf_ID'] = df['#IID'] + '_' + df['#IID']
    st.session_state['covariates'] = cov_df.drop('#IID', axis=1).columns.tolist()
    st.session_state['phenotypes'] = pheno_df.drop('#IID', axis=1).columns.tolist()
    st.session_state['pheno_df'] = df

    # Initialisation
    st.session_state['pos_ls'] = dict(finemapping=[], highest=[], manual=[], loaded=[])
    if not os.path.isdir( os.path.join(st.session_state['folder'], '.temp')):
        os.makedirs(os.path.join(st.session_state['folder'], '.temp'))
    st.session_state['pos_file'] = os.path.join(st.session_state['folder'], '.temp', f'list_pos_{st.session_state["chr"]}.temp')

    if not os.path.exists(st.session_state['pos_file']):
        open(st.session_state['pos_file'], 'w').close()
        st.session_state['selected_points'] = []
    else:
        try:
            st.session_state['pos_ls']['loaded'] = pd.read_table(st.session_state['pos_file'], header=None)[0].tolist()
            st.session_state['selected_points'] = st.session_state['pos_ls']['loaded']
            update_selection(mode='a')

        except pd.errors.EmptyDataError:
            st.session_state['selected_points'] = []
        
    st.session_state['selection_n'] = 0


def selection_options():
    # SIDEBAR options for the selected points
    option_cols = st.columns(4)
    if option_cols[0].button('Reset all', use_container_width=True):
        confirmation_box(container=st)
    if option_cols[1].button('Append selected', use_container_width=True):
        update_selection(mode='a')
    if option_cols[2].button('Replace selected', use_container_width=True):
        update_selection(mode='w')
    if option_cols[3].button('Save list of selected positions', use_container_width=True):
        update_selection(mode='s')


def confirmation_box(container=st):
    confirmation_container = container.empty()
    with confirmation_container.container(border=True):
        container.error("Are you sure?")
        col_conf = container.columns(2)
        yes_btn = col_conf[0].button("Yes", on_click=update_selection, kwargs=dict(mode='rm'))
        cancel_btn = col_conf[1].button('Cancel')


def update_selection(mode='a'):
    list_pos = [el for poslist in st.session_state['pos_ls'].values() for el in poslist]
    if mode == 'w':
        st.session_state['selected_points'] = list_pos
        st.session_state['gwas_df']['selected'] = False
        for key, val in st.session_state['pos_ls'].items():
             st.session_state['gwas_df'].loc[st.session_state['gwas_df']['pos'].isin(val), 'selected'] =  key
    elif mode == 'a':
        st.session_state['selected_points'] = list(set(list_pos + st.session_state['selected_points']))
        for key, val in st.session_state['pos_ls'].items():
             st.session_state['gwas_df'].loc[st.session_state['gwas_df']['pos'].isin(val), 'selected'] =  key
    elif mode == 'rm':    
        st.session_state['selected_points'] = []
        st.session_state['pos_ls'] = dict(finemapping=[], highest=[], manual=[], loaded=[])
        st.session_state['gwas_df']['selected'] = False
        st.session_state['action_submit'] == False
    elif mode == 's':
        pd.Series(st.session_state['selected_points']).drop_duplicates().to_csv(st.session_state['pos_file'], 
            sep='\t', mode='w', index=False, header=False)


def plot_gwas_association(df, color='selected', selection=False, key=''):
    n=1
    color_dict={"False":'#949191', "highest":"red", "finemapping":"orange", "manual":"blue", "loaded":'purple'}
    df['selected'] = df['selected'].astype(str)
    fig = px.scatter(data_frame=df, x='pos', y='log_p', color=color, color_discrete_map=color_dict, opacity=0.5, title='GWAS ' + st.session_state['gwas_phenotype'])
    while selection:
        selected_points = plotly_events(fig, click_event=True, hover_event=False, select_event=True, key=key + str(n))
        if selected_points:
            selected_points = pd.DataFrame(selected_points)
            st.session_state['pos_ls']['manual'] = selected_points['x'].tolist()
            df.loc[df['pos'].isin(st.session_state['pos_ls']['manual']), 'selected'] = 'manual'
            selected_points = None
            n += 1                   
            fig = px.scatter(data_frame=df, x='pos', y='log_p', color=color, color_discrete_map=color_dict, opacity=0.5, title='GWAS ' + st.session_state['gwas_phenotype'])
    st.session_state['gwas_df'] = df
    st.session_state['gwas_fig'] = fig
    return fig


def tab_function():
    tabs = st.tabs([s.center(20,"\u2001") for s in ['Finemapping', 'Highest', 'Manual selection']])
    with tabs[0]:
        # Finemapping
        finemapping = st.slider('finemapping', min_value=0, max_value=10000000, value=0, step=100000)
        if finemapping > 0:
            lead_snp(st.session_state['gwas_df'], finemapping=finemapping)  # Cached   

    with tabs[1]:
        high_cols = st.columns([2,1,1,1])
        num = high_cols[0].slider('Number of SNPs', min_value=0, max_value=10, step=1, value=0)
        with high_cols[1]:
            th_selection = st_btn_select(('Permutation', 'Bonferroni correction', 'Custom'))
            if th_selection == 'Permutation':
                thresh = 0
            elif th_selection == 'Bonferroni correction':
                thresh = 0.05/st.session_state['gwas_df']['pos']
            else:
                thresh = st.number_input('Threshold', min_value=0.0, max_value=1.0, step=0.05)
        
        if num > 0:
            select_highest(st.session_state['gwas_df'], thresh, num)

    with tabs[2]:
        data_edit(st.session_state['gwas_df'])


@st.cache_data()
def lead_snp(df, finemapping=0):
    if finemapping > 0: 
        df['distance_lead'] = None
        df['lead'] = 'False'
        lead = df.loc[df['lrt_p'].idxmin(), 'pos']
        while df[df['distance_lead'].isna()].shape[0] !=0:
            try:
                df.loc[df['pos'] == lead, 'lead'] = 'True'
                index_lead = (df['pos'] >= lead - finemapping ) & (df['pos'] <= lead + finemapping) & (df['distance_lead'].isna())
                df.loc[index_lead, 'distance_lead'] = df.loc[index_lead, 'pos'] - lead
                min_pval = df[df['distance_lead'].isna()].idxmin()['lrt_p']
                lead = df.loc[min_pval, 'pos']
            except:
                break
    # else:
    #     df.loc[df['pos'] == df.idxmin()['lrt_p'], 'lead'] = 'True'
    #     df['distance_lead'] = df['pos'] - df.loc[df.idxmin()['lrt_p'], 'pos']  # A SNP before the lead will be negative
        st.session_state['pos_ls']['finemapping'] = df[df['lead'] == 'True']['pos'].tolist()
        df.loc[df['pos'].isin(st.session_state['pos_ls']['finemapping']), 'selected'] = 'finemapping'
        # st.session_state['gwas_df'] = df.drop('lead', axis=1)


@st.cache_data()
def select_highest(df, thresh, num):
    if num >0:
        highpos_ls = df[df['lrt_p'] > thresh].sort_values('lrt_p').head(num)['pos'].tolist()
        st.session_state['pos_ls']['highest'] = highpos_ls
        df.loc[df['pos'].isin(st.session_state['pos_ls']['highest']), 'selected'] = 'highest'
        # st.session_state['gwas_df'] = df


def data_edit(df):
    selected_df = df[df['selected'] != 'False']
    selected_df['Selection'] = True
    edited_df = st.expander('Refinement of selected positions').data_editor(
    selected_df[["chr", "pos", "lrt_p", "log_p",  "selected", "Selection"]].dropna(),
    column_config={
        "chr": "Chromosome",
        "pos": "Position",
        "ref": "Reference Allele",
        "alt": "Alternative Allele",
        "lrt_p": "P-value",
        "log_p": "-log(P-value)",
        "selected": "Method"
    },
    disabled=["pos", "chr", "lrt_p", "log_p", "ref", "alt"],
    hide_index=True,
    use_container_width=True
)

    st.session_state['selected_points'] =  edited_df[edited_df['Selection'] == True]['pos'].tolist()
    df.loc[(df['pos'].isin(st.session_state['pos_ls'].values()) == False), 'selected'] = 'False'
    # st.session_state['gwas_df'] = df

@st.cache_data()
def get_vcf(vcf_file, chr, pos_ls):
    # Get VCF 
    #TODO: Optimize the function
    try:
        from pysam import VariantFile
        reader = VariantFile(vcf_file)
        vcf_df = pd.DataFrame(index=pos_ls)
        for pos in pos_ls:
            for rec in reader.fetch(str(chr), int(pos)-1, int(pos)):
                info_pos = {'pos':rec.pos, 'ref':rec.ref, 'alt':rec.alts[0]}
                info_pos.update({rec.samples[sample].name:str(rec.samples[sample]['GT']) for _, sample in enumerate(rec.samples)})
                vcf_df = pd.concat([vcf_df, pd.DataFrame(info_pos, index=[pos])])
    except ModuleNotFoundError:
        import vcf
        reader = vcf.Reader(filename=vcf_file)
        for pos in pos_ls:
            for rec in reader.fetch(str(chr), int(pos)-1 , int(pos)):
                info_pos = {'pos':rec.POS, 'ref':rec.REF, 'alt':rec.ALT[0]}
                info_pos.update({rec.samples[sample].sample:str(rec.samples[sample]['GT']) for sample, _ in enumerate(rec.samples)})
                vcf_df = pd.concat([vcf_df, pd.DataFrame(info_pos, index=[pos])])
    try:
        vcf_df.drop_duplicates(inplace=True)
        vcf_df.dropna(inplace=True, how='all')
        vcf_df = vcf_df.melt(id_vars=['pos', 'ref', 'alt'], value_name='genotype')
        vcf_df['id'] = vcf_df['variable'].str.rsplit('lane1', n=1, expand=True)[1]
        return vcf_df
    except:
        return pd.DataFrame()


def collect_data():
    status = st.status("Processing...", state="running", expanded=True)
    status.write("Downloading VCF...")
    vcf_df = get_vcf(st.session_state['vcf_file'], st.session_state['chr'], st.session_state['selected_points'])
    if vcf_df.empty:
        status.update(label="No data found!", state="error", expanded=True)
    else:
        status.write("VCF retrieved!")
        status.update(label="VCF retrieved!", state="complete", expanded=False)
        merge = pd.merge(vcf_df, st.session_state['pheno_df'], left_on='variable', right_on='#IID')
        merge.drop(['variable', '#IID'], axis=1, inplace=True)
        pheno_cols = st.session_state['pheno_df'].drop('#IID', axis=1).columns
        id_cols = merge.drop(st.session_state['phenotypes'] , axis=1).columns
        merge = merge.melt(id_vars=id_cols, value_vars=pheno_cols)
        if merge.empty:
            status.update(label="Compilation failed", state="error", expanded=True)
        else: 
            status.update(label="Combination of data complete", state="complete", expanded=False)
            st.session_state['vcf_df'] = merge

def locus_correlation(df):
    # Transpose the DataFrame so that samples are in rows and loci in columns
    transposed_df = df.transpose()

    # Calculate pairwise squared correlation between loci
    correlation_df = transposed_df.corr().pow(2)

    return correlation_df

def make_boxplots(params=dict(x=None, facet_col=None, facet_row=None, color=None)):
            x_list = params.pop('x')
            replicates_n = st.session_state['vcf_df'].groupby(['genotype', 'variable', params['facet_col'], params['facet_row']])['id'].count().reset_index()
            replicates_n.rename(columns=dict(id='n'), inplace=True)
            # st.write(replicates_n)
            df = pd.merge(replicates_n, st.session_state['vcf_df'], on=['genotype', 'variable', params['facet_col'], params['facet_row']])
            # st.write(df, df.shape, st.session_state['vcf_df'].shape)
            df.sort_values(['pos', 'genotype', 'strain'], inplace=True)
            fig1 = px.box(df[df['variable'].isin(x_list)], x='variable', y='value', hover_data=df.drop('vcf_ID', axis=1).columns,facet_row_spacing=0.01, **params)
            if params['facet_row'] != None:
                n_row = df.loc[df['variable'].isin(x_list), params['facet_row']].nunique()
                fig1.update_layout(height=200 * int(n_row))
            st.plotly_chart(fig1, use_container_width=True)
            st.session_state['box_fig'] = fig1
            return fig1



def save_page(name_file='page_gwasbox_graph.html'):
    list_folders = name_file.split('/')
    # for n in range(len(list_folders)):
    if not os.path.exists(os.path.dirname(name_file)):
        os.makedirs(os.path.dirname(name_file))

    try:
        fig = plot_gwas_association(st.session_state['gwas_df'], color='selected', selection=False)       
        with open(name_file, 'w') as f:
            f.write(st.session_state['box_fig'].to_html(full_html=False, include_plotlyjs='cdn'))
            f.write(fig.to_html(full_html=False, include_plotlyjs=False))
    except Exception as err:
        print('could not save\n' + str(err))


# def save_png(name_file='page_gwasbox_graph.html'):
#     import pdfkit
#     pdfkit.html_to_pdf('http://localhost:8501/', 'out2.pdf')  
#     pdfkit.from_url(0, 'out.pdf', False)
#     # pdfkit.from_file('p_graph.html', 'out.pdf')  
#     # name_document = name_file
#     # try:
#     #     fig = plot_gwas_association(st.session_state['gwas_df'], color='selected', selection=False)
#     #     with open('test.png', 'w') as f:
#     #         f.write(st.session_state['box_fig'].write_image())
#     #         f.write(fig.write_image())
        
#     #     with open(name_document, 'w') as f:
#     #         f.write(st.session_state['box_fig'].to_html(full_html=False, include_plotlyjs='cdn'))
#     #         f.write(fig.to_html(full_html=False, include_plotlyjs='cdn'))
#     # except:
#     #     pass


def main(folder):
    st.title('GWAS association study')

    if 'pos_ls' not in st.session_state:
        st.session_state['pos_ls'] = dict(finemapping=[], highest=[], manual=[], loaded=[])
    if 'submitted' not in st.session_state:
        st.info('Submit the data you want to analyse on the left-side sidebar.')
        st.session_state['submitted'] = False

    # INITIAL SIDEBAR FORM
    with st.sidebar.form("file_selection", border = False):
        st.markdown('## Selection of the files: ')
        st.session_state['chr'] = st.slider(label='Chromosome', min_value=1, max_value=24, step=1)
        phenotype_ls = {f.split('/')[-2] for f in glob.glob(os.path.join(folder, '*/*.gwas*'))}
        st.session_state['gwas_phenotype']  = st.selectbox('GWAS', sorted(phenotype_ls), index=len(phenotype_ls)-1)
        st.session_state['vcf_file'] = st.selectbox('VCF', glob.glob(os.path.join(folder, '*vcf.gz*')), format_func=lambda x: x.split('/')[-1])
        select_pheno = st.selectbox('Phenotypes measures', glob.glob(os.path.join(folder, '*phenotype*')), format_func=lambda x: x.split('/')[-1])
        select_cov = st.selectbox('Covariates', glob.glob(os.path.join(folder, '*covariates*')), format_func=lambda x: x.split('/')[-1])
        
        submitted = st.form_submit_button("Submit", use_container_width=True)
        if submitted:
            with st.spinner('Initialization...'):
                initialization(folder, select_pheno, select_cov)
                # initialization(folder, select_chr, select_gwas, select_vcf, select_pheno, select_cov)
    
    if st.session_state['submitted'] :
        tab_function()
        selection_options()  # SNPs SELECTION TABs
        st.session_state['gwas_df']['type'] =  st.session_state['gwas_df']['selected'].apply(lambda x:str(type(x)))
        with st.empty():
            try:
                st.session_state['gwas_fig'] = plot_gwas_association(st.session_state['gwas_df'], selection=True, key='selection' + str(st.session_state['selection_n']))
            except:
                pass
        
        with st.form('action_form', border=False):
            action_list = ['Boxpot Genotype x Phenotype', 'Linkage desiquilibrium', 'Region']
            action_cols = st.columns([5,1])
            st.session_state['select_action'] = action_cols[0].selectbox('Select what you want to do :', action_list, label_visibility='collapsed')
            action_submit = action_cols[1].form_submit_button('Go')
            if action_submit:
                st.session_state['action_submit'] = True 
        
        # Action
        if not 'action_submit' in  st.session_state:
            st.session_state['action_submit']= False
        if not 'select_action' in  st.session_state:
            st.session_state['select_action']= None
        if not 'plot_submit' in st.session_state:
            st.session_state['plot_submit'] = False

        if st.session_state['action_submit'] == True:
            if st.session_state['select_action'] == action_list[0]:
                with st.form('Parameters for the boxplots: '):
                    st.header('Boxplot parameters: ')
                    params_col = st.columns(2)
                    x = st.multiselect('X-axis', st.session_state['phenotypes'], default=[col for col in st.session_state['phenotypes'] if 'loop' in col])

                    # hue = params_col[0].selectbox('Color parameter', st.session_state['pheno_df'].columns.tolist()  + [None])# index=st.session_state['pheno_df'].columns.tolist().index('genotype'))
                    facet_row = params_col[0].selectbox('Split row parameter',  st.session_state['covariates'] + [None])
                    facet_col = params_col[1].selectbox('Split column parameter', ['pos'] + [None])
                    dict_params = dict(x=x, facet_col=facet_col, facet_row=facet_row, color='genotype')
                    if st.checkbox('Display all points'): dict_params.update(dict(points='all'))
                    exp = st.expander('Filtering parameters')
                    # with exp:
                    #     cross_filter = st.multiselect('Cross selection: ', st.session_state['pheno_df']['strain'].unique(), st.session_state['pheno_df']['strain'].unique())
                    plot_submit = st.form_submit_button("Plot")
                    if plot_submit:
                        st.session_state['plot_submit'] =  True
                    #     pos_filter = st.multiselect('Position selection: ', st.session_state['selected_points'], st.session_state['selected_points'])


                if st.session_state['plot_submit'] == True:
                    collect_data()
                    st.session_state['box_fig'] = make_boxplots(params=dict_params)
                    with st.container(border=True):
                        save_cols = st.columns([2,1,1])
                        default_name = os.path.join(st.session_state['folder'], 'output', f"page_gwasbox_graph_{st.session_state['chr']}.html")
                        name_file = save_cols[0].text_input('Enter file path to save', value=default_name, label_visibility='collapsed')
                        saving = save_cols[1].button('**Save HTML page**', use_container_width=True, on_click=save_page, kwargs=dict(name_file=name_file))
                        # saving_png = save_cols[2].button('**Save PNG plots**', use_container_width=True, on_click=save_png, kwargs=dict(name_file=name_file))
                        if saving and os.path.exists(name_file):
                            st.success(f'The file {name_file} has been saved successfully!')
                        elif saving and (os.path.exists(name_file) ==False):
                            st.error(f'There has been an issue while saving file {name_file}')
            
            elif st.session_state['select_action'] == action_list[1]:
                vcf_df = get_vcf(st.session_state['vcf_file'], st.session_state['chr'], st.session_state['selected_points'])
                vcf_df = vcf_df[vcf_df['genotype'] != '(None, None)']
                st.write(vcf_df['genotype'].unique())
                vcf_df['gen_dum'] = vcf_df['genotype'].str[1].astype(int) + vcf_df['genotype'].str[4].astype(int)
                new_df = vcf_df.groupby(['pos', 'variable'])['gen_dum'].first().unstack()
                corr = locus_correlation(new_df)
                st.write(corr.head())
                st.write(new_df.head())
                # st.write(corr.min().min(), corr.max().max())
                labels = [str(c) for c in corr.columns]
                corr.columns = labels
                corr.index = [str(c) for c in corr.index]
                corr_fig = px.imshow(corr, zmin=0, zmax=1, x=corr.columns, y=corr.columns)
                st.plotly_chart(corr_fig, use_container_width=True)
            
   
     
if __name__ == '__main__':
    assert os.path.isdir(sys.argv[1]), 'Folder cannot be found'
    main(sys.argv[1])


# TODO: IGV embedding ?
# TODO: calculation of distance from lead SNP AND linkqge desiquilirium
# TODO: info about - gene -LOF - ...
# TODO: Option to save in PNG not HTML
# TODO: Include fitering ?