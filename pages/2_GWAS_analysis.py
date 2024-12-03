import sys
import os
import numpy as np
import pandas as pd
import socket
import subprocess
from io import BytesIO
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
from streamlit_utils import confirmation_box, filter_dataframe
pio.templates.default = "simple_white"
st.set_page_config(
    page_title="MIKK Toolbox",
    page_icon="🐟",
    layout="wide",
)

#IDEA: Save specific loci in wanted file + Load the selection you want
#TODO: remove general functions eg. get_vcf 


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
    elif mode == 's':
        pd.Series(st.session_state['selected_points']).drop_duplicates().to_csv(st.session_state['list_pos_filename'], 
            sep='\t', mode='w', index=False, header=False)


def plot_gwas_association(df, color='selected', selection=False, key='', cont=st):
    df = st.session_state['gwas_df'].copy()
    n=1
    color_dict={"False":'#949191', "highest":"red", "finemapping":"orange", "manual":"blue", "loaded":'purple'}
    df['selected'] = df['selected'].astype(str)
    fig = px.scatter(data_frame=df, x='pos', y='log_p', color=color, color_discrete_map=color_dict, opacity=0.5, title='GWAS ' + st.session_state['gwas_phenotype'])
    with cont:
        while selection:
            selected_points = plotly_events(fig, click_event=True, hover_event=False, select_event=True, key=key + str(n))
            if selected_points:
                selected_points = pd.DataFrame(selected_points)
                st.session_state['pos_ls']['manual'] = selected_points['x'].tolist()
                selected_points = None
                n += 1                   
                fig = px.scatter(data_frame=df, x='pos', y='log_p', color=color, color_discrete_map=color_dict, opacity=0.5, title='GWAS ' + st.session_state['gwas_phenotype'])
        st.session_state['gwas_fig'] = fig
    return fig


def tab_function():
    st.markdown('## Select the SNPs to analyze')
    tabs = st.tabs([s.center(20,"\u2001") for s in ['Finemapping', 'Highest', 'Load list']])
    with tabs[0]:
        # Finemapping
        fmap_form = st.form('finemapping_form', border=False)
        fmap_cols = fmap_form.columns([1,3,2])
        fmap_cols[0].caption('Finemapping limit (bp)')
        finemapping = fmap_cols[1].slider('finemapping', min_value=0, max_value=10000000, value=0, step=100000, label_visibility='collapsed')
        fmap_submit = fmap_cols[2].form_submit_button('Run Finemapping', use_container_width=True)
        if fmap_submit:
             st.session_state['pos_ls']['finemapping'] = lead_snp(st.session_state['gwas_df'], finemapping=finemapping)  # Cached   

    with tabs[1]:
        high_form = st.form('select_high_form', border=False)
        high_cols = high_form.columns([1,2,1,2,1])
        high_cols[0].caption('Number of SNPs to select')
        num = high_cols[1].slider('Number of SNPs', min_value=0, max_value=10, step=1, value=0, label_visibility='collapsed')
        high_cols[2].caption(' Threshold')
        thresh = high_cols[3].number_input('Threshold', min_value=0.0, max_value=1.0, step=0.05, label_visibility='collapsed')
        high_submit = high_cols[4].form_submit_button('Select highest', use_container_width=True)
        if high_submit:
            st.session_state['pos_ls']['highest'] = select_highest(st.session_state['gwas_df'], thresh, num)

    with tabs[2]:
        load_form = st.form('load_form', border=False)
        load_cols = load_form.columns([2,1])
        load_cols[0].file_uploader("Choose a txt file", accept_multiple_files=False, key='load_file', label_visibility='collapsed')
        load_submit = load_cols[1].form_submit_button('Load SNPs from file', use_container_width=True)
        if (st.session_state['load_file'] is not None) & load_submit:
            st.session_state['pos_ls']['loaded'] = pd.read_table(st.session_state['load_file'], header=None)[0].tolist()

    option_cols = st.columns(4)
    if option_cols[0].button('Reset all', use_container_width=True):
        confirmation_box(update_selection, click_kwargs=dict(mode='rm'), container=st)
    if option_cols[1].button('Append selected', use_container_width=True):
        update_selection(mode='a')
    if option_cols[2].button('Replace selected', use_container_width=True):
        update_selection(mode='w')
    if option_cols[3].button('Save list of selected positions', use_container_width=True):
        with st.form('Save list', border=False):
            save_param_cols = st.columns([3,1])
            save_param_cols[0].text_input('Saving path', value=os.path.join(st.session_state['output_dir'], 'list_snp.txt'), label_visibility='collapsed', key='list_pos_filename')
            save_param_cols[1].form_submit_button('Save', use_container_width=True, on_click=update_selection, kwargs=dict(mode='s'))
    

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
    return df[df['lead'] == 'True']['pos'].tolist()



@st.cache_data()
def select_highest(df, thresh, num):
    if num >0:
        highpos_ls = df[df['lrt_p'] > thresh].sort_values('lrt_p').head(num)['pos'].tolist()
        return highpos_ls


def data_edit(df, cont=st):
    if df.shape[0] >10:
        st.session_state['filtered_gwas'] = filter_dataframe(df, key='data_edit', cont=cont)
    else:
        df['Selection'] = True
        st.session_state['filtered_gwas']  = st.data_editor(
            df[["chr", "pos", "lrt_p", "log_p",  "selected", "Selection"]].dropna(),
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
    st.session_state['filtered_gwas'] = st.session_state['filtered_gwas'][st.session_state['filtered_gwas']['Selection'] == True]
    st.session_state['pos_ls']['manual'] =   st.session_state['filtered_gwas']['pos'].tolist()


@st.cache_data()
def get_vcf(vcf_file, chr, pos_ls):
    # Get VCF 
    #TODO: Optimize the function or put in utils
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
        merge = pd.merge(merge, st.session_state['gwas_df'][st.session_state['gwas_df']['pos'].isin(st.session_state['selected_points'])][['pos', 'ref', 'alt', 'log_p', 'lrt_p', 'lrt_df', 'lrt_chisq']], on=['pos', 'ref', 'alt'])
        merge.drop(['variable', '#IID'], axis=1, inplace=True)
        pheno_cols = st.session_state['pheno_df'].drop('#IID', axis=1).columns
        id_cols = merge.drop(st.session_state['phenotypes'] , axis=1).columns
        melt = merge.melt(id_vars=id_cols, value_vars=pheno_cols)
        if merge.empty:
            status.update(label="Compilation failed", state="error", expanded=True)
        else: 
            status.update(label="Combination of data complete", state="complete", expanded=False)
            melt['strain'] = melt['strain'].astype('category')
            st.session_state['vcf_df'] = melt


def locus_correlation(df):
    # Transpose the DataFrame so that samples are in rows and loci in columns
    transposed_df = df.transpose()

    # Calculate pairwise squared correlation between loci
    correlation_df = transposed_df.corr().pow(2)

    return correlation_df


def make_boxplots(vcf_df, params=dict(x=None, facet_col=None, facet_row=None, color=None)):
            x_list = params.pop('x')
            if (params['facet_row'] != None) & (params['polymorphic'] ==  True): 
                row_gen = vcf_df[vcf_df['genotype'] != '(None, None)'].groupby([params['facet_row']])['genotype'].nunique().reset_index()
                row_pol = row_gen[row_gen['genotype'] >= 3][params['facet_row']]
                vcf_df = vcf_df[vcf_df[params['facet_row']].isin(row_pol)]
            params.pop('polymorphic')
            replicates_n = vcf_df.groupby(['genotype', 'variable', params['facet_col'], params['facet_row']])['id'].count().reset_index()
            replicates_n.rename(columns=dict(id='n'), inplace=True)
            # st.write(replicates_n)
            df = pd.merge(replicates_n, vcf_df, on=['genotype', 'variable', params['facet_col'], params['facet_row']])
            # st.write(df, df.shape, st.session_state['vcf_df'].shape)
            df.sort_values(['pos', 'genotype', 'strain'], inplace=True)
            fig1 = px.box(df[df['variable'].isin(x_list)], x='variable', y='value', hover_data=df.drop('vcf_ID', axis=1).columns,facet_row_spacing=0.01, **params)
            if params['facet_row'] != None:
                n_row = df.loc[df['variable'].isin(x_list), params['facet_row']].nunique()
                fig1.update_layout(height=200 * max(1, int(n_row)))
            st.plotly_chart(fig1, use_container_width=True)
            st.session_state['box_fig'] = fig1
            return fig1


def save_page(name_file='page_gwasbox_graph.html'):
    with st.form('Html page saving', border=True):
        save_cols = st.columns([2,1,1])
        default_name = os.path.join(st.session_state['flexlmm_dir'], 'output', f"page_gwasbox_graph_{st.session_state['chr']}.html")
        name_file = save_cols[0].text_input('Enter file path to save', value=default_name, label_visibility='collapsed')
        saving = save_cols[1].form_submit_button('**Save HTML page**', use_container_width=True)
    if saving :
        if not os.path.exists(os.path.dirname(name_file)):
            os.mkdir(os.path.dirname(name_file))

        fig = plot_gwas_association(st.session_state['gwas_df'], color='selected', key='saving', cont=st.container())
               
        try:
            with open(name_file, 'w') as f:
                from markdown import markdown
                f.write(markdown(''' <style>
                                    body {background-color: white;}
                                    h1   {color: #1c205a; font-family: Arial, Helvetica, sans-serif;}
                                    h2   {color: #414452; font-family: Arial, Helvetica, sans-serif;}
                                    h3   {color: #414452; font-family: Arial, Helvetica, sans-serif;}
                                    p    {color: black; font-family: Arial, Helvetica, sans-serif;}
                                    </style>
                                 ''')) 
                styles = [
                dict(selector="tr:hover",
                            props=[("background", "#f4f4f4")]),
                dict(selector="th", props=[("color", "#fff"),
                                        ("border", "1px solid #eee"),
                                        ("padding", "12px 30px"),
                                        ("border-collapse", "collapse"),
                                        ("background", "#343838"),
                                        ("text-transform", "uppercase"),
                                        ("font-size", "14px"),
                                        ("font-family" , 'Arial'),
                                        ]),
                dict(selector="td", props=[("color", "#999"),
                                        ("border", "1px solid #eee"),
                                        ("padding", "12px 30px"),
                                        ("border-collapse", "collapse"),
                                        ("font-size", "12px"),
                                        ("font-family" , 'Arial'),
                                        ]),
                dict(selector="table", props=[
                                                ("font-family" , 'Arial'),
                                                ("margin" , "25px auto"),
                                                ("border-collapse" , "collapse"),
                                                ("border" , "1px solid #eee"),
                                                ("border-bottom" , "2px solid #343838"),                                    
                                                ]),
                dict(selector="caption", props=[("caption-side", "bottom")])]
                f.write(markdown('# Report analysis'))
                f.write(markdown('### GWAS'))

                gwas_out = st.session_state['gwas_df'][st.session_state['gwas_df']['pos'].isin(st.session_state['selected_points'])].drop('beta', axis=1)
                f.write(gwas_out.reset_index(drop=True).style.set_table_styles(styles).background_gradient().to_html(index=False, justify='center', border=1))
                f.write(fig.to_html(full_html=False, include_plotlyjs='cdn'))
                if 'box_fig' in st.session_state:
                    f.write(markdown('### Association Boxplots'))
                    f.write(st.session_state['box_fig'].to_html(full_html=False, include_plotlyjs='cdn'))
                if 'LD_fig' in st.session_state:
                    f.write(markdown('### Linkage desiquilibrium'))
                    f.write(st.session_state['LD_fig'].to_html(full_html=False, include_plotlyjs='cdn'))
            if (os.path.exists(name_file) == True):
                st.success(f'The file {name_file} has been saved successfully!')
            elif (os.path.exists(name_file) == False):
                st.error(f'There has been an issue while saving file {name_file}')
        except Exception as err:
            st.error('could not save' + str(err))


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

# @st.cache_data()
def initialization():

    # Initiation of the state session
    st.session_state['submit_gwas'] = True

    # Reading tables: 
    gwas_file = os.path.join(st.session_state['flexlmm_dir'], 'gwas', st.session_state['gwas_phenotype'] , f'input_{st.session_state["chr"]}_{st.session_state["gwas_phenotype"] }.gwas.tsv.gz')
    # st.write(gwas_file)
    df_gwas = pd.read_table(gwas_file)
    df_gwas['log_p'] = -np.log10(df_gwas['lrt_p'])
    df_gwas['selected'] = 'False'
    st.session_state['gwas_df'] = df_gwas

    pheno_df = pd.read_table(st.session_state['phenotype_file'])
    st.session_state['phenotypes'] = pheno_df.drop('#IID', axis=1).columns.tolist()
    cov_df = pd.read_table(st.session_state['covariate_file'])
    df = pd.merge(pheno_df, cov_df, on='#IID')
    st.session_state['covariates'] = cov_df.drop('#IID', axis=1).columns.tolist()
    df['vcf_ID'] = df['#IID'] + '_' + df['#IID']
    st.session_state['pheno_df'] = df

    # Initialisation
    st.session_state['pos_ls'] = dict(finemapping=[], highest=[], manual=[], loaded=[])  
    st.session_state['selected_points']  = []    
    st.session_state['selection_n'] = 0


def main():
    st.title('GWAS association study')
    if 'submit_gwas' not in st.session_state:
        st.info('Submit the data you want to analyse on the left-side sidebar.')
        st.session_state['submit_gwas'] = False

    if len(set(['vcf_file', 'flexlmm_dir', 'phenotype_file', 'covariate_file', 'output_dir']) - set(st.session_state.keys())) > 0: 
        with st.expander('Files parameters').form('General params', border=False):
            st.header('Enter paths manually')
            params = {}
            for var in ['output_dir', 'vcf_file', 'flexlmm_dir', 'phenotype_file','covariate_file']:
                if (var not in st.session_state) or (st.session_state[var] == ''):
                    params[var] = st.text_input(var.replace('_', ' ').replace('dir', 'directory').title(), value='')
            submit_params = st.form_submit_button('Submit', use_container_width=True)     
            if submit_params:
                st.session_state.update(params)

    if 'pos_ls' not in st.session_state:
        st.session_state['pos_ls'] = dict(finemapping=[], highest=[], manual=[], loaded=[])
   
    if 'flexlmm_dir' in st.session_state:
        with st.sidebar.form("params_selection", border=False):
            st.markdown('## Selection of the parameters: ')
            st.session_state['chr'] = st.slider(label='Chromosome', min_value=1, max_value=24, step=1)
            phenotype_ls = {f.split('/')[-2] for f in glob.glob(os.path.join(st.session_state['flexlmm_dir'], 'gwas/*/*.gwas*'))}
            st.session_state['gwas_phenotype']  = st.selectbox('GWAS', sorted(phenotype_ls), index=len(phenotype_ls)-1)
            submitted = st.form_submit_button("Submit", use_container_width=True)
            if submitted:
                with st.spinner('Initialization...'):
                    initialization()

        # with st.sidebar.form("saving_params", border=True):
        #     col_sav = st.columns([2,1])
        #     new_path = col_sav[0].text_input('Path for list of position to be saved', value=os.path.join(os.getcwd(), f'list_pos.temp'), label_visibility='collapsed')
        #     change_save_path = col_sav[1].form_submit_button("Change path", use_container_width=True)
        #     if change_save_path:
        #         st.session_state['list_pos_filename'] = new_path

    
    if st.session_state['submit_gwas'] :
        try:
            st.image(os.path.join(st.session_state['flexlmm_dir'], 'plots', 'manhattan', st.session_state['gwas_phenotype'] + '.png' ), use_column_width=True)
        except:
            st.warning('Could not find full GWAS plot')
        
        st.sidebar.write('Positions', st.session_state['pos_ls'])
        st.sidebar.write('Selected points', st.session_state['selected_points'])
        
        tab_function()
        placeholder_gwas = st.empty()
        try:
            st.session_state['gwas_fig'] = plot_gwas_association(st.session_state['gwas_df'], selection=True, key='selection' + str(st.session_state['selection_n']), cont=placeholder_gwas)
        except:
            pass
        
        if 'vcf_df' not in st.session_state:
            st.session_state['vcf_df'] = pd.DataFrame()
        st.button('Download genotype data', type="primary", on_click=collect_data, use_container_width=True)

        if 'filtered_vcf' not in st.session_state:
            st.session_state['filtered_vcf'] = pd.DataFrame()

        if not st.session_state['vcf_df'].empty:
            st.session_state['boxplot_submit'] = False
            plot_tabs = st.tabs(['Boxpot Genotype x Phenotype', 'Linkage desiquilibrium', 'Get the datas'])
            with plot_tabs[0]:
                cont_box = st.container(border=True)
                cont_box.header('Boxplot parameters: ')
                filter_exp = cont_box.expander('Further Filtering')
                filter_dataframe(st.session_state['vcf_df'], key='boxplot_df', cont=filter_exp, session_state_var='filtered_vcf')
                with cont_box.form('Parameters for the boxplots: ', border=False):
                    params_col = st.columns(2)
                    x = st.multiselect('X-axis', st.session_state['phenotypes'], default=st.session_state['gwas_phenotype'])
                    # hue = params_col[0].selectbox('Color parameter', st.session_state['pheno_df'].columns.tolist()  + [None])# index=st.session_state['pheno_df'].columns.tolist().index('genotype'))
                    facet_col = params_col[0].selectbox('Split column parameter', ['pos'] + [None])
                    facet_row = params_col[1].selectbox('Split row parameter',  st.session_state['covariates'] + [None])
                    dict_params = dict(x=x, facet_col=facet_col, facet_row=facet_row, color='genotype')
                    if params_col[1].checkbox('Display all points'): dict_params.update(dict(points='all'))
                    if params_col[0].checkbox('Display only polymorphic rows (work with facet_row)'): 
                        dict_params.update(dict(polymorphic=True))
                    else: 
                        dict_params.update(dict(polymorphic=False))
                    boxplot_submit = st.form_submit_button("Plot")
                    if boxplot_submit:
                        st.session_state['boxplot_submit'] =  True

                if st.session_state['boxplot_submit'] == True:
                    if not st.session_state['filtered_vcf'].empty:
                        st.session_state['box_fig'] = make_boxplots(vcf_df=st.session_state['filtered_vcf'], params=dict_params)
                    else:
                        st.session_state['box_fig'] = make_boxplots(vcf_df=st.session_state['vcf_df'], params=dict_params)

            
            with plot_tabs[1]:
                vcf_df = st.session_state['vcf_df'][st.session_state['vcf_df']['genotype'] != '(None, None)'].copy()
                vcf_df = vcf_df[vcf_df['variable'] == st.session_state['gwas_phenotype']]
                vcf_df['gen_dum'] = vcf_df['genotype'].str[1].astype(int) + vcf_df['genotype'].str[4].astype(int)
                new_df = vcf_df.groupby(['pos', 'id'])['gen_dum'].first().unstack()
                corr = locus_correlation(new_df)
                # st.write(corr)
                labels = [str(c) for c in corr.columns]
                corr.columns = labels
                corr.index = [str(c) for c in corr.index]
                corr_fig = px.imshow(corr, zmin=0, zmax=1, x=corr.columns, y=corr.columns)
                st.session_state['LD_fig']  = corr_fig
                st.plotly_chart(corr_fig, use_container_width=True)
            
            with plot_tabs[2]:
                vcf_df = st.session_state['vcf_df'][st.session_state['vcf_df']['genotype'] != '(None, None)'].copy()
                vcf_df = vcf_df[vcf_df['variable'] == st.session_state['gwas_phenotype']]
                vcf_df['gen_dum'] = vcf_df['genotype'].str[1].astype(int) + vcf_df['genotype'].str[4].astype(int)
                new_df = vcf_df.groupby(['pos', 'id'])['gen_dum'].first().unstack()
                st.write('### GWAS data:')
                st.dataframe(vcf_df.filter(regex='pos|ref|alt|(^l)|(variable)').drop_duplicates(), hide_index=True, use_container_width=True)
                st.write('### Genotype data:')
                st.dataframe(new_df.T, use_container_width=True)


            
            save_page()

            
     
if __name__ == '__main__':
    main()


# TODO: IGV embedding ?
# TODO: Option to save in PNG not HTML
# TODO: Include fitering ?