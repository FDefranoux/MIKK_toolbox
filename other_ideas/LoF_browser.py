import pandas as pd
import sqlite3
import streamlit as st
import seaborn as sns
import matplotlib.pyplot as plt
import plotly.express as px
import sys
import socket
host = socket.gethostname()
if 'Fanny' in host:
    PATH_UTILS = '/home/fanny/Work/EBI/Utils'
    ABS_PATH = '/home/fanny/Work/EBI/covid_nanopore'
else:
    PATH_UTILS = '/nfs/research/birney/users/fanny/Utils'
    ABS_PATH = '/nfs/research/birney/users/fanny/covid_nanopore'
sys.path.insert(0, PATH_UTILS)
from utils import *
# from utils_plots import *
from goatools import obo_parser
from PIL import Image
import plotly.graph_objects as graph_obj

def static_image_plotly(picture):
    fig = graph_obj.Figure()

    # Constants
    img_width, img_height = picture.size
    # img_width = 1600
    # img_height = 900
    scale_factor = 0.5

    # Add invisible scatter trace.
    # This trace is added to help the autoresize logic work.
    # fig.add_trace(
    #     graph_obj.Scatter(
    #         x=[0, img_width * scale_factor],
    #         y=[0, img_height * scale_factor],
    #         mode="markers",
    #         marker_opacity=0
    #     )
    # )

    # Configure axes
    fig.update_xaxes(
        visible=False,
        range=[0, img_width * scale_factor]
    )

    fig.update_yaxes(
        visible=False,
        range=[0, img_height * scale_factor],
        # the scaleanchor attribute ensures that the aspect ratio stays constant
        scaleanchor="x"
    )

    # Add image
    fig.add_layout_image(
        dict(
            x=0,
            sizex=img_width * scale_factor,
            y=img_height * scale_factor,
            sizey=img_height * scale_factor,
            xref="x",
            yref="y",
            opacity=1.0,
            # layer="below",
            # sizing="stretch",
            source=picture)
    )

    # Configure other layout
    # fig.update_layout(
    #     width=img_width * scale_factor,
    #     height=img_height * scale_factor,
    #     margin={"l": 0, "r": 0, "t": 0, "b": 0},
    # )
    return fig

@st.cache
def open_dataframes(genome_file, go_file):
    genome_file='/home/fanny/Work/EBI/Indigene/high_results.txt'
    go_file='/home/fanny/Work/EBI/Indigene/datafiles/Medaka_Gene_Database3.csv'
    st.session_state['db_name']  = 'datafiles/medaka_genome_part.db'
    st.session_state['df'] = pd.read_csv(genome_file)
    st.session_state['go_df'] = pd.read_csv(go_file)
    # df = pd.read_csv(genome_file)
    # go_df = pd.read_csv(go_file)
    go_obo_file = '/home/fanny/Work/EBI/Indigene/datafiles/GO_obo/go-basic.obo'
    st.session_state['go_dict'] = obo_parser.GODag(go_obo_file)


# blou = go_df.loc[go_df['name_1006'].astype(str).str.contains('temperature', regex=False)].copy()
# blou.nunique()
# blou['goslim_goa_description'].unique()
# blou.loc[blou['goslim_goa_description'].astype(str).str.contains('temperature', regex=False)]
# blou['gp_name'] = blou['go_id'].apply(lambda x: get_go_info(x, go, var='name'))
# blou[blou['goslim_goa_description'] == blou['name_1006']]

def form_for_filtering(df, go_df):
    st.sidebar.title("Filtering options")
    dict_path = go_df[['go_id', 'name_1006']].drop_duplicates().sort_values('name_1006').set_index(
        'name_1006').to_dict()['go_id']
    if st.sidebar.selectbox('', ['Plot pathways', 'Filter dataset']) == 'Plot pathways':
        pathway_plots(go_df, dict_path)
    else:
        with st.sidebar.form("Filtering options"):
            df['variant'] = df['variant'].str.split(',', expand=True)[0]
            df['line'] = df['line'].str.replace('-', '_')

            line_ls = st.multiselect('Line selection:', df['line'].unique(), default=[])
            impact_ls = st.multiselect('Impact selection:', df['impact'].unique(), default=[])
            variant_ls = st.multiselect('Variants selection:', df['variant'].unique(), default=[])
            list_path_new = st.multiselect('Pathway selection', dict_path.keys())
            list_id = [dict_path[path] for path in list_path_new]
            if st.checkbox('Include all children'):
                for id in list_id:
                    list_id += st.session_state['go_dict'][id].get_all_children()
            # Submit button:
            submitted = st.form_submit_button("Filter")
            if submitted:
                dataset_filtering('df',{ 'variant': variant_ls, 'line': line_ls, 'impact': impact_ls})
                dataset_filtering('go_df', dict(go_id=list_id))
                gene_ls = list(st.session_state['filt_df']['gene'].unique())
                # st.write(list(st.session_state['filt_go_df']['go_id'].unique()) == list_id)
                blou = pd.read_sql(f"SELECT ensembl_gene_id AS gene, go_id, uniprot_gn_symbol FROM medaka_genes WHERE ensembl_gene_id IN ({str(gene_ls)[1:-1]}) AND go_id IN ({str(list_id)[1:-1]})", con=sqlite3.connect(st.session_state['db_name']))
                merge_filtered_df(blou)

            else:
                # st.session_state['filt_df'] = st.session_state['df']
                # st.session_state['filt_go_df'] = st.session_state['go_df']
                st.session_state['final_df'] = pd.DataFrame()
        # st.write(list(st.session_state['filt_go_df']['go_id'].unique()) == list_id)
        # st.write(blou['gene'].unique())
        # st.write(blou['go_id'].unique())
        # st.dataframe(st.session_state['final_df'])
        go_ids = list(st.session_state['filt_go_df']['go_id'].unique())
        go_ids = list(st.session_state['filt_go_df']['go_id'].unique())
        gene_df = pd.read_sql(f"SELECT ensembl_gene_id AS gene, go_id, uniprot_gn_symbol FROM medaka_genes WHERE go_id IN ({str(go_ids)[1:-1]})", con=sqlite3.connect(st.session_state['db_name']))
        # gene_df.drop_duplicates().to_csv('Temp_heartrate_gene_go_df.csv', index=False)
        # st.dataframe(st.session_state['filt_go_df'][['go_id', 'name_1006']].drop_duplicates())


def merge_filtered_df(df):
    st.session_state['final_df'] = pd.merge(st.session_state['filt_df'], df, on='gene')
    st.session_state['final_df'] = pd.merge(st.session_state['final_df'], st.session_state['filt_go_df'], on='go_id')
    st.write('Filtered')
# def find_pathway_interest():
#     with st.sidebar.form("Pathway selection Filtering pathways"):
#         text_filt = st.text_input('Filtering pattern')
#         list_path = go_df.loc[go_df['definition_1006'].astype(str).str.contains(text_filt, regex=st.checkbox('Regex'))]['name_1006'].unique()
#         list_path_new = st.select('Pathway selection', list_path)
#

def pathway_plots(go_df, dict_path):
    st.sidebar.write(" ## Filtering pathways")
    # text_filt = st.sidebar.text_input('Filtering pattern')
    path_df = go_df[['go_id', 'name_1006']].drop_duplicates().sort_values('name_1006')
    # dict_path = go_df.loc[go_df['name_1006'].astype(str).str.contains(
    #     text_filt, regex=st.sidebar.checkbox('regex'))][['go_id', 'name_1006']].drop_duplicates().sort_values('name_1006').set_index(
    #     'name_1006').to_dict()['go_id']
    with st.sidebar.form("Plotting pathways"):
        paths = st.multiselect('Plot these pathway', path_df['name_1006'].sort_values().unique())
        cols = st.columns(2)
        cols[1] = st.checkbox('Include parents')
        cols[0] = st.checkbox('Include children')
        if st.form_submit_button("Plot"):
            dict_path = path_df.set_index('name_1006').to_dict()['go_id']
            st.session_state['go_dict'].draw_lineage([st.session_state['go_dict'][dict_path[name]] for name in paths],
            draw_parents=cols[1], draw_children=cols[0], dpi=200, output='test_go.png')
    st.title('Graphical view of pathway')
    picture = Image.open('test_go.png')
    fig = static_image_plotly(picture)
    st.plotly_chart(fig, use_container_width=True)


def pathway_filtering(go_df):
    st.sidebar.title("Filtering pathways")
    text_filt = st.sidebar.text_input('Filtering pattern')
    dict_path = go_df.loc[go_df['name_1006'].astype(str).str.contains(
        text_filt, regex=st.sidebar.checkbox('regex'))][['go_id', 'name_1006']].drop_duplicates().sort_values('name_1006').set_index(
        'name_1006').to_dict()['go_id']
    if st.sidebar.checkbox('Graphical view of the pathways'):
        pathway_plots(go_df, dict_path)

    with st.sidebar.form("Filtering pathways"):
        list_path_new = st.multiselect('Pathway selection', dict_path.keys())
        list_id = [dict_path[path] for path in list_path_new]
        st.write(list_id)
        if st.checkbox('Include all children'):
            list_id += [st.session_state['go_dict'][id].get_all_children() for id in list_id]#
        submit= st.form_submit_button("Filter pathways")
        if submit:
            dataset_filtering('go_df', dict(go_id=list_id))
        else:
            st.session_state['filt_go_df'] = st.session_state['go_df']


def dataset_filtering(state_session_key, dict_select):
    df = st.session_state[state_session_key].copy()
    filter_mask = df.iloc[:, 0] == df.iloc[:, 0]
    for key, val in dict_select.items():
        if val:
            filter_mask &= st.session_state[state_session_key][key].isin(list(val))
    # if variant_ls:
    #     filter_mask &= df['variant'].isin(variant_ls)
    # if line_ls:
    #     filter_mask &= df['line'].isin(line_ls)
    st.session_state[f'filt_{state_session_key}'] = st.session_state[state_session_key].loc[filter_mask]


@st.cache(suppress_st_warning=True, allow_output_mutation=True)
def plots_count_tot(df):
    # TODO: Add the total counts when you will have them
    st.warning('Fixed dataframe, no filtering possible at the moment')
    count = pd.read_csv('effect_prediction_analysis/impact_counts_results.txt')
    count.rename(columns={'COUNT(DISTINCT snp)': 'count'}, inplace=True)
    count = count.pivot(index='line', columns = 'impact')['count']
    count.index = count.index.str.replace('-', '_')

    fig_count = px.bar(count, barmode='stack')
    tot_count = px.pie(count.mean().reset_index(), names='impact', values=0)

    return fig_count, tot_count


@st.cache(suppress_st_warning=True, allow_output_mutation=True)
def plots_count(df):
    st.warning('The filtered dataframe is currently only on HIGH impact variants')
    ### PLOTS ####
    var = 'variant'   #st.sidebar.radio('variants', ['all', 'first', 'last'])
    if var == 'all':
        df_variant = st.session_state['filt_df'].explode(column='variant_ls').groupby(['line', 'variant']).size().unstack().sort_values('frameshift_variant')
        var = 'variant_ls'
    else:
        df_variant = st.session_state['filt_df'].groupby(['line', var]).size().unstack().sort_values('frameshift_variant')

    fig_line = px.bar(df_variant, barmode='stack')
    fig_tot = px.pie(df_variant.mean().reset_index(), names=var, values=0)
    return fig_line,  fig_tot


def dict_fig_plots(fig_ls, title):
    fig_dict = {f'{title}_perline': fig_ls[0], f'{title}_tot': fig_ls[1]}
    for title, plot in fig_dict.items():
        try:
            st.plotly_chart(plot)
            if st.button(f'Save {title}'):
                save_plot_report(title.replace(' ', '_'), plot)
        except Exception as err:
            st.error(f'Maybe wrong file type:\n{err}')


def Merge_Gene_Info(df, db):
    # Merge with the medaka_gene table
    gene_ls = str(list(df['gene'].unique()))[1:-1]
    genes = pd.read_sql(f"SELECT DISTINCT ensembl_gene_id, description, uniprot_gn_symbol, go_id FROM medaka_genes WHERE ensembl_gene_id IN ({gene_ls})", con=db)
    genes.rename(columns={'ensembl_gene_id':'gene'}, inplace=True)
    df = pd.merge(df, genes, on='gene')
    return df


def Merge_GO_Info(df):
    from goatools import obo_parser
    go_obo_file = '/home/fanny/Work/EBI/Indigene/datafiles/GO_obo/go-basic.obo'
    go = obo_parser.GODag(go_obo_file)
    gene_ls = str(list(df['gene'].unique()))[1:-1]
    go_df = pd.read_sql(f"SELECT DISTINCT go_id FROM medaka_genes WHERE ensembl_gene_id IN ({gene_ls})", con=db)
    for var in ['name', 'depth', 'level', 'is_obsolete']:
        go_df[var] = go_df['go_id'].apply(lambda x: get_go_info(x, go, var=var))
    for var in ['_parents', 'get_all_children']:
        go_df['Num' + var] = go_df['go_id'].apply(lambda x: get_go_info(x, go, var=var)).str.len()
    go_df.dropna(thresh=4, inplace=True)
    go_df[go_df['Numget_all_children']< 15].sort_values('Numget_all_children')
    len(go['GO:0002027']._parents)
    go_df[go_df['go_id'].isin(go['GO:0008016'].get_all_children())]

    go_df[go_df['name'].str.contains('ethanol')]
    dir(go['GO:0070815'])
    getattr(go['GO:0070815'], '_parents')
    go['GO:0070815'].get_all_children()
    # 24 obsolete GOID
    pd.DataFrame(dict(go))
    pd.DataFrame(go)
    go_df
    help(go.draw_lineage)
    go.draw_lineage([go['GO:0009266'], go['GO:0002027'], go['GO:0045471']], draw_parents=True, output='test_go.png')
    help(go.draw_lineage)
    from IPython.display import Image
    Image('GO_lineage.png')


def get_go_info(x, go, var='name'):
    if 'get_all' in var:
        try:
            return getattr(go.get(x), var)()
        except:
            return None
    else:
        try:
            return getattr(go.get(x), var)
        except:
            return None


def cluster_analysis(df):
    # PATHWAY
    path = st.text_input('Enter pathway')
    cross = st.checkbox('Cross selection')
    if st.sidebar.button('Run'):
        if (not lines) :
            lines = df['line'].unique()
        if cross:
            assert len(lines) == 2, 'Cannot perform cross variant analysis if more than 2 lines'
            df = df.loc[df['line'].isin(lines)]

        blou = df[(df['go_name'].astype(str).str.contains(path))
            & (df['line'].isin(lines))].groupby(['line',
            'go_name']).size().unstack().dropna(how='all').fillna(0)

        dict_colors = {'KM': 'viridis', 'HC': 'Set2'}
        blou['KM'] = Component_analysis(blou, type='blou').KMeans_clustering(n_clust=5)
        blou['HC'] = Component_analysis(blou, type='blou').hierarchical_clustering()
        clustermap = px.imshow(blou.sort_values('HC').T, labels={'y': False})
        st.plotly_chart(clustermap)
        clust = myclustermap(blou.sort_values('HC'), dict_color=dict_colors,
            row_cluster=True, col_cluster=True, orientation='horizontal')
        st.pyplot(clust)

def main():
    st.set_page_config(
        page_title="Variant MIKK panel analyses",
        page_icon="/home/fanny/Pictures/bioinfo.gif",
        layout="wide", initial_sidebar_state="expanded",
        menu_items={'About': "# MIKK\n ## Variant analysis"})


    # read dataframes
    open_dataframes(genome_file='/home/fanny/Work/EBI/Indigene/high_results.txt',
                    go_file='/home/fanny/Work/EBI/Indigene/datafiles/Medaka_Gene_Database3.csv')

    # Filtering of the dataframes
    form_for_filtering(st.session_state['df'], st.session_state['go_df'])

    # TODO: charge bettina rank or dataset
    # TODO: perform cluster/correlation analysis
    # dataset = st.file_uploader("Choose a CSV file", type=['csv', 'txt', 'xls', 'xlsx'])
    # dict_table_read = dict(csv=pd.read_csv, txt=pd.read_table, xls=pd.read_excel, xlsx=pd.read_excel)
    # rank_df = dict_table_read[dataset.name.split('.')[-1]](dataset)
    # with st.form('Data association'):
    # rank_df=pd.read_excel('/home/fanny/Work/EBI/Indigene/MIKK-Panel_Heart_Rate_Data_Rank.xlsx')
    # if rank_df.filter(regex='[Ll]ine').shape[1] != 1:
    #     line_col = st.selectbox('Line ID column', rank_df.columns)
    # else:
    #     line_col = rank_df.filter(regex='line').columns[0]
    # rank_df.rename(columns={line_col:'line'}, inplace=True)
    # rank_df['line'] = rank_df['line'].str.replace('-', '_')
    # # cols_ls = st.multiselect('Columns to mix (must include the line id)', rank_df.columns, default='line')
    # cols_ls = ['response_rank', 'temperature_combined_rank']
    # st.write(rank_df.filter(regex='line'))
    # st.write(line_col)
    # rank_df=rank_df[set(cols_ls+['line'])]
    # # st.write(rank_df.head())
    # # var = st.selectbox('Variable to analyze', st.session_state['final_df'].columns)
    # var='gene'
    # # if st.form_submit_button('merge'):
    # # st.write(st.session_state['final_df'].groupby(['line', var]).nunique())
    # roux = pd.merge(st.session_state['final_df'].groupby(['line', 'go_id', 'gene'])['snp'].count().reset_index().fillna(0), rank_df, on='line')
    # roux.to_csv('Results_merge_LOF.csv', index=False)
    # st.dataframe(roux)

    # # Verification printing
    with st.expander('Verification prints'):
    #     cols = st.columns(2)
    #     cols[0].dataframe(st.session_state['final_df'].drop(['goslim_goa_description', 'goslim_goa_accession', 'alt_num', 'alt', 'allele2'], axis=1).drop_duplicates())
        st.dataframe(st.session_state['final_df']['gene'].unique())

    st.dataframe(st.session_state['final_df'].drop(['goslim_goa_description', 'goslim_goa_accession', 'alt_num', 'alt', 'allele2'], axis=1).drop_duplicates())
        # cols[1].dataframe(roux)
    #     cols[0].write(st.session_state['final_df']['line'].unique())
    #     cols[1].write(st.session_state['final_df']['variant'].unique())
    #     cols[2].write(st.session_state['final_df']['name_1006'].unique())
    #     cols[3].write(st.session_state['final_df']['go_id'].unique())
    #     cols[4].write(st.session_state['final_df']['gene'].unique())

    # fig = px.scatter(roux, x='response_rank', y='snp', color='gene', hover_data=['line', 'response_rank', 'gene', 'go_id'])
    # st.plotly_chart(fig)
    # fig = px.scatter(roux, x='temperature_combined_rank', y='snp', color='gene', hover_data=['line', 'temperature_combined_rank', 'gene', 'go_id'])
    # st.plotly_chart(fig)
    #
    # doux = pd.merge(st.session_state['final_df'].groupby(['line', 'go_id'])['gene'].count().reset_index().fillna(0), rank_df, on='line')
    # fig = px.scatter(doux.sort_values('response_rank'), x='response_rank', y='gene', color='go_id', hover_data=['line', 'response_rank', 'gene', 'go_id'])
    # st.plotly_chart(fig)
    # fig = px.scatter(doux.sort_values('temperature_combined_rank'), x='temperature_combined_rank', y='gene', color='go_id', hover_data=['line', 'temperature_combined_rank', 'gene', 'go_id'])
    # st.plotly_chart(fig)
    # gnou = pd.merge(st.session_state['final_df'].groupby(['line', 'go_id'])['gene'].count().unstack().reset_index().fillna(0), rank_df, on='line')
    # lut = dict(zip(gnou['response_rank'].unique(), sns.color_palette("viridis")))
    # row_colors = gnou['response_rank'].map(lut)
    # gnoux = gnou.set_index('line').drop(['temperature_combined_rank'], axis=1).copy()
    # gnoux['col'] = gnou['response_rank'].map(lut)
    # clust = sns.clustermap(gnoux, row_colors=row_colors)
    # st.pyplot(clust)

def cross():
    cross = pd.read_csv('/home/fanny/Work/EBI/Indigene/high_impact_gene_homref_included.txt')
    cross['alt_num'].unique()
    cross.groupby(['chr', 'gene', 'line']).mean()['allele1'].unstack()
    cross.groupby(['gene', 'line']).mean()['allele1'].unstack().std(axis=1)
    cross.groupby(['gene', 'line']).mean()['allele1'].unstack()[cross.groupby(['gene', 'line']).mean()['allele1'].unstack().std(axis=1) != 0]
    cross.groupby(['gene', 'line'])['allele1'].count()[cross.groupby(['gene', 'line'])['allele1'].count() == 1].unstack()


def blou():
### OPENING OF DF
    # Genome
    df = pd.read_csv('/home/fanny/Work/EBI/Indigene/grep_tempheart_genes.txt', header=None, names=['line', 'allele1', 'allele2', 'alt_num', 'snp', 'alt', 'gene', 'type', 'variant', 'impact'])
    df['variant'] = df['variant'].str.split(',', expand=True)[0]
    df['line'] = df['line'].str.replace('-', '_')
    #   GO
    go_df = pd.read_csv('/home/fanny/Work/EBI/Indigene/Temp_heartrate_gene_go_df.csv')
    df=pd.merge(go_df, df, on='gene')
    go_obo_file = '/home/fanny/Work/EBI/Indigene/datafiles/GO_obo/go-basic.obo'
    go = obo_parser.GODag(go_obo_file)
    df['go_name'] = df['go_id'].apply(lambda x: get_go_info(x, go))
    count = df.groupby(['go_id', 'gene', 'line', 'impact'])['snp'].count().unstack().reset_index()
    count_all = df.groupby(['go_id', 'gene', 'line'])['snp'].count().reset_index()

    #RNA expression
    samples = pd.read_csv('/home/fanny/Work/EBI/Indigene/Heart.samplesheet.csv')
    samples['sample'] = samples['sample'].str.split('_', n=1, expand=True)[0].str.replace('-', '.')
    samples['line'] = samples['fastq_1'].str.rsplit('/', n=1, expand=True)[1].str.split('-', expand=True)[2]
    dict_samples = samples.set_index('sample')['line'].to_dict()
    rna = pd.read_table('/home/fanny/Work/EBI/Indigene/salmon.merged.gene_tpm.tsv')
    samples['season'] = samples['sample'].str.split('.', n=2, expand=True)[1]
    rna.rename(columns=dict_samples, inplace=True)
    rna_df = pd.melt(rna[rna['gene_id'].isin(df['gene'].unique())], id_vars=['gene_name', 'gene_id'], var_name='line_id',
    value_name='rna_exp')
    rna_df.rename(columns={'line_id':'line', 'gene_id':'gene'}, inplace=True)
    rna_df = pd.merge(rna_df, samples[['line', 'season']], on='line')
    rna_df = rna_df[rna_df['season'] == 'S']
    rna_count = pd.merge(rna_df, count, on=['gene', 'line'], how='outer').drop(['season', 'gene_name'], axis=1)
    rna_count = pd.merge(rank_df, rna_count, on='line').fillna(0)
    rna_count.head(2)
    rna_count['TOT'] = rna_count[['HIGH', 'MODERATE', 'MODIFIER', 'LOW']].sum(axis=1)
    rna_count.corr()
    [['HIGH', 'rna_exp']].corr()
    sns.clustermap(rna_count.groupby(['response_rank', 'gene'])['HIGH'].sum().unstack().fillna(0), cmap='Blues', row_cluster=False)

    # RANK TEMP
    rank_df=pd.read_excel('/home/fanny/Work/EBI/Indigene/MIKK-Panel_Heart_Rate_Data_Rank.xlsx')
    rank_df.rename(columns={'line_id':'line'}, inplace=True)
    rank_df['line'] = rank_df['line'].str.replace('-', '_')
    rank_df = rank_df[['line', 'response_rank', 'temperature_combined_rank']]


    # ANALYSIS FOR RNA / RANK COMPARISON
    rna.head(2)
    rank_df.head(2)
    rna_df.head(2)

    rna_rank = pd.merge(rna_df, rank_df, on='line')
    rna_rank=pd.merge(go_df, rna_rank, on='gene')
    go_obo_file = '/home/fanny/Work/EBI/Indigene/datafiles/GO_obo/go-basic.obo'
    go = obo_parser.GODag(go_obo_file)
    rna_rank['go_name'] = rna_rank['go_id'].apply(lambda x: get_go_info(x, go))
    sns.catplot(data=rna_rank[rna_rank['season'] == 'S'].sort_values('response_rank'), kind='bar', x='gene', y='rna_exp', hue='response_rank', col='go_name', col_wrap=4, sharex=False, sharey=False, palette='viridis')
    sns.catplot(data=rna_rank.sort_values('response_rank'), kind='bar', x='gene', y='rna_exp', hue='response_rank', col='go_name', row='season', sharex=False, sharey=False, palette='viridis')
    sns.catplot(data=rna_rank.sort_values('temperature_combined_rank'), kind='bar', x='gene', y='rna_exp', hue='temperature_combined_rank', col='go_name', col_wrap=4, sharex=False, sharey=False, palette='flare')


    # ANALYSIS FOR VARIANT COUNTS
    df.drop(['allele1', 'allele2', 'alt_num', 'alt', 'type'], axis=1, inplace=True)
    df.head(2)

    df['go_name'] = df['go_id'].apply(lambda x: get_go_info(x, go))

    px.bar(df.groupby(['line','go_name', 'impact']).nunique().reset_index(),
                x='line', y='gene',  color='impact',
                barmode='stack', facet_col="go_name", facet_col_wrap=4, width=2000, height=2000)
    px.bar(df.groupby(['line','go_name', 'variant', 'impact']).nunique().reset_index(),
                x='line', y='gene',  color='variant',
                barmode='stack', facet_col="go_name", facet_col_wrap=4, width=2000, height=2000)

    blou = df.groupby(['line', 'gene', 'go_name', 'go_id', 'uniprot_gn_symbol'])['impact'].agg(lambda x: set(x)).reset_index()
    blou.loc[blou['impact'].astype(str).str.contains('HIGH'), 'impact'] = 'HIGH'
    blou.loc[blou['impact'].astype(str).str.contains('MODERATE'), 'impact'] = 'MODERATE'
    blou.loc[blou['impact'].astype(str).str.contains('LOW'), 'impact'] = 'LOW'
    blou.loc[blou['impact'].astype(str).str.contains('MODIFIER'), 'impact'] = 'MODIFIER'
    blou.groupby(['gene', 'impact', 'line']).count().tail(50)


    blou = pd.merge(rna_df.groupby(['line', 'gene']).mean(), blou, on=['line', 'gene'])
    blou.head(2)
    fig = px.bar( blou.groupby(['line', 'gene', 'go_name'])['rna_exp'].mean().reset_index(), y='rna_exp', x='gene', color='line', height=2000, width=1000, barmode='group', facet_col='go_name', facet_col_wrap=3)
    fig = px.box(blou, x='impact', y='rna_exp', height=1000, width=1000, facet_col='gene', facet_col_wrap=5, facet_row_spacing=0.01)
    fig.update_xaxes(matches=None)
    fig.update_yaxes(matches=None)

    rank_df = rank_df.set_index('line_id').filter(regex='rank').reset_index()
    rank_df['line'] = rank_df['line_id'].str.replace('-', '_')
    blou
    rank = pd.merge(blou, rank_df, left_on='line', right_on='line').drop('line_id', axis=1)
    rank['line'].unique()
    lut = dict(zip(rank['temperature_combined_rank'].unique(), sns.color_palette("viridis", 60)))
    row_colors = rank.set_index('line')['temperature_combined_rank'].drop_duplicates().map(lut)
    sns.clustermap(rank[rank['impact'] == 'LOW'].groupby(['line', 'go_name'])['gene'].nunique().unstack().fillna(0), row_colors=row_colors, row_cluster=False)



    g = sns.catplot(data=rna_df, kind='bar', y='line_id', x='rna_exp', col='gene_id', hue='line_id', dodge=False, sharex=False)
    rna_df[rna_df['gene_name'] == 'trpm7']
    df.shape
    merge = pd.merge(rna_df, df, on=['line', 'gene'])
    merge.head()




    df['line'].nunique()
    lut = dict(zip(df.set_index('line')['temperature_combined_rank'].astype(int).unique(), sns.color_palette("viridis", 60)))
    blou = df.groupby(['line', 'temperature_combined_rank', 'response_rank', 'gene'])['snp'].sum().unstack().fillna(0).reset_index()
    row_colors = pd.DataFrame()
    row_colors['temperature_combined_rank'] = merge.set_index('line')['temperature_combined_rank'].map(lut)
    row_colors['response_rank'] = merge.set_index('line')['response_rank'].map(lut)
    sns.clustermap(merge.drop([ 'snp', 'go_id', 'gene_name'], axis=1).groupby(['line','temperature_combined_rank', 'response_rank', 'gene']).mean()['rna_exp'].unstack().dropna(how='all', axis=1).fillna(0).reset_index().set_index('line').sort_values('response_rank').drop(['temperature_combined_rank', 'response_rank', 'ENSORLG00000027716', 'ENSORLG00000022886'], axis=1),row_cluster=False)
    sns.clustermap(rna.drop(['gene_name'], axis=1).groupby('gene_id').mean(), figsize=(5,15), row_colors=row_colors, row_cluster=False)


    # st.dataframe(st.session_state['filt_go_df'])
    # st.dataframe(st.session_state['filt_go_df'])

    # Analysis part
    # st.write('# Analysis and plots')
    # dict_part = {'Variants counts': plots_count, 'Total impacts': plots_count_tot}
    # part = st.selectbox('Type of analysis', dict_part.keys())
    # figs = dict_part[part](st.session_state['final_df'])
    # dict_fig_plots(figs, title=part)

    # go_df[go_df['name'].str.contains('temperature')]
    # go_ls = go['GO:0009266'].get_all_children()
    # go_df[go_df['go_id'].isin(go_ls | {'GO:0009266'})]

    # Check problems with res files
    # df[df['variant'].isna()][['snp', 'alt']].drop_duplicates()
    # list_snp = df[df['variant'].isna()][['snp', 'alt']].drop_duplicates()['snp'].tolist()
    # list_snp = str(list_snp)[1:-1]
    # blou = pd.read_sql(f'SELECT DISTINCT id, gene, alt, variant, impact FROM vep WHERE id LIKE "1_19930_A/%"', con=db)
    # blou.head(10)

    # # Associating information
    # db_name = 'datafiles/medaka_genome_part.db'
    # db = sqlite3.connect(db_name)
    # df = Merge_Gene_Info(df, db)

    # #### WAY TO GET THE GO UNIQUE ID PER GO_ID
    # df['go_name'] = df['go_id'].apply(lambda x: get_go_info(x))
    # df = Merge_GO_Info(df)

#     # CLUSTERMAPS
#     lut = dict(zip(set(blou['cross']), sns.color_palette('husl', len(set(blou['cross'])))))
#     row_colors = blou['cross'].map(lut)
#     row_col = pd.concat([row_col, row_colors], axis=1)
#     sns.clustermap(blou[(blou['go_name'].str.contains('heart')) & blou['line'].isin(['72-2', '79-2', '15-1', '139-4', '55-2', '62-2', '68-1', '22-1', 'iCab-F25'])].groupby(['go_name', 'line']).size().unstack().fillna(0), figsize=(20,15), row_cluster=False, cmap="vlag", row_colors=row_colors)
#     sns.clustermap(blou[(blou['go_name'].str.contains('ethanol'))].groupby(['line', 'go_name']).size().unstack().fillna(0), figsize=(10,25), row_cluster=False, cmap="vlag")

## ## Query to the obo databases
# # Import the OBO parser from GOATools
# import wget
# import os

if __name__ == '__main__':
    main()
