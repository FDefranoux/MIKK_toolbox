import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import argparse
import os
import sys
import re
import streamlit as st
import plotly.graph_objects as go
import plotly.express as px

import plotly.graph_objects as go
import plotly.express as px
import pysam
st.set_page_config(
    page_title="MIKK Toolbox",
    page_icon="🐟",
    layout="wide",
)

def find_region_name(x):
    try:
        name = x.split(';')[0]
        name = re.sub('.*(ID|id).*:', '', name)
    except IndexError:
        name = re.sub('.*Parent', '', x)
    except:
        print('Error name to recognized: ', x)
        name=''
    return name


@st.cache_data()
def get_info(info_file, chr, start, end):
    df = pd.DataFrame()
    try:
        tbx = pysam.TabixFile(info_file)
        for n, rec in enumerate(tbx.fetch(str(chr), start, end)):
            info_pos = rec.split('\t')
            df[n] = info_pos
    except:
        pass
    df = df.T
    if df.shape[1] != 0:
        df.columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
        df['name'] = df['attributes'].apply(find_region_name)
        df['type_name'] = df['type'] + '_' + df['name']
        return df 
    else:
        return pd.DataFrame()

@st.cache_data()
def genomic_info_plot2(info_df, mark_ls=[], legend_type='color'):
    # Recuperation of the genomic info
    # dict_color = {type: color for color, type in zip(px.colors.qualitative.Set1, info_df['type'].unique())}
    st.write(info_df.head())
    dict_height = {name_type: n + 10 for n, name_type in enumerate(info_df['type'].unique())}
    info_df['height'] = info_df['type'].map(dict_height)

    # Modify info_df to plot 
    info_df['coord-x'] = (info_df['start']).astype(str) + '_' + (info_df['start']).astype(str) + '_' + (info_df['end']).astype(str) + '_' + (info_df['end']).astype(str) + '_' + (info_df['start']).astype(str) + '_ '
    info_df['coord-y'] = info_df['height'].astype(str) + '_' + (info_df['height'] +1).astype(str) + '_' + (info_df['height']+1).astype(str) + '_' + info_df['height'].astype(str) + '_' + info_df['height'].astype(str) + '_ '

    # Create the figure
    fig = go.Figure()

    # Plots the genomic regions
    info_df.sort_values(['type', 'start', 'end'], inplace=True)
    for row in info_df.index:
        
        fig.add_trace(
            go.Scatter(
                x=info_df.loc[row, 'coord-x'].split('_'),
                y=info_df.loc[row, 'coord-y'].split('_'),
                # legendgroup=dict_color[info_df.loc[row, 'type']],
                # legend='legend1',
                # mode='lines',
                # line=dict(color=dict_color[info_df.loc[row, 'type']], width=1),
                name=info_df.loc[row, 'type_name'],
                text=info_df.loc[row, 'name'],
                # mode="text",
                fill="toself",
                textposition="middle center",
                # hovertemplate='<b>%{text}</b>', 
                opacity=0.5,
                # hoveron="fills",
                customdata = info_df.loc[row].tolist(),
                # hovertemplate='Name:%{customdata[9]}<br><b>CHR:%{customdata[0]} <br><b>Start:%{customdata[3]}</b><br>End: %{customdata[4]}</b><br>Other: %{customdata[8]} ',
                showlegend=True
            )
            )

    # Add positional marker as blue vertical lines
    for mark in mark_ls:
        if mark != None:
            fig.add_shape(
                go.layout.Shape(
                    type='line',
                    name='Marker position',
                    x0=mark,
                    x1=mark,
                    y0=info_df['height'].min(),
                    y1=info_df['height'].max() +1,
                    line=dict(color='red', width=2),
                    showlegend=True
                )
            )

    # Get a legend by color
    if legend_type == 'color':
        fig.update_yaxes(visible=False)
        fig.update_layout(legend=dict(traceorder='reversed', x=1, y=0.5))

    return fig

@st.cache_data()
def get_df(file, chr, start, end):
    import subprocess 
    res = subprocess.Popen(f"zcat {file} | awk -F_ '$1 == {chr} && $2 >= {start} && $2 <= {end}'", stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    res_out = res.stdout.read().decode('utf8')
    res_out = res_out.split('\n')
    res_out = [row.split('\t') for row in res_out if row != '']
    return pd.DataFrame(res_out)

@st.cache_data()
def load_vep(vep_file, chr, start, end):
    import plotly.figure_factory as ff
    df = get_df(vep_file, chr, start, end)
    st.write(df)
    if df.iloc[:, 0].dtypes == 'str':
        df[['chr', 'pos']] = df.iloc[:, 0].str.split('_', expand=True, n=3)[[0,1]].astype(int)
    df['impact'] = df.iloc[:, 13].str[7:].str.split(';', n=1,expand=True)[0]
    df['gene_transc'] = df.iloc[:, 3] + '_' + df.iloc[:, 4] + '_' + df.iloc[:, 6]
    impact_ls = ['HIGH', 'MODERATE', 'LOW', 'MODIFIER']
    group_labels = ['HIGH', 'MODERATE', 'LOW', 'MODIFIER']
    
    data_ls = []
    rug_text = []
    for imp in impact_ls:
        sub_df = df[df['impact'] == imp]
        if sub_df.shape[0] == 0:
            group_labels.remove(imp)
        else:
            data_ls.append(sub_df['pos'].tolist())
            rug_text.append(sub_df['gene_transc'].tolist()) 
    
    colors_dict = {'HIGH':'#e62c11', 'MODERATE':'#e6a311', 'LOW':'#dce611', 'MODIFIER':'#11e6b6'}
    colors = [col for name,col in colors_dict.items() if name in group_labels]
    fig = ff.create_distplot(data_ls, group_labels=group_labels, colors=colors, show_curve=False, show_hist=False, rug_text=rug_text)
    st.write(df.head())
    return df, fig

def vep_filters():
    vep_df = st.session_state['vep_df'].drop([0,1,2,7,8,9,10,11,12,13,'chr','gene_transc'], axis=1)
    vep_df.rename(columns={3:'gene', 4:'transcript', 5:'type', 6:'variant'}, inplace=True)
    st.session_state['vep_cols'] = st.multiselect('Columns', vep_df.columns)
    ls_filt = vep_df.iloc[:, 0] == vep_df.iloc[:, 0]
    with st.form('vep_filt', border=False):
        for col in st.session_state['vep_cols']:
            vals = st.multiselect(f'Filtering criteria for {col}', vep_df[col].unique())
            ls_filt &= vep_df[col].isin(vals)
        submit_filt = st.form_submit_button('Submit filters')
    if submit_filt:
        st.write(vep_df[ls_filt].drop_duplicates(), use_container_width=True)


def initialisation(chr, start, end, marker, vep):
    
        st.session_state['chr'] = chr
        st.session_state['start'] = start
        st.session_state['end'] = end
        st.session_state['marker'] = marker
        st.session_state['load_vep'] = vep


def save_page(fig, name_file='region_graph.html'):
    with st.form('Html page saving', border=True):
        save_cols = st.columns([2,1,1])
        default_name = os.path.join(os.path.join(st.session_state['output_dir'], 'Region_plot_' + name_file + ".html"))
        name_file = save_cols[0].text_input('Enter file path to save', value=default_name, label_visibility='collapsed')
        saving = save_cols[1].form_submit_button('**Save Figure**', use_container_width=True)
    if saving :
        err=''
        try:
            # fig.to_html(full_html=True, include_plotlyjs='cdn')
            fig.write_html(name_file)

        except Exception as err:
            pass
        if (os.path.exists(name_file) == True):
            st.success(f'The file {name_file} has been saved successfully!')
        elif (os.path.exists(name_file) == False):
            st.error(f'There has been an issue while saving file {name_file}')
            if err:
                st.error(str(err))


def main():
    # st.header('Region information')
    if 'submit_info' not in st.session_state:
        st.info('Submit the data you want to analyse on the left-side sidebar.')
        st.session_state['submit_info'] = False

    if len(set(['info_file', 'vep_file']) - set(st.session_state.keys())) > 0: 
        with st.expander('Files parameters').form('General params', border=False):
            st.header('Enter paths manually')
            params = {}
            for var in ['output_dir', 'info_file']:
                if (var not in st.session_state) or (st.session_state[var] == ''):
                    st.write(var in st.session_state)
                    params[var] = st.text_input(var.replace('_', ' ').replace('dir', 'directory').title(), value='')
            submit_params = st.form_submit_button('Submit', use_container_width=True)
            if submit_params:
                st.session_state.update(params)    
    with st.sidebar.form('Genomic regions'):
        chr = st.number_input('Chromosome', min_value=1, max_value=24, step=1)
        start = st.number_input('Start', min_value=1, max_value=100000000, step=10000)
        end = st.number_input('End', min_value=1, max_value=100000000, step=10000)
        marker = st.number_input('Maker position', min_value=0, max_value=100000000, value=None)
        vep = st.checkbox('Load VEP information')
        submit_info = st.form_submit_button('Submit')
    if submit_info:
        initialisation(chr, start, end, marker, vep)
        st.session_state['submit_info'] = True

    if st.session_state['submit_info']:
        from plotly.subplots import make_subplots
        fig = make_subplots(rows=2, cols=1, shared_xaxes=True, vertical_spacing=0.02)
        st.session_state['info_df'] = get_info(st.session_state['info_file'] , st.session_state['chr'], st.session_state['start'], st.session_state['end'])
        if st.session_state['info_df'].shape[0] != 0:
            st.session_state['fig1'] = genomic_info_plot2(st.session_state['info_df'], legend_type='color', mark_ls=[st.session_state['marker']])
            for i in st.session_state['fig1'].data :
                fig.add_trace(i, row=1, col=1)
        
        if ('vep_file' in st.session_state) and st.session_state['load_vep'] == True:
            title= f"{st.session_state['chr']}:{st.session_state['info_df']['start'].min()}-{st.session_state['info_df']['end'].max()}"
            st.header(f"Genomic plot for region {title}")
            st.session_state['vep_df'], st.session_state['fig2'] = load_vep(st.session_state['vep_file'], st.session_state['chr'], st.session_state['info_df']['start'].min(),st.session_state['info_df']['end'].max())
            for i in st.session_state['fig2'].data :    
                fig.add_trace(i, row=2, col=1)

        st.plotly_chart(fig, use_container_width=True)
        save_page(fig, name_file=title)
        
        with st.expander('More information: '):
            st.write(st.session_state['info_df'], use_container_width=True)
            if ('vep_file' in st.session_state) and st.session_state['load_vep'] == True:
                vep_filters()


if __name__ == '__main__':
    main()