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
st.set_page_config(layout="wide")

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


# @st.cache_data()
def get_info(info_file, chr, start, end):
    print(chr, start, end)
    df = pd.DataFrame()
    try:
        tbx = pysam.TabixFile(info_file)
        for n, rec in enumerate(tbx.fetch(str(chr), start, end)):
            info_pos = rec.split('\t')
            df[n] = info_pos
    except:
        pass
    df = df.T
    df.columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    df['name'] = df['attributes'].apply(find_region_name)
    df['type_name'] = df['type'] + '_' + df['name']
    return df 


def genomic_info_plot2(info_df, mark_ls=[], legend_type='color'):
    # Recuperation of the genomic info
    dict_color = {type: color for color, type in zip(px.colors.qualitative.Set1, info_df['type'].unique())}
    dict_height = {name_type: n + 10 for n, name_type in enumerate(info_df['type'].unique())}
    info_df['height'] = info_df['type'].map(dict_height)

    # Modify info_df to plot 
    info_df['coord-x'] = (info_df['start']).astype(str) + '_' + (info_df['start']).astype(str) + '_' + (info_df['end']).astype(str) + '_' + (info_df['end']).astype(str) + '_' + (info_df['start']).astype(str) + '_ '
    info_df['coord-y'] = info_df['height'].astype(str) + '_' + (info_df['height'] +1).astype(str) + '_' + (info_df['height']+1).astype(str) + '_' + info_df['height'].astype(str) + '_' + info_df['height'].astype(str) + '_ '

    # Create the figure
    fig = go.Figure()

    # Plots the genomic regions
    place = 'top'
    info_df.sort_values(['type', 'start', 'end'], inplace=True)
    for row in info_df.index:
        
        fig.add_trace(
            go.Scatter(
                x=info_df.loc[row, 'coord-x'].split('_'),
                y=info_df.loc[row, 'coord-y'].split('_'),
                # legendgroup=dict_color[info_df.loc[row, 'type']],
                # legend='legend1',
                # mode='lines',
                line=dict(color=dict_color[info_df.loc[row, 'type']], width=1),
                name=info_df.loc[row, 'type_name'],
                text=info_df.loc[row, 'name'],
                # mode="text",
                fill="toself",
                textposition="middle center",
                hovertemplate='<b>%{text}</b>', 
                opacity=0.5,
                hoveron="fills",
                customdata = info_df.loc[row].tolist(),
                # hovertemplate='Name:%{customdata[9]}<br><b>CHR:%{customdata[0]} <br><b>Start:%{customdata[3]}</b><br>End: %{customdata[4]}</b><br>Other: %{customdata[8]} ',
                showlegend=True
            )
            )

    # Add positional marker as blue vertical lines
    for mark in mark_ls:
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


def initialisation(chr, start, end, file_info, marker):
        st.session_state['chr'] = chr
        st.session_state['start'] = start
        st.session_state['end'] = end
        st.session_state['file_info'] = file_info
        if marker > 0:
            st.session_state['marker'] = marker
        else:
            st.session_state['marker'] = None


def main():
    with st.sidebar.form('Genomic regions'):
        chr = st.number_input('Chromosome', min_value=1, max_value=24, step=1)
        start = st.number_input('Start', min_value=1, max_value=100000000, step=10000)
        end = st.number_input('End', min_value=1, max_value=100000000, step=10000)
        file_info = st.text_input('Path', value='/home/fanny/Documents/Work/HeartRate_Tox/finemapping/Oryzias_latipes-GCA_002234675.1-2022_04-genes_clean.sorted.colored.gff3.gz') 
        marker = st.number_input('Maker position', min_value=0, max_value=100000000)
        submit = st.form_submit_button('Submit')
    if submit:
        initialisation(chr, start, end, file_info, marker)
        st.session_state['submit'] = True
        try:
            st.sidebar.write(st.session_state['chr'], st.session_state['start'], st.session_state['end'])
        except:
            pass
    if not 'submit' in st.session_state:
        st.session_state['submit'] = False
    if st.session_state['submit']:
        st.session_state['info_df'] = get_info(st.session_state['file_info'] , st.session_state['chr'], st.session_state['start'], st.session_state['end'])
    fig2 = genomic_info_plot2(st.session_state['info_df'], legend_type='color', mark_ls=[st.session_state['marker']])
    st.plotly_chart(fig2, use_container_width=True)
    with st.expander('More information: '):
        st.write(st.session_state['info_df'], use_container_width=True)


if __name__ == '__main__':
    main()