import os
import streamlit as st
from io import StringIO
import yaml 
st.set_page_config(
    page_title="MIKK Toolbox",
    page_icon="🐟",
    layout="wide",
)

def save_yaml(dict_params, cont=st):
    try:
        with open(st.session_state['yml_filename'], 'w') as file:
            yaml.dump(dict_params, file)
        cont.success(st.session_state['yml_filename'] + ' saved successfully !')
    except Exception as e:
        cont.exception(e)


def main():
    st.header("Choose a YAML file")
    yaml_file = st.file_uploader("Choose a YAML file", key='yaml_file', label_visibility="collapsed")
    if yaml_file is not None:
        for key in st.session_state.keys(): del st.session_state[key]
        params = yaml.safe_load(yaml_file)
        st.session_state.update(params)
        st.success('Parameters loaded successfully !')
    st.header('OR')
    cont = st.container(border=True)
    with cont.form('General params', border=False):
        st.header('Enter paths manually')
        params = {}
        for var in ['output_dir', 'vcf_file', 'fasta_dir', 'sample_file', 'flexlmm_dir', 'phenotype_file', 'covariate_file', 'info_file', 'vep_file']:
            params[var] = st.text_input(var.replace('_', ' ').replace('dir', 'directory').title(), value='')
        submit_params = st.form_submit_button('Submit', use_container_width=True)
    if submit_params:
        for key in st.session_state.keys(): del st.session_state[key]
        st.session_state.update({key:item for key, item in params.items() if item})
        st.success('Parameters updated successfully !')
    
        with cont.form('Save config', border=False):
            save_param_cols = st.columns([3,1])
            save_param_cols[0].text_input('Saving path', value=os.path.join(os.getcwd(), 'config.yml'), label_visibility='collapsed', key='yml_filename')
            save_param_cols[1].form_submit_button('Save', use_container_width=True, on_click=save_yaml, args=(params, cont))
    
    # st.write(st.session_state)

if __name__ == '__main__':
    main()
