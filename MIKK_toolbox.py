import streamlit as st
import argparse 
import yaml 
import os

st.set_page_config(
    page_title="MIKK Toolbox",
    page_icon="🐟",
    layout="wide",
)


def main(yaml_file=''):

    st.header('Welcome to the MIKK panel Toolbox')

    # Setting sesssion state
    if yaml_file:
        with open(yaml_file, 'r') as y:
            params = yaml.safe_load(''.join(y.readlines()))
            st.session_state.update(params)
        # st.write(st.session_state)
        st.success('YAML info loaded successfully')
    else:
        with st.container(border=True):
            st.page_link("pages/4_Modify_parameters.py", label="**Upload your parameters here**", icon="➡️", use_container_width=True)
        
    if ('output_dir' in st.session_state):
        if os.path.isdir(st.session_state['output_dir']) == False:
            os.mkdir(st.session_state['output_dir'])

    with open('README.md', 'r') as f:
        readme = f.readlines()
    st.markdown(''.join(readme))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
                        prog='MIKK Toolbox',
                        description='Description: ')
    parser.add_argument('-y', '--yaml_file', required=False, help="YAML file containing the arguments")     
    args = parser.parse_args()
    main(**vars(args))
    # main()