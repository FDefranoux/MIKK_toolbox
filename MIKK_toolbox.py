import streamlit as st

st.set_page_config(
    page_title="MIKK Toolbox",
    page_icon="🐟",
    layout="wide",
)


def main():
    with open('README.md', 'r') as f:
        readme = f.readlines()
    st.markdown(''.join(readme))

if __name__ == '__main__':
    main()