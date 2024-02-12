import streamlit as st
import pandas as pd

def confirmation_box(onClick, click_kwargs=dict(), container=st, ):
    confirmation_container = container.empty()
    with confirmation_container.container(border=True):
        container.error("Are you sure?")
        col_conf = container.columns(2)
        col_conf[0].button("Yes", on_click=onClick, kwargs=click_kwargs)
        col_conf[1].button('Cancel')



def filter_dataframe(df, key='filtering', cont=st, session_state_var='filtered_df'):
    from pandas.api.types import (
        is_categorical_dtype,
        is_datetime64_any_dtype,
        is_numeric_dtype,
        is_object_dtype,
    )
    df = df.copy()

    # Try to convert datetimes into a standard format (datetime, no timezone)
    for col in df.columns:
        if is_object_dtype(df[col]):
            try:
                df[col] = pd.to_datetime(df[col], format='mixed')
            except Exception:
                pass

        if is_datetime64_any_dtype(df[col]):
            df[col] = df[col].dt.tz_localize(None)

    with cont:
        to_filter_columns = cont.multiselect("Filter dataframe on", df.columns)
        with st.form(key, border=False):
            ls_filter = df.iloc[:, 0] == df.iloc[:, 0]
            for column in to_filter_columns:
                left, right = st.columns((1, 20))
                # Treat columns with < 10 unique values as categorical
                if is_categorical_dtype(df[column]) or df[column].nunique() < 10:
                    user_cat_input = right.multiselect(
                        f"Values for {column}",
                        df[column].unique(),
                        default=list(df[column].unique()),
                    )
                    ls_filter &= df[column].isin(user_cat_input)
                    # df = df[df[column].isin(user_cat_input)]
                elif is_numeric_dtype(df[column]):
                    _min = float(df[column].min())
                    _max = float(df[column].max())
                    step = (_max - _min) / 100
                    user_num_input = right.slider(
                        f"Values for {column}",
                        min_value=_min,
                        max_value=_max,
                        value=(_min, _max),
                        step=step,
                    )
                    ls_filter &= df[column].between(*user_num_input)
                    # df = df[df[column].between(*user_num_input)]
                elif is_datetime64_any_dtype(df[column]):
                    user_date_input = right.date_input(
                        f"Values for {column}",
                        value=(
                            df[column].min(),
                            df[column].max(),
                        ),
                    )
                    if len(user_date_input) == 2:
                        user_date_input = tuple(map(pd.to_datetime, user_date_input))
                        start_date, end_date = user_date_input
                        # df = df.loc[df[column].between(start_date, end_date)]
                        ls_filter &= df[column].between(start_date, end_date)
                else:
                    user_text_input = right.text_input(
                        f"Substring or regex in {column}",
                    )
                    if user_text_input:
                        ls_filter &= df[column].astype(str).str.contains(user_text_input)
                        # df = df[df[column].astype(str).str.contains(user_text_input)]
            filter_submit = st.form_submit_button('Filter')
    if filter_submit:
        if session_state_var:
            try:
                st.session_state[session_state_var] = df[ls_filter]
            except:
                st.session_state[session_state_var] = df
                
        else:
            try:
                return df[ls_filter]
            except:
                return df
