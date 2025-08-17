import streamlit as st


st.session_state["name"] = 'anonymous'

st.write("# Welcome to HEtools web!")
st.write('Hello %s.' %(st.session_state["name"]))
st.image("pages/image.png", use_column_width=False,\
         width=500, \
         caption="In Chinese, the word for “rabbit” is pronounced the same as “tools” in English. Interestingly, this rabbit is actually my pet.")    

st.sidebar.page_link("pages/RapidGBM.py", label="🛰️ RapidGBM")
