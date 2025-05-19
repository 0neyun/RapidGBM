
import streamlit as st
#import ollama
import matplotlib.pyplot as plt
import numpy as np
import os

#try:
    #if st.session_state["authentication_status"]:
        #pass
#except:
    #st.switch_page("login.py")

def mkdir(path):
    folder = os.path.exists(path)
    if folder:
        pass
    else:
        os.makedirs(path)
        
user_path = st.session_state["name"]+'_Amati'
mkdir(user_path)

botton0 = st.sidebar.button('Back mune')
if botton0:
    st.switch_page("menu.py")    
    
st.title('Amati relation with Redshift 🐰')
st.markdown("""
- **Ep,z - Eiso,gamma (Amati)** relation: Intrinsic spectra and energetics of BeppoSAX Gamma–Ray Bursts with known redshifts [DOI: 10.1051/0004-6361:20020700](https://www.aanda.org/articles/aa/abs/2002/28/aah3486/aah3486.html)
- Data from：**Minaev & Pozanenko (2020a, 2020b)**  
   The Ep,i–Eiso correlation: type I gamma-ray bursts and the new classification method  
    [DOI:10.1093/mnras/stz3611](https://academic.oup.com/mnras/article/492/2/1919/5687348)  
   GRB 200415A: Magnetar Giant Flare or Short Gamma-Ray Burst?  
    [DOI:10.1134/S1063773720090042](https://link.springer.com/article/10.1134/S1063773720090042)
""")


st.sidebar.markdown(r"## Set $E_p$ and $S_{\gamma}$")
str_Ep = st.sidebar.text_input(r'$E_p$ [keV]', '100')
str_Fluence = st.sidebar.text_input(r'$S_{\gamma}$ (ingore k-correlation) [erg/cm$^2$]', '1e-5')

st.sidebar.markdown(r"### If You known redshift")
str_z = st.sidebar.text_input(r'z', '1')

st.sidebar.markdown(r"### If You Don't known redshift")
plot_z_line = st.sidebar.checkbox('Plot the curve as a function of redshift?')
z_range = st.sidebar.text_input(r'interpolation range', '0.001, 9')

botton1 = st.sidebar.button('Plot')

def linesplit_2(line):
    line = line.split(',')
    return float(line[0]), float(line[1])
    
import sys
sys.path.append('./base_func')
from amati_plot import plot_Digram,plot_line,plot_point,plot_text

#@st.cache_data(max_entries=5)
def plot_amati(Ep, Fluence, z, z_range):
    plot_Digram()
    #if Ep is not None and Fluence is not None:
    Ep = np.float64(Ep)
    Fluence = np.float64(Fluence)

    if plot_z_line:
        z1, z2 = linesplit_2(z_range)
        z_list = np.logspace(np.log10(z1),np.log10(z2), 300)
        #z_list = np.linspace((z1),(z2), 300)
        plot_line(z_list=z_list,Ep=Ep,fluence=Fluence,color='k')
        z_list2 = np.logspace(np.log10(z1),np.log10(z2), 6)
        #z_list2 = np.linspace((z1),(z2), 6)
        plot_text(z_list=z_list2,Ep=Ep,fluence=Fluence,color='k')
        
    #if z is not None:
    z = np.float64(z)
    plot_point(z=z,Ep=Ep,fluence=Fluence, color='k')

    plt.savefig(f'{user_path}/_amati.png')
        
    st.image(f'{user_path}/_amati.png', use_column_width=True)
if botton1:
    plot_amati(Ep=str_Ep, Fluence=str_Fluence, z=str_z, z_range=z_range)
    
