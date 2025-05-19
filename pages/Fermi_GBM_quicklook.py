import streamlit as st
import glob, os
import matplotlib.pyplot as plt
import sys
sys.path.append('./base_func')
from quicklook_base_func import * 

#try:
    #if st.session_state["authentication_status"]:
        #pass
#except:
    #st.switch_page("login.py")

user_path = st.session_state["name"]+'_quicklook'
mkdir(user_path)
File_for_Fit =  f'{user_path}/0_File_for_Fit'
Result_path = f'{user_path}/1_Result_path'

back_btn  = st.sidebar.button("Back menu", key="back_menu")
clear_btn = st.sidebar.button("🗑 CLEAR USER FOLDER",key="clear_all",)

if back_btn:
    preserve = {'authentication_status', 'name'}
    for k in list(st.session_state.keys()):
        if k not in preserve:
            del st.session_state[k]    
    st.switch_page("menu.py")

import shutil
if clear_btn:
    for entry in os.listdir(user_path):
        entry_path = os.path.join(user_path, entry)
        try:
            if os.path.isfile(entry_path) or os.path.islink(entry_path):
                os.remove(entry_path)        
            elif os.path.isdir(entry_path):
                shutil.rmtree(entry_path)     
        except Exception as e:
            st.sidebar.error(f"Error deleting {entry}: {e}")
    st.sidebar.success("All files and folders under user_path have been deleted.")

mkdir(File_for_Fit)
mkdir(Result_path)

# Title the app
st.title('RapidGBM 🐰')

st.markdown("""
 * Use the menu at left to select data
 * Your plots will appear below
""")


######################################################################################################
check_pos_botton = False
plot_lc_botton = False

######################################################################################################
#part 1
st.sidebar.header("Set time and location")
str_t0 = st.sidebar.text_input('UTC Time', '2024-06-17T12:19:13.00')
str_coordinates = st.sidebar.text_input('Input source coordinates (RA, DEC in degree)', '285.03,-22.56')
#pre_orbit = st.sidebar.text_input('Input the number of orbits before', '30')
pre_orbit = 30
check_spechist_pre_orbit = st.sidebar.checkbox(
    'Use the earlier poshist file (30 orbits prior)?',
    value=False
)
check_pos_botton = st.sidebar.button('check location')
st.sidebar.divider()
######################################################################################################
#part 2
# 1. Multi-select detectors
st.sidebar.header("Set plot for lightcurves")
det_select = st.sidebar.multiselect(
    label="### Which detector(s)?",
    options=det_list,
    default=[],
)

# 2. NaI energy range in two columns
with st.sidebar.expander("NaI Energy Band (default: 8-900 keV)", expanded=False):
        nai_col1, nai_col2 = st.columns(2)
        nai_e_min = nai_col1.number_input(
            "Min", min_value=0.0, value=8.0, step=1.0,
            format="%.1f", key="nai_e_min"
        )
        nai_e_max = nai_col2.number_input(
            "Max", min_value=0.0, value=900.0, step=10.0,
            format="%.1f", key="nai_e_max"
        )
        energy_range1 = (nai_e_min, nai_e_max)


# 3. BGO energy range in two columns
with st.sidebar.expander("BGO Energy Band (default: 200-40000 keV)", expanded=False):
    bgo_col1, bgo_col2 = st.columns(2)
    bgo_e_min = bgo_col1.number_input(
        "Min", min_value=0.0, value=200.0, step=10.0,
        format="%.1f", key="bgo_e_min"
    )
    bgo_e_max = bgo_col2.number_input(
        "Max", min_value=0.0, value=40000.0, step=100.0,
        format="%.1f", key="bgo_e_max"
    )
    energy_range2 = (bgo_e_min, bgo_e_max)

# 4. Time bin size
st.sidebar.markdown("### Set Lightcurves Interval")
tbin = st.sidebar.number_input(
    label="Time bin size (s)",
    min_value=1e-4,
    max_value=100.0,
    value=2.0,
    step=0.01,
    format="%.3f",
)

# 5. Time range
t_col1, t_col2 = st.sidebar.columns(2)
t_start = t_col1.number_input(
    label="Time start (s)",
    value=-100.0,
    step=1.,
    format="%.1f",
)
t_end = t_col2.number_input(
    label="Time end (s)",
    value=500.0,
    step=1.,
    format="%.1f",
)
lc_time_range = (t_start, t_end)
if t_end <= t_start:
    st.sidebar.error("Error: End time must be greater than start time.")
    
plot_lc_botton = st.sidebar.button('plot lightcurves')
Build_Files_botton = False

col1, col2 = st.sidebar.columns([4, 1])
col1.markdown("### Subtract background & Build files")
subtract_bg = col2.checkbox("", value=False, key="subtract_bg")

if subtract_bg:
    b1_col1, b1_col2 = st.sidebar.columns(2)
    b1_start = b1_col1.number_input(
        "Bkg1 Start", 
        value=-100.0, 
        key="b1_start", 
        format="%.1f", 
        step=1.0
    )
    b1_end = b1_col2.number_input(
        "Bkg1 End", 
        value=-50.0, 
        key="b1_end", 
        format="%.1f", 
        step=1.0
    )
    if b1_end <= b1_start:
        st.sidebar.error("Error: Bkg1 end must be > start.")
    background1 = (b1_start, b1_end)

    b2_col1, b2_col2 = st.sidebar.columns(2)
    b2_start = b2_col1.number_input(
        "Bkg2 Start", 
        value=450.0, 
        key="b2_start", 
        format="%.1f", 
        step=1.0
    )
    b2_end = b2_col2.number_input(
        "Bkg2 End", 
        value=500.0, 
        key="b2_end", 
        format="%.1f", 
        step=1.0
    )
    if b2_end <= b2_start:
        st.sidebar.error("Error: Bkg2 end must be > start.")
    background2 = (b2_start, b2_end)
    
    snr_enabled = st.sidebar.checkbox(
        "plot SNR", 
        value=True, 
        key="snr_enabled",
    )
    snr_threshold = st.sidebar.number_input(
            "SNR threshold",
            min_value=0.0,
            value=5.0,
            step=1.,
            format="%.1f",
            key="snr_threshold",
        )
    
    gen_spec = st.sidebar.checkbox(
        "Generate spectrum", 
        value=False, 
        key="gen_spectrum",
    )    
    spectrum_time_range = None
    if gen_spec:
        st.sidebar.markdown("### Spectrum Time Range (s)")
        spec_col1, spec_col2 = st.sidebar.columns(2)
        spec_start = spec_col1.number_input(
            "Spec Start", 
            #value=-0.3,
            value=252.,
            key="spec_start", 
            format="%.2f", 
            step=0.1,
        )
        spec_end = spec_col2.number_input(
            "Spec End", 
            #value=.3,
            value=288.,
            key="spec_end", 
            format="%.2f", 
            step=0.1,
        )
        if spec_end <= spec_start:
            st.sidebar.error("Error: Spectrum end time must be greater than start time.")
        spectrum_time_range = (spec_start, spec_end)
        
        poshist_files = [f for f in os.listdir(user_path) if 'glg_poshist_all' in f]
        default = st.session_state.get('poshist_filename', None)
        if default in poshist_files:
            default_index = poshist_files.index(default)
        else:
            default_index = 0
        selected_poshist = st.sidebar.selectbox(
            label="Select poshist file",
            options=poshist_files,
            index=default_index,
            key="selected_poshist"
        )


        try:
            value_temp = st.session_state.met_time_reCal
        except:
            value_temp = None
        
        RSP_time = st.sidebar.number_input(
            "RSP_MET",
            value=value_temp,
            key="RSP_time", 
            format="%.2f",        
        )
        Build_Files_botton = st.sidebar.button('Build Files (PHA/BAK/RSP) for Fit')
    else:
        RSP_time = None
        selected_poshist = None
    
else:
    background1 = None
    background2 = None
    snr_enabled = True
    snr_threshold = None
    spectrum_time_range = None
    gen_spec = False
    RSP_time = None
    selected_poshist = None

st.sidebar.divider()

st.sidebar.header("Spectral Analysis")
st.sidebar.markdown("### Upload Files")

uploaded_files = st.sidebar.file_uploader(
    label="Upload files",
    type=None,                   
    accept_multiple_files=True
)

if uploaded_files:
    for uploaded_file in uploaded_files:
        save_path = os.path.join(File_for_Fit, uploaded_file.name)
        with open(save_path, "wb") as f:
            f.write(uploaded_file.getbuffer())
    st.sidebar.success(f"Saved {len(uploaded_files)} files.")

if 'spectra' not in st.session_state:
    st.session_state['spectra'] = []

if st.sidebar.button("Add dataset"):
    st.session_state['spectra'].append({'pha': None, 'bak': None, 'rsp': None})

pha_paths = sorted(glob(os.path.join(File_for_Fit, "*.pha*"))+glob(os.path.join(File_for_Fit, "*.pi*")))
bak_paths = sorted(glob(os.path.join(File_for_Fit, "*.bak")))
rsp_paths = sorted(glob(os.path.join(File_for_Fit, "*.rsp")))

pha_files = [os.path.basename(p) for p in pha_paths] + ['None']
bak_files = [os.path.basename(p) for p in bak_paths] + ['None']
rsp_files = [os.path.basename(p) for p in rsp_paths] + ['None']

for idx, cfg in enumerate(st.session_state['spectra']):
    st.sidebar.markdown(f"**Dataset {idx+1}**")
    if pha_files:
        cfg['pha'] = st.sidebar.selectbox(
            f"PHA file #{idx+1}",
            options=pha_files,
            index=pha_files.index(cfg['pha']) if cfg['pha'] in pha_files else 0,
            key=f"pha_{idx}"
        )
    else:
        st.sidebar.warning("No .pha files found")

    if bak_files:
        cfg['bak'] = st.sidebar.selectbox(
            f"BAK file #{idx+1}",
            options=bak_files,
            index=bak_files.index(cfg['bak']) if cfg['bak'] in bak_files else 0,
            key=f"bak_{idx}"
        )
    else:
        st.sidebar.warning("No .bak files found")

    if rsp_files:
        cfg['rsp'] = st.sidebar.selectbox(
            f"RSP file #{idx+1}",
            options=rsp_files,
            index=rsp_files.index(cfg['rsp']) if cfg['rsp'] in rsp_files else 0,
            key=f"rsp_{idx}"
        )
    else:
        st.sidebar.warning("No .rsp files found")
        
        
    cfg['energy_range'] = st.sidebar.text_input(
        label=f"Energy range #{idx+1} (keV)",
        value=cfg.get('energy_range', '**-8.0 900.0-**'),
        key=f"energy_{idx}", 
        help='use xspec command'
    )
    
    cfg['statistic'] = st.sidebar.text_input(
        label=f"statistic #{idx+1}",
        value=cfg.get('statistic', 'pgstat'),
        key=f"statistic_{idx}", 
        help='use xspec command'
    )        
    
    if st.sidebar.button(f"Delete Dataset {idx+1}", key=f"del_{idx}"):
        st.session_state['spectra'].pop(idx)
        st.experimental_rerun()

        
all_src = []
for cfg in st.session_state['spectra']:
    src = {
        'det': cfg['pha'],
        'pha_file': cfg['pha'],
        'bak_file': cfg['bak'],
        'rsp_file': cfg['rsp'],
        'arf_file': None,
        'ign_range': cfg['energy_range'],
        'statistic': cfg['statistic'],
    }
    all_src.append(src)
    
st.sidebar.markdown("### Model Settings")

model_name = st.sidebar.text_input(
label="Model Name",
value=st.session_state.get('model_name', "cutoff"),
key="model_name",
help="use XSPEC model, e.g. po, cutoff, grbm")

model_params = st.sidebar.text_input(
    label="Model Parameters",
    value=st.session_state.get('model_params', None),
    key="model_params",
    help="Enter parameters as a XSPEC string (e.g. 1:1 0.01 0 0 2 2;2:30 0.01 10 10 200 500, use ';' split)"
)

if model_params:
    params_set = parse_intervals(model_params)
else:
    params_set = None
    
fit_button = st.sidebar.button('Fit')
######################################################################################################
main_func = quick_GBM(user_path=user_path, fit_path=File_for_Fit, result_path=Result_path,
                      time_str=str_t0,
                      coordinates_str=str_coordinates,
                      pre_orbit=pre_orbit,
                      check_spechist_pre_orbit=check_spechist_pre_orbit,
                      det_select=det_select,
                      tbin=tbin,
                      lc_time_range=lc_time_range,
                      energy_range_NaI=energy_range1,
                      energy_range_BGO=energy_range2, 
                      subtract_bg=subtract_bg,
                      background1=background1, 
                      background2=background2,
                      snr_enabled=snr_enabled,
                      snr_threshold=snr_threshold,
                      gen_spec=gen_spec,
                      spectrum_time_range=spectrum_time_range,
                      Build_Files_botton=Build_Files_botton,
                      selected_poshist=selected_poshist,
                      RSP_time=RSP_time,
                      all_src=all_src, 
                      model_name=model_name,
                      params_set=params_set, 
                      )
        

tab1, tab2, tab3 = st.tabs(["Position map", "Lightcurves & File produce", "Spectrum", ])
with tab1:
    #-- Position map
    if 'position_run' not in st.session_state:
        st.session_state.position_run = False
    
    if check_pos_botton:
        main_func.check_position()
    
    if st.session_state.position_run:
        st.subheader('Position map')
        st.write(f'Poshist file: {st.session_state.poshist_filename}')
        st.write(f'MET: {st.session_state.met_time_reCal:.2f}')
        st.write('Visible:', st.session_state.visible)
        st.write('Nearest 2 NaI detectors and their angles to the source:')
        for det, ang in st.session_state.nearest_three:
            st.write(f"• {det}: {ang:.1f}°")    
        st.image(st.session_state.mappng, use_column_width=True)
        st.image(st.session_state.orbitpng, use_column_width=True)
        

with tab2:
    if 'plotlc_run' not in st.session_state:
        st.session_state.plotlc_run = False
        st.session_state.plotlc_figs = []
    
    if Build_Files_botton:
        plot_lc_botton = True
        
    if plot_lc_botton:
        main_func.plot_LC_build_file()

    if st.session_state.plotlc_run:
        st.subheader('Lightcurves')
        for fig in st.session_state.plotlc_figs:
            st.pyplot(fig)


with tab3:
    if st.session_state['spectra']:
        st.markdown("## Current Spectral Datasets")
        for i, cfg in enumerate(st.session_state['spectra'], start=1):
            st.write(
                f"- **Dataset {i}:** "
                f"PHA=`{cfg['pha']}`, "
                f"BAK=`{cfg['bak']}`, "
                f"RSP=`{cfg['rsp']}`, "
                f"Energy range: `{cfg['energy_range']}`,"
                f"statistic: `{cfg['statistic']}`,"
            )    
    
    if fit_button:
        main_func.spectrum_fit()
    
    st.subheader("Fit Results")
    txt_file = f"{model_name}_fitresult.txt"
    txt_path = os.path.join(Result_path, txt_file)
    if os.path.exists(txt_path):
        with open(txt_path) as f:
            st.text(f.read())
    else:
        st.warning(f"No result file found: {txt_file}")
        
    st.subheader("Spectrum plots")
    for suffix in ["ldata", "euf", "eeuf", "corner"]:
        img_file = f"{model_name}_{suffix}.png"
        img_path = os.path.join(Result_path, img_file)
        if os.path.exists(img_path):
            st.image(img_path, caption=img_file, use_column_width=True)
        else:
            st.warning(f"Missing plot: {img_file}")
    
    files = sorted(os.listdir(Result_path))
    
    if not files:
        st.sidebar.warning("No files available for download.")
        selected_file = None
    else:
        options = [None] + files
        selected_file = st.sidebar.selectbox(
            "Select a file to download",
            options=options,
            index=0,
            help="Choose a file from the directory (or leave blank)"
        )
    
    if selected_file:
        file_path = os.path.join(Result_path, selected_file)
        with open(file_path, "rb") as f:
            data = f.read()
        st.sidebar.download_button(
            label="Download",
            data=data,
            file_name=selected_file,
            mime="application/octet-stream"
        )