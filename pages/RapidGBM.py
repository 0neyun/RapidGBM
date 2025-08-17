import streamlit as st
import glob, os
import matplotlib.pyplot as plt
import sys
sys.path.append('./base_func')
from RapidGBM_base_func import * 

#try:
    #if st.session_state["authentication_status"]:
        #pass
#except:
    #st.switch_page("login.py")
    
if "name" not in st.session_state:
    st.session_state["name"] = "anonymous"
user_path = 'User_fold/'+ st.session_state["name"]+'_RapidGBM'
mkdir(user_path)
File_for_Fit =  f'{user_path}/0_File_for_Fit'
Result_path = f'{user_path}/1_Result_path'

back_btn  = st.sidebar.button("Back menu", key="back_menu")
if back_btn:
    preserve = {'authentication_status', 'name'}
    for k in list(st.session_state.keys()):
        if k not in preserve:
            del st.session_state[k]    
    st.switch_page("menu.py")
    
  
if st.sidebar.button("🗑 CLEAR CACHE and USER FOLDER"):
    import shutil  
    st.cache_data.clear()
    for entry in os.listdir(user_path):
        entry_path = os.path.join(user_path, entry)
        try:
            if os.path.isfile(entry_path) or os.path.islink(entry_path):
                os.remove(entry_path)        
            elif os.path.isdir(entry_path):
                shutil.rmtree(entry_path)     
        except Exception as e:
            st.sidebar.error(f"Error deleting {entry}: {e}")
    #st.sidebar.success("All files and folders under user_path have been deleted.")    
    st.sidebar.success("cleared.")

#clear_btn = st.sidebar.button("🗑 CLEAR USER FOLDER",key="clear_all",)

#
#if clear_btn:


mkdir(File_for_Fit)
mkdir(Result_path)

st.title('RapidGBM')
# Use columns for side-by-side layout
col1, col2 = st.columns([2, 4])  # Adjust width ratio as needed

with col1:
    st.image("pages/rapidgbm_logo.png", width=200)

with col2:
    st.write(f"Hi {st.session_state['name']}.")
    st.markdown("""
    
    
 
    - Use the menu at left to input parameters  
    - Your plots will appear below  
    - Code available: [https://github.com/0neyun/RapidGBM](https://github.com/0neyun/RapidGBM)
    """)
    
######################################################################################################
check_pos_botton = False
plot_lc_botton = False

######################################################################################################
# ===== Form 1: Set time and location =====
with st.sidebar.form("form_time_location", clear_on_submit=False):
    st.header("1. Set time and location")
    
    _time_default = st.session_state.get('time_str', '2024-06-17T12:19:13.00')
    _coord_default = st.session_state.get('coordinates_str', '285.03,-22.56')

    str_t0_input = st.text_input('UTC Time', _time_default, key="time_str")
    str_coordinates_input = st.text_input('Input source coordinates (RA, DEC in degree)',
                                          _coord_default, key="coordinates_str")

    pre_orbit = 30
    check_spechist_pre_orbit_input = st.checkbox(
        'Use the earlier poshist file (30 orbits prior) if the current one is missing',
        value=st.session_state.get('check_spechist_pre_orbit', False),
        key="check_spechist_pre_orbit",
    )
    
    # check position
    check_pos_botton = st.form_submit_button('check position', use_container_width=True)


str_t0 = st.session_state["time_str"]
str_coordinates = st.session_state["coordinates_str"]
check_spechist_pre_orbit = st.session_state["check_spechist_pre_orbit"]

######################################################################################################
# ===== Form 2: Set plot for lightcurves =====
with st.sidebar.form("form_lightcurves", clear_on_submit=False):
    st.header("2. Set plot for lightcurves")

    # 1) det
    det_select = st.multiselect(
        label="### Which detector(s)?",
        options=det_list,
        default=st.session_state.get('det_select', []),
        key="det_select",
    )

    # 2) NaI 
    with st.expander("NaI Energy Band (default: 8-900 keV)", expanded=False):
        nai_col1, nai_col2 = st.columns(2)
        nai_e_min = nai_col1.number_input(
            "Min", min_value=0.0, value=float(st.session_state.get('nai_e_min', 8.0)),
            step=1.0, format="%.1f", key="nai_e_min"
        )
        nai_e_max = nai_col2.number_input(
            "Max", min_value=0.0, value=float(st.session_state.get('nai_e_max', 900.0)),
            step=10.0, format="%.1f", key="nai_e_max"
        )
    energy_range1 = (st.session_state["nai_e_min"], st.session_state["nai_e_max"])

    # 3) BGO 
    with st.expander("BGO Energy Band (default: 200-40000 keV)", expanded=False):
        bgo_col1, bgo_col2 = st.columns(2)
        bgo_e_min = bgo_col1.number_input(
            "Min", min_value=0.0, value=float(st.session_state.get('bgo_e_min', 200.0)),
            step=10.0, format="%.1f", key="bgo_e_min"
        )
        bgo_e_max = bgo_col2.number_input(
            "Max", min_value=0.0, value=float(st.session_state.get('bgo_e_max', 40000.0)),
            step=100.0, format="%.1f", key="bgo_e_max"
        )
    energy_range2 = (st.session_state["bgo_e_min"], st.session_state["bgo_e_max"])

    # 4) LC set
    st.markdown("### Set Lightcurves Interval")
    tbin = st.number_input(
        label="Time bin size (s)",
        min_value=1e-4, max_value=100.0,
        value=float(st.session_state.get('tbin', 2.0)),
        step=0.01, format="%.3f", key="tbin"
    )
    t_col1, t_col2 = st.columns(2)
    t_start = t_col1.number_input(
        label="Time start (s)",
        value=float(st.session_state.get('t_start', -120.0)),
        step=1., format="%.1f", key="t_start"
    )
    t_end = t_col2.number_input(
        label="Time end (s)",
        value=float(st.session_state.get('t_end', 470.0)),
        step=1., format="%.1f", key="t_end"
    )
    if t_end <= t_start:
        st.error("Error: End time must be greater than start time.")
    lc_time_range = (st.session_state["t_start"], st.session_state["t_end"])

    # ===== 在 sub_bk 之前增加一个“画光变”按钮（label 唯一）=====
    _do_plot_early = st.form_submit_button('plot')

    # 5) sub_bk
    col1_b, col2_b = st.columns([4, 1])
    col1_b.markdown("### Subtract background & Build files")
    subtract_bg = col2_b.checkbox("", value=st.session_state.get("subtract_bg", False), key="subtract_bg")

    background1 = None
    background2 = None
    snr_enabled = True
    snr_threshold = None
    spectrum_time_range = None
    gen_spec = False
    RSP_time = None
    selected_poshist = None

    # bk1
    b1_col1, b1_col2 = st.columns(2)
    b1_start = b1_col1.number_input(
        "Bkg1 Start", value=float(st.session_state.get("b1_start", -100.0)),
        key="b1_start", format="%.1f", step=1.0
    )
    b1_end = b1_col2.number_input(
        "Bkg1 End", value=float(st.session_state.get("b1_end", -50.0)),
        key="b1_end", format="%.1f", step=1.0
    )
    if b1_end <= b1_start:
        st.error("Error: Bkg1 end must be > start.")
    background1 = (st.session_state["b1_start"], st.session_state["b1_end"])

    # bk2
    b2_col1, b2_col2 = st.columns(2)
    b2_start = b2_col1.number_input(
        "Bkg2 Start", value=float(st.session_state.get("b2_start", 400.0)),
        key="b2_start", format="%.1f", step=1.0
    )
    b2_end = b2_col2.number_input(
        "Bkg2 End", value=float(st.session_state.get("b2_end", 450.0)),
        key="b2_end", format="%.1f", step=1.0
    )
    if b2_end <= b2_start:
        st.error("Error: Bkg2 end must be > start.")
    background2 = (st.session_state["b2_start"], st.session_state["b2_end"])

    # Polynomial order (degree) — default 1
    poly_order = st.number_input(
        "Polynomial order",
        min_value=1,
        max_value=5,
        value=st.session_state.get("poly_order", 1),
        step=1,
        key="poly_order",
        help="Set the polynomial order for background fitting (1–5)."
    )

    # SNR
    snr_enabled = True
    snr_threshold = st.number_input(
        "SNR threshold", min_value=0.0,
        value=float(st.session_state.get("snr_threshold", 3.0)),
        step=1., format="%.1f", key="snr_threshold"
    )

    st.markdown("### Spectrum Time Range (s)")
    spec_col1, spec_col2 = st.columns(2)
    spec_start = spec_col1.number_input(
        "Spec Start", value=float(st.session_state.get("spec_start", 227.)),
        key="spec_start", format="%.2f", step=0.1
    )
    spec_end = spec_col2.number_input(
        "Spec End", value=float(st.session_state.get("spec_end", 347.)),
        key="spec_end", format="%.2f", step=0.1
    )
    if spec_end <= spec_start:
        st.error("Error: Spectrum end time must be greater than start time.")
    spectrum_time_range = (st.session_state["spec_start"], st.session_state["spec_end"])

    poshist_files = [f for f in os.listdir(user_path) if 'glg_poshist_all' in f]
    default = st.session_state.get('poshist_filename', None)
    if default in poshist_files:
        default_index = poshist_files.index(default)
    else:
        default_index = 0 if poshist_files else 0
    selected_poshist = st.selectbox(
        label="Select poshist file",
        options=poshist_files if poshist_files else ["(no poshist files)"],
        index=default_index,
        key="selected_poshist"
    )

    value_temp = st.session_state.get("met_time_reCal", None)
    RSP_time = st.number_input(
        "RSP_MET",
        value=float(value_temp) if value_temp is not None else 0.0,
        key="RSP_time", format="%.2f"
    )

    gen_spec = st.checkbox(
        "Generate spectrum",
        value=st.session_state.get("gen_spectrum", False),
        key="gen_spectrum"
    )

    # ===== 底部提交按钮（label 保证唯一）=====
    _do_plot = st.form_submit_button('plot all')
    _do_build = st.form_submit_button('Build Files (PHA/BAK/RSP) for Fit')

# 表单外聚合动作
plot_lc_botton = _do_plot_early or _do_plot or _do_build
Build_Files_botton = _do_build

subtract_bg = st.session_state.get("subtract_bg", False)
if subtract_bg:
    background1 = (st.session_state["b1_start"], st.session_state["b1_end"])
    background2 = (st.session_state["b2_start"], st.session_state["b2_end"])
else:
    background1 = None
    background2 = None

snr_enabled = st.session_state.get("snr_enabled", True)
snr_threshold = st.session_state.get("snr_threshold", None)
gen_spec = st.session_state.get("gen_spectrum", False)
spectrum_time_range = (st.session_state.get("spec_start"), st.session_state.get("spec_end")) if gen_spec else None
selected_poshist = st.session_state.get("selected_poshist", None) if gen_spec else None
RSP_time = st.session_state.get("RSP_time", None) if gen_spec else None


######################################################################################################
# ===== Form 3: Set Spectral Analysis =====
st.sidebar.header("3. Spectral Analysis")

with st.sidebar.form("spectral_analysis_form", clear_on_submit=False):
    st.markdown("### Upload Files")

    # Upload
    uploaded_files = st.file_uploader(
        label="",
        type=None,                   
        accept_multiple_files=True
    )

    # save file
    if uploaded_files:
        for uploaded_file in uploaded_files:
            save_path = os.path.join(File_for_Fit, uploaded_file.name)
            with open(save_path, "wb") as f:
                f.write(uploaded_file.getbuffer())
        st.sidebar.success(f"Saved {len(uploaded_files)} files.")

    # spec list
    if 'spectra' not in st.session_state:
        st.session_state['spectra'] = []

    # add
    add_dataset_button = st.form_submit_button("Add dataset")
    if add_dataset_button:
        st.session_state['spectra'].append({'pha': None, 'bak': None, 'rsp': None})

    pha_paths = sorted(glob(os.path.join(File_for_Fit, "*.pha*"))+glob(os.path.join(File_for_Fit, "*.pi*")))
    bak_paths = sorted(glob(os.path.join(File_for_Fit, "*.bak")))
    rsp_paths = sorted(glob(os.path.join(File_for_Fit, "*.rsp")))

    pha_files = [os.path.basename(p) for p in pha_paths] + ['None']
    bak_files = [os.path.basename(p) for p in bak_paths] + ['None']
    rsp_files = [os.path.basename(p) for p in rsp_paths] + ['None']

    for idx, cfg in enumerate(st.session_state['spectra']):
        st.markdown(f"**Dataset {idx+1}**")
        if pha_files:
            cfg['pha'] = st.selectbox(
                f"PHA file #{idx+1}",
                options=pha_files,
                index=pha_files.index(cfg['pha']) if cfg['pha'] in pha_files else 0,
                key=f"pha_{idx}"
            )
        else:
            st.sidebar.warning("No .pha files found")

        if bak_files:
            cfg['bak'] = st.selectbox(
                f"BAK file #{idx+1}",
                options=bak_files,
                index=bak_files.index(cfg['bak']) if cfg['bak'] in bak_files else 0,
                key=f"bak_{idx}"
            )
        else:
            st.sidebar.warning("No .bak files found")

        if rsp_files:
            cfg['rsp'] = st.selectbox(
                f"RSP file #{idx+1}",
                options=rsp_files,
                index=rsp_files.index(cfg['rsp']) if cfg['rsp'] in rsp_files else 0,
                key=f"rsp_{idx}"
            )
        else:
            st.sidebar.warning("No .rsp files found")

        # Energy range input
        cfg['energy_range'] = st.text_input(
            label=f"Energy range #{idx+1} (keV)",
            value=cfg.get('energy_range', '**-8.0 900.0-**'),
            key=f"energy_{idx}", 
            help='use xspec command'
        )
        
        # Statistic input
        cfg['statistic'] = st.text_input(
            label=f"statistic #{idx+1}",
            value=cfg.get('statistic', 'pgstat'),
            key=f"statistic_{idx}", 
            help='use xspec command'
        )        

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

    # ==== model ====
    st.markdown("### Model Settings")
    model_name = st.text_input(
        label="Model Name",
        value=st.session_state.get('model_name', "cutoff"),
        key="model_name",
        help="use XSPEC model, e.g. po, cutoff, grbm"
    )

    model_params = st.text_input(
        label="Model Parameters",
        value=st.session_state.get('model_params', None),
        key="model_params",
        help="Enter parameters as a XSPEC string (e.g. 1:1 0.01 0 0 2 2;2:30 0.01 10 10 200 500, use ';' split)"
    )

    if model_params:
        params_set = parse_intervals(model_params)
    else:
        params_set = None

    # ==== fit ====
    fit_button = st.form_submit_button('Fit')
    
    # ==== delect ====
    st.markdown("### Delete Dataset")
    
    # Default selection is the first dataset or 'None'
    delete_idx = st.selectbox(
        "Select Dataset to delete", 
        options=[f"Dataset {i+1}" for i in range(len(st.session_state['spectra']))] + ["None"],
        index=0,  # 默认选中第一个数据集
        key="delete_idx"
    )

    if delete_idx != "None":
        delete_button = st.form_submit_button("Delete Dataset")
        if delete_button:
            idx_to_delete = int(delete_idx.split()[-1]) - 1 
            del st.session_state['spectra'][idx_to_delete]  
            st.experimental_rerun()   


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
                      poly_order=poly_order
                      )
        

tab1, tab2, tab3 = st.tabs(["1.Position map", "2.Lightcurves & File produce", "3.Spectrum", ])
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
        st.write('Nearest three NaI detectors and their angles to the source:')
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
                f"statistic: `{cfg['statistic']}`"
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
    
    #files = sorted(os.listdir(Result_path))
    
    #if not files:
        #st.sidebar.warning("No files available for download.")
        #selected_file = None
    #else:
        #options = [None] + files
        #selected_file = st.sidebar.selectbox(
            #"Select a file to download",
            #options=options,
            #index=0,
            #help="Choose a file from the directory (or leave blank)"
        #)
    
    #if selected_file:
        #file_path = os.path.join(Result_path, selected_file)
        #with open(file_path, "rb") as f:
            #data = f.read()
        #st.sidebar.download_button(
            #label="Download",
            #data=data,
            #file_name=selected_file,
            #mime="application/octet-stream"
        #)
