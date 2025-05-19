import glob, os
import math

#import __init__
from astropy.io import fits as pf
import streamlit as st
from gbm.plot import EarthPlot
from gbm.data import TTE
from gbm.binning.unbinned import bin_by_time
from gbm.background import BackgroundFitter
from gbm.background.binned import Polynomial
from gbm.plot import SkyPlot
from gbm.data import PosHist
from gbm.time import Met
import matplotlib.pyplot as plt

from requests_html import HTMLSession
import wget
import numpy as np
from glob import glob


def linesplit_2(line):
    line = line.split(',')
    return float(line[0]), float(line[1])

def mkdir(path):
    folder = os.path.exists(path)
    if folder:
        pass
    else:
        os.makedirs(path)
        
def file_exists(dir_path: str, filename: str) -> bool:
    return os.path.isfile(os.path.join(dir_path, filename))
    
def read_poshist(pos_file, verbose = True):
    '''
    Extract Quaternions, Position, Time & Geo Coordinates from file.
    Poshist files for days prior to March 2009 either have the spacecraft lat &
    lon set to zero or the fields are missing altogheter. This should be caught
    by the try except block in place.
    '''
    dtorad=180./math.acos(-1.)
    data=pf.getdata(pos_file,ext=1)
    nt=np.size(data)
    sc_time=data.SCLK_UTC
    sc_quat=np.zeros((nt,4),float)
    sc_pos=np.zeros((nt,3),float)
    sc_coords=np.zeros((nt,2),float)
    try:
        sc_coords[:,0]=data.SC_LON
        sc_coords[:,1]=data.SC_LAT
    except:
        if verbose:
            mes = ''
            mes += '*** No geographical coordinates available '
            mes += 'for this file: %s' %pos_file
            print (mes)
            
    sc_quat[:,0]=data.QSJ_1
    sc_quat[:,1]=data.QSJ_2
    sc_quat[:,2]=data.QSJ_3
    sc_quat[:,3]=data.QSJ_4
    sc_pos[:,0]=data.POS_X
    sc_pos[:,1]=data.POS_Y
    sc_pos[:,2]=data.POS_Z
    return sc_time,sc_pos,sc_quat,sc_coords
    
def calc_period(sc_pos):
    '''
    Calculate orbital period of Fermi using the position of the spacecraft and
    assuming circular motion.   
    '''
    G = 6.67428e-11    # m^3 kg^-1 s^-2
    M = 5.9722e24      # kg Mass Earth  
    r = (np.sum(sc_pos**2.,1))**(1/2.)
    r_avg = np.average(r)
    r_cubed = (r_avg)**3.
    factor = r_cubed/(G*M)
    period = 2. * np.pi * np.sqrt(factor)
    return period

@st.cache_data(max_entries=5)
def request_ftp(yyyy, mm, dd):
    url = 'https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/daily/%s/%s/%s/current/' %(yyyy, mm, dd)
    session = HTMLSession()
    r = session.get(url)
    file_list = np.array(np.sort(list(r.html.absolute_links)))[5:]
    return url, file_list

def download_with_progress(url, dest_path): 
    filename = os.path.basename(dest_path)
    
    progress_text = st.empty()
    progress_bar  = st.empty()

    def wget_bar(current, total, width=0):
        if total:
            frac = current / total
            progress_bar.progress(min(frac, 1.0))
            progress_text.text(f"Downloading {filename}: {current}/{total} bytes ({frac*100:.1f}%)")
        else:
            progress_text.text(f"Downloading {filename}: {current} bytes")

    progress_text.text(f"Starting download of {filename} ...")
    wget.download(url, out=dest_path, bar=wget_bar)

    progress_text.text(f"Finished downloading {filename}")
    progress_bar.empty()

def check_ttefile(det_list, user_path, file_list, ymd_h):
    tte_fns = []
    for det_name in det_list:
        links = [f for f in file_list if f"_{det_name}_{ymd_h}" in f]
        if not links:
            st.warning(f"No TTE file found for {det_name}")
            request_ftp.clear()
            continue
        url = links[0]
        fname = url.split('/')[-1]
        out_path = os.path.join(user_path, fname)
        
        if not os.path.exists(out_path):
            #if download_button:
            st.info(f"Downloading {fname} ...")
            download_with_progress(url, out_path)
            st.success(f"{fname} downloaded successfully.")
            #else:
                #st.warning(f"{fname} not found. Click the download button to fetch it.")
                #continue
        else:
            st.info(f"{fname} already exists, skipping download.")

        tte_fns.append(fname)

    return np.array(tte_fns)

def parse_intervals(s: str):
    result = []
    for segment in s.split(';'):
        seg = segment.strip()
        if not seg:
            continue
        key_str, val_str = seg.split(':', 1)
        key = int(key_str)
        val = val_str.strip()
        result.append({ key: val })
    return result

def get_lc(ttefile_path,tbin,t1,t2,E1,E2):
    tte = TTE.open(ttefile_path)
    phaii = tte.to_phaii(bin_by_time, tbin, time_range=(t1, t2))
    lc_data = phaii.to_lightcurve(energy_range=(E1,E2))
    y = lc_data.counts
    yerr = lc_data.count_uncertainty
    x = lc_data.centroids
    xerr = lc_data.centroids-lc_data.hi_edges
    return x[:-1],xerr[:-1],y[:-1],yerr[:-1],phaii,tte
    
det_list = ['n0', 'n1', 'n2', 'n3', 'n4', 'n5', 'n6', 'n7', 'n8', 'n9', 'na', 'nb', 'b0', 'b1']
det_list_NaI = ['n0', 'n1', 'n2', 'n3', 'n4', 'n5', 'n6', 'n7', 'n8', 'n9', 'na', 'nb']
det_list_BGO = ['b0','b1']
class quick_GBM:
    def __init__(self,user_path, fit_path, result_path, 
                 time_str,
                 coordinates_str,
                 pre_orbit,
                 check_spechist_pre_orbit,
                 det_select,
                 tbin,
                 lc_time_range,
                 energy_range_NaI, 
                 energy_range_BGO, 
                 subtract_bg,
                 background1,
                 background2,
                 snr_enabled,
                 snr_threshold,
                 gen_spec,
                 spectrum_time_range, 
                 Build_Files_botton,
                 selected_poshist,
                 RSP_time,
                 all_src,
                 model_name,
                 params_set,
                 ):
        #set path
        self.user_path = user_path
        self.fit_path = fit_path
        self.result_path = result_path
        
        #-- set time and coordinates
        ###################################################################################
        met_time = Met.from_iso(time_str).met
        ymd = Met.from_iso(time_str).ymd
        ymd_h = Met.from_iso(time_str).ymd_h
        self.yyyy = '20%s' %(ymd[0:2])
        self.mm = ymd[2:4]
        self.dd = ymd[4:6]
        
        default_orbit_period = 5737.70910239
        met_time2 = Met(Met.from_iso(time_str).met - float(pre_orbit) *default_orbit_period)
        ymd2 = met_time2.ymd
        #print (ymd2)
        self.yyyy2 = '20%s' %(ymd2[0:2])
        self.mm2 = ymd2[2:4]
        self.dd2 = ymd2[4:6]        
        
        met_time3 = Met(Met.from_iso(time_str).met - float(pre_orbit/2) *default_orbit_period)
        ymd3 = met_time3.ymd
        #print (ymd3)
        self.yyyy3 = '20%s' %(ymd3[0:2])
        self.mm3 = ymd3[2:4]
        self.dd3 = ymd3[4:6]            
        
        RA,DEC = linesplit_2(coordinates_str)

        self.met_time = met_time
        #self.met_orbit = met_time2
        self.ra = RA
        self.dec = DEC
        ###################################################################################
        #set params for part1

        self.pre_orbit = float(pre_orbit)
        self.check_spechist_pre_orbit = check_spechist_pre_orbit
        ###################################################################################
        #set params for part2
        self.det_select = det_select
        self.tbin = tbin
        self.ymd_h = ymd_h
        self.lc_time_range = lc_time_range
        self.energy_range_NaI = energy_range_NaI
        self.energy_range_BGO = energy_range_BGO
        self.subtract_bg = subtract_bg
        self.background1 = background1
        self.background2 = background2
        self.snr_enabled = snr_enabled
        self.snr_threshold = snr_threshold
        self.gen_spec = gen_spec
        self.spectrum_time_range = spectrum_time_range
        self.Build_Files_botton = Build_Files_botton
        self.selected_poshist = selected_poshist
        self.RSP_time = RSP_time
        self.model_name = model_name
        self.params_set = params_set
        self.all_src = all_src
        
    #@st.cache_data(max_entries=5)
    def check_position(self, ):
        if self.check_spechist_pre_orbit:
            st.warning('Notice: Using the earlier poshist file, the pointing may have deviations!')
            self.ftp_url_orbit, self.ftp_file_list_orbit = request_ftp(self.yyyy2, self.mm2, self.dd2)
            self.poshist_link_orbit = [f for f in self.ftp_file_list_orbit if 'glg_poshist_all_' in f]
            self.poshist_filename_orbit = self.poshist_link_orbit[0].split('/')[-1]        
            poshist_filename = self.poshist_filename_orbit
            poshist_link = self.poshist_link_orbit
            
        else:
            self.ftp_url, self.ftp_file_list = request_ftp(self.yyyy, self.mm, self.dd)
            self.poshist_link = [f for f in self.ftp_file_list if 'glg_poshist_all_' in f]
            self.poshist_filename = self.poshist_link[0].split('/')[-1]
            poshist_filename = self.poshist_filename
            poshist_link = self.poshist_link
        
        if poshist_link:
            if os.path.exists(self.user_path+'/'+poshist_filename):
                pass
            else:
                download_with_progress(poshist_link[0], self.user_path+'/'+poshist_filename)
            
            if self.check_spechist_pre_orbit:
                ftp_url_orbit3, ftp_file_list_orbit3 = request_ftp(self.yyyy3, self.mm3, self.dd3)
                poshist_link_orbit3 = [f for f in ftp_file_list_orbit3 if 'glg_poshist_all_' in f]
                poshist_filename_orbit3 = poshist_link_orbit3[0].split('/')[-1]        
                poshist_filename3 = poshist_filename_orbit3
                poshist_link3 = poshist_link_orbit3
                if os.path.exists(self.user_path+'/'+poshist_filename3):
                    pass
                else:                
                    download_with_progress(poshist_link3[0], self.user_path+'/'+poshist_filename3)
                
                _, sc_pos1, __, ___ = read_poshist(self.user_path+'/'+poshist_filename)
                _, sc_pos2, __, ___ = read_poshist(self.user_path+'/'+poshist_filename3)
                
                new_period = calc_period(np.concatenate((sc_pos1, sc_pos2), axis=0))
                #new_period = calc_period(sc_pos1)
                #print ('new_period: ', new_period)
                met_time_reCal = self.met_time - (self.pre_orbit) *new_period
        
            else:
                met_time_reCal = self.met_time
                quick_GBM
            poshist = PosHist.open(self.user_path+'/'+poshist_filename)
            skyplot = SkyPlot()
            skyplot.add_poshist(poshist, trigtime=met_time_reCal)
            skyplot.plot_point(self.ra, self.dec, color='r', s=200, alpha=1, marker='*')
            mappng = f'{self.user_path}/{self.ra}_{self.dec}_{met_time_reCal}_position_map_{self.check_spechist_pre_orbit}.png'
            plt.savefig(mappng)
            plt.close()
            
            t0 = float(met_time_reCal)
            #lat = poshist.get_latitude(float(t0))
            #lon = poshist.get_longitude(float(t0)) # note: East longitude
            #alt = poshist.get_altitude(float(t0))    
            
            earthplot = EarthPlot()
            earthplot.add_poshist(poshist, trigtime=t0, time_range=(t0-1200, t0+1200))
            orbitpng = f'{self.user_path}/{self.ra}_{self.dec}_{met_time_reCal}_location_map_{self.check_spechist_pre_orbit}.png'
            plt.savefig(orbitpng)
            plt.close()    
            
            visible = poshist.location_visible(self.ra, self.dec, met_time_reCal)
            
            angle_list = []
            for det in det_list_NaI:
                angle = poshist.detector_angle(self.ra, self.dec, det, met_time_reCal)
                angle_list.append((det, angle))
        
            nearest_three = sorted(angle_list, key=lambda x: x[1])[:3]
            
            st.session_state.position_run = True
            st.session_state.visible = visible[0]
            st.session_state.nearest_three = nearest_three
            st.session_state.mappng = mappng
            st.session_state.orbitpng = orbitpng
            st.session_state.osv = self.check_spechist_pre_orbit
            st.session_state.poshist_filename = poshist_filename
            st.session_state.met_time_reCal = met_time_reCal            
            
        else:
            st.warning('The required poshist file not exist.')
            st.session_state.position_run = False
            request_ftp.clear()
            

    #@st.cache_data
    def plot_LC_build_file(self, ):
        if not self.det_select:
            st.warning('No detector been selected')
        else:
            self.ftp_url, self.ftp_file_list = request_ftp(self.yyyy, self.mm, self.dd)
            self.ftp_url_orbit, self.ftp_file_list_orbit = request_ftp(self.yyyy2, self.mm2, self.dd2)
            self.tte_links = [f for f in self.ftp_file_list if self.ymd_h in f]
            if not self.tte_links:
                st.warning('The required TTE file not exist.')        
            else:
                t1, t2 = (self.lc_time_range)
                progress_text = st.empty()
                progress_bar  = st.empty()
                
                tte_filename_list = check_ttefile(self.det_select,
                                                  user_path=self.user_path,
                                                  file_list=self.ftp_file_list,
                                                  ymd_h=self.ymd_h)
                #if self.subtract_bg:
                    #t1_ = self.background1[0]
                    #t2_ = self.background2[1]
                #else:
                t1_ = t1
                t2_ = t2        
                
                figs = []
                for fn in tte_filename_list:
                    if ('_b0_' in fn) or ('_b1_' in fn):
                        E1,E2 = self.energy_range_BGO
                    else:
                        E1,E2 = self.energy_range_NaI
        
                    fig = plt.figure(figsize=(12, 5))
                    x,xerr,y,yerr,_,tte = get_lc(self.user_path+'/'+fn, self.tbin, t1+self.met_time, t2+self.met_time, E1, E2)
                    x = x - self.met_time
                    plt.axvline(x=0,ls="-",c="red")
                    plt.errorbar(x=x, y=y, yerr=yerr, fmt='.', color='k', markersize=0.1)
                    plt.step(x, y, where='mid', color='k', label=f'{E1}-{E2} keV')
                    if self.subtract_bg:
                        _,__,___,____,phaii,_____ = get_lc(self.user_path+'/'+fn, 2.048, t1_+self.met_time, t2_+self.met_time, E1, E2)
                        bak_interval = [(self.background1[0] + self.met_time, self.background1[1] + self.met_time),\
                                        (self.background2[0] + self.met_time, self.background2[1] + self.met_time)]
                        backfitter = BackgroundFitter.from_phaii(phaii, Polynomial, time_ranges=bak_interval)
                        backfitter.fit(order=1)
                        bkgd = backfitter.interpolate_bins(phaii.data.tstart, phaii.data.tstop)    
                        lc_bkgd = bkgd.integrate_energy(*(E1, E2))        
                        bak_xx = lc_bkgd.time_centroids
                        bak_yy = lc_bkgd.rates*(x[1]-x[0])
                        plt.plot(bak_xx-self.met_time,bak_yy, lw=3, label='background')
                        plt.axvspan(
                            self.background1[0],      
                            self.background1[1],        
                            color='blue',
                            alpha=0.2,
                        )
                        plt.axvspan(
                            self.background2[0],      
                            self.background2[1],        
                            color='blue',
                            alpha=0.2,
                        )                         
                        
                        if self.snr_enabled:
                            snr_ = bak_yy + self.snr_threshold * bak_yy ** 0.5
                            plt.plot(bak_xx-self.met_time,snr_, label=f'SNR={self.snr_threshold}', lw=3)
                        if self.gen_spec:
                            src_t1, src_t2 = self.spectrum_time_range
                            plt.axvspan(
                                src_t1,      
                                src_t2,        
                                color='green',
                                alpha=0.2,
                            )
                            plt.axvline(
                                src_t1,
                                color='green',
                                linestyle='--',
                                linewidth=1
                            )
                            plt.axvline(
                                src_t2,
                                color='green',
                                linestyle='--',
                                linewidth=1
                            )
                            if self.Build_Files_botton:
                                gen_file_path = f'{self.user_path}/{fn[:10]}_{src_t1:.2f}_{src_t2:.2f}_datafile'
                                mkdir(gen_file_path)
                                
                                if not file_exists(gen_file_path, self.selected_poshist):
                                    os.system(f'cp {self.user_path}/{self.selected_poshist} {gen_file_path}/{self.selected_poshist}')
                                
                                pha_name = f'{fn[:10]}_{src_t1:.2f}_{src_t2:.2f}.pha'
                                if not file_exists(gen_file_path, pha_name):
                                    srcpha = tte.to_pha(time_ranges=(src_t1+self.met_time, src_t2+self.met_time))
                                    srcpha.write(f'{gen_file_path}/', filename=pha_name, backfile='none')
                                    st.success(f"Written {pha_name}")
                                    if not file_exists(self.fit_path, pha_name):
                                        os.system(f'cp {gen_file_path}/{pha_name} {self.fit_path}/{pha_name}')
                                else:
                                    st.info(f"{pha_name} already exists, skipping PHA write")
                                
                                bak_name = f'{fn[:10]}_{src_t1:.2f}_{src_t2:.2f}.bak'
                                if not file_exists(gen_file_path, bak_name):
                                    bak = bkgd.to_bak(time_range=(src_t1+self.met_time, src_t2+self.met_time))
                                    bak.write(f'{gen_file_path}/', filename=bak_name)
                                    st.success(f"Written {bak_name}")
                                    if not file_exists(self.fit_path, bak_name):
                                        os.system(f'cp {gen_file_path}/{bak_name} {self.fit_path}/{bak_name}')                    
                                else:
                                    st.info(f"{bak_name} already exists, skipping BAK write")                
                                
                                #'''Bulid RSP'''
                                det = fn[8:10]
                                
                                if self.check_spechist_pre_orbit:
                                    st.warning('Notice: Using the earlier poshist file, the pointing may have deviations!')
                                    cspec_links = [f for f in self.ftp_file_list_orbit if f'glg_cspec_{det}' in f]
                                else:
                                    cspec_links = [f for f in self.ftp_file_list if f'glg_cspec_{det}' in f]
                                
                                url = cspec_links[0]
                                cspec_name = url.split('/')[-1]
                                if not file_exists(self.user_path, cspec_name):
                                    wget.download(url, self.user_path+'/')
        
                                if not file_exists(gen_file_path, cspec_name):
                                    os.system(f'cp {self.user_path}/{cspec_name} {gen_file_path}/{cspec_name}')                
                                
                                rsp_name = f"{fn[:10]}_{src_t1:.2f}_{src_t2:.2f}_{self.selected_poshist[16:22]}.rsp"
                                if not file_exists(gen_file_path, rsp_name):
                                    rsp_t1 = self.RSP_time + src_t1 - 1
                                    rsp_t2 = self.RSP_time + src_t2 + 1
                                    cmd = f"SA_GBM_RSP_Gen.pl -Ccspec -R{self.ra} -D{self.dec} -S{rsp_t1} -E{rsp_t2} ./{gen_file_path}"
                                    os.system(cmd)
                                    pattern = f"./{gen_file_path}/glg_cspec_{det}*.rsp*"
                                    matches = glob(pattern)
                                    if not matches:
                                        st.error(f"No RSP file was build")
                                    else:
                                        old_path = matches[0]
                                        new_path = os.path.join(gen_file_path, rsp_name)
                                        os.rename(old_path, new_path)
                                        st.success(f"{rsp_name} was build")
                                        if not file_exists(self.fit_path, rsp_name):
                                            os.system(f'cp {gen_file_path}/{rsp_name} {self.fit_path}/{rsp_name}')                         
                                else:
                                    st.info(f"{rsp_name} already exists, skipping RSP build")
                            
                            
                    plt.title(fn)
                    plt.xlim(x[0], x[-1])
                    #plt.ylim(np.min(y)/1.1, np.max(y)*1.1)
                    plt.xlabel('Time since input UTC [s]')
                    plt.ylabel('Counts/bin')
                    plt.legend(loc=1)
                    figs.append(fig)
                    plt.close(fig)            
                
                st.session_state.plotlc_figs = figs
                st.session_state.plotlc_run = True
                progress_text.empty()
                progress_bar.empty()        

    
    
    def spectrum_fit(self, ):
        import PyXspec_fit_MCMC
        model_set = {'name': self.model_name, 'params_set': self.params_set}

        PyXspec_fit_MCMC.Fit(self.fit_path, self.result_path, self.all_src, model_set, HPD=False)
