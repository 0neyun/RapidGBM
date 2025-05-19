import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
import bilby
import os
import streamlit as st

def model(x,a,b):
    return a+b*(x)

def PL(x,log_norm,index1):
    norm = 10**log_norm
    return norm * (x)**(index1)

#RMSE
def rmse(y_real, y_model):
    return np.sqrt(((y_real-y_model) **2).mean())

def point2_func(x1,x2,y1,y2,x):
    y = (x-x1)/(x2-x1)*(y2-y1)+y1
    return y

c = ['#1f77b4', '#ff7f0e', 'k', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf','#EE82EE']
def plot_Digram(rerun=False):
    plt.figure(figsize=(10,8))

    current_path = os.path.abspath(os.path.dirname(__file__))

    '''plot line'''
    line = np.loadtxt(f'{current_path}/fit_line.txt', usecols=(0, 1), unpack=True)
    
    #plt.plot(line[0][0:2],line[1][0:2],'--',alpha=1,c=c[0])
    #plt.plot(line[0][2:4],line[1][2:4],alpha=1,c=c[0])
    #plt.plot(line[0][4:6],line[1][4:6],'--',alpha=1,c=c[0])
    
    #plt.plot(line[0][6:8],line[1][6:8],'--',alpha=1,c='grey')
    #plt.plot(line[0][8:10],line[1][8:10],alpha=1,c='grey')
    #plt.plot(line[0][10:12],line[1][10:12],'--',alpha=1,c='grey')        
    
    alpha_line=0.5
    x_list = np.linspace(42,57,100)
    y_list = point2_func(x1=np.log10(line[0][0]), x2=np.log10(line[0][1]), y1=np.log10(line[1][0]), y2=np.log10(line[1][1]), x=x_list)
    plt.plot(10**x_list,10**y_list,'--',alpha=alpha_line,c='b')
    y_list = point2_func(x1=np.log10(line[0][2]), x2=np.log10(line[0][3]), y1=np.log10(line[1][2]), y2=np.log10(line[1][3]), x=x_list)
    plt.plot(10**x_list,10**y_list,'-',alpha=alpha_line,c='b')
    y_list = point2_func(x1=np.log10(line[0][4]), x2=np.log10(line[0][5]), y1=np.log10(line[1][4]), y2=np.log10(line[1][5]), x=x_list)
    plt.plot(10**x_list,10**y_list,'--',alpha=alpha_line,c='b')
    
    y_list = point2_func(x1=np.log10(line[0][6]), x2=np.log10(line[0][7]), y1=np.log10(line[1][6]), y2=np.log10(line[1][7]), x=x_list)
    plt.plot(10**x_list,10**y_list,'--',alpha=alpha_line,c='r')
    y_list = point2_func(x1=np.log10(line[0][8]), x2=np.log10(line[0][9]), y1=np.log10(line[1][8]), y2=np.log10(line[1][9]), x=x_list)
    plt.plot(10**x_list,10**y_list,'-',alpha=alpha_line,c='r')
    y_list = point2_func(x1=np.log10(line[0][10]), x2=np.log10(line[0][11]), y1=np.log10(line[1][10]), y2=np.log10(line[1][11]), x=x_list)
    plt.plot(10**x_list,10**y_list,'--',alpha=alpha_line,c='r')
    
    
    line = np.loadtxt(f'{current_path}/SGR_line.txt',usecols=(0,1),unpack=True)
    y_list = point2_func(x1=np.log10(line[0][0]), x2=np.log10(line[0][1]), y1=np.log10(line[1][0]), y2=np.log10(line[1][1]), x=x_list)
    plt.plot(10**x_list,10**y_list,'--',alpha=alpha_line,c='m')
    y_list = point2_func(x1=np.log10(line[0][2]), x2=np.log10(line[0][3]), y1=np.log10(line[1][2]), y2=np.log10(line[1][3]), x=x_list)
    plt.plot(10**x_list,10**y_list,'-',alpha=alpha_line,c='m')
    y_list = point2_func(x1=np.log10(line[0][4]), x2=np.log10(line[0][5]), y1=np.log10(line[1][4]), y2=np.log10(line[1][5]), x=x_list)
    plt.plot(10**x_list,10**y_list,'--',alpha=alpha_line,c='m')
    
    #plt.show()
    
    #l_data = np.loadtxt('./Amati/LGRB.txt',usecols=(0,1,2,3,4,5), unpack=True)
    #s_data = np.loadtxt('./Amati/SGRB.txt',usecols=(0,1,2,3,4,5), unpack=True)
    s_data = np.loadtxt(f'{current_path}/1.txt',usecols=(0,1,2,3,4,5,6,7), unpack=True)
    l_data = np.loadtxt(f'{current_path}/2.txt',usecols=(0,1,2,3,4,5,6,7), unpack=True)

    Eiso = s_data[2]*1e51
    Eiso_lo = s_data[3]*1e51
    Eiso_hi = s_data[4]*1e51
    
    Ep = s_data[5]
    Ep_lo = s_data[6]
    Ep_hi = s_data[7]

    plt.errorbar(x=Eiso, y=Ep, yerr=(Ep_lo,Ep_hi), xerr=(Eiso_lo,Eiso_hi), fmt='.',capsize=3, alpha=0.15,label='Tpye I GRB', c='b')
    

    ''' for long_burst '''
    Eiso = l_data[2]*1e51
    Eiso_lo = l_data[3]*1e51
    Eiso_hi = l_data[4]*1e51
    Ep = l_data[5]
    Ep_lo = l_data[6]
    Ep_hi = l_data[7]
    
    plt.errorbar(x=Eiso, y=Ep, yerr=(Ep_lo,Ep_hi), xerr=(Eiso_lo,Eiso_hi), fmt='.',capsize=3, alpha=0.15,label='Tpye II GRB',color='r')
    
    
    '''for SGR'''
    sgr_data = np.loadtxt(f'{current_path}/SGR_point.txt',usecols=(0,1,2,3), unpack=True)
    plt.errorbar(x=sgr_data[0], y=sgr_data[1], yerr=(sgr_data[3],sgr_data[2]), fmt='s',capsize=3, alpha=0.15,label='SGR GF',color='m')
        
    
    plt.ylim(1e0,4e4)
    plt.xlim(2e42,3e56)
    plt.legend(fontsize=13)
    plt.ylabel(r'E$_{p,z}$ [keV]',fontsize=25)
    plt.xlabel(r'E$_{\gamma,iso}$ [erg]',fontsize=25)
    plt.tick_params(axis='both',labelsize=20)    
    plt.legend()
    plt.loglog()

from f_CosmoCalc import _cal
#from scipy import interpolate
import math
def plot_line(z_list,Ep,fluence,color,label=None):
    Epz,Eiso = [],[]
    
    z = z_list

    for i in range(len(z)):
        DL_sq = (_cal(z[i], H0=69.6, WM=0.286, WV=0.714)*3.08568e24)**2
        Epz = np.append(Epz, Ep*(1+z[i]))
        Eiso = np.append(Eiso, 4*math.pi*DL_sq*fluence/(1+z[i]))

    plt.plot(Eiso,Epz,color=color,label=label, linewidth=3)
    
def plot_text(z_list,Ep,fluence,color,label=None):
    Epz,Eiso = [],[]
    
    z = z_list

    for i in range(len(z)):
        DL_sq = (_cal(z[i], H0=69.6, WM=0.286, WV=0.714)*3.08568e24)**2
        Epz = Ep*(1+z[i])
        Eiso = 4*math.pi*DL_sq*fluence/(1+z[i])
        plt.plot(Eiso, Epz-0.1*Epz, '^', markersize=6, color=color)
        plt.text(Eiso, Epz-0.3*Epz, f'z={z[i]:.3f}', fontsize=12, color=color)
        
def plot_point(z,Ep,fluence,color,label=None):
    
    DL_sq = (_cal(z, H0=69.6, WM=0.286, WV=0.714)*3.08568e24)**2
    Epz = Ep*(1+z)
    Eiso = 4*math.pi*DL_sq*fluence/(1+z)

    plt.plot(Eiso,Epz,'*',color=color,markersize=20,label=label)

