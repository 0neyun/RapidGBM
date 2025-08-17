'''
HEtools Spectral fit module demo

wangyun@pmo.ac.cn
2024.9
'''

import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
os.environ["OMP_NUM_THREADS"] = '1'

import numpy as np
import bilby
from astropy.io import fits
import xspec
import matplotlib.pyplot as plt
import corner
from astropy.io import fits
from scipy.stats import gaussian_kde
from scipy.optimize import *
from scipy.integrate import quad


def mkdir(path):
    if os.path.exists(path):
        pass
    else:
        os.makedirs(path)  
        
def find_interval2(chain,conf=0.6826):
    KDE = gaussian_kde(chain)
    minv = np.sort(chain)[0]
    maxv = np.sort(chain)[-1]
    lx = np.linspace(minv,maxv,10000)
    ly = KDE(lx)
    loc = np.where(ly==max(ly))
    best_gues = lx[loc]
    
    def func_KDE(x):
        return -KDE(x)
    minimum = fmin(func_KDE, best_gues)
    bestv = float(minimum)
    print ('HPD_value =',bestv)
    def find_ins(ins, func):
        def f(x):
            pdf = func(x)
            return pdf if pdf > ins else 0
        prob = quad(f, minv, maxv)[0]
        return prob-conf
    def find_root(func, ins, low, hig):
        return ridder(lambda x: func(x)-ins, low, hig)
    
    ins = ridder(find_ins, 0, KDE(bestv), KDE)
    hpd_low = find_root(KDE, ins, minv, bestv)
    hpd_hig = find_root(KDE, ins, bestv, maxv)
    print ('KDE_integrate =',KDE.integrate_box_1d(low=hpd_low, high=hpd_hig))
    return hpd_low,bestv,hpd_hig

def find_interval3(chain,conf=0.6826):
    KDE = gaussian_kde(chain)
    minv = np.sort(chain)[0]
    maxv = np.sort(chain)[-1]    
    lx = np.linspace(minv,maxv,10000)
    ly = KDE(lx)
    loc = np.where(ly==max(ly))
    best_gues = lx[loc]
    
    def func_KDE(x):
        return -KDE(x)
    minimum = fmin(func_KDE, best_gues)
    bestv = float(minimum)
    print ('HPD_value =',bestv)
    
    step = (maxv-minv)/100000
    i,j = bestv,bestv
    while minv < i and i < maxv and minv < j and j < maxv and KDE.integrate_box_1d(low=i, high=j) < conf:
        if KDE(i) >= KDE(j):
            i = i - step
        else:
            j = j + step
    hpd_low = i
    hpd_hig = j
    print ('KDE_integrate =',KDE.integrate_box_1d(low=hpd_low, high=hpd_hig))
    return hpd_low,bestv,hpd_hig

def plot_corner_for_chain(chain_file,savepath,HPD=False):
    temp = fits.open(chain_file)
    labels = temp['CHAIN'].data.names[:-1]
    para_chain = []
    para_meadin = []
    for i in range(len(labels)):
        para_chain.append(temp['CHAIN'].data[labels[i]])
        para_meadin.append(np.median(temp['CHAIN'].data[labels[i]]))
    sample = np.reshape(para_chain,(len(labels),len(temp['CHAIN'].data[labels[i]]))).T
    if HPD is not True:
        fig = corner.corner(sample,labels=labels,quantiles=[0.16, 0.5, 0.84],\
                            truths=para_meadin,show_titles=True,smooth=True)
        
    else:
        #get parameter boundary
        para_bound = []
        for i in range(len(labels)):
            try:
                para_bound = np.append(para_bound, find_interval2(chain=temp['CHAIN'].data[labels[i]]))
            except:
                para_bound = np.append(para_bound, find_interval3(chain=temp['CHAIN'].data[labels[i]]))    
        para_bound = np.reshape(para_bound, (len(labels),3))
        
        fig = corner.corner(sample,labels=labels,smooth=True)
        corner.overplot_lines(fig, para_bound[:,1])
        corner.overplot_points(fig, np.array(para_bound[:,1])[None],marker='s')
        axes = np.array(fig.axes).reshape(len(labels),len(labels))
        
        mcmc_result = open('%s_result.txt' %chain_file[:-5], 'w')
        for i in range(len(labels)):
            ax = axes[i,i]
            ax.axvline(x=para_bound[i][0],linestyle='--')
            ax.axvline(x=para_bound[i][2],linestyle='--')
            ax.set_title(r'%s = $%4.2f_{-%4.2f}^{+%4.2f}$' %(labels[i],para_bound[i][1],para_bound[i][1]-para_bound[i][0],para_bound[i][2]-para_bound[i][1]))
            mcmc_result.write('%s %s %s %s\n' %(str(labels[i]),str(para_bound[i][1]),str(para_bound[i][1]-para_bound[i][0]),str(para_bound[i][2]-para_bound[i][1])))
        mcmc_result.close()
    plt.savefig(savepath)
    plt.close()
    
    
xspec.Xset.allowPrompting = False
#xspec.Xset.chatter = 0
xspec.Xset.abund = "wilm"
xspec.AllModels.lmod("grblocalmod",'/home/wangyun/grblocalmodels')

def Fit(data_path, result_path, tot_src, model_set, HPD=False):
    xspec.AllData.clear()
    xspec.AllModels.clear()
    xspec.AllChains.clear()    
    for i in range(len(tot_src)):
        if '.pi' in tot_src[i]['pha_file']:
            os.chdir(data_path)
            s1 = xspec.Spectrum(dataFile=tot_src[i]['pha_file'],)
            os.chdir('../../../')
        else:
            s1 = xspec.Spectrum(dataFile=data_path+'/'+tot_src[i]['pha_file'],
                                backFile=data_path+'/'+tot_src[i]['bak_file'],
                                respFile=data_path+'/'+tot_src[i]['rsp_file']+'{1}',
                                arfFile=tot_src[i]['arf_file'])
        s1.ignore(tot_src[i]['ign_range'])
        xspec.Fit.statMethod = '%s %s' %(str(i+1), tot_src[i]['statistic'])    
    xspec.AllData.ignore("bad")
    a = xspec.Model(model_set['name'])
    
    if model_set['params_set'] is not None:
        for i in range(len(model_set['params_set'])):
            a.setPars(model_set['params_set'][i])
    
    
    xspec.Fit.renorm()
    xspec.Fit.query= "yes"
    xspec.Fit.perform()
    
    try:
        chain_exists = os.path.exists("%s/%s_chain.fits" %(result_path,model_set['name']))
        if chain_exists:
            os.remove("%s/%s_chain.fits" %(result_path,model_set['name']))
        c1 = xspec.Chain("%s/%s_chain.fits" %(result_path,model_set['name']),\
                         runLength=10000,burn=2000, walkers=20)
        c1.run()
        print ('mcmc done.')
    except Exception as e:
        print (f'MCMC ERROR: {e}')
    
    else:
        plot_corner_for_chain(chain_file="%s/%s_chain.fits" %(result_path,model_set['name']),\
                              savepath="%s/%s_corner.png" %(result_path,model_set['name']),HPD=HPD)
    
    fit_result = open('%s/%s_fitresult.txt' %(result_path,model_set['name']), 'w+')
    
    model_para = xspec.AllModels(1)
    for i in range(xspec.AllModels(1).nParameters):
        i = i + 1
        fit_result.write('%s: %s +/- %s    |' %(str(model_para(i).name), str(model_para(i).values[0]), str(model_para(i).sigma)))
        xspec.Fit.error(f"1.0 {str(i)}")
        par_temp = xspec.AllModels(1)(i)
        fit_result.write('err: '+str(par_temp.error)+'\n')
    

    #print ('k = ',xspec.AllModels(1).nParameters)
    #print ('k = ',xspec.Fit.nVarPars)
    fit_result.write('stat/dof: %s/%s = %s\n' %(str(xspec.Fit.statistic),str(xspec.Fit.dof),str(xspec.Fit.statistic/xspec.Fit.dof)))
    chiq = xspec.Fit.statistic
    k = xspec.Fit.nVarPars
    N = xspec.Fit.dof + k
    bic = chiq + np.log(N)*k
    aic = chiq + 2*k
    fit_result.write('BIC: %s\n' %str(bic))
    fit_result.write('AIC: %s\n' %str(aic))
    fit_result.close()
    

    c = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf', \
         "#9467bd","#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf","#aec7e8", "#ffbb78","#dbdb8d", "#9edae5"]
    
    xspec.Plot.xAxis='KeV'
    xspec.Plot.yLog='true'
    xspec.Plot.setRebin(minSig=3,maxBins=5)    
    
    xspec.Plot("ldata")
    Energy_bound_list = []
    for i in range(len(tot_src)):
        plt.errorbar(x=xspec.Plot.x(i+1),xerr=xspec.Plot.xErr(i+1),y=xspec.Plot.y(i+1),yerr=xspec.Plot.yErr(i+1),fmt='.', color=c[i], label=tot_src[i]['det'])
        plt.plot(xspec.Plot.x(i+1),xspec.Plot.model(i+1), color=c[i])
        Energy_bound_list = np.append(Energy_bound_list, np.array(xspec.Plot.x(i+1))[0])
        Energy_bound_list = np.append(Energy_bound_list, np.array(xspec.Plot.x(i+1))[-1])
    
    E_min = np.min(Energy_bound_list)    
    E_max = np.max(Energy_bound_list)   
    plt.loglog()
    plt.legend(loc=0)
    plt.xlim(E_min,E_max)
    plt.ylim(1e-7,)
    plt.xlabel('Energy (keV)')
    plt.ylabel('counts/s/keV')    
    plt.savefig('%s/%s_ldata.png' %(result_path,model_set['name']))
    plt.close()
    
    xspec.Plot("euf")
    for i in range(len(tot_src)):
        plt.errorbar(x=xspec.Plot.x(i+1),xerr=xspec.Plot.xErr(i+1),y=xspec.Plot.y(i+1),yerr=xspec.Plot.yErr(i+1),fmt='.', color=c[i], label=tot_src[i]['det'])
        plt.plot(xspec.Plot.x(i+1),xspec.Plot.model(i+1), color=c[i])
    plt.loglog()
    plt.legend(loc=0)
    plt.xlim(E_min,E_max)
    plt.ylim(1e-5,)
    plt.xlabel('Energy (keV)')
    plt.ylabel(r'keV (Photons/cm${^2}$/s/keV)')    
    plt.savefig('%s/%s_euf.png' %(result_path,model_set['name']))
    plt.close()    
    
    xspec.Plot("eeuf")
    for i in range(len(tot_src)):
        plt.errorbar(x=xspec.Plot.x(i+1),xerr=xspec.Plot.xErr(i+1),y=xspec.Plot.y(i+1),yerr=xspec.Plot.yErr(i+1),fmt='.', color=c[i], label=tot_src[i]['det'])
        plt.plot(xspec.Plot.x(i+1),xspec.Plot.model(i+1), color=c[i])
    plt.loglog()
    plt.legend(loc=0)
    plt.xlim(E_min,E_max)
    plt.ylim(1e-1,)
    plt.xlabel('Energy (keV)')
    plt.ylabel(r'keV${^2}$ (Photons/cm${^2}$/s/keV)')
    plt.savefig('%s/%s_eeuf.png'%(result_path,model_set['name']))
    plt.close()    
    
    
    xspec.AllData.clear()
    xspec.AllModels.clear()
    xspec.AllChains.clear()
    
       
       
       
