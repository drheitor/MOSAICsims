#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 10:35:40 2020

@author: Heitor
"""


import numpy as np
from specutils import Spectrum1D
from astropy import units as u
from astropy.visualization import quantity_support
from specutils.fitting import fit_generic_continuum
from specutils import SpectralRegion
from specutils.analysis import equivalent_width
from specutils.analysis import fwhm
import matplotlib.pyplot as plt
import seaborn as sns
import os


sns.set_style("white")
sns.set_context("paper", font_scale=2.0, rc={"lines.linewidth": 1.5})
#------------------------------------------
#spec7  s4500g1.00m1.0z-1.00t2.0a0.40_3000_9000.spec.txt
#spec8  s4500g1.00m1.0z-1.50t2.0a0.40_3000_9000.spec.txt
#spec9  s4500g1.00m1.0z-2.00t2.0a0.40_3000_9000.spec.txt
#spec10 s4500g1.00m1.0z-2.50t2.0a0.40_3000_9000.spec.txt
#spec11 s4500g1.00m1.0z-3.00t2.0a0.40_3000_9000.spec.txt
#spec12 s4500g1.00m1.0z-4.00t2.0a0.40_3000_9000.spec.txt

def EW(specname,name):

    lamb, flux= np.genfromtxt(specname, skip_header=1, unpack=True)

    flux = flux * u.Unit('J cm-2 s-1 AA-1') 
    #flux = flux * u.Unit('erg cm-2 s-1 AA-1') 
    lamb= lamb * u.AA 
    spec = Spectrum1D(spectral_axis=lamb, flux=flux) 
    # normalization is not so good 
    cont_norm_spec = spec / fit_generic_continuum(spec)(spec.spectral_axis) 

    print('-----------'+name+'------------')
#line A
    EWa = equivalent_width(cont_norm_spec, regions=SpectralRegion(8493*u.AA, 8502*u.AA))
    #FWHMa = fwhm(cont_norm_spec, regions=SpectralRegion(8493*u.AA, 8502*u.AA))
    print('EW A line: '+str(EWa))
#line B
    EWb = equivalent_width(cont_norm_spec, regions=SpectralRegion(8533*u.AA, 8551*u.AA))
    print('EW B line: '+str(EWb))
#line C
    EWc = equivalent_width(cont_norm_spec, regions=SpectralRegion(8655*u.AA, 8670*u.AA))
    print('EW C line: '+str(EWc))
#open log file 
    
    #nonlinear to metal-poor
    V_VHB = -2.0
    
    EWbc= (EWb+EWc)
    EWbc= float(EWbc/(1. * u.AA))
    
    EWp = (EWbc)**(-1.5) 
    
    #nonlinear to metal-poor
    #Wl = float(EWb / (1. * u.AA)) + float(EWc / (1. * u.AA)) + (0.64 * V_VHB)
    #FeH= -2.81 + 0.44*Wl
    # FeH constants to V-VHB
    
    a=-2.87
    b=0.195
    c=0.458
    d=-0.913
    e=0.0155
    
    #float all
    
    FeH = a + b * V_VHB + c * EWbc + d * EWp + e * EWbc * V_VHB 
    
    
    
    print('[Fe/H]: '+str(FeH))

    #change relampled spectrum to noise spectrum 
    LOG = open('./EWs/EWfile-'+name+'.txt', 'w')
    #LOG = open('./EWs/EWfileRE-'+name+'.txt', 'w')
    LOG.write('Log file of '+ name +' \n \n')
    LOG.write('Input Spectrum:   '+ specname +' \n \n')
    LOG.write('EW A line:             '+ str(EWa) +' \n') 
    LOG.write('EW B line:             '+ str(EWb) +' \n') 
    LOG.write('EW C line:             '+ str(EWc) +' \n') 
    LOG.write('[Fe/H]_CaT:             '+ str(FeH) +' \n') 

    
    f1 = plt.figure(figsize=(16,9))

    ax = f1.add_subplot(111)
    ax.plot(cont_norm_spec.spectral_axis, cont_norm_spec.flux)
    ax.set_xlim([8480,8690])
    ax.set_ylabel('Flux (J cm-2 s-1 AA-1)')
    ax.set_xlabel('Wavelength ( $\AA$ )')
    ax.axvspan(8498-float(EWa / (2. * u.AA))  , 8498+float(EWa / (2. * u.AA))   , alpha=0.2, color='red')
    ax.axvspan(8542-float(EWb / (2. * u.AA))  , 8542+float(EWb / (2. * u.AA))   , alpha=0.2, color='red')
    ax.axvspan(8662-float(EWc / (2. * u.AA))  , 8662+float(EWc / (2. * u.AA))   , alpha=0.2, color='red')
    
    #change relampled spectrum to noise spectrum 
    plt.savefig('./EWs/EW-figs/EW'+name+'.pdf')
    #plt.savefig('./EWs/EW-figs/EWRE'+name+'.pdf')
    
#------------------------------------------
   
    
try:    
    os.system("mkdir EWs")
except:
    print('...')
    
try:    
    os.system("mkdir EWs/EW-figs")
except:
    print('...')


#spec'+dirr+'.txt'

#specname='specI_20_spec1_1800-6_figs.txt'

os.system("ls -d *norm_1800-6_figs > listEW")

names = np.genfromtxt('listEW', skip_header=1, unpack=True, dtype=str)



for name in names:
    specname = name+'/spec'+name+'.txt'
    #specname = name+'/specRESAMPLED'+name+'.txt'
    EW(specname, name)



#check names e best way to do it and see the diferences between resampled and spec















































