#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 11:35:43 2020

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



#spec7  s4500g1.00m1.0z-1.00t2.0a0.40_3000_9000.spec.txt
#spec8  s4500g1.00m1.0z-1.50t2.0a0.40_3000_9000.spec.txt
#spec9  s4500g1.00m1.0z-2.00t2.0a0.40_3000_9000.spec.txt
#spec10 s4500g1.00m1.0z-2.50t2.0a0.40_3000_9000.spec.txt
#spec11 s4500g1.00m1.0z-3.00t2.0a0.40_3000_9000.spec.txt
#spec12 s4500g1.00m1.0z-4.00t2.0a0.40_3000_9000.spec.txt

#EWfileRE-I_25_spec3nPSF_norm_1800-6_figs.txt
#EWfile-I_20_spec7_norm_1800-6_figs.txt
#EWfile-I_22_spec9nPSF_norm_1800-6_figs.txt

sns.set_style("white")
sns.set_context("paper", font_scale=2.0, rc={"lines.linewidth": 1.5})


def readEW(spec):


    #EWfile-I_20_spec11_1800-6_figs.txt
    specname_20='./EWs/EWfile-I_20_'+spec+'nPSF_norm_1800-6_figs.txt'
    specname_21='./EWs/EWfile-I_21_'+spec+'nPSF_norm_1800-6_figs.txt'
    specname_22='./EWs/EWfile-I_22_'+spec+'nPSF_norm_1800-6_figs.txt'
    specname_23='./EWs/EWfile-I_23_'+spec+'nPSF_norm_1800-6_figs.txt'
    specname_24='./EWs/EWfile-I_24_'+spec+'nPSF_norm_1800-6_figs.txt'
    specname_25='./EWs/EWfile-I_25_'+spec+'nPSF_norm_1800-6_figs.txt'


    arg_20, num_20= np.genfromtxt(specname_20, skip_header=4, unpack=True, delimiter = ':')
    arg_21, num_21= np.genfromtxt(specname_21, skip_header=4, unpack=True, delimiter = ':')
    arg_22, num_22= np.genfromtxt(specname_22, skip_header=4, unpack=True, delimiter = ':')
    arg_23, num_23= np.genfromtxt(specname_23, skip_header=4, unpack=True, delimiter = ':')
    arg_24, num_24= np.genfromtxt(specname_24, skip_header=4, unpack=True, delimiter = ':')
    arg_25, num_25= np.genfromtxt(specname_25, skip_header=4, unpack=True, delimiter = ':')

    feh_20 = num_20[3]
    feh_21 = num_21[3]
    feh_22 = num_22[3]
    feh_23 = num_23[3]
    feh_24 = num_24[3]
    feh_25 = num_25[3]
              
    Fehlist=[feh_20,feh_21,feh_22,feh_23,feh_24, feh_25]
    
    return Fehlist

def readEWRE(spec):

    # RE means RESAMPLED WHICH IS THE SPECTRUM  WITHOUT NOISE
    #EWfileRE-I_24_spec9_1800-6_figs.txt
    specname_20='./EWs/EWfileRE-I_20_'+spec+'nPSF_norm_1800-6_figs.txt'
    specname_21='./EWs/EWfileRE-I_21_'+spec+'nPSF_norm_1800-6_figs.txt'
    specname_22='./EWs/EWfileRE-I_22_'+spec+'nPSF_norm_1800-6_figs.txt'
    specname_23='./EWs/EWfileRE-I_23_'+spec+'nPSF_norm_1800-6_figs.txt'
    specname_24='./EWs/EWfileRE-I_24_'+spec+'nPSF_norm_1800-6_figs.txt'
    specname_25='./EWs/EWfileRE-I_25_'+spec+'nPSF_norm_1800-6_figs.txt'


    arg_20, num_20= np.genfromtxt(specname_20, skip_header=4, unpack=True, delimiter = ':')
    arg_21, num_21= np.genfromtxt(specname_21, skip_header=4, unpack=True, delimiter = ':')
    arg_22, num_22= np.genfromtxt(specname_22, skip_header=4, unpack=True, delimiter = ':')
    arg_23, num_23= np.genfromtxt(specname_23, skip_header=4, unpack=True, delimiter = ':')
    arg_24, num_24= np.genfromtxt(specname_24, skip_header=4, unpack=True, delimiter = ':')
    arg_25, num_25= np.genfromtxt(specname_25, skip_header=4, unpack=True, delimiter = ':')


    feh_20 = num_20[3]
    feh_21 = num_21[3]
    feh_22 = num_22[3]
    feh_23 = num_23[3]
    feh_24 = num_24[3]
    feh_25 = num_25[3]
              
    Fehlist=[feh_20,feh_21,feh_22,feh_23,feh_24, feh_25]
    
    return Fehlist


# return delta ferro e fazer uma funcao com um for para usar todas as RGBs
    
#spec='spec11'
#RealFeH = - 3.00

specs  = ['spec7', 'spec8','spec9','spec10','spec11','spec12']

specs_colours   = sns.color_palette("viridis", n_colors=6)
#specs_colours   = ['firebrick', 'red','yellow','green','steelblue','blue']
RealFeHs = [-1.0, -1.5,-2.0,-2.5,-3.0,-4.0]

I = [20,21,22,23,24,25]


f0 = plt.figure(figsize=(10,6))
ax = f0.add_subplot(111)

n=0
for spec in specs:
    Feh_cat = readEW(spec) 
    Feh_catRE = readEWRE(spec) 
    
    Delta=[fehcat_i - fehcatre_i for fehcat_i, fehcatre_i in zip(Feh_cat, Feh_catRE)]
    
    #for fe in Feh_cat:
       # h = fe - Feh_catRE[n]
        #Delta.append(h)
             
    ax.scatter(I,Delta, color=specs_colours[n], label='[Fe/H]='+str(RealFeHs[n]), s=470, alpha=0.8)
    ax.plot([19,25], [0,0], linestyle='dashed', linewidth=0.5, color='black', alpha=0.5)
    
    n=n+1


ax.set_ylim([-1.0,1.0])
ax.set_xlim([19.5,24.5])
ax.set_title('$\Delta$[Fe/H] vs. I')
ax.set_xlabel('I')
ax.set_ylabel('$\Delta$[Fe/H]')
ax.legend(loc=3, ncol=2)
ax.text(20, 0.8 , 'V$_{HB}$= -2.0 ')

plt.show()





