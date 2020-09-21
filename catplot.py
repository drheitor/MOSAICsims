#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 12:38:16 2020

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
#---------------------------------------

#spec7  s4500g1.00m1.0z-1.00t2.0a0.40_3000_9000.spec.txt
#spec8  s4500g1.00m1.0z-1.50t2.0a0.40_3000_9000.spec.txt
#spec9  s4500g1.00m1.0z-2.00t2.0a0.40_3000_9000.spec.txt
#spec10 s4500g1.00m1.0z-2.50t2.0a0.40_3000_9000.spec.txt
#spec11 s4500g1.00m1.0z-3.00t2.0a0.40_3000_9000.spec.txt
#spec12 s4500g1.00m1.0z-4.00t2.0a0.40_3000_9000.spec.txt


def norm(wl,fl):
    
    fl = fl * u.Unit('J cm-2 s-1 AA-1') 
    wl = wl * u.AA 
    spec = Spectrum1D(spectral_axis=wl, flux=fl) 
    fln = spec / fit_generic_continuum(spec)(spec.spectral_axis) 
    
    fln = [float(i) for i in fln.flux ]
    fln = np.array(fln)
    return fln


#spec='spec5'
spec='spec11'

specname_20='./allspecs/specI_20_'+spec+'_1800-6_figs.txt'
specname_21='./allspecs/specI_21_'+spec+'_1800-6_figs.txt'
specname_22='./allspecs/specI_22_'+spec+'_1800-6_figs.txt'
specname_23='./allspecs/specI_23_'+spec+'_1800-6_figs.txt'
specname_24='./allspecs/specI_24_'+spec+'_1800-6_figs.txt'



wl_20, fl_20= np.genfromtxt(specname_20, skip_header=1, unpack=True)
wl_21, fl_21= np.genfromtxt(specname_21, skip_header=1, unpack=True)
wl_22, fl_22= np.genfromtxt(specname_22, skip_header=1, unpack=True)
wl_23, fl_23= np.genfromtxt(specname_23, skip_header=1, unpack=True)
wl_24, fl_24= np.genfromtxt(specname_24, skip_header=1, unpack=True)






f0 = plt.figure(figsize=(7,10))

ax = f0.add_subplot(111)

ax.plot(wl_20[4:],norm(wl_20,fl_20)[4:]/1.5, color='black', label='I=20')
ax.plot(wl_21[4:],norm(wl_21,fl_21)[4:]/2.4, color='gray', label='I=21')
ax.plot(wl_22[4:],norm(wl_22,fl_22)[4:]/4, color='purple', label='I=22')
ax.plot(wl_23[4:],norm(wl_23,fl_23)[4:]/7, color='firebrick', label='I=23')
ax.plot(wl_24[4:],norm(wl_24,fl_24)[4:]/18, color='blue', label='I=24')


#spec5
#snr=[187,92,56,22,10]

#spec11
snr=[198,101,53,22,10]

ax.text(8400, norm(wl_20,fl_20)[4]/1.7 -0.05, 'S/N: '+str(snr[0]) )
ax.text(8400, norm(wl_21,fl_21)[4]/2.8-0.05, 'S/N: '+str(snr[1]),color= 'gray')
ax.text(8400, norm(wl_22,fl_22)[4]/4  -0.1, 'S/N: '+str(snr[2]),color= 'purple')
ax.text(8400, norm(wl_23,fl_23)[4]/7  -0.05, 'S/N: '+str(snr[3]),color= 'firebrick')
ax.text(8400, norm(wl_24,fl_24)[4]/18 -0.05, 'S/N: '+str(snr[4]),color= 'blue')


ax.text(8714, norm(wl_20,fl_20)[4]/1.7 -0.05, 'I=20  ')
ax.text(8714, norm(wl_21,fl_21)[4]/2.8-0.05, 'I=21  ',color= 'gray')
ax.text(8714, norm(wl_22,fl_22)[4]/4  -0.1, 'I=22  ',color= 'purple')
ax.text(8714, norm(wl_23,fl_23)[4]/7  -0.05, 'I=23  ',color= 'firebrick')
ax.text(8714, norm(wl_24,fl_24)[4]/18 -0.05, 'I=24  ',color= 'blue')


#spec11
ax.set_title('Teff=4500 logg=1.00 [Fe/H]=-3.00')
#spec5
#ax.set_title('Teff=6000 logg=4.00 [Fe/H]=-3.00')

ax.set_xlabel('Wavelength ( $\AA$ )')
ax.set_ylabel('Flux')
#ax.set_xlim([min(wl_20),max(wl_20)])
#ax.set_ylim([min(fl_24)-10,max(fl_20)])



#ax.legend(loc=1, ncol=1)
plt.tight_layout()
#plt.savefig('test_1.pdf')
























































