#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  6 11:36:49 2020

@author: Heitor
"""

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from specutils import Spectrum1D
import astropy.units as u
from astropy.wcs import WCS
from spectral_cube import SpectralCube
from matplotlib.colors import LogNorm
from subprocess import call
import scipy 
import pyfiglet
import sys
import seaborn as sns


sns.set_style("white")
sns.set_context("paper", font_scale=2.0, rc={"lines.linewidth": 1.5})

#------------------------------------------------------------------------

def snr(a,axis,ddof):
        
    a = np.asanyarray(a)
    m = a.mean(axis)
    sd = a.std(axis=axis, ddof=ddof)
    
    return np.where(sd == 0, 0, m/sd)


# =====================================================================================

def DER_SNR(flux):
   
# =====================================================================================
   """
   DESCRIPTION This function computes the signal to noise ratio DER_SNR following the
               definition set forth by the Spectral Container Working Group of ST-ECF,
	       MAST and CADC. 

               signal = median(flux)      
               noise  = 1.482602 / sqrt(6) median(abs(2 flux_i - flux_i-2 - flux_i+2))
	       snr    = signal / noise
               values with padded zeros are skipped

   USAGE       snr = DER_SNR(flux)
   PARAMETERS  none
   INPUT       flux (the computation is unit independent)
   OUTPUT      the estimated signal-to-noise ratio [dimensionless]
   USES        numpy      
   NOTES       The DER_SNR algorithm is an unbiased estimator describing the spectrum 
	       as a whole as long as
               * the noise is uncorrelated in wavelength bins spaced two pixels apart
               * the noise is Normal distributed
               * for large wavelength regions, the signal over the scale of 5 or
	         more pixels can be approximated by a straight line
 
               For most spectra, these conditions are met.

   REFERENCES  * ST-ECF Newsletter, Issue #42:
               www.spacetelescope.org/about/further_information/newsletters/html/newsletter_42.html
               * Software:
	       www.stecf.org/software/ASTROsoft/DER_SNR/
   AUTHOR      Felix Stoehr, ST-ECF
               24.05.2007, fst, initial import
               01.01.2007, fst, added more help text
               28.04.2010, fst, return value is a float now instead of a numpy.float64
   """
   from numpy import array, where, median, abs 

   flux = array(flux)

   # Values that are exactly zero (padded) are skipped
   flux = array(flux[where(flux != 0.0)])
   n    = len(flux)      

   # For spectra shorter than this, no value can be returned
   if (n>4):
      signal = median(flux)

      noise  = 0.6052697 * median(abs(2.0 * flux[2:n-2] - flux[0:n-4] - flux[4:n]))

      return float(signal / noise)  

   else:

      return 0.0

# end DER_SNR -------------------------------------------------------------------------




#------------------------------------------------------------------------

#file names
print('\n')
print("*******************************************************\n")


ascii_banner = pyfiglet.figlet_format("MOSAICsim PLOTS")
print(ascii_banner)

print("*******************************************************\n")
print('\n')
print('\n')
print('READING ALL FILES ... \n')


#------------------------------------------------------------------------

#--------------------


filename = sys.argv[1]
#filename='I_23_spec11_1800-6'



Original_spec = sys.argv[2]
#Original_spec='./Inspec/p6000g4.00m0.0z-1.00t2.0a0.40_3000_9000.spec.txt'
#Original_spec='p6000g4.00m0.0z-1.00t2.0a0.40_3000_9000.spec'
#Original_spec='flux.norm.nulbad.0.200'
#Original_spec='flat_spec850nm.ascii'


#--------------------

dirr=filename + '_figs'

# mkdir folder
try:    
    os.system("mkdir " + dirr)
except:
    print('...')
    
    
try:    
    os.system("mkdir allspecs")
except:
    print('...')
    
    
    
#--------------------------------------------- 
#output figs names
OUTPUTname=[dirr+'/spec_output.pdf',dirr+'/spec_resample.pdf',dirr+'/fibermask.pdf'
            ,dirr+'/psf.pdf', dirr+'/specRED-Horizon.pdf',dirr+'/SPEC-RED.pdf',dirr+'/SPEC-vert-input.pdf', dirr+'/SKY.pdf', 
            dirr+'/Spintg.pdf', dirr+'/MASKEDregions.pdf' , dirr+'/Redured-Telluric.pdf', 
            dirr+'/Spintg-TelluricFINAL.pdf', dirr+'/Telluric.pdf', dirr+'/PSFprofile.pdf'   ]
#---------------------------------------------
#fits names

Spec_resamp = filename+'/'+filename+'_resampled_template_spec.fits'

Spec = filename+'/'+filename+'_spintg.fits'

Im_maskfiber = filename+'/'+filename+'_mask_fiber_in_aperture.fits'

Im_psf = filename+'/'+filename+'_psf_resampled.fits'

cube = filename+'/'+filename+'_cube_seeing.fits'

reduced = filename+'/'+filename+'_reduced.fits'

reduced_te = filename+'/'+filename+'_reduced_telluric.fits'

SKY = filename+'/'+filename+'_sky.fits'

spec_telluric = filename+'/'+filename+'_spintg_telluric.fits' 

Telluric = filename+'/'+filename+'_telluric.fits'


#---------------------------------------------
#read spec

wl,flux = np.genfromtxt(Original_spec, skip_header=1, unpack=True)


#---------------------------------------------

# SPEC READING 

print("------------------------------\n")

print('OPENNING RESAMPLED SPEC... \n')

print("------------------------------\n")


fits_spec = Spec


hdul = fits.open(fits_spec)

hdul.info()

print(hdul[0].header)

Header1 = open(dirr+'/HeaderResampledSpec-'+dirr+'.txt', 'w')

Header1.write('Resempled spectra Header \n \n')

Header1.write(str(hdul[0].header)+'\n')


dataSPEC = hdul[0].data

initial=hdul[0].header['CRVAL1']
pace=hdul[0].header['CDELT1']
times=hdul[0].header['NAXIS1']

uni=hdul[0].header['BUNIT']

mean = np.mean(dataSPEC)

#     RECOVER wavelen

wave=[]
t=0
while t < times:
    wave.append(initial + (pace*t))
    #print(initial + (pace*t))
    t=t+1
    
    
#----------

f1 = plt.figure(figsize=(10,7))

ax11 = f1.add_subplot(211)

ax11.plot(wave,dataSPEC, color='black', label='spintg Spec')
ax11.set_xlabel('Wavelength ( $\AA$ )')
ax11.set_ylabel(str(uni))
ax11.set_xlim([min(wave),max(wave)])

ax11.legend(loc=4, ncol=1)



ax12 = f1.add_subplot(212, sharex=ax11)

fluxs=flux/np.mean(flux)*mean*1.7
ax12.plot(wl,fluxs, color='gray', lw=2.0 , label='Input Spec')
ax12.plot(wave, dataSPEC, color='black', label='spintg Spec')

ax12.set_ylim([min(dataSPEC)*0.5,max(dataSPEC)*1.1])
ax12.set_xlim([min(wave),max(wave)])
ax12.set_xlabel('Wavelength ( $\AA$ )')
ax12.set_ylabel(str(uni))

ax12.legend(loc=4, ncol=1)

plt.tight_layout()

plt.savefig(OUTPUTname[0])

hdul.close()

print("------------------------------\n")


#---------------------------------------------

# SPEC Resampled template READING 

print('OPENNING SPINTG SPEC... \n')

print("------------------------------\n")



fits_specr = Spec_resamp

hdu = fits.open(fits_specr)

hdu.info()
print(hdu[0].header)

Header2 = open(dirr+'/HeaderSPINTGSpec-'+dirr+'.txt', 'w')

Header2.write('SPINTG spec Header \n \n')

Header2.write(str(hdu[0].header)+'\n')


data_resamp = hdu[0].data


initial=0
pace=1
times=hdu[0].header['NAXIS1']


mean = np.mean(dataSPEC)

#     RECOVER wavelen

wave_resamp=[]
t=0
while t < times:
    wave_resamp.append(initial + (pace*t))
    #print(initial + (pace*t))
    t=t+1


#---------------------------------------------

# Resampled spec 

f2 = plt.figure(figsize=(10,7))

ax21 = f2.add_subplot(111)
ax21.set_title('Resampled spec')

ax21.plot(wave_resamp,data_resamp, color='blue', label='Resampled Spec' )
ax21.set_xlabel('Pix')
ax21.set_ylabel(str(uni))
ax21.legend(loc=4, ncol=1)

plt.tight_layout()

plt.savefig(OUTPUTname[1])

hdu.close()

print("------------------------------\n")

#---------------------------------------------
#IMAGE MASK FIBER

print('OPENNING MASK FIBER... \n')

image_data = fits.getdata(Im_maskfiber)
print(type(image_data))
print(image_data.shape)

Header3 = open(dirr+'/HeaderMASKFIBER-'+dirr+'.txt', 'w')

Header3.write('MASK FIBER Header \n \n')

Header3.write(str(type(image_data))+'\n')
Header3.write(str(image_data.shape)+'\n')


f3 = plt.figure(figsize=(10,7))
ax31 = f3.add_subplot(111)


plt.imshow(image_data, cmap='gray')
plt.colorbar()

plt.savefig(OUTPUTname[2])

print("------------------------------\n")


#---------------------------------------------

#IMAGE PSF

print('OPENNING IMAGE PSF... \n')

image_data = fits.getdata(Im_psf)
print(type(image_data))
print(image_data.shape)

Header4 = open(dirr+'/HeaderPSF-'+dirr+'.txt', 'w')

Header4.write('PSF Header \n \n')

Header4.write(str(type(image_data))+'\n')
Header4.write(str(image_data.shape)+'\n')

x = np.linspace(1, 617, num=617)

f4 = plt.figure(figsize=(10,7))
ax41 = f4.add_subplot(111)


plt.imshow(image_data, cmap='gray',norm=LogNorm())
plt.colorbar()

plt.savefig(OUTPUTname[3])


f13 = plt.figure(figsize=(10,7))
ax131 = f13.add_subplot(111)
ax131.plot(x, image_data[309][1:])

ax131.set_title('PSF')
ax131.set_xlabel('Pix')
ax131.set_ylabel('Flux')

plt.savefig(OUTPUTname[13])


#---------------------------------------------
print("------------------------------\n")
print('OPENNING REDUCED SPEC... \n')
print("------------------------------\n")

# Horizontal Projection Reduced spec

fits_red = reduced

hdur = fits.open(fits_red)

hdur.info()

print(hdur[0].header)


Header5 = open(dirr+'/HeaderREDUCEDSPEC-'+dirr+'.txt', 'w')

Header5.write('Reduced Spec Header \n \n')

Header5.write(str(hdur[0].header)+'\n')


data_reduced = hdur[0].data


lenn=len(hdur[0].data[1])

a=np.linspace(1,lenn, num=lenn)


#horizontal cut
f5 = plt.figure(figsize=(10,7))
ax51 = f5.add_subplot(111)

ax51.plot(a, hdur[0].data[1])
ax51.set_title('Horizontal Projection')
ax21.set_xlabel('Pix')
plt.savefig(OUTPUTname[4])



#---------------------------------------------


#Vertical mean plot

lennb=len(hdur[0].data)

vertical=[]
n=0
while n < lennb:
    meanvert=np.mean(hdur[0].data[n])
    vertical.append(meanvert)
    n=n+1
    
    
b=np.linspace(1,lennb, num=lennb)
    


hdur.close()



#---------------------------------------------
print("------------------------------\n")
print('OPENNING REDUCED TELLURIC SPEC... \n')
print("------------------------------\n")

# Horizontal Projection Reduced spec

fits_red_te = reduced_te

hdurt = fits.open(fits_red_te)

hdurt.info()

print(hdurt[0].header)


Header8 = open(dirr+'/HeaderREDUCEDSPECTELLURIC-'+dirr+'.txt', 'w')

Header8.write('Reduced_telluric Spec Header \n \n')

Header8.write(str(hdurt[0].header)+'\n')


data_reduced_te = hdurt[0].data


#---------------------------------------------


#Vertical mean plot

lennbt=len(hdurt[0].data)

vertical_t=[]
n=0
while n < lennbt:
    meanvert_t=np.mean(hdurt[0].data[n])
    vertical_t.append(meanvert_t)
    n=n+1
    
    
b_t=np.linspace(1,lennbt, num=lennbt)
    



#---------------------------------------------
print("------------------------------\n")
print('OPENNING REDUCED SKY... \n')
print("------------------------------\n")

# SKY

fits_sky= SKY

hdus = fits.open(fits_sky)

hdus.info()

print(hdus[0].header)


Header6 = open(dirr+'/HeaderSKY-'+dirr+'.txt', 'w')

Header6.write('Reduced Spec Header \n \n')

Header6.write(str(hdus[0].header)+'\n')


data_sky = hdus[0].data



f8 = plt.figure(figsize=(10,7))
ax81 = f8.add_subplot(111)

ax81.plot(wave, data_sky, linewidth= 0.9)
ax81.set_title('Sky spectrum')
ax81.set_xlabel('Wavelength ( $\AA$ )')
plt.savefig(OUTPUTname[7])


print("------------------------------\n")
print('OPENNING SPINTG TELLURICS ... \n')
print("------------------------------\n")


fits_stelu= spec_telluric

hdust = fits.open(fits_stelu)

hdust.info()

print(hdust[0].header)


Header9 = open(dirr+'/HeaderSpintgTelluric-'+dirr+'.txt', 'w')

Header9.write('Spintg Telluric Header \n \n')

Header9.write(str(hdust[0].header)+'\n')


data_stelu = hdust[0].data





print("------------------------------\n")
print('OPENNING TELLURICS        ... \n')
print("------------------------------\n")


fits_telu= Telluric

hdutelu = fits.open(fits_telu)

hdust.info()

print(hdutelu[0].header)


Header10 = open(dirr+'/HeaderTelluric-'+dirr+'.txt', 'w')

Header10.write('Telluric Header \n \n')

Header10.write(str(hdutelu[0].header)+'\n')


data_telu = hdutelu[0].data



f12 = plt.figure(figsize=(12,7))
ax121 = f12.add_subplot(111)

ax121.plot(wave, data_telu, linewidth= 0.9)
ax121.set_title('Telluric')
ax121.set_xlabel('Wavelength ( $\AA$ )')
plt.savefig(OUTPUTname[12])




#---------------------------------------------


    
# plot mask I band 
wave = np.array(wave)

#line A B C 

mask_CATa = (wave > 8493) & (wave < 8502)

mask_CATb = (wave > 8532) & (wave < 8550)

mask_CATc = (wave > 8655) & (wave < 8670)

mask_tails = (wave < 8404)


mask_CAT = [any(tup) for tup in zip(mask_CATa, mask_CATb, mask_CATc, mask_tails)]
Imask_CAT = np.invert(mask_CAT)


# sky mask
SKYE = np.array(data_sky)
   
mask_sky = (SKYE > 59.0) 

Imask_sky = np.invert(mask_sky)

mask_ALL = [any(tup) for tup in zip(mask_CAT, mask_sky)]

Imask_ALL = np.invert(mask_ALL)

vertical = np.array(vertical)
vertical_t = np.array(vertical_t)
dataSPEC = np.array(dataSPEC)
data_stelu = np.array(data_stelu)


print("------------------------------\n")
print('        Signal to Noise       \n')
print("------------------------------\n")



print('Vertical S/N: ' + str(int(DER_SNR(vertical))))

print('Vertical Masked S/N: ' + str( int(DER_SNR(vertical[Imask_ALL]))))

print('Vertical Telluric Masked S/N: ' + str( int(DER_SNR(vertical_t[Imask_ALL]))))

print('Spingtg S/N: ' + str(int(DER_SNR(dataSPEC))))

print('Spingtg Masked S/N: ' + str(int(DER_SNR(dataSPEC[Imask_ALL]))))


#PLOT MASK

f10 = plt.figure(figsize=(16,9))

ax101 = f10.add_subplot(111)

ax101.set_title('Masked regions')

#sky
ax101.plot(wave[mask_sky], data_sky[mask_sky])

ax101.plot(wave,dataSPEC, color='black', label='Spintg Spec')
ax101.set_xlabel('Wavelength ( $\AA$ )')
ax101.set_ylabel(str(uni))
ax101.set_xlim([min(wave),max(wave)])
#ax101.text(min(wave), min(vertical), 'S/N: '+ str(int(DER_SNR(dataSPEC[Imask_ALL]))), fontsize=20, color='purple') 

ax101.axvspan(min(wave[mask_CATa]), max(wave[mask_CATa]), alpha=0.2, color='red')
ax101.axvspan(min(wave[mask_CATb]), max(wave[mask_CATb]), alpha=0.2, color='red')
ax101.axvspan(min(wave[mask_CATc]), max(wave[mask_CATc]), alpha=0.2, color='red')


#y2 = np.zeros(len(dataSPEC))
#plt.fill_between(wave[Imask_sky], dataSPEC[Imask_sky], y2[Imask_sky], alpha= 0.3)


ax101.legend(loc=4, ncol=1)

plt.savefig(OUTPUTname[9])




#---------------------------------------------

# FINAL SPEC

f6 = plt.figure(figsize=(10,7))
ax61 = f6.add_subplot(111)
ax61.set_title('Vertical Projection')

ax61.plot(wave, vertical, color='purple', lw=2.0, label='Reduced Spectrum')
ax61.set_xlabel('Wavelength ( $\AA$ )')
ax61.set_ylabel('FLUX')


ax61.legend(loc=4, ncol=1)
plt.savefig(OUTPUTname[5])



f7 = plt.figure(figsize=(10,7))

ax71 = f7.add_subplot(111)
ax71.set_title('Vertical Projection')

mean=np.mean(vertical)
fluxv=flux/np.mean(flux)*mean*1.7

ax71.plot(wave, vertical, color='purple', lw=2.0, label='Reduced Spectrum')
ax71.plot(wl,fluxv, color='gray', lw=2.0 , label='Input Spec')

ax71.set_ylim([min(vertical)*0.5,max(vertical)*1.1])
ax71.set_xlim([min(wave),max(wave)])
ax71.set_xlabel('Wavelength ( $\AA$ )')
ax71.set_ylabel('FLUX')
ax71.text(min(wave), min(vertical), 'S/N: '+str(int(DER_SNR(vertical))) , fontsize=20, color='purple') 

ax71.legend(loc=4, ncol=1)

plt.savefig(OUTPUTname[6])




f9 = plt.figure(figsize=(16,9))

ax91 = f9.add_subplot(111)

ax91.plot(wave,dataSPEC, color='black', label='Spintg Spec')
ax91.set_xlabel('Wavelength ( $\AA$ )')
ax91.set_ylabel(str(uni))
ax91.set_xlim([min(wave),max(wave)])
ax91.text(min(wave), min(vertical), 'S/N: '+ str(int(DER_SNR(dataSPEC[Imask_ALL]))) , color='purple') 

ax91.legend(loc=4, ncol=1)

plt.savefig(OUTPUTname[8])


#-------------------

#FINAL
data_stelu= np.array(data_stelu)

f12 = plt.figure(figsize=(12,7))
ax121 = f12.add_subplot(111)

ax121.plot(wave[5:], data_stelu[5:], linewidth= 0.9, color='black')
ax121.set_title('Spintg Telluric')
ax121.set_xlabel('Wavelength ( $\AA$ )')

ax121.text(min(wave[5:]), min(data_stelu[5:]), 'S/N: '+str(int(DER_SNR(data_stelu[Imask_ALL]))) , fontsize=20, color='purple') 


plt.savefig(OUTPUTname[11])



#---------------------------------------------


#  SAVE spec


file = open(dirr+'/spec'+dirr+'.txt', 'w')

file.write('# Wave Flux \n')
n=0   
while n < lennb:
    WW  = wave[n]
    flu = data_stelu[n]
    file.write('%7.3f %7.3f\n'%(WW,flu))
    n=n+1
  

fil = open(dirr+'/specRESAMPLED'+dirr+'.txt', 'w')

fil.write('# Wave Flux \n')
n=0   
while n < lennb:
    WWr  = wave[n]
    flur = data_resamp[n]
    fil.write('%7.3f %7.3f\n'%(WWr,flur))
    n=n+1  
      

    
    
file = open('allspecs/spec'+dirr+'.txt', 'w')

file.write('# Wave Flux \n')
n=0   
while n < lennb:
    WW  = wave[n]
    flu = data_stelu[n]
    file.write('%7.3f %7.3f\n'%(WW,flu))
    n=n+1
    

fil = open('allspecs/specRESAMPLED'+dirr+'.txt', 'w')

fil.write('# Wave Flux \n')
n=0   
while n < lennb:
    WWr  = wave[n]
    flur = data_resamp[n]
    fil.write('%7.3f %7.3f\n'%(WWr,flur))
    n=n+1  
    
 
#---------------------------------------------

#LOGFILE
    
    
#creating a mask SN small regions, snr all spec , snr combined 
# name of inputspec 
#name of psf
# exp time ...
    
#open log file 
LOG = open(dirr+'/LOGfile-'+dirr+'.txt', 'w')

LOG.write(' Log file of '+ dirr +' \n \n')

LOG.write('PSF:              '+ str(hdul[0].header['PSF']) +' \n')
LOG.write('Input Spectrum:   '+ str(Original_spec) +' \n \n')


# read all this from teh resampled spec header

print(hdul[0].header)

LOG.write('Band:             '+ str(dirr[0])+' \n') 
LOG.write('Mag:              '+ str(hdul[0].header['MAGN']) +' \n') 
LOG.write('R:                '+ str(hdul[0].header['R']) +' \n') 
LOG.write('dit:              '+ str(hdul[0].header['DIT']) +' \n') 
LOG.write('ndit:             '+ str(hdul[0].header['NDIT']) +' \n') 


LOG.write('SpingtgS/N:       ' + str(int(DER_SNR(dataSPEC))) +' \n')

LOG.write('SpingtgTelS/N:    ' + str(int(DER_SNR(data_stelu))) +' \n')

LOG.write('VerticalTelS/N:   ' + str(int(DER_SNR(vertical_t[Imask_ALL]))) +' \n')

LOG.write('SpingtgTelMaskS/N:' + str(int(DER_SNR(data_stelu[Imask_ALL]))) +' \n')

LOG.write('SpingtgMaskedS/N: ' + str(int(DER_SNR(dataSPEC[Imask_ALL]))) +' \n')




print('LOG file saved \n')


    

print("------------------------------\n")

print('SAVING spec in Spec'+filename+'.txt \n')

print("------------------------------\n")

print('DONE!!!\n \n')




#---------------------------------------------

