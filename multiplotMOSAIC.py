#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 14:46:23 2020

@author: Heitor
"""

import numpy as np 
from astropy.io import fits
import os 
#import time
#import datetime


# read header to use the spectrum name before run each step


#list of directories ( I_mag_spec1...) format
os.system('ls -d *nPSF_norm_1800-6 > list')


#list of spectra format: spec1 name.spec.txt
NAMES = np.genfromtxt('list', dtype=str, unpack=True)

tags, specs = np.genfromtxt('./inspec/lista', dtype=str, unpack=True)



for name in NAMES:
    n=0
    
    for tag in tags:
        if tag == name[5:10]:
            
            command = 'python3 MOSAICfits.py '+ name +' '+ './inspec/'+str(specs[n])
            print('****'+str(specs[n]+'****'))
            print(tag+ '\n')
            print(name[5:10]+ '\n')
            print(command+ '\n')
            os.system(command)
        n=n+1
    print(name+' DONE!!!')
    


print('*------------DONE------------*')
    
    





















