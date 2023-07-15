# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 13:49:54 2023

@author: dansi
"""

import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d 
import scipy.signal
from mpl_toolkits.mplot3d import Axes3D
import h5py
from nexusformat.nexus import *
from matplotlib import cm
from scipy.optimize import curve_fit
import sys
import json
from scipy import interpolate
import matplotlib.colors as colors
import pystyle
from pystyle import Colors, Colorate

def read_HDF(filename):
    f = h5py.File(filename, 'r')

    #Print the full path inside the file 
    for key in f.keys():
        print(key) #Names of the root level object names in HDF5 file - can be 
            #groups or datasets.
        print(type(f[key])) # get the object type: usually group or dataset
        print(f[key])

    #Get the HDF5 group; key needs to be a group name from above
    group = f[key]

    #Checkout what keys are inside that group.
    for key in group.keys():
        print('####################################')
        print(key)
        print(type(group[key]))
        print(group[key])

    i_dataset = group['instrument']

    for key in i_dataset.keys():
        print('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')
        print(key) #Names of the root level object names in HDF5 file - can be groups or datasets.
        print(type(i_dataset[key])) # get the object type: usually group or dataset
        print(i_dataset[key])

    d_dataset = i_dataset['detector']

    for key in d_dataset.keys():
        print('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
        print(key) #Names of the root level object names in HDF5 file - can be groups or datasets.
        print(type(d_dataset[key])) # get the object type: usually group or dataset
        print(d_dataset[key])

    # This assumes group[some_key_inside_the_group] is a dataset,
    # and returns a np.array:

    dataset = group['data']

    for key in dataset.keys():
        print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        print(key) #Names of the root level object names in HDF5 file - can be groups or datasets.
        print(type(dataset[key])) # get the object type: usually group or dataset
        print(dataset[key])
    return dataset['data']

File_Num = '376149'

filename = 'C:/Users/dansi/OneDrive/Documents/I13/pco1-376149.hdf'

hdf = read_HDF(filename)

camera1 = hdf[0,1810:2185,1100:1300]

fig, ax = plt.subplots()
c = ax.pcolormesh( camera1, cmap='rainbow')
cbar=fig.colorbar(c,orientation='vertical').set_label(label='Counts',size=14)
plt.show()


#Read in data from the nexus file 
nxs=nxload('C:/Users/dansi/OneDrive/Documents/I13/376149.nxs')
#Print the full file to see what we've got
print(nxs.tree)