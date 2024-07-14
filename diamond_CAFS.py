## Preliminary analysis code for the DLS CAFS experiment.

# Data will come in as MDA files.
# Need to be able to extract the frame data from the MDAs as well as the log info for I0, and motor positions.

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
#from pystyle import Colors, Colorate


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

import os, fnmatch
def find(pattern, path):
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))
    return result




input_dir = "data/"

#test = find('pco1-*', input_dir)
#print(test)
test_file = input_dir + "387780.nxs"


#input_dir = 'C:/Users/pauld/CAFS-DLS-Analysis/I13/'
# Test_scan_Ni = 376112 -> 376142
# Test_scan_Ni_3 = 376143 -> 376153
#fnum1 = 376175
#fnum2 = 376176

#file_nums = np.arange(fnum1, fnum2+1, 2, dtype=int)

#frames_pre_energy=10
#Kedge = 8770
#energy_step_array = [(-50, -10, 1)]
#print(energy_step_array[0][0])
#energy_steps = np.arange(energy_step_array[0][0], energy_step_array[0][1] + energy_step_array[0][2], energy_step_array[0][2], dtype=float)
i = 0
I = []
#for files in test[-30:-1]:
    #filename = 'pco1-' + str(fnum) + '.hdf'
    #FULL_FILE = input_dir + filename
f = read_HDF(test_file) 

DF = f[0,0:500,0:500]
DF_mean = np.mean(DF)
print(DF_mean)

c1 = f[0,:,:] - DF_mean
c1[c1<0] = 0
c2 = f[1,:,:] - DF_mean
c2[c2<0] = 0
c3 = f[2,:,:] - DF_mean
c3[c3<0] = 0
c_all = np.sum((c1 + c2 + c3)/3)
print(c_all)

#fig, ax = plt.subplots()
#c = ax.pcolormesh(c_all, cmap='rainbow')
#cbar=fig.colorbar(c,orientation='vertical').set_label(label='Counts',size=14)
#plt.show()


plt.figure()
plt.imshow(c1)
plt.show()

DFs = []
#print(f1)
#start for loop
#for E in energy_steps:
    #import frame data
    
# For alignment: scipy.signal.correlate2d(in1, in2, mode='full', boundary='fill', fillvalue=0)
# Might just stick to a CoM alignment algorithm for now