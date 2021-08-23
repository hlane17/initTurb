#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import h5py
from glob import glob
import numpy as np
from os.path import isdir
import load_from_snapshot
from sys import argv
filenames = argv[1:]
# initialize lists to store all the stuff we will want in the final data file

for filename in filenames:
    f = h5py.File(filename, 'r')
    print(f)
    n = 29.9 * np.array(f["PartType0"]["Density"])
        
    logn= np.log10(n)
    bin_counts = np.histogram(logn,bins=np.linspace(-2,6,100))[0]
    
    ids = (n>10)
    density = np.array(f["PartType0"]["Density"])[ids]
        
    time = f["Header"].attrs["Time"]       
        
    f.close() 


    np.savetxt("PDFOne_" + filename + ".dat", np.c_[bin_counts], 
            header = "#(0) time (1) bincounts"
)

