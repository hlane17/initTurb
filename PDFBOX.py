#!/usr/bin/env python
# coding: utf-8

# In[3]:


import h5py
from glob import glob
import numpy as np
from os.path import isdir
import load_from_snapshot
from sys import argv
import multiprocessing 
from multiprocessing import Pool
sims_dir = "/scratch3/03532/mgrudic/stirring/box_rerun/M2e4_R10_S0_T1_B0.01_Res271_n2_sol0.5_42_BOX/"
filenames = argv[1:]
# initialize lists to store all the stuff we will want in the final data file

nproc = multiprocessing.cpu_count()
for filename in filenames:
    snaps = sorted(glob(sims_dir + filename+ "/snapshot*.hdf5"))
    def Function(snap):
        f = h5py.File(snap, 'r')
        print(f)
        n = 29.9 * np.array(f["PartType0"]["Density"])
        
        logn= np.log10(n)
        bin_counts = np.histogram(logn,bins=np.linspace(-2,6,100))[0]
        
        ids = (n>10)
        density = np.array(f["PartType0"]["Density"])[ids]
        
        time = f["Header"].attrs["Time"]       
        
        f.close()
        return time, bin_counts
 

    data = Pool(nproc).map(Function,snaps)
    data = [np.concatenate([d[0]], d[1]]) for d in data]
    np.savetxt("PDFBOX_" + filename + ".dat", data, 
                header = "#(0) time (1) bincounts"
    )


