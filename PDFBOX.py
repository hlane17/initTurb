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
#filenames = ["a_1120_1_10_2e7_y_1_0.01", "a_800_1_10_2e7_y_0.5_0.01", "a_560_1_10_2e7_y_1_0.01", "a_800_1_10_2e7_y_1_0.01", "a_800_0.5_10_2e7_y_1_0.01", "a_800_1_10_2e7_y_1_0.1", "a_800_1_10_1e5_y_1_0.01", "a_800_1_10_2e7_y_1_1", "a_800_1_10_2e6_y_1_0.01", "a_800_1_10_2e7_y_2_0.01", "a_800_2_10_2e7_y_1_0.01",]    
filenames = argv[1:]
# initialize lists to store all the stuff we will want in the final data file

nproc = multiprocessing.cpu_count()
for filename in filenames:
    snaps = sorted(glob(sims_dir + filename+ "/snapshot*.hdf5"))
    def Function(snap):
        rhoList = []
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

    np.savetxt("PDFBOX_" + filename + ".dat", data, 
                header = "#(0) time (1) bincounts"
    )


