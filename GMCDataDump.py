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
sims_dir = "/scratch3/08056/hlane17/GMCTurb/"
#filenames = ["a_1120_1_10_2e7_y_1_0.01", "a_800_1_10_2e7_y_0.5_0.01", "a_560_1_10_2e7_y_1_0.01", "a_800_1_10_2e7_y_1_0.01", "a_800_0.5_10_2e7_y_1_0.01", "a_800_1_10_2e7_y_1_0.1", "a_800_1_10_1e5_y_1_0.01", "a_800_1_10_2e7_y_1_1", "a_800_1_10_2e6_y_1_0.01", "a_800_1_10_2e7_y_2_0.01", "a_800_2_10_2e7_y_1_0.01",]    
filenames = argv[1:]
# initialize lists to store all the stuff we will want in the final data file

nproc = multiprocessing.cpu_count()
magConstant = 5.88 * 10**24
for filename in filenames:
    snaps = sorted(glob(sims_dir+filename+"/output/snapshot*.hdf5"))
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

        masses = np.array(f["PartType0"]["Masses"])[ids]
        massDensity10 = np.sum(masses)
        
        velocities = np.array(f["PartType0"]["Velocities"])[ids]
        
        kineticEnergy = np.sum(0.5 * (masses[:,None]) * (velocities ** 2)) #calculates kinetic energy of the system
        
        Volumes = masses / density
        
        magneticEnergy = np.sum(Volumes*(((np.array(f["PartType0"]["MagneticField"])[:,0][ids])**2 + (np.array(f["PartType0"]["MagneticField"])[:,1][ids])**2 + (np.array(f["PartType0"]["MagneticField"])[:,2][ids])**2) * magConstant))
        gravPotEnergy = 0.5 * np.sum(np.array(f["PartType0"]["Potential"])[ids] * masses)
        
        positionsX = np.array(f["PartType0"]["Coordinates"])[:,0][ids]
        positionsY = np.array(f["PartType0"]["Coordinates"])[:,1][ids]
        positionsZ = np.array(f["PartType0"]["Coordinates"])[:,2][ids]

        distances = np.array(f["PartType0"]["Coordinates"])[ids]
        cmassX = np.sum(masses*positionsX)/np.sum(masses)
        cmassY = np.sum(masses*positionsY)/np.sum(masses)
        cmassZ = np.sum(masses*positionsZ)/np.sum(masses)

        distX = distances[:,0]-cmassX
        distY = distances[:,1]-cmassY
        distZ = distances[:,2]-cmassZ
        distance = np.sqrt((distX**2) + (distY**2) + (distZ**2))
        
        rmsDistCom = (np.sqrt(np.sum((distance)**2)/len(distance)))
        
        dx = np.array(f["PartType0"]["Coordinates"])[ids] - np.array([cmassX, cmassY, cmassZ]) # vector from COM
        distances = np.sqrt((dx*dx).sum(1)) # computes the distances sqrt(dX^2 + dY^2 + dZ^2)
        
        medianDistCom = np.median(distances)
        
        f.close()
        return time, massDensity10, kineticEnergy, magneticEnergy, gravPotEnergy, rmsDistCom, medianDistCom, bin_counts
 

    data = Pool(nproc).map(Function,snaps)
    data1 = [d[:-1] for d in data]
    time = [d[0] for d in data1]
    np.savetxt(sims_dir + filename + "/GMC_" + filename + ".dat", data1, 
                header = "#(0) time (1) mDensity10 (2) kinetic energy (3) magnetic energy (4) gravitational potential energy (5) rmsDistCom (6) medianDistCom"
    )
    bincounts = [d[-1] for d in data]
    np.savetxt(sims_dir + filename + "/PDF_" + filename + ".dat", np.c_[time, bincounts], 
                header = "#(0) time (1) bincounts"
    )

