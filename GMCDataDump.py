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
sims_dir = "/scratch1/08056/hlane17/GMCTurb/"
#filenames = ["a_1120_1_10_2e7_y_1_0.01", "a_800_1_10_2e7_y_0.5_0.01", "a_560_1_10_2e7_y_1_0.01", "a_800_1_10_2e7_y_1_0.01", "a_800_0.5_10_2e7_y_1_0.01", "a_800_1_10_2e7_y_1_0.1", "a_800_1_10_1e5_y_1_0.01", "a_800_1_10_2e7_y_1_1", "a_800_1_10_2e6_y_1_0.01", "a_800_1_10_2e7_y_2_0.01", "a_800_2_10_2e7_y_1_0.01",]    
filenames = argv[1:]
# initialize lists to store all the stuff we will want in the final data file

nproc = multiprocessing.cpu_count()
magConstant = 5.88 * 10**24
for filename in filenames:
    snaps = sorted(glob(sims_dir+filename+"/output/snapshot*.hdf5"))
    i=0
    def Function(snap):

        f = h5py.File(snap, 'r')
        print(i)
        i+=1
        n = 29.9 * np.array(f["PartType0"]["Density"])
        ids = (n>10)
        density = np.array(f["PartType0"]["Density"])[ids]
        time = f["Header"].attrs["Time"] 
        masses = np.array(f["PartType0"]["Masses"])[ids]
        
        massDensity10 = np.sum(masses)
        
        velocities = np.array(f["PartType0"]["Velocities"])[ids]
        
        kineticEnergy = np.sum(0.5 * (masses[:,None]) * (velocities ** 2)) #calculates kinetic energy of the system
        
        Volumes = masses / density
        
        magneticEnergy = np.sum(Volumes*(((np.array(f["PartType0"]["MagneticField"])[:,0][ids])**2 + (np.array(f["PartType0"]["MagneticField"])[:,1][ids])**2 + (np.array(f["PartType0"]["MagneticField"])[:,2][ids])**2) * magConstant))
        
        
        positionsX = np.array(f["PartType0"]["Coordinates"])[:,0][ids]
        positionsY = np.array(f["PartType0"]["Coordinates"])[:,1][ids]
        positionsZ = np.array(f["PartType0"]["Coordinates"])[:,2][ids]

        distances = np.array(f["PartType0"]["Coordinates"])[ids]
        cmassX = np.sum(masses*positionsX)/np.sum(masses)
        massY = np.sum(masses*positionsY)/np.sum(masses)
        cmassZ = np.sum(masses*positionsZ)/np.sum(masses)

        distX = distances[:,0]-cmassX
        distY = distances[:,1]-cmassY
        distZ = distances[:,2]-cmassZ
        distance = np.sqrt((distX**2) + (distY**2) + (distZ**2))
        
        rmsDistCOM = (np.sqrt(np.sum((distance)**2)/len(distance)))
        
        dx = np.array(f["PartType0"]["Coordinates"])[ids] - np.array([cmassX, cmassY, cmassZ]) # vector from COM
        distances = np.sqrt((dx*dx).sum(1)) # computes the distances sqrt(dX^2 + dY^2 + dZ^2)
        
        medianDistCOM = np.median(distances))
        
        f.close()
        return time, massDensity10, kineticEnergy, magneticEnergy, rmsDistCom, medianDistCom
 

    data = Pool(nproc).map(Function,snaps)
    np.savetxt(sims_dir + filename + "/GMC_" + filename + ".dat", data, 
                header = "#(0) time (1) mDensity10 (2) kinetic energy (3) magnetic energy (4) rmsDistCOM (5) medianDistCOM"
    )

