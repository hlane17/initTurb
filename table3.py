#!/usr/bin/env python
# coding: utf-8

# In[ ]:


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
snapshots = argv[1:]
# initialize lists to store all the stuff we will want in the final data file

magConstant = 5.88 * 10**24
cms = 200
for snapshot in snapshots:
    f = h5py.File(snapshot, 'r')
    print(f)
       
    n = 29.9 * np.array(f["PartType0"]["Density"])      
    ids = (n>10)
    masses = np.array(f["PartType0"]["Masses"])[ids]
    density = np.array(f["PartType0"]["Density"])[ids]

    Volumes = masses / density
    velGrad = np.array(f["PartType0"]["Velocities"])[ids]    
    MachRMS = np.sqrt(np.sum(velGrad**2)/len(velGrad))/cms
        
    magneticEnergy = np.sum(Volumes*(((np.array(f["PartType0"]["MagneticField"])[:,0][ids])**2 + (np.array(f["PartType0"]["MagneticField"])[:,1][ids])**2 + (np.array(f["PartType0"]["MagneticField"])[:,2][ids])**2) * magConstant))
    gravPotEnergy = 0.5 * np.sum(np.array(f["PartType0"]["Potential"])[ids] * masses)
    mus = 0.42*np.sqrt(np.abs(gravPotEnergy)/magneticEnergy)
    time = f["Header"].attrs["Time"]       
        
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
        
        #rmsDistCom = (np.sqrt(np.sum((distance)**2)/len(distance)))
        
    dx = np.array(f["PartType0"]["Coordinates"])[ids] - np.array([cmassX, cmassY, cmassZ]) # vector from COM
    distancesFinal = np.sqrt((dx*dx).sum(1)) # computes the distances sqrt(dX^2 + dY^2 + dZ^2)
        
    medianDistCOM = np.median(distancesFinal)
        
    velocities = np.array(f["PartType0"]["Velocities"])[ids]
    kineticEnergy = np.sum(0.5 * (masses[:,None]) * (velocities ** 2)) #calculates kinetic energy of the system
          
    alpha = 2*kineticEnergy/(np.abs(gravPotEnergy))
        
    f.close()
 

    np.savetxt(sims_dir + "/" + snapshot + ".dat", np.c_[time, mus, MachRMS, medianDistCOM, alpha], 
                header = "#(0) time (1) mus (2) MachRMS (3) R50 (4) alpha "
    )


