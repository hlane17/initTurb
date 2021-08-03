#!/usr/bin/env python
# coding: utf-8

# In[3]:


import h5py
from glob import glob
import numpy as np
from os.path import isdir
import load_from_snapshot
import os
from os import system, mkdir, chdir
import multiprocessing 
from multiprocessing import Pool
# In[8]:
magConstant = 5.88 * 10**24

nproc = multiprocessing.cpu_count()

sims_dir = "/scratch3/08056/hlane17/GMCTurb/"
filename = "a_1120_1_10_2e7_y_1_0.01" #["a_1120_1_10_2e7_y_1_0.01", "a_800_1_10_2e7_y_0.5_0.01", "a_560_1_10_2e7_y_1_0.01", "a_800_1_10_2e7_y_1_0.01", "a_800_1_10_2e7_y_1_0.1", "a_800_1_10_1e5_y_1_0.01", "a_800_1_10_2e7_y_1_1", "a_800_1_10_2e6_y_1_0.01", "a_800_1_10_2e7_y_2_0.01", "a_800_2_10_2e7_y_1_0.01",]    

def DataDump(sim):
    for dir in glob(sims_dir+ filename + "/output"):
        sim = sim.split("/")[-1]
        time = []
        massDensity10 = []
        KEs = []
        MEs = []
        GEs = []
        rmsDistCOM = []
        medianDistCOM = []
        for i in range(751):
            ext='00'+str(i);
            if (i>=10): ext='0'+str(i)                                     #This resolves naming issues
            if (i>=100): ext=str(i)
            current_snap = sorted(glob(dir+"/snapshot*.hdf5"))[i]            
            print("snap_" + str(i))
            f = h5py.File(current_snap, "r")  #opens file
            n = 29.9 * np.array(f["PartType0"]["Density"])
            ids = (n>10)
            density = np.array(f["PartType0"]["Density"])[ids]
            time.append(load_from_snapshot.load_from_snapshot("Time",0,dir,i))
            masses = np.array(f["PartType0"]["Masses"])[ids]
            massDensity10.append(np.sum(masses))
            velocities = np.array(f(["PartType0"]["Velocities"]))[ids]
            kineticEnergy = np.sum(0.5 * (masses[:,None]) * (velocities ** 2)) #calculates kinetic energy of the system
            KEs.append(kineticEnergy)    
            Volumes = masses * density
            MagneticEnergy = np.sum(Volumes*(((np.array(f["PartType0"]["MagneticField"])[:,0][ids])**2 + (np.array(f["PartType0"]["MagneticField"])[:,1][ids])**2 + (np.array(f["PartType0"]["MagneticField"])[:,2][ids])**2) * magConstant))
            MEs.append(MagneticEnergy)
        
            GPE = 0.5 * np.sum(np.array(f(["PartType0"]["Potential"]))[ids] * masses)
            GEs.append(GPE)
        
            positionsX = np.array(f["PartType0"]["Coordinates"])[:,0][ids]
            positionsY = np.array(f["PartType0"]["Coordinates"])[:,1][ids]
            positionsZ = np.array(f["PartType0"]["Coordinates"])[:,2][ids]

            distances = (f["PartType0"]["Coordinates"])[ids]

            cmassX = np.sum(masses*positionsX)/np.sum(masses)
            cmassY = np.sum(masses*positionsY)/np.sum(masses)
            cmassZ = np.sum(masses*positionsZ)/np.sum(masses)
    
            distX = distances[:,0]-cmassX
            distY = distances[:,1]-cmassY
            distZ = distances[:,2]-cmassZ
            distance = np.sqrt((distX**2) + (distY**2) + (distZ**2))
            rmsDistCOM.append(np.sqrt(np.sum((distance)**2)/len(distance)))
            dx = np.array(f["PartType0"]["Coordinates"])[ids] - np.array([cmassX, cmassY, cmassZ]) # vector from COM
            distances = np.sqrt((dx*dx).sum(1)) # computes the distances sqrt(dX^2 + dY^2 + dZ^2)
            medianDistCOM.append(np.median(distances))
            f.close()
 
        np.savetxt("GMC_" + filename + ".dat", np.c_[time, massDensity10, KEs, MEs, GEs, rmsDistCOM, medianDistCOM], 
                    header = "#(0) time (1) mDensity10 (2) kinetic energy (3) magnetic energy (4) gravitational energy (5) rmsDistCOM (6) medianDistCOM"
        )
dirs = sims_dir + filename
Pool(nproc).map(DataDump, dirs)

