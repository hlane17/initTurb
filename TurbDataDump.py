#!/usr/bin/env python
# coding: utf-8

# In[1]:


import h5py
from glob import glob
import numpy as np
from os.path import isdir
import load_from_snapshot
import os


# In[17]:


fileName = "initTurb_stirring_0.25_0.5_10" 
datafolder = "/scratch1/08056/hlane17/initTurb/" + fileName + "/output"
#datafolder = "/home/hlane/scpData/initTurb_stirring_0.25_0.5"
snapNum = 751
vRMSFinal = []
timeFinal = []
centerMassX = []
centerMassY = []
centerMassZ = []
rmsDist = []
medianDist = []
for i in range(snapNum): #This is for taking the FINAL Mach and Fraction#.
    ext='00'+str(i);
    if (i>=10): ext='0'+str(i)                                     #This resolves naming issues
    if (i>=100): ext=str(i)
    f = h5py.File(datafolder + "/snapshot_" + ext + ".hdf5", "r")  #opens file
    ids = np.array(f["PartType0"]["ParticleIDs"])
    n = 29.9 * np.array(f["PartType0"]["Density"]) # density in protons cm^-3
    filter = (n>100)
    velGrad = np.array(f["PartType0"]["Velocities"])
    
    vRMSFinal.append(np.sqrt(np.sum(velGrad**2)/len(velGrad)))

    
    time = load_from_snapshot.load_from_snapshot("Time",0,datafolder,i)
    timeFinal.append(time)

    mass = np.array(f["PartType0"]["Masses"])[filter]

    positionsX = np.array(f["PartType0"]["Coordinates"])[:,0][filter]
    positionsY = np.array(f["PartType0"]["Coordinates"])[:,1][filter]
    positionsZ = np.array(f["PartType0"]["Coordinates"])[:,2][filter]
    centerMassX.append(np.sum(mass*positionsX)/np.sum(mass))
    centerMassY.append(np.sum(mass*positionsY)/np.sum(mass))
    centerMassZ.append(np.sum(mass*positionsZ)/np.sum(mass))
    distances = np.array(f["PartType0"]["Coordinates"])[filter]

    #print(len(np.array(f["PartType0"]["Masses"])))
    cmassX = np.sum(mass*positionsX)/np.sum(mass)
    cmassY = np.sum(mass*positionsY)/np.sum(mass)
    cmassZ = np.sum(mass*positionsZ)/np.sum(mass)
    
    
    distList = []
    distX = distances[:,0]-cmassX
    distY = distances[:,1]-cmassY
    distZ = distances[:,2]-cmassZ
    distance = np.sqrt((distX**2) + (distY**2) + (distZ**2))
    #medianDist.append(np.median(distance))
    rmsDist.append(np.sqrt(np.sum((distance)**2)/len(distance)))
    dx = np.array(f["PartType0"]["Coordinates"])[filter] - np.array([cmassX, cmassY, cmassZ]) # vector from COM
    distances = np.sqrt((dx*dx).sum(1)) # computes the distances sqrt(dX^2 + dY^2 + dZ^2)
    medianDist.append(np.median(distances))


np.savetxt(str(fileName) + ".dat", np.c_[vRMSFinal, timeFinal, centerMassX, centerMassY, centerMassZ, rmsDist, medianDist,], 
           header = "#(0) vRMS (1) time (2) center mass X (3) center mass Y (4) center mass Z (5) rms distance (6) median distance"
)


# In[ ]:




