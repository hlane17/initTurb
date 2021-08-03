#!/usr/bin/env python

# Computes gas self-gravity potential for input snapshots
# USAGE: ComputePotential.py /path/to/snapshots*.hdf5
# WARNING: this will irreversibly alter the snapshot you run it on if a potential field is already present - make sure you want this!
# WARNING: This runs in parallel by default - 14 python threads x 4 openMP threads. Don't run on the head node!

from sys import argv
from pytreegrav import Potential
import h5py
import numpy as np
from multiprocessing import cpu_count, Pool
from numba import set_num_threads 
from time import time


cpus = cpu_count()
threads_per_file = 4 # number of threads to run pytreegrav on
num_tasks = cpus//threads_per_file # number of python multiprocessing threads/simultaneous files
set_num_threads(threads_per_file)
def RecomputePotential(f):
    with h5py.File(f,'r') as F:
        print(f) 

        rho = np.array(F["PartType0"]["Density"])
        n = 29.9 * rho
        cut = (n > 10)         # in this case we need to do a cut because we're only looking at stuff with n>10

        m = np.array(F["PartType0"]["Masses"])[cut]
        x = np.array(F["PartType0"]["Coordinates"])[cut]
        h = np.array(F["PartType0"]["SmoothingLength"])[cut]
        t = time()
        phi = Potential(x,m,h,G=4300.7,theta=1.,parallel=True)
        print(time() - t)
    with h5py.File(f,'r+') as F:
        if not "Potential" in F["PartType0"].keys():
            phi0 = np.zeros_like(np.array(F["PartType0"]["Masses"]))
            phi0[cut] = phi
            F["PartType0"].create_dataset("Potential", data=phi0)
        else: 
            phi0 = np.array(F["PartType0"]["Potential"])
            phi0[cut] = phi
            F["PartType0/Potential"].write_direct(phi0)
        
        F.close()
#[RecomputePotential(a) for a in argv[1:]]
Pool(num_tasks).map(RecomputePotential,argv[1:])
