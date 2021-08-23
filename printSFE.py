#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import h5py
from glob import glob
import numpy as np
from os.path import isdir
import load_from_snapshot
import os
from os import system, mkdir, chdir


# In[ ]:


mtot = 2e4
snapSFE = argv[1:]
for snapshot in snapSFE:
    f = h5py.File(snapshot, 'r')
    starmassTot = np.sum(np.array(f["PartType5"]["Masses"]))
    SFE = starmassTot/mtot
    print(SFE)

