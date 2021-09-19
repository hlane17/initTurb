#!/usr/bin/env python
from os import system, chdir, getcwd
from os.path import realpath
from glob import glob
from sys import argv 
from os import popen
cwd = getcwd()

num_resubmits = 3

for f in argv[1:]:
    if "Config.sh" in f or "submit.sh" in f: continue
    rpath = realpath(f)
    chdir(realpath("".join(f.split("/")[:-1])))
    
    for i in range(num_resubmits):
        if i==0: 
            command = "sbatch "+ rpath
        else: 
            command = "sbatch --dependency=afterok:%s "%last_job_id + rpath
        print(command)
        msg = popen(command).read()
        last_job_id = msg.split()[-1]
        print(last_job_id)
        system('''sed -i "s/.txt 0/.txt 1/g" ''' + rpath)
        system('''sed -i "s/.txt 2/.txt 1/g" ''' + rpath)
#        print()
    chdir(cwd)
