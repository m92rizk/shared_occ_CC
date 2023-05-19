#!/usr/bin/python3.8
# 92mrizk@gmail.com
# slurm parameters
#SBATCH -p low
#SBATCH -n 10
#SBATCH -t 4:0:0
#SBATCH -J OCC_CC_stats

import os, sys, re, shutil, random, time, copy
import numpy as np
from pathlib import Path
from optparse import OptionParser

mydir = os.getcwd()

parser = OptionParser()

parser.add_option("-l", "--ligand", dest="lig")
parser.add_option("-p", "--model", dest="pdb")
parser.add_option("-d", "--directory", dest="directory", default="best_solution")
parser.add_option("-f", "--hkl", dest="fhkl", default="merged.hkl", help="hkl file, usually merged.hkl or scaled.hkl from HCA, sometimes XDS_ASCII.HKL from full rotation.")

(options, args) = parser.parse_args()


if options.lig is None:
    lig = input ("Enter ligand (or amino acid) name: ")
else:
    lig = options.lig
LIG = lig.upper()


if options.pdb is None:
    pdb = input ("Enter pdb name: ")
else:
    pdb = options.pdb
pdbname = Path(pdb).stem


with open('occ_rscc_'+str(LIG)+'-'+str(options.directory)+'.log', 'w') as outd:
    outd.write('occ and RSCC in '+str(options.directory)+' \ndirs\n\t\tOCC\t\tRSCC\n')

#                     CREATING TMP DIRECTORY

rand=int(round(random.random(), 6)*100000)
tmp="tmpfiles-"+str(rand)
try:
    os.mkdir(tmp)
except:
    print(str(tmp)+" already exists .. weird!")

#                        GETTING THE PDB
if os.path.isfile(pdb):
    print('============================================\n getting the pdb')
    try:
        shutil.copyfile(str(pdb), str(pdbname)+'.pdb')
    except shutil.SameFileError:
        print('pdb already here')
else:
    print('============================================\n pdb does not exist, will download it')
    try:
        os.system('wget http://www.rcsb.org/pdb/files/'+str(pdbname)+'.pdb >/dev/null')
        os.system('fixpdb.sh '+str(pdbname)+' '+str(LIG))
#        shutil.move(pdb, str(tmp)+'/')
    except:
        print("pdb name is not recognizable .. quitting")
        quit
        
os.chdir(tmp)

for folder in os.listdir(mydir):
    if folder.startswith(str(options.directory)):
        for hkl in os.listdir(mydir+'/'+folder):
            if hkl == options.fhkl:
                shutil.copyfile('../'+str(folder)+'/'+str(hkl), hkl)
                os.system('convert_mtz_hkl.py -f '+str(hkl))
                os.system('sbatch --wait get_occ_parallel.py -l '+str(LIG)+' -p ../'+str(pdbname)+'.pdb --overwrite')
                shutil.copy('details_occ_'+str(LIG)+'.txt', '../'+str(folder)+'/details_occ_'+str(LIG)+'-'+str(folder)+'.txt')
                os.system('get_rscc.py -l '+str(LIG)+' -p ../'+str(pdbname)+'.pdb')
                file = open('details_occ_'+str(LIG)+'.txt', "r")
                for line in file:
                    if re.search('accurate', line):
                        occ = line.split()[5]
                file = open('RSCC_'+str(pdbname)+'_'+str(LIG)+'.txt', "r")
                for line in file:
                        if (re.search(str(LIG), line[5:8])):
                            rscc = line[30:36]
                with open(str(mydir)+'/occ_rscc_'+str(LIG)+'-'+str(options.directory)+'.log', 'a') as outd:
                    outd.write(str(folder)+'\t\t'+str(occ)+'\t\t'+str(rscc)+'\n')
                if (os.getcwd().endswith(str(rand))):
                    for dele in os.listdir():
                        os.remove(dele)
