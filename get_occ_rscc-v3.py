#!/usr/bin/python3.8
# 92mrizk@gmail.com
#SBATCH -p low
#SBATCH -n 10
#SBATCH -t 4:0:0
#SBATCH -J OCC_CC_stats


import os, sys, re, shutil, random, time, copy
import numpy as np
from pathlib import Path
from optparse import OptionParser

mydir = os.popen("pwd").read().strip()

parser = OptionParser()
parser.add_option("-l", "--ligand", dest="lig")
parser.add_option("-p", "--model", dest="pdb")
parser.add_option("-c", "--cif", dest="cif", help="the cif file")
parser.add_option("-d", "--directory", dest="directory", default=".", help="be default this will run here. But it can take a directory, or a bunch of directories starting with the same string; for example Group_1 Group_2 .., it is enough to add '-d Group'")
parser.add_option("-r", "--software", dest="ref_software", default="refmac", help="by default it runs refmac5, but you can choose phenix '-r phenix'")
parser.add_option("-f", "--hkl", dest="fhkl", default="merged.hkl", help="hkl file, usually merged.hkl or scaled.hkl from HCA, sometimes XDS_ASCII.HKL from full rotation.")
parser.add_option("", "--slurm", action="store_true", dest="slurm", default=False)
(options, args) = parser.parse_args()

ref_sof = options.ref_software

if options.lig is None:
    lig = input ("Enter ligand or residue name: ")
else:
    lig = options.lig
LIG = lig.upper()
if options.slurm:
    slurm = "--slurm"
else:
    slurm = ""
if options.cif:
    cifpath = os.path.abspath(options.cif)
    cif = ' --cif '+str(cifpath)
else:
    cif = " "

if options.pdb is None:
    pdb = input ("Enter pdb name: ")
else:
    pdb = options.pdb 
pdbname = Path(pdb).stem


with open('occ_rscc_'+str(LIG)+'-'+str(options.directory)+'_'+str(ref_sof)+'.log', 'w') as outd:
    outd.write('occ and RSCC in '+str(options.directory)+' \ndirs\t\tOCC\t\tRSCC\n')

                    # CREATING TMP DIRECTORY
rand=int(round(random.random(), 6)*100000)
tmp="tmpfiles-"+str(rand)
try:
    os.mkdir(tmp)
except:
    print(str(tmp)+" already exists .. weird!")

                    #    GETTING THE PDB
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
    except:
        print("pdb name is not recognizable .. quitting")
        quit
    
os.chdir(tmp)

def run_prog(folder, input_dir):
    for hkl in os.listdir(folder):
        folder_name=os.path.basename(folder)
        if hkl == options.fhkl:
            shutil.copyfile(str(folder)+'/'+str(hkl), hkl)
            os.system('convert_hkl2mtz.py -f '+str(hkl))
            os.system('sbatch --wait rocc_phe-ref.py  \
                -l '+str(LIG)+' \
                -p '+str(mydir)+'/'+str(pdbname)+'.pdb \
                '+str(cif)+'  \
                --'+str(ref_sof)+'  \
                '+str(slurm))
            shutil.copy('details_occ_'+str(LIG)+'.txt', +str(folder)+'/details_occ_'+str(LIG)+'-'+str(folder_name)+'.txt')
            os.system('get_rscc-v2.py -l '+str(LIG)+' -p '+str(mydir)+'/'+str(pdbname)+'.pdb -r '+str(ref_sof))
            file = open('details_occ_'+str(LIG)+'_'+str(tmp)+'_'+str(ref_sof)+'.txt', "r")
            for line in file:
                if re.search('accurate', line):
                    occ = line.split()[5]
            shutil.copy('RSCC_'+str(pdbname)+'_'+str(LIG)+'_'+str(tmp)+'.txt', str(folder)+'/RSCC_'+str(pdbname)+'_'+str(LIG)+'_'+str(folder_name)+'.txt')
            file = open('RSCC_'+str(pdbname)+'_'+str(LIG)+'_'+str(tmp)+'.txt', "r")
            for line in file:
                    if (re.search(str(LIG), line[5:8])):
                        rscc = line[30:36]
            with open(str(mydir)+'/occ_rscc_'+str(LIG)+'_'+str(ref_sof)+str(input_dir)+'.log', 'a') as outd:
                outd.write(str(folder)+'\t\t'+str(occ)+'\t\t'+str(rscc)+'\n')
            if (os.getcwd().endswith(str(rand))):
                for dele in os.listdir():
                    os.remove(dele)
                    
if options.directory == "." or options.directory == "./":
    input_dir = ""
    run_prog(mydir,input_dir)
else:
    for folder in os.listdir(mydir):
        if folder.startswith(str(options.directory)):
            target_dir=mydir+'/'+str(folder)
            input_dir="-"+str(options.directory)
            run_prog(target_dir, input_dir)
            