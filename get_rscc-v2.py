#!/usr/bin/python3.8

import os, sys, re, shutil,copy
from pathlib import Path
from optparse import OptionParser

mydir = os.getcwd()

parser = OptionParser()
parser.add_option("-l", "--ligand", dest="lig")
parser.add_option("-p", "--pdb", dest="pdb")
parser.add_option("-m", "--mtz", dest="mtz")
parser.add_option("-d", "--directory", dest="directory", default=".")
parser.add_option("-r", "--software", dest="ref_software", default="refmac")
(options, args) = parser.parse_args()

ref_sof = options.ref_software

if options.lig is None:
    lig = input ("Enter ligand or residue name: ")
else:
    lig = options.lig
LIG = lig.upper()

if options.pdb:
    pdb = options.pdb    
    pdbname = Path(pdb).stem

mtz = options.mtz
folder = options.directory
if folder == "." or folder == "./" :
    folder_name = os.path.basename(os.getcwd()).strip()
    here = True
else:
    here = False
    folder_stripped = folder
    if folder.endswith('/'):
        folder_stripped = folder.strip('/')
    folder_name = os.path.basename(folder_stripped)
    
def rscc(folder,folder_name,clpdb,clmtz):
    print('Calculating real space correlation coefficient of '+str(LIG)+' in '+str(folder)+'/'+str(clpdb))
    os.system('phenix.real_space_correlation '+str(folder)+'/'+str(clpdb)+' '+str(folder)+'/'+str(clmtz)+' detail=residue > '+str(folder)+'/RSCC_'+str(LIG)+'_'+str(folder_name)+'.txt')
    file = open(str(folder)+'/RSCC_'+str(LIG)+'_'+str(folder_name)+'.txt', "r")
    for line in file:
        if (re.search(str(LIG), line[5:8])):
            rscc = line[30:36]
            print ('The RSCC of '+str(LIG)+' in '+str(clpdb)+' is: '+rscc)


if here:
    if mtz:  # the given mtz and pdb are taken exactly as input
        clpdb = pdb
        clmtz = mtz
    else:   # will look for the output of rocc_phe-ref.py
        for file in os.listdir(mydir):
            if file.startswith('close_'+str(LIG)+'_'+str(ref_sof)+'_'):
                if file.endswith('.pdb'):
                    clpdb = file
                if file.endswith('.mtz'):
                    clmtz = file        
    rscc('.', folder_name, clpdb, clmtz)
else:
    for malaf in os.listdir(mydir):        # will look for pdb and mtz in each folder starting with the given directory
        if malaf.startswith(folder_name):
            for file in os.listdir(mydir+"/"+malaf):
                if options.mtz:    # will look for pdb and mtz as given, like pattern
                    if file.startswith(pdb) and file.endswith('.pdb'):
                        clpdb = file
                    if file.startswith(mtz) and file.endswith('.mtz'):
                        clmtz = file
                else:     # will look for the output of rocc_phe-ref.py in each folder
                    if file.startswith('close_'+str(LIG)+'_'+str(ref_sof)):
                        if file.endswith('.pdb'):
                            clpdb = file
                        if file.endswith('.mtz'):
                            clmtz = file
            rscc(malaf, malaf, clpdb, clmtz)
