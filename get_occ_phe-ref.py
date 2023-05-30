#!/usr/bin/python3.8
# 92mrizk@gmail.com
#SBATCH -p low
#SBATCH -N 1
#SBATCH -n 10
#SBATCH -t 0:20:0
#SBATCH -J rzk_OCC
#SBATCH --array=1,3,5,7,9

import os, sys, re, shutil, copy
import numpy as np
from pathlib import Path
from optparse import OptionParser

mydir = os.popen("pwd").read().strip()

parser = OptionParser()
parser.add_option("-l", "--ligand", dest="lig", help="the 3 letters name of the residue")
parser.add_option("-p", "--pdb", dest="pdb", help="the pdb file, {model.pdb}")
parser.add_option("-c", "--cif", dest="cif", help="the cif file")
parser.add_option("-n", "--ncyc", dest="n", default=5)
parser.add_option("-o", "--ncyc_occ", dest="occn", default=10, help="number of occupancy cycles")
parser.add_option("-d", "--directory", dest="directory", default=".", help="leaving this empty will take current directory")
parser.add_option("-f", "--hkl", dest="hkl", default="", help="the reflection hkl file")
#parser.add_option("", "--hkl_pattern", dest="hkl_pattern", default="", help="the reflection hkl file name pattern, by default {merged.hkl}")
parser.add_option("-m", "--mtz", dest="mtz", default="junk_xdsconv.mtz", help="mtz file {MTZfile.mtz}")
#parser.add_option("", "--mtz_pattern", dest="mtz_pattern", default="junk_xdsconv.mtz", help="mtz file name pattern")
parser.add_option("", "--output_result", dest="output_result")
parser.add_option("", "--output_dir", dest="output_dir")
parser.add_option("", "--phenix", action="store_true", dest="phenix", default=False)
parser.add_option("", "--refmac", action="store_true", dest="refmac", default=False)
parser.add_option("", "--slurm", action="store_true", dest="slurm", default=False)
parser.add_option("", "--overwrite", action="store_true", dest="overwrite", default=False)
parser.add_option("", "--fast", action="store_true", dest="fast", default=False)

(options, args) = parser.parse_args()

if options.slurm:
    arguments = []
    arg_length = len(sys.argv)
    for x in range(0, arg_length):
        if sys.argv[x] != "--slurm":
            arguments.append(sys.argv[x])
    command_line = ' '.join(map(str, arguments))
    print("running this: sbatch "+str(command_line))
    os.system("sbatch  "+command_line)
    quit()


if options.fast:
    occupancies = [5]
else:
    occupancies = [1,5,9]
    
output_dir = options.output_dir
if options.refmac and options.phenix:
    print("either refmac or phenix at a time, quitting")
    quit()
elif options.refmac:
    ref_sof = "refmac"
elif options.phenix:
    ref_sof = "phenix"
else:
    print("either refmac or phenix should be selected, quitting")
    quit()
    
if options.lig is None or options.pdb is None:
    print("the pdb and ligand are needed! \nquitting.")
    quit()
LIG = options.lig.upper()
if options.cif:
    cifpath = os.path.abspath(options.cif)
    cif = ' libin '+str(cifpath)
else:
    cif = " "
    cifpath = " "
if str(options.hkl) != "" :
    hkl="-f "+str(options.hkl)
else:
    hkl=""
jobid = os.getenv('SLURM_ARRAY_TASK_ID')
if options.slurm:
    print ("\n\nTHE JOBID IS: "+str(jobid)+" \n\n")
    if jobid == None:
        quit()

pdb = options.pdb    
pdbname = Path(pdb).stem

    
if options.hkl:
    hkl = options.hkl
else:
    hkl = ""

MTZ = options.mtz
if hkl == "" and MTZ == "junk_xdsconv.mtz":
    print("neither a mtz nor a hkl is provided, quitting")
    quit()

file = open(pdb)
for line in file:
    if (re.search(str(LIG), line) and re.search("HETATM", line)) or (re.search(str(LIG), line) and re.search("ATOM", line)):
        res=line[22:26]
        ch=line[21]
        break

def modify_occ(jobid, occ):
    with open('f_'+str(jobid)+'.eff', 'w') as outeff:
        outeff.write('modify {\noccupancies {\n\tatom_selection=resname '+str(LIG)+'\n\tset='+str(occ)+'\n}\n}\noutput {\n\tsuffix = "_modified_'+str(jobid)+'"\n}')
    os.system('phenix.pdbtools '+str(pdb)+' f_'+str(jobid)+'.eff | grep "Setting occupancies to"')
    print("Changed occupancy to "+str(occ))
def get_best_occ(occ,output_dir):
    OCC=[]
    print("Done refining, xyzout: "+str(LIG)+"_"+str(ref_sof)+"_"+str(occ)+".pdb")
    file = open(str(folder)+'/'+str(LIG)+'_'+str(ref_sof)+'_'+str(occ)+'.pdb', "r")
    length = 0
    for line in file:
        if (re.search(str(LIG), line) and re.search("HETATM", line)) or (re.search(str(LIG), line) and re.search("ATOM", line)):
            val = float(line[56:60]) #refined occ
            OCC.append(val)
    print("individual atomic occ: "+str(OCC))
    avg = round(sum(OCC)/len(OCC), 2)
    delta = abs(avg - float(occ))
    close_occ = float(avg)
    print("for "+str(LIG)+"_"+str(ref_sof)+"_"+str(occ)+".pdb: the average is "+str(avg))
    with open(details_output, 'a') as outf:
        outf.write(str(occ)+'\t\t'+str(avg)+'\t\t\t'+str(delta)+'\n')
    #modifying the accurate occupancy from each sub-run
    with open(details_output, 'r') as dfile:
        data = dfile.readlines()
    with open(occlist, 'r') as ofile:
        lst = ofile.readlines()
    arr = [] # to store the stats in new array without spaces and tabs
    for i in range (3,len(data)):
        arr.append(data[i].split())
    acc_occ = sorted(arr,key=lambda x: x[2])[0][1]
    cl_occ = sorted(arr,key=lambda x: x[2])[0][0]
    for clfile in os.listdir(str(folder)):
        if clfile.startswith('close_'+str(LIG)+'_'+str(ref_sof)+'_'):
            os.remove(str(folder)+'/'+clfile)
    for clfile in os.listdir(str(output_dir)):
        if clfile.startswith('close_'+str(LIG)+'_'+str(ref_sof)+'_'):
            os.remove(str(output_dir)+'/'+clfile)            
    for file in os.listdir(str(folder)):
        if file.endswith(str(cl_occ)+'.mtz'):
            shutil.copyfile(str(folder)+'/'+str(LIG)+'_'+str(ref_sof)+'_'+str(cl_occ)+'.mtz', str(folder)+'/close_'+str(LIG)+'_'+str(ref_sof)+'_'+str(cl_occ)+'.mtz')
        if file.endswith(str(cl_occ)+'.pdb'):
            shutil.copyfile(str(folder)+'/'+str(LIG)+'_'+str(ref_sof)+'_'+str(cl_occ)+'.pdb', str(folder)+'/close_'+str(LIG)+'_'+str(ref_sof)+'_'+str(cl_occ)+'.pdb')
    try:
        shutil.copyfile(str(folder)+'/close_'+str(LIG)+'_'+str(ref_sof)+'_'+str(cl_occ)+'.pdb', str(output_dir)+'/close_'+str(LIG)+'_'+str(ref_sof)+'_'+str(cl_occ)+'.pdb')
        shutil.copyfile(str(folder)+'/close_'+str(LIG)+'_'+str(ref_sof)+'_'+str(cl_occ)+'.mtz', str(output_dir)+'/close_'+str(LIG)+'_'+str(ref_sof)+'_'+str(cl_occ)+'.mtz')
    except shutil.SameFileError:
        pass
        
    data[1] = 'Final average occupancy (more accurate): '+str(acc_occ)+'\n'
    with open(details_output, 'w') as file:
        file.writelines( data )
    edited = False
    for index, line in enumerate(lst):
        if line.split()[0] == str(folder_name):
            lst[index] = str(folder_name)+'\t\t'+str(acc_occ)+'\n'
            edited = True
            break
    if not edited:
        lst.append(str(folder_name)+'\t\t'+str(acc_occ)+'\n')
    with open(occlist, 'w') as outlist:
        outlist.writelines( lst )
def refmac_refine(folder, output_dir, jobid, occ):
    os.system('refmac5 \
    hklin '+str(folder)+'/'+str(options.mtz)+' \
    xyzin '+str(pdbname)+'_modified_'+str(jobid)+'.pdb \
    '+str(cif)+' \
    xyzout '+str(folder)+'/'+str(LIG)+'_refmac_'+str(occ)+'.pdb \
    hklout '+str(folder)+'/'+str(LIG)+'_refmac_'+str(occ)+'.mtz <<eor >> '+str(folder)+'/OUT_refmac_'+str(jobid)+'.log \n \
    ncyc '+str(options.n)+'\n\
    occupancy group id 1 chain '+str(ch)+' residue '+str(res)+'\n\
    occupancy refine ncycle '+str(options.occn)+'\n\
    occupancy refine\n\
    eor')
    #shutil.copyfile(str(folder)+'/'+str(LIG)+'_refmac_'+str(occ)+'.mtz', str(output_dir)+'/'+str(LIG)+'_refmac_'+str(occ)+'.mtz')
    #shutil.copyfile(str(folder)+'/'+str(LIG)+'_refmac_'+str(occ)+'.pdb', str(output_dir)+'/'+str(LIG)+'_refmac_'+str(occ)+'.pdb')
def phenix_refine(folder, output_dir, jobid, occ):
    if os.path.isfile('params'+str(LIG)+'.eff') :
        print('params'+str(LIG)+' file already exists')
    else:
        os.system('phenix.refine --show_defaults=0 > params'+str(LIG)+'.eff')
        os.system('sed -zi "s/individual = None/individual = resname '+str(LIG)+'/2"/" params'+str(LIG)+'.eff')
    os.system('phenix.refine \
    '+str(folder)+'/'+str(options.mtz)+' \
    '+str(pdbname)+'_modified_'+str(jobid)+'.pdb \
    '+str(cifpath)+' \
    params'+str(LIG)+'.eff \
    main.number_of_macro_cycles='+str(options.n)+' \
    main.nproc=16 \
    write_geo_file=False \
    write_def_file=False \
    write_model_cif_file=False \
    base_output_dir="folder" \
    --overwrite >> '+str(folder)+'/OUT_phenix_'+str(jobid)+'.log')
    #shutil.copyfile(str(folder)+'/'+str(pdbname)+'_modified_'+str(jobid)+'_refine_001.mtz', str(output_dir)+'/'+str(LIG)+'_phenix_'+str(occ)+'.mtz')
    #shutil.copyfile(str(folder)+'/'+str(pdbname)+'_modified_'+str(jobid)+'_refine_001.pdb', str(output_dir)+'/'+str(LIG)+'_phenix_'+str(occ)+'.pdb')
def refine_occupancy(folder, folder_name, output_dir):
    OCC = []
    if options.overwrite:
        try:
            print('overwriting '+str(folder)+'/junk_xdsconv.mtz'+' ..')
            os.remove(str(folder)+'/junk_xdsconv.mtz')
        except FileNotFoundError:
            pass
    if os.path.isfile(str(folder)+'/'+str(MTZ)):
        print('mtz file already exists')
    else:    
        try:
            os.system('convert_hkl2mtz.py -d '+str(folder)+' '+str(hkl))
        except:
            print("====================================\n\
            something wrong happened while converting hkl into mtz! maybe try running it manually using: \n\
            convert_hkl2mtz.py\n\
            ====================================\n")
    
    if options.slurm:
        occ = "0."+str(jobid)
        modify_occ(jobid, occ)
        #refmac refinement
        if options.refmac:
            refmac_refine(folder,output_dir,jobid,occ)
        #phenix refinement
        if options.phenix:
            phenix_refine(folder,output_dir,jobid,occ)
        get_best_occ(occ,output_dir)
    else:
        for job in occupancies:
            occ = "0."+str(job)
            modify_occ(job, occ)
            #refmac refinement
            if options.refmac:
                refmac_refine(folder,output_dir,job,occ)
            #phenix refinement
            if options.phenix:
                phenix_refine(folder,output_dir,job,occ)
        get_best_occ(occ,output_dir)
    shutil.copy
    
    
    
    
    
    
    
    
    
    
    
folder = options.directory
if folder == "." or folder == "./" :
    folder_name = os.path.basename(os.getcwd()).strip()
elif folder.endswith('/'):
    folder_name = os.path.basename(folder).strip('/')
else:
    folder_name = os.path.basename(folder)

occlist = 'occ'+str(LIG)+'_list_'+str(folder_name)+'_'+str(ref_sof)+'.txt'
if options.output_result:
    details_output = options.output_result
else:
    details_output = str(folder)+'/details_occ_'+str(LIG)+'_'+str(ref_sof)+str(folder_name)+'.txt'

with open(occlist, 'w') as outlist:
   outlist.write('Directory\tOccupancy Refined\n')
with open(details_output, 'w') as outf:
    outf.write(str(ref_sof)+'\t occupancy refinement\nFinal average occupancy (more accurate): 0\n')
    outf.write('starting_occ\tavg_refined_occ\t\tdifference\n')

if folder == "." or folder == "./" or folder.endswith("/"):
    folder = folder.strip('/')
    if output_dir:
        outdir = output_dir
    else:
        outdir = mydir
    refine_occupancy(folder, folder_name, outdir)
else:
    for x in sorted(os.listdir(mydir)):
        if x.startswith(str(folder)):
            xdir = mydir+'/'+x
            with open(details_output, 'a') as outf:
                outf.write('*** Working Dir: '+str(x)+' ***\n')
            refine_occupancy(x, folder_name, xdir)

