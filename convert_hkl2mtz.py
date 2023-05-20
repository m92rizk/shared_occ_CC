#!/usr/bin/python
# rizk@esrf.fr
# converts hkl to mtz and adds a freeR flag to 5%
import os, sys, re, shutil, copy
import numpy as np
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-d", "--directory", dest="directory", default=".", help="by default it takes current current directory")
parser.add_option("-f", "--filename", dest="hklfile", default="merged.hkl", help="looking for merged.hkl file by default")
parser.add_option("-u", "--unit_cell", dest="unit", default="", help="unit cell parameters, by default it will be guessed by the program")
parser.add_option("-s", "--space_group", dest="sg", default="", help="spacegroup, by default it will be guessed by the program")
parser.add_option("", "--friedel_false", action="store_false", dest="friedel", default=True)
(options, args) = parser.parse_args()

mydir = os.getcwd()

if options.unit == "" :
    unit=""
else:
    unit="UNIT_CELL_CONSTANTS=    "+str(options.unit)+"\n"

if options.sg == "" :
    sg=""
else:
    sg="SPACE_GROUP_NUMBER=    "+str(options.sg)+"\n"

if options.friedel:
    Fr_law = "TRUE"
    labelIN = "E1=FP E2=SIGFP E3=FreeRflag"
    labelOUT = "E1=FP E2=SIGFP E3=FreeRflag"
else:
    Fr_law = "FALSE"
    labelIN = "E1=FP E2=SIGFP E3=DANO E4=SIGDANO E5=ISYM E6=FreeRflag"
    labelOUT = "E1=FP E2=SIGFP E3=DANO_sulf E4=SIGDANO_sulf E5=ISYM_sulf E6=FreeRflag"


def convert():
    if os.path.isfile(str(options.hklfile)):
        hklfile = str(options.hklfile)
        print("Found "+str(hklfile))
    else:
        print("Couldn't find any hkl file")
        return 0
    with open('XDSCONV.INP', 'w') as outxc:
        outxc.write("INPUT_FILE="+str(hklfile)+"\n"+str(unit)+str(sg)+"INCLUDE_RESOLUTION_RANGE=50 1\nOUTPUT_FILE=temp.hkl  CCP4\nFRIEDEL'S_LAW="+str(Fr_law)+"\nGENERATE_FRACTION_OF_TEST_REFLECTIONS=0.05")
    os.system('xdsconv >/dev/null')
    os.system('f2mtz HKLOUT temp.mtz<F2MTZ.INP >/dev/null')
    os.system('cad HKLIN1 temp.mtz HKLOUT junk_xdsconv.mtz<<EOF >/dev/null\nLABIN FILE 1 '+str(labelIN)+'\nLABOUT FILE 1 '+str(labelOUT)+'\nEND\nEOF ')
    

if options.directory == "." :
    convert()
    print("done converting hkl to mtz here")
else:
    for folder in os.listdir(mydir):
        if folder.startswith(str(options.directory)):
            os.chdir(str(mydir)+'/'+str(folder))
            convert()
            print("done converting hkl to mtz in "+str(folder))
    os.chdir(str(mydir))

print("you may now run get_phases_v4.py")