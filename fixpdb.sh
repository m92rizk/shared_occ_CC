#!/bin/bash


if [ -z "${1}" ] ;
then
    read -p "Enter pdb name: " pdb
else
    pdb="${1}"
fi
if [ -z "${2}" ] ;
then
    read -p "Enter residue name: " res
else
    res="${2}"
fi

A=`grep  "${res^^}" $pdb.pdb | grep -E 'HETATM'\|'ATOM' | cut  -b 22-22 | head -1`

sed -i "s/${res^^} $A/${res^^} Z/g" $pdb.pdb

grep -v ANISOU $pdb.pdb > renamethis.pdb
mv renamethis.pdb $pdb.pdb
