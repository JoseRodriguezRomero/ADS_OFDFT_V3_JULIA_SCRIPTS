#!/bin/bash

element1=$1
element2=$2

rm "VibrationalData/scan_${element1}${element2}"*
python MakeVibrationalInputData.py $element1 $element2
        
psi4 "VibrationalData/scan_${element1}${element2}.in"
