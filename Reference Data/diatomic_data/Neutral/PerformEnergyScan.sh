#!/bin/bash

sh PerformVibrationalCalc.sh $1 $2

element1=$1
element2=$2

rm "EnergyScans/scan_${element1}${element2}_"*
python MakeScanInputData.py $element1 $element2
        
echo "scanning $1 + $2"

for ((i=25; i>0; i--)); do
    echo $i
    psi4 "EnergyScans/scan_${1}${2}_${i}.in"
done

for ((i=25; i<51; i++)); do
    echo $i
    psi4 "EnergyScans/scan_${element1}${element2}_${i}.in"
done

python ParseEnergyScan.py $1 $2
