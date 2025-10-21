#!/bin/bash

rm -rf Basis\ Set
mkdir Basis\ Set

for el in H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar
do
    mkdir Basis\ Set/${el}
    cp ../Reference\ Data/electron_density/Psi4\ Files/${el}/*.txt Basis\ Set/${el}/
done
