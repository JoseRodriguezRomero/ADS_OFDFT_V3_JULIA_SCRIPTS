#!/bin/bash

for el in H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar
do
    cd $el
    for axis in x y z
    do
        cd $axis
        rm *.cube
        rm *.log
        rm *.xyz
        rm *.dat
        rm *.out
        psi4 $el.in
        cd ..
    done
    
    echo "done with $el"
    cd ..
done
