#!/bin/bash

elements=("H" "He" "Li" "Be" "B" "C" "N" "O" "F" "Ne" \
          "Na" "Mg" "Al" "Si" "P" "S" "Cl" "Ar")
num_elements=${#elements[@]}

multiplicities=("singlet" "doublet" "triplet" "quartet"\
                "quintet" "sextet")
num_multiplicities=${#multiplicities[@]}

# Loop through the indices and print the elements
for ((i = 0; i < num_elements; i++)); do
    python MakeInputFile.py $i
    
    for ((j = 0; j < num_multiplicities; j++)); do
        # echo "${multiplicities[j]}_${elements[i]}.in"
        psi4 "Psi4 Files/${multiplicities[j]}_${elements[i]}.in"
    done
done
