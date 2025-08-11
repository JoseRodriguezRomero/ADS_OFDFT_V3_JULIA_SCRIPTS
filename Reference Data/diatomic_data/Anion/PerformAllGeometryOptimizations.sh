#!/bin/bash

elements=(H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar)

for i in {1..18}; do
    for (( j = i; j <= 18; j++ )); do
        sum=$((i+j+1))
        sum_rem=$((sum%2))
        e1="${elements[i-1]}"
        e2="${elements[j-1]}"
        
        echo "Optimizing $e1 + $e2"
        python MakeUnoptimizedInitialInput.py $e1 $e2 $sum_rem
        
        if ((sum_rem == 0)); then
            psi4 "GeometryOptimizations/singlet_$e1$e2.in"
            psi4 "GeometryOptimizations/triplet_$e1$e2.in"
            psi4 "GeometryOptimizations/quintet_$e1$e2.in"
        fi
        
        if ((sum_rem == 1)); then
            psi4 "GeometryOptimizations/doublet_$e1$e2.in"
            psi4 "GeometryOptimizations/quartet_$e1$e2.in"
            psi4 "GeometryOptimizations/sextet_$e1$e2.in"
        fi
  done
done
