import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt


def Norm3D(x,y,z):
    return ((x**2.0) + (y**2.0) + (z**2.0))**0.5


def ReadOutputFile(file_name):
    # Reads the Psi4 output file of the diatomic system and returns the bond
    # lenght, energy and a boolean value of whether the output file exists
    # and/or is valid.
    is_valid = False
    total_energy = 0.0
    
    if not Path(file_name).exists():
        return is_valid, total_energy
    
    fileID = open(file_name, "r")
    lines = fileID.readlines()
    
    alpha_homo = float("nan");
    alpha_lumo = float("nan");
    
    beta_homo = float("nan");
    beta_lumo = float("nan");
    
    for i in range(len(lines)):
        if "*** Psi4 exiting successfully." in lines[i]:
            is_valid = True
    
        if "Total Energy =" in lines[i]:
            line_splitted = lines[i].split()
            total_energy = float(line_splitted[-1])
            
        if "Alpha Virtual:" in lines[i]:
            line_splitted = lines[i-2].split()
            if len(line_splitted) > 0:
                alpha_homo = float(line_splitted[-1])
            
            line_splitted = lines[i+2].split()
            alpha_lumo = float(line_splitted[1])
            
        if "Beta Virtual:" in lines[i]:
            line_splitted = lines[i-2].split()
            if len(line_splitted) > 0:
                beta_homo = float(line_splitted[-1])
                
            line_splitted = lines[i+2].split()
            beta_lumo = float(line_splitted[1])
            
    fileID.close()
    
    if np.isnan(beta_homo):
        beta_homo = alpha_homo
        beta_lumo = alpha_lumo
    
    homo = max(alpha_homo,beta_homo)
    lumo = min(alpha_lumo,beta_lumo)
    
    return is_valid, total_energy, homo, lumo


elements = ["H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar"]
multiplicities = ["singlet","doublet","triplet","quartet","quintet","sextet"]
all_energies = []

fileID_1 = open("elem_total_energy.txt", "w")
fileID_2 = open("elem_orb_energies.txt", "w")
for element in elements:
    curr_total_energy = 1000.0;
    curr_homo_energy = 0.0;
    curr_lumo_energy = 0.0;
    which_mult = ""
    for multiplicity in multiplicities:
        file_name = multiplicity + "_" + element + ".out"
        file_name = "Psi4 Files/" + file_name
        is_valid, total_energy, homo, lumo = ReadOutputFile(file_name)
        
        if is_valid:
            if curr_total_energy > total_energy:
                which_mult = multiplicity
                curr_total_energy = total_energy
                curr_homo_energy = homo
                curr_lumo_energy = lumo
        
    
    print([element, which_mult])
    all_energies.append(curr_total_energy)
    fileID_1.write('{:.8e}'.format(curr_total_energy) + "\n")
    
    fileID_2.write('{:20.8e}'.format(curr_homo_energy))
    fileID_2.write('{:20.8e}'.format(curr_lumo_energy))
    fileID_2.write("\n")

fileID_1.close()
fileID_2.close()

plt.plot(np.linspace(1,18,18), all_energies,'.-')
plt.show()
