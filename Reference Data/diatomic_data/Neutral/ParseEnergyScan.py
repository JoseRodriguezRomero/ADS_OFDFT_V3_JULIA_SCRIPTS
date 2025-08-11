import sys
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt


def Norm3D(x,y,z):
    return ((x**2.0) + (y**2.0) + (z**2.0))**0.5


def ReadOrbEnergies(lines, i):
    # Auxiliary function to read table blocks with orbital energies in a
    # Psi4 output file. Can be used for both Alpha and Beta orbitals,
    # occupied or unnocupied.
    orb_energies = []
            
    j = i + 1
    while True:
        j = j + 1
        line_splitted = lines[j].split()
                
        if len(line_splitted) < 2:
            break
                
        for k in range(int(len(line_splitted)/2)):
            orb_energies.append(float(line_splitted[2*k+1]))
            
    return orb_energies
    

def ReadOutputFile(file_name):
    # Reads the Psi4 output file of the diatomic system and returns the bond
    # lenght, energy and a boolean value of whether the output file exists
    # and/or is valid.
    bond_length = 0.0
    energy = 0.0
    is_valid = False
    q1 = 0.0;
    q2 = 0.0;
    
    alpha_occ_orb_energies = []
    alpha_ucc_orb_energies = []
    
    beta_occ_orb_energies = []
    beta_ucc_orb_energies = []
    
    if not Path(file_name).exists():
        return is_valid, bond_length, energy, q1 ,q2
    
    fileID = open(file_name, "r")
    lines = fileID.readlines()
    
    for i in range(len(lines)):
        if "*** Psi4 exiting successfully." in lines[i]:
            is_valid = True
    
        if "Center              X                  Y" in lines[i]:
            at1_line = lines[i+2].split()
            at2_line = lines[i+3].split()
            
            dx = float(at1_line[1]) - float(at2_line[1])
            dy = float(at1_line[2]) - float(at2_line[2])
            dz = float(at1_line[3]) - float(at2_line[3])
            bond_length = Norm3D(dx, dy, dz)
            
        if "Total Energy =" in lines[i]:
            energy = float(lines[i].split()[-1])
            
        if "Mulliken Charges: (a.u.)" in lines[i]:
            at1_line = lines[i+2].split()
            at2_line = lines[i+3].split()
            q1 = float(at1_line[-1])
            q2 = float(at2_line[-1])
            
        if "MBIS Charges: (a.u.)" in lines[i]:
            at1_line = lines[i+2].split()
            at2_line = lines[i+3].split()
            q1 = float(at1_line[-1])
            q2 = float(at2_line[-1])
                
        if "Alpha Occupied:" in lines[i]:
            alpha_occ_orb_energies = ReadOrbEnergies(lines,i)
            
        if "Alpha Virtual:" in lines[i]:
            alpha_ucc_orb_energies = ReadOrbEnergies(lines,i)
            
        if "Beta Occupied:" in lines[i]:
            beta_occ_orb_energies = ReadOrbEnergies(lines,i)
            
        if "Beta Virtual:" in lines[i]:
            beta_ucc_orb_energies = ReadOrbEnergies(lines,i)
                    
    fileID.close()
    
    if (len(beta_occ_orb_energies) == 0):
        e_homo = max(alpha_occ_orb_energies)
    else:
        e_homo = max([max(alpha_occ_orb_energies), max(beta_occ_orb_energies)])
    
    if (len(beta_ucc_orb_energies) == 0):
        e_lumo = min(alpha_ucc_orb_energies)
    else:
        e_lumo = min([min(alpha_ucc_orb_energies), min(beta_ucc_orb_energies)])
    
    return is_valid, bond_length, energy, q1, q2, e_homo, e_lumo

elements = ["H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na",
    "Mg", "Al", "Si", "P", "S", "Cl", "Ar"]

elem1 = sys.argv[1]
elem2 = sys.argv[2]
        
Q1 = []
Q2 = []
energies = []
distances = []
all_e_homo = []
all_e_lumo = []

for i in range(50):
    file_name = "EnergyScans/scan_" + elem1 + elem2 + "_" + str(i+1) + ".out"
    is_valid, bond_length, energy, q1, q2, e_homo, e_lumo = ReadOutputFile(file_name)
    
    if is_valid:
        energies.append(energy)
        distances.append(bond_length)
        Q1.append(q1)
        Q2.append(q2)
        all_e_homo.append(e_homo)
        all_e_lumo.append(e_lumo)

fileID = open("ParsedData/scan_" + elem1  + elem2 + ".txt", "w")

# Print the column labels
fileID.write("{:>22}".format("separation"))
fileID.write("{:>22}".format("total energy"))
fileID.write("{:>22}".format("partial charge"))
fileID.write("{:>22}".format("partial charge"))
fileID.write("{:>22}".format("HOMO"))
fileID.write("{:>22}".format("LUMO"))
fileID.write("\n")

# Print the column units
fileID.write("{:>22}".format("[Angstrom]"))
fileID.write("{:>22}".format("[Hartree]"))
fileID.write("{:>22}".format("[" + elem1 + " atom]"))
fileID.write("{:>22}".format("[" + elem2 + " atom]"))
fileID.write("{:>22}".format("[Hartree]"))
fileID.write("{:>22}".format("[Hartree]"))
fileID.write("\n")

for i in range(len(energies)):
    fileID.write("{:>22.6f}".format(distances[i]))
    fileID.write("{:>22.6f}".format(energies[i]))
    fileID.write("{:>22.6f}".format(Q1[i]))
    fileID.write("{:>22.6f}".format(Q2[i]))
    fileID.write("{:>22.6f}".format(all_e_homo[i]))
    fileID.write("{:>22.6f}".format(all_e_lumo[i]))
    fileID.write("\n")
fileID.close()

