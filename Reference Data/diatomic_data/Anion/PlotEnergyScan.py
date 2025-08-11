import sys
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt


def Norm3D(x,y,z):
    return ((x**2.0) + (y**2.0) + (z**2.0))**0.5


def ReadOutputFile(file_name):
    # Reads the Psi4 output file of the diatomic system and returns the bond
    # lenght, energy and a boolean value of whether the output file exists
    # and/or is valid.
    bond_length = 0.0
    energy = 0.0
    is_valid = False
    q1 = 0.0;
    q2 = 0.0;
    
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
            
    
    fileID.close()
    return is_valid, bond_length, energy, q1, q2

Q1 = []
Q2 = []
energies = []
distances = []
for i in range(50):
    file_name = "EnergyScans/scan_" + sys.argv[1] + sys.argv[2] + "_" + str(i+1) + ".out"
    is_valid, bond_length, energy, q1, q2 = ReadOutputFile(file_name)
    
    if is_valid:
        Q1.append(q1)
        Q2.append(q2)
        energies.append(energy)
        distances.append(bond_length)

min_e = min(energies)
for i in range(len(energies)):
    energies[i] = (energies[i] - min_e) * (distances[i] ** 0.0)


num_coeffs = 2
coeffs = [0.0, -1.0, -2.5]

M = np.zeros((len(distances), len(coeffs)))
Y = np.zeros(len(distances))

for i in range(len(distances)):
    Y[i] = energies[i]
    for j in range(len(coeffs)):
        M[i,j] = distances[i] ** coeffs[j]
    
plt.subplot(2, 1, 1)
plt.plot(distances, energies, '.-')
plt.title(sys.argv[1] + " + " + sys.argv[2] + " (anion)")
plt.ylabel("Energy Difference [Hartree]")

plt.subplot(2, 1, 2)

plt.plot(distances, Q1, '.-')
plt.plot(distances, Q2, '.-')
plt.legend([sys.argv[1], sys.argv[2]])

plt.ylabel("Partial Charge")
plt.xlabel("Interatomic Separation [Angstrom]")

plt.tight_layout()

plt.show()
