import sys
from pathlib import Path


def Norm3D(x,y,z):
    return ((x**2.0) + (y**2.0) + (z**2.0))**0.5


def ReadOutputFile(file_name):
    # Reads the Psi4 output file of the diatomic system and returns the bond
    # lenght, energy and a boolean value of whether the output file exists
    # and/or is valid.
    bond_length = 0.0
    energy = 0.0
    is_valid = False
    
    if not Path(file_name).exists():
        return is_valid, bond_length, energy
    
    fileID = open(file_name, "r")
    lines = fileID.readlines()
    
    for i in range(len(lines)):
        if "Psi4 exiting successfully." in lines[i]:
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
    
    fileID.close()
    return is_valid, bond_length, energy


def MakeInputFile(elem1, elem2, bond_length, multiplicity, file_name):
    # Creates a Psi4 input file of a point-energy calculation of the diatomic
    # system with the interatomic distance being bond_length in angstrom.
    fileID = open(file_name, "w")
    fileID.write("molecule { \n")
    fileID.write("0 " + str(multiplicity) + "\n")
    
    format_string = "{element:} {x:18.6f} {y:18.6f} {z:18.6f}\n"
    fileID.write(format_string.format(element=elem2, x=0.0, y=0.0, z=0.0))
    fileID.write(format_string.format(element=elem1, x=0.0, y=0.0, z=bond_length))
    fileID.write("}\n\n")
        
    fileID.write("set { \n")
    fileID.write("basis              cc-pvtz\n")
    fileID.write("reference          uks\n")
    fileID.write("maxiter            500\n")
    fileID.write("DAMPING_PERCENTAGE 1\n")
    fileID.write("E_CONVERGENCE      1.0E-8\n")
    fileID.write("D_CONVERGENCE      1.0E-8\n")
    fileID.write("DIIS_MAX_VECS      400\n")
    fileID.write("MBIS_MAXITER       8000\n")
    
    fileID.write("guess              SAP\n")
    fileID.write("}\n\n")
    fileID.write("E, wfn = optimize('b3lyp', return_wfn=True)\n")
    fileID.write("frequencies('b3lyp', wfn=wfn, return_wfn=True)\n")
    fileID.write("properties('b3lyp','MULLIKEN_CHARGES', wfn=wfn)\n")
    
    fileID.close()
    return


output_files = [
    "GeometryOptimizations/singlet_" + sys.argv[1] + sys.argv[2] + ".out",
    "GeometryOptimizations/doublet_" + sys.argv[1] + sys.argv[2] + ".out",
    "GeometryOptimizations/triplet_" + sys.argv[1] + sys.argv[2] + ".out",
    "GeometryOptimizations/quartet_" + sys.argv[1] + sys.argv[2] + ".out",
    "GeometryOptimizations/quintet_" + sys.argv[1] + sys.argv[2] + ".out",
    "GeometryOptimizations/sextet_" + sys.argv[1] + sys.argv[2] + ".out"
]

file_is_valid, bond_length, old_energy = ReadOutputFile(output_files[0])
old_energy = old_energy + 1.0

for i in range(len(output_files)):
    print(output_files[i])
    multiplicity = i + 1
    output_file = output_files[i]
    file_is_valid, bond_length, energy = ReadOutputFile(output_file)
    
    if file_is_valid and (energy < old_energy):
        old_energy = energy
        file_name = "VibrationalData/scan_" + sys.argv[1] + sys.argv[2] + ".in"
        MakeInputFile(sys.argv[1], sys.argv[2], bond_length, multiplicity, file_name)
