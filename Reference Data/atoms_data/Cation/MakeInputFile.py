import sys


def MakeInputFile(element, multiplicity, file_name):
    # Writes a Psi4 input file for a monoatomic system with the specified
    # element and number of unpaired electrons.
    fileID = open(file_name, "w")
    fileID.write("molecule { \n")
    fileID.write("1 " + str(multiplicity) + "\n")
    
    format_string = "{element:} {x:18.6f} {y:18.6f} {z:18.6f}\n"
    fileID.write(format_string.format(element=element, x=0.0, y=0.0, z=0.0))
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
    fileID.write("energy('b3lyp')\n")
    fileID.close()
    
    return


elements = ["H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar"]
    
element = elements[int(sys.argv[1])]
multiplicities = ["singlet","doublet","triplet","quartet","quintet","sextet"]

for i in range(len(multiplicities)):
    file_name = multiplicities[i] + "_" + element + ".in"
    file_name = "Psi4 Files/" + file_name
    MakeInputFile(element, i+1, file_name)

