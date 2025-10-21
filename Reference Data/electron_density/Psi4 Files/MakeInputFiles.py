def make_input_file(element_name, multiplicity, which_orbitals, axis):
    file_name = element_name + "/" + axis + "/" + element_name + ".in"
    file_id = open(file_name, "w")
    
    file_id.write("molecule " + element_name +  " {\n")
    file_id.write("0 " + str(multiplicity) + " \n")
    file_id.write(element_name + " 0.0 0.0 0.0 \n")
    file_id.write("}\n\n")
    
    file_id.write("set reference uks\n")
    file_id.write("set basis aug-cc-pvtz\n")
    file_id.write("set cubeprop_tasks ['orbitals']\n")
    
    if axis == "x":
        file_id.write("set cubic_grid_spacing [0.01, 0.0, 0.0] \n")
        file_id.write("set cubic_grid_overage [20, 0, 0] \n")

    if axis == "y":
        file_id.write("set cubic_grid_spacing [0.0, 0.01, 0.0] \n")
        file_id.write("set cubic_grid_overage [0, 20, 0] \n")
        
    if axis == "z":
        file_id.write("set cubic_grid_spacing [0.0, 0.0, 0.01] \n")
        file_id.write("set cubic_grid_overage [0, 0, 20] \n")
    
    file_id.write("set cubeprop_orbitals [")
    for orbital in which_orbitals:
        if orbital == which_orbitals[-1]:
            file_id.write(str(orbital))
        else:
            file_id.write(str(orbital)+",")
    file_id.write("]\n\n")
    
    file_id.write("E, wfn = energy('b3lyp', return_wfn=True)\n")
    file_id.write("cubeprop(wfn)\n\n")
    file_id.close()
    

elements = [
    "H", "He",
    "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar"
]

multiplicities = [
    2, 1,
    2, 1, 2, 3, 4, 3, 2, 1,
    2, 1, 2, 3, 4, 3, 2, 1,
]

for i in range(18):
    N = i + 1
    n_paired_orbitals = int((N - multiplicities[i] + 1) / 2)
    n_alpha_orbitals = n_paired_orbitals + multiplicities[i] - 1;
    n_beta_orbitals = n_paired_orbitals
    
    which_orbitals = []
    for n in range(n_alpha_orbitals):
        which_orbitals.append(n+1)
        
    for n in range(n_beta_orbitals):
        which_orbitals.append(-n-1)
    
    element = elements[i]
    multiplicity = multiplicities[i]
    make_input_file(element, multiplicity, which_orbitals, "x")
    make_input_file(element, multiplicity, which_orbitals, "y")
    make_input_file(element, multiplicity, which_orbitals, "z")
