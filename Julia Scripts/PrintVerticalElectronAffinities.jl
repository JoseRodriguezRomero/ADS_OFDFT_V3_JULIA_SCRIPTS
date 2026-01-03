include("FitCoeffs_General.jl")

function print_table()
    println("\\begin{table}[b]");
    println("\\begin{ruledtabular}");
    println("\\begin{tabular}{lll}");
    println("& \\multicolumn{1}{c}{KS-DFT} & ");
    println("  \\multicolumn{1}{c}{ADS OF-DFT} \\\\ \\hline \\\\[-0.275cm]");
    
    atomic_numbers = [1,6,7,8];

    hartree_to_ev = 27.2114;
    dft_neutral_energies, dft_cation_energies, dft_anion_energies = 
        get_reference_atom_total_energy();

    for atomic_number in atomic_numbers
        elem_label = get_element_symbol(atomic_number);
        print("    \$\\text{");
        print(elem_label);

        at_neutral = make_monoatomic_system(atomic_number,0);
        at_anion = make_monoatomic_system(atomic_number,-1);

        dft_e_neutral = dft_neutral_energies[atomic_number];
        dft_e_anion = dft_anion_energies[atomic_number];

        model_vea = total_energy(at_neutral) - total_energy(at_anion);
        dft_vea = dft_e_neutral - dft_e_anion;

        model_vea *= hartree_to_ev;
        dft_vea *= hartree_to_ev;

        print(@sprintf "}\$   & ");
        print(@sprintf "\\multicolumn{1}{r}{%8.4lf} & " dft_vea);
        print(@sprintf "\n        ");
        print(@sprintf "\\multicolumn{1}{r}{%8.4lf} \\\\[0.1cm]" model_vea);
        print("\n");
    end

    for atomic_number in atomic_numbers
        elem_label = get_element_symbol(atomic_number);
        print("    \$\\text{");
        print(elem_label);

        at_cation = make_monoatomic_system(atomic_number,1);
        at_neutral = make_monoatomic_system(atomic_number,0);

        dft_e_cation = dft_cation_energies[atomic_number];
        dft_e_neutral = dft_neutral_energies[atomic_number];

        model_vea = total_energy(at_cation) - total_energy(at_neutral);
        dft_vea =  dft_e_cation - dft_e_neutral;

        model_vea *= hartree_to_ev;
        dft_vea *= hartree_to_ev;

        print(@sprintf "}^+\$ & ");
        print(@sprintf "\\multicolumn{1}{r}{%8.4lf} & " dft_vea);
        print(@sprintf "\n        ");
        print(@sprintf "\\multicolumn{1}{r}{%8.4lf} \\\\[0.1cm]" model_vea);
        print("\n");
    end

    println("\\end{tabular}");
    println("\\end{ruledtabular}");
    println("\\end{table}");      
end

print_table();
