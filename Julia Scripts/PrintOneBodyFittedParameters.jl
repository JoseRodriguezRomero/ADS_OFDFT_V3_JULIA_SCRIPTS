using Printf;

include("FitCoeffs_General.jl")

function print_to_latex(number::Float64, digits::Int = 4)
    if number == 0
        return "\$0\$";
    end
    exp = floor(Int, log10(abs(number)));
    mant = round(number / 10.0^exp; digits = digits);
    mant_str = @sprintf("%.*f", digits, mant);

    ret_str = "\$ $(mant_str) \\times 10^{ $(exp) } \$";

    if number >= 0
        ret_str = "\\hphantom{\$-\$}"*ret_str;
    end

    return ret_str;
end

function print_tot_e_xc()
    atomic_numbers = [1,6,7,8];

    println("\\begin{tabular}{cllll}");
    println("   \\hline \\\\[-0.35cm] \\hline \\\\[-0.35cm]");
    println("   \$Z\$ & ");
    println("       \\multicolumn{1}{c}{\$A_\\text{1B} \\left( Z \\right)\$} & ");
    println("       \\multicolumn{1}{c}{\$B_\\text{1B} \\left( Z \\right)\$} & ");
    println("       \\multicolumn{1}{c}{\$C_\\text{1B} \\left( Z \\right)\$} & ");
    println("       \\multicolumn{1}{c}{\$D_\\text{1B} \\left( Z \\right)\$} ");
    println("       \\\\[0.1cm] \\hline \\\\[-0.35cm] \\hline \\\\[-0.35cm]");
    
    for Z in atomic_numbers
        sim = make_monoatomic_system(Z,0);
        println(@sprintf "        %d & " Z);

        xc_a_1b = sim.tot_e_coeffs.xc_coeffs.xc_a_1b[Z];
        xc_b_1b = sim.tot_e_coeffs.xc_coeffs.xc_b_1b[Z];
        xc_c_1b = sim.tot_e_coeffs.xc_coeffs.xc_c_1b[Z];
        xc_d_1b = sim.tot_e_coeffs.xc_coeffs.xc_d_1b[Z];

        println(@sprintf "        %s & " print_to_latex(xc_a_1b));
        println(@sprintf "        %s & " print_to_latex(xc_b_1b));
        println(@sprintf "        %s & " print_to_latex(xc_c_1b));
        println(@sprintf "        %s \\\\" print_to_latex(xc_d_1b));
    end

    println("   \\\\[-0.35cm]");
    println("   \\hline \\\\[-0.35cm] \\hline \\\\[-0.35cm]");
    println("\\end{tabular}")
end

function print_tot_e_ke()
    atomic_numbers = [1,6,7,8];

    println("\\begin{tabular}{cll}");
    println("   \\hline \\\\[-0.35cm] \\hline \\\\[-0.35cm]");
    println("   \$Z\$ & ");
    println("       \\multicolumn{1}{c}{\$E_\\text{1B} \\left( Z \\right)\$} & ");
    println("       \\multicolumn{1}{c}{\$F_\\text{1B} \\left( Z \\right)\$} ");
    println("       \\\\[0.1cm] \\hline \\\\[-0.35cm] \\hline \\\\[-0.35cm]");
    
    for Z in atomic_numbers
        sim = make_monoatomic_system(Z,0);
        println(@sprintf "        %d & " Z);

        ke_e_1b = sim.tot_e_coeffs.ke_coeffs.ke_e_1b[Z];
        ke_f_1b = sim.tot_e_coeffs.ke_coeffs.ke_f_1b[Z];

        println(@sprintf "        %s & " print_to_latex(ke_e_1b));
        println(@sprintf "        %s \\\\" print_to_latex(ke_f_1b));
    end

    println("   \\\\[-0.35cm]");
    println("   \\hline \\\\[-0.35cm] \\hline \\\\[-0.35cm]");
    println("\\end{tabular}")
end

function print_pol_e_xc()
    atomic_numbers = [1,6,7,8];

    println("\\begin{tabular}{cllll}");
    println("   \\hline \\\\[-0.35cm] \\hline \\\\[-0.35cm]");
    println("   \$Z\$ & ");
    println("       \\multicolumn{1}{c}{\$A_\\text{1B} \\left( Z \\right)\$} & ");
    println("       \\multicolumn{1}{c}{\$B_\\text{1B} \\left( Z \\right)\$} & ");
    println("       \\multicolumn{1}{c}{\$C_\\text{1B} \\left( Z \\right)\$} & ");
    println("       \\multicolumn{1}{c}{\$D_\\text{1B} \\left( Z \\right)\$} ");
    println("       \\\\[0.1cm] \\hline \\\\[-0.35cm] \\hline \\\\[-0.35cm]");
    
    for Z in atomic_numbers
        sim = make_monoatomic_system(Z,0);
        println(@sprintf "        %d & " Z);

        xc_a_1b = sim.pol_e_coeffs.xc_coeffs.xc_a_1b[Z];
        xc_b_1b = sim.pol_e_coeffs.xc_coeffs.xc_b_1b[Z];
        xc_c_1b = sim.pol_e_coeffs.xc_coeffs.xc_c_1b[Z];
        xc_d_1b = sim.pol_e_coeffs.xc_coeffs.xc_d_1b[Z];

        println(@sprintf "        %s & " print_to_latex(xc_a_1b));
        println(@sprintf "        %s & " print_to_latex(xc_b_1b));
        println(@sprintf "        %s & " print_to_latex(xc_c_1b));
        println(@sprintf "        %s \\\\" print_to_latex(xc_d_1b));
    end

    println("   \\\\[-0.35cm]");
    println("   \\hline \\\\[-0.35cm] \\hline \\\\[-0.35cm]");
    println("\\end{tabular}")
end

function print_pol_e_ke()
    atomic_numbers = [1,6,7,8];

    println("\\begin{tabular}{cll}");
    println("   \\hline \\\\[-0.35cm] \\hline \\\\[-0.35cm]");
    println("   \$Z\$ & ");
    println("       \\multicolumn{1}{c}{\$E_\\text{1B} \\left( Z \\right)\$} & ");
    println("       \\multicolumn{1}{c}{\$F_\\text{1B} \\left( Z \\right)\$} ");
    println("       \\\\[0.1cm] \\hline \\\\[-0.35cm] \\hline \\\\[-0.35cm]");
    
    for Z in atomic_numbers
        sim = make_monoatomic_system(Z,0);
        println(@sprintf "        %d & " Z);

        ke_e_1b = sim.pol_e_coeffs.ke_coeffs.ke_e_1b[Z];
        ke_f_1b = sim.pol_e_coeffs.ke_coeffs.ke_f_1b[Z];

        println(@sprintf "        %s & " print_to_latex(ke_e_1b));
        println(@sprintf "        %s \\\\" print_to_latex(ke_f_1b));
    end

    println("   \\\\[-0.35cm]");
    println("   \\hline \\\\[-0.35cm] \\hline \\\\[-0.35cm]");
    println("\\end{tabular}")
end
