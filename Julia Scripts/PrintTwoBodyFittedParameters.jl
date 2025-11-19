using Printf;

include("FitCoeffs_General.jl")

function print_to_latex(number::Float64, digits::Int = 4)
    if abs(number) < 1.0E-8
        return "\\hphantom{\$-\$}\$ 0.0000 \\times 10^{ 0 } \$";
    end

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
    atomic_pairs = [
        [1,1],[6,6],[7,7],[8,8],
        [1,6],[1,8],[6,8],[7,8]
    ];

    println("\\begin{tabular}{ccllll}");
    println("   \\hline \\\\[-0.35cm] \\hline \\\\[-0.35cm]");
    println("   \$Z_1\$ & \$Z_2\$ & ");
    println("       \\multicolumn{1}{c}{\$A_\\text{1B} \\left( Z_1, Z_2 \\right)\$} & ");
    println("       \\multicolumn{1}{c}{\$B_\\text{1B} \\left( Z_1, Z_2 \\right)\$} & ");
    println("       \\multicolumn{1}{c}{\$C_\\text{1B} \\left( Z_1, Z_2 \\right)\$} & ");
    println("       \\multicolumn{1}{c}{\$D_\\text{1B} \\left( Z_1, Z_2 \\right)\$} ");
    println("       \\\\[0.1cm] \\hline \\\\[-0.35cm] \\hline \\\\[-0.35cm]");
    
    for atomic_pair in atomic_pairs
        Z1 = atomic_pair[1];
        Z2 = atomic_pair[2];
        sim = make_diatomic_system(Z1,Z2,2.0,0);
        println(@sprintf "        %d & %d & " Z1 Z2);

        xc_a_2b = sim.tot_e_coeffs.xc_coeffs.xc_a_2b[(Z1,Z2)];
        xc_b_2b = sim.tot_e_coeffs.xc_coeffs.xc_b_2b[(Z1,Z2)];
        xc_c_2b = sim.tot_e_coeffs.xc_coeffs.xc_c_2b[(Z1,Z2)];
        xc_d_2b = sim.tot_e_coeffs.xc_coeffs.xc_d_2b[(Z1,Z2)];

        println(@sprintf "        %s & " print_to_latex(xc_a_2b));
        println(@sprintf "        %s & " print_to_latex(xc_b_2b));
        println(@sprintf "        %s & " print_to_latex(xc_c_2b));
        println(@sprintf "        %s \\\\" print_to_latex(xc_d_2b));
    end

    println("   \\\\[-0.35cm]");
    println("   \\hline \\\\[-0.35cm] \\hline \\\\[-0.35cm]");
    println("\\end{tabular}")
end

function print_pol_e_xc()
    atomic_pairs = [
        [1,1],[6,6],[7,7],[8,8],
        [1,6],[1,8],[6,8],[7,8]
    ];

    println("\\begin{tabular}{ccllll}");
    println("   \\hline \\\\[-0.35cm] \\hline \\\\[-0.35cm]");
    println("   \$Z_1\$ & \$Z_2\$ & ");
    println("       \\multicolumn{1}{c}{\$A_\\text{1B} \\left( Z_1, Z_2 \\right)\$} & ");
    println("       \\multicolumn{1}{c}{\$B_\\text{1B} \\left( Z_1, Z_2 \\right)\$} & ");
    println("       \\multicolumn{1}{c}{\$C_\\text{1B} \\left( Z_1, Z_2 \\right)\$} & ");
    println("       \\multicolumn{1}{c}{\$D_\\text{1B} \\left( Z_1, Z_2 \\right)\$} ");
    println("       \\\\[0.1cm] \\hline \\\\[-0.35cm] \\hline \\\\[-0.35cm]");
    
    for atomic_pair in atomic_pairs
        Z1 = atomic_pair[1];
        Z2 = atomic_pair[2];
        sim = make_diatomic_system(Z1,Z2,2.0,0);
        println(@sprintf "        %d & %d & " Z1 Z2);

        xc_a_2b = sim.pol_e_coeffs.xc_coeffs.xc_a_2b[(Z1,Z2)];
        xc_b_2b = sim.pol_e_coeffs.xc_coeffs.xc_b_2b[(Z1,Z2)];
        xc_c_2b = sim.pol_e_coeffs.xc_coeffs.xc_c_2b[(Z1,Z2)];
        xc_d_2b = sim.pol_e_coeffs.xc_coeffs.xc_d_2b[(Z1,Z2)];

        println(@sprintf "        %s & " print_to_latex(xc_a_2b));
        println(@sprintf "        %s & " print_to_latex(xc_b_2b));
        println(@sprintf "        %s & " print_to_latex(xc_c_2b));
        println(@sprintf "        %s \\\\" print_to_latex(xc_d_2b));
    end

    println("   \\\\[-0.35cm]");
    println("   \\hline \\\\[-0.35cm] \\hline \\\\[-0.35cm]");
    println("\\end{tabular}")
end

function print_tot_e_xc()
    atomic_pairs = [
        [1,1],[6,6],[7,7],[8,8],
        [1,6],[1,8],[6,8],[7,8]
    ];

    println("\\begin{tabular}{ccllll}");
    println("   \\hline \\\\[-0.35cm] \\hline \\\\[-0.35cm]");
    println("   \$Z_1\$ & \$Z_2\$ & ");
    println("       \\multicolumn{1}{c}{\$A_\\text{1B} \\left( Z_1, Z_2 \\right)\$} & ");
    println("       \\multicolumn{1}{c}{\$B_\\text{1B} \\left( Z_1, Z_2 \\right)\$} & ");
    println("       \\multicolumn{1}{c}{\$C_\\text{1B} \\left( Z_1, Z_2 \\right)\$} & ");
    println("       \\multicolumn{1}{c}{\$D_\\text{1B} \\left( Z_1, Z_2 \\right)\$} ");
    println("       \\\\[0.1cm] \\hline \\\\[-0.35cm] \\hline \\\\[-0.35cm]");
    
    for atomic_pair in atomic_pairs
        Z1 = atomic_pair[1];
        Z2 = atomic_pair[2];
        sim = make_diatomic_system(Z1,Z2,2.0,0);
        println(@sprintf "        %d & %d & " Z1 Z2);

        xc_a_2b = sim.tot_e_coeffs.xc_coeffs.xc_a_2b[(Z1,Z2)];
        xc_b_2b = sim.tot_e_coeffs.xc_coeffs.xc_b_2b[(Z1,Z2)];
        xc_c_2b = sim.tot_e_coeffs.xc_coeffs.xc_c_2b[(Z1,Z2)];
        xc_d_2b = sim.tot_e_coeffs.xc_coeffs.xc_d_2b[(Z1,Z2)];

        println(@sprintf "        %s & " print_to_latex(xc_a_2b));
        println(@sprintf "        %s & " print_to_latex(xc_b_2b));
        println(@sprintf "        %s & " print_to_latex(xc_c_2b));
        println(@sprintf "        %s \\\\" print_to_latex(xc_d_2b));
    end

    println("   \\\\[-0.35cm]");
    println("   \\hline \\\\[-0.35cm] \\hline \\\\[-0.35cm]");
    println("\\end{tabular}")
end

function print_morse_params()
    atomic_pairs = [
        [1,1],[6,6],[7,7],[8,8],
        [1,6],[1,8],[6,8],[7,8]
    ];

    println("\\begin{tabular}{cclll}");
    println("   \\hline \\\\[-0.35cm] \\hline \\\\[-0.35cm]");
    println("   \$Z_1\$ & \$Z_2\$ & ");
    println("       \\multicolumn{1}{c}{\$M^{(A)}_\\text{1B} \\left( Z_1, Z_2 \\right)\$} & ");
    println("       \\multicolumn{1}{c}{\$M^{(B)}_\\text{1B} \\left( Z_1, Z_2 \\right)\$} & ");
    println("       \\multicolumn{1}{c}{\$M^{(C)}_\\text{1B} \\left( Z_1, Z_2 \\right)\$} ");
    println("       \\\\[0.1cm] \\hline \\\\[-0.35cm] \\hline \\\\[-0.35cm]");
    
    for atomic_pair in atomic_pairs
        Z1 = atomic_pair[1];
        Z2 = atomic_pair[2];
        sim = make_diatomic_system(Z1,Z2,2.0,0);
        println(@sprintf "        %d & %d & " Z1 Z2);

        A = sim.tot_e_coeffs.non_polarizable_coeffs.depth[(Z1,Z2)];
        B = sim.tot_e_coeffs.non_polarizable_coeffs.stiffness_parameter[(Z1,Z2)];
        C = sim.tot_e_coeffs.non_polarizable_coeffs.equilibrium_distance[(Z1,Z2)];

        println(@sprintf "        %s & " print_to_latex(A));
        println(@sprintf "        %s & " print_to_latex(B));
        println(@sprintf "        %s \\\\" print_to_latex(C));
    end

    println("   \\\\[-0.35cm]");
    println("   \\hline \\\\[-0.35cm] \\hline \\\\[-0.35cm]");
    println("\\end{tabular}")
end
