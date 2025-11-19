using Printf;

include("FitCoeffs_General.jl")

function print_to_latex(number::Float64, digits::Int = 6)
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

function print_tot_e_xc(atomic_numer::Number)
    atomic_pairs = [
        [1,1],[6,6],[7,7],[8,8],
        [1,6],[1,8],[6,8],[7,8]
    ];

    println("\\begin{tabular}{ccll}");
    println("   \\hline \\\\[-0.35cm] \\hline \\\\[-0.35cm]");
    println("   Shell & \$n\$ & ");
    println("       \\multicolumn{1}{c}{\$\\Lambda_n^\\text{Shell} \\left( Z_m \\right)\$} & ");
    println("       \\multicolumn{1}{c}{\$C_n^\\text{Shell} \\left( Z_m \\right)\$} ");
    println("       \\\\[0.1cm] \\hline \\\\[-0.35cm] \\hline \\\\[-0.35cm]");

    sim = make_monoatomic_system(atomic_numer,0);
    electron_shells = sim.system.molecules[1].atoms[1].electron_cloud_shells;

    for i in eachindex(electron_shells)
        electron_shell = electron_shells[i];
        C = electron_shell.basis_function_amplitude;
        Λ = electron_shell.basis_function_decay;

        for j in eachindex(C)
            shell_name = "";
            if i == 1
                shell_name = "K";
            elseif i == 2
                shell_name = "L";
            else
                shell_name = "M";
            end

            if j == 1
                println(@sprintf "        %s & %d & %s & " shell_name j print_to_latex(Λ[j]));
            else
                println(@sprintf "          & %d & %s & " j print_to_latex(Λ[j]));
            end
            

            if i != length(electron_shells) && j == length(C)
                println(@sprintf "                %s \\\\[0.2cm]" print_to_latex(C[j]));
            else
                println(@sprintf "                %s \\\\" print_to_latex(C[j]));
            end
        end
    end

    println("   \\\\[-0.35cm]");
    println("   \\hline \\\\[-0.35cm] \\hline \\\\[-0.35cm]");
    println("\\end{tabular}")
end
