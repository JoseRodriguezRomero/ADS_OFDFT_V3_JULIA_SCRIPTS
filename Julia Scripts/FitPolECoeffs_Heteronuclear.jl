using Optim, Plots;
using Printf, LinearAlgebra;
using LaTeXStrings, Latexify, Measures;
using Base.Threads;

include("FitHomonuclearCoeffs_General.jl")

which_atomic_numbers = [[1,6],[1,8],[6,8],[7,8]];
for atomic_numbers in which_atomic_numbers
    Z1 = atomic_numbers[1];
    Z2 = atomic_numbers[2];

    neutral_data, cation_data, anion_data = read_all_sanitized_data(Z1,Z2);
    all_data = vcat(neutral_data, cation_data, anion_data);

    num_2b_coeffs = 4;
    num_vars = 3*num_2b_coeffs;
    aux_X = zeros(Float64,num_vars);

    n_threads = Threads.nthreads();
    simulation = Vector{SimulationSystem}();
    resize!(simulation,n_threads);

    for thread_id in 1:n_threads
        simulation[thread_id] = make_system_from_parsed_file(all_data[1]);
    end
    
    needs_casting = true;
    function cost_func(aux_X::Vector)
        aux_type = typeof(aux_X[1]);

        if needs_casting
            needs_casting = false;
            for thread_id in 1:n_threads
                cast_pol_e_coeffs_to_type!(simulation[thread_id],aux_type);
            end
        end
        
        ret_val = zeros(aux_type,n_threads);
        @threads for thread_id in 1:n_threads
            set_fitted_pol_e_coeffs!(simulation[thread_id],Z1,Z2,aux_X);
            cast_type = simulation[thread_id].cast_types.cast_to_xc_coeff_type;

            for i in thread_id:n_threads:(length(all_data))
                set_diatomic_system_to_parsed_file!(
                    simulation[thread_id],all_data[i]);

                aux_m, aux_y = polarization_matrix_problem(
                    simulation[thread_id],cast_type);

                ζ1 = simulation[thread_id].system.molecules[1].cloud_data[1,6];
                ζ2 = simulation[thread_id].system.molecules[1].cloud_data[4,6];
                μ = simulation[thread_id].system.chemical_potential;

                aux_x = zeros(aux_type,3);
                aux_x[1] = ζ1;
                aux_x[2] = ζ2;
                aux_x[3] = μ;

                aux_diff = (aux_m * aux_x) - aux_y;
                ret_val[thread_id] += dot(aux_diff,aux_diff);
            end
        end

        return sum(ret_val);
    end

    sol = Optim.optimize(cost_func, aux_X, LBFGS(), autodiff=:forward,
        Optim.Options(show_trace=true,iterations=8000));
    aux_X = Optim.minimizer(sol);

    simulation = make_system_from_parsed_file(all_data[1]);
    set_fitted_pol_e_coeffs!(simulation,Z1,Z2,aux_X);

    save_fitted_coeffs(simulation);
end

function test_tesult_chemical_potential(Z1::Integer, Z2::Integer)
    neutral_data, cation_data, anion_data = read_all_sanitized_data(Z1,Z2);
    all_data = vcat(neutral_data, cation_data, anion_data);

    model_chem_μ = zeros(Float64,length(all_data));
    dft_chem_μ = zeros(Float64,length(all_data));

    function get_data(data::Vector{ParsedFile})
        model_chem_μ = zeros(Float64,length(data));
        dft_chem_μ = zeros(Float64,length(data));

        n_threads = Threads.nthreads();
        @threads for thread_id in 1:n_threads
            simulation = make_system_from_parsed_file(data[1]);
            for i in thread_id:n_threads:length(data)
                set_diatomic_system_to_parsed_file!(simulation,data[i]);
                polarize_molecules!(simulation);
                
                model_μ = simulation.system.chemical_potential;
                dft_μ = (data[i].HOMO_energy + data[i].LUMO_energy)/2;
                
                model_chem_μ[i] = model_μ;
                dft_chem_μ[i] = dft_μ
            end
        end

        return model_chem_μ, dft_chem_μ;
    end

    neutral_model_chem_μ, neutral_dft_chem_μ = get_data(neutral_data);
    cation_model_chem_μ, cation_dft_chem_μ = get_data(cation_data);
    anion_model_chem_μ, anion_dft_chem_μ = get_data(anion_data);

    model_chem_μ = 
        vcat(neutral_model_chem_μ,cation_model_chem_μ,anion_model_chem_μ);
    dft_chem_μ = 
        vcat(neutral_dft_chem_μ,cation_dft_chem_μ,anion_dft_chem_μ);

    mean_dft_chem_μ = sum(dft_chem_μ) / length(dft_chem_μ);
    SS_res = sum((model_chem_μ - dft_chem_μ).^2.0);
    SS_tot = sum((dft_chem_μ .- mean_dft_chem_μ).^2.0);
    R² = 1.0 - SS_res / SS_tot;

    min_x = minimum([minimum(model_chem_μ), minimum(dft_chem_μ)]);
    max_x = maximum([maximum(model_chem_μ), maximum(dft_chem_μ)]);

    elem_symbol_1 = get_element_symbol(Z1);
    elem_symbol_2 = get_element_symbol(Z2);

    neutral_label = elem_symbol_1*" + "*elem_symbol_2;
    cation_label = "("*elem_symbol_1*" + "*elem_symbol_2*")⁺";
    anion_label = "("*elem_symbol_1*" + "*elem_symbol_2*")⁻";

    Pts = [-5,5];
    p = plot(Pts,Pts,label=false);
    scatter!(neutral_dft_chem_μ, neutral_model_chem_μ, label=neutral_label);
    scatter!(cation_dft_chem_μ, cation_model_chem_μ, label=cation_label);
    scatter!(anion_dft_chem_μ, anion_model_chem_μ, label=anion_label);
    plot!(framestyle = :box);
    plot!(legend=:topleft);

    return p, R²;
end

function compare_chemical_potential()
    y_label_all = L"$\tilde{\mu} \quad \mathrm{(This \ Work)}$";
    x_label_all = L"$\tilde{\mu} \quad (\mathrm{KS{-}DFT})$";

    R²_rel_pos_x = (5.0/8.0);
    R²_rel_pos_y = (1.0/8.0);

    # HC Plots
    pHC, R² = test_tesult_chemical_potential(1,6);
    aux_lims = [-1.0,0.6];
    aux_ticks = -1.0:0.4:0.6;
    plot!(xlims=aux_lims);
    plot!(ylims=aux_lims);
    plot!(xticks=aux_ticks);
    plot!(yticks=aux_ticks);
    plot!(ylabel=y_label_all);
    plot!(xlabel=x_label_all);
    plot!(left_margin=6mm, bottom_margin=8mm);

    l_x_pos = aux_lims[1] + R²_rel_pos_x*(aux_lims[2] - aux_lims[1]);
    l_y_pos = aux_lims[1] + R²_rel_pos_y*(aux_lims[2] - aux_lims[1]);
    R² = round(R²,digits=5);
    annotate!(l_x_pos, l_y_pos, text("R² = "*string(R²), :center, 10));

    # HO Plots
    pHO, R² = test_tesult_chemical_potential(1,8);
    aux_lims = [-1.0,0.6];
    aux_ticks = -1.0:0.4:0.6;
    plot!(xlims=aux_lims);
    plot!(ylims=aux_lims);
    plot!(xticks=aux_ticks);
    plot!(yticks=aux_ticks);
    plot!(xlabel=x_label_all);

    l_x_pos = aux_lims[1] + R²_rel_pos_x*(aux_lims[2] - aux_lims[1]);
    l_y_pos = aux_lims[1] + R²_rel_pos_y*(aux_lims[2] - aux_lims[1]);
    R² = round(R²,digits=5);
    annotate!(l_x_pos, l_y_pos, text("R² = "*string(R²), :center, 10));

    # CO Plots
    pCO, R² = test_tesult_chemical_potential(6,8);
    aux_lims = [-1.0,0.6];
    aux_ticks = -1.0:0.4:0.6;
    plot!(xlims=aux_lims);
    plot!(ylims=aux_lims);
    plot!(xticks=aux_ticks);
    plot!(yticks=aux_ticks);
    plot!(xlabel=x_label_all);

    l_x_pos = aux_lims[1] + R²_rel_pos_x*(aux_lims[2] - aux_lims[1]);
    l_y_pos = aux_lims[1] + R²_rel_pos_y*(aux_lims[2] - aux_lims[1]);
    R² = round(R²,digits=5);
    annotate!(l_x_pos, l_y_pos, text("R² = "*string(R²), :center, 10));

    # NO Plots
    pNO, R² = test_tesult_chemical_potential(7,8);
    aux_lims = [-1.0,0.6];
    aux_ticks = -1.0:0.4:0.6;
    plot!(xlims=aux_lims);
    plot!(ylims=aux_lims);
    plot!(xticks=aux_ticks);
    plot!(yticks=aux_ticks);
    plot!(xlabel=x_label_all);

    l_x_pos = aux_lims[1] + R²_rel_pos_x*(aux_lims[2] - aux_lims[1]);
    l_y_pos = aux_lims[1] + R²_rel_pos_y*(aux_lims[2] - aux_lims[1]);
    R² = round(R²,digits=5);
    annotate!(l_x_pos, l_y_pos, text("R² = "*string(R²), :center, 10));

    p = plot(pHC,pHO, pCO, pNO, layout=(1,4), size = (1100,235));
    savefig("Figures/HeteronuclearFitComps.pdf");

    return p;
end

compare_chemical_potential()
