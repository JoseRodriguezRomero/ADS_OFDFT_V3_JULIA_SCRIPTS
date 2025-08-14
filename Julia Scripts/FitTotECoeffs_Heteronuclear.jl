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
    
    neutral_at1 = make_atom_system(Z1,0);
    neutral_at2 = make_atom_system(Z2,0);

    dft_neutral_at_e, _, _ = get_reference_atom_total_energy();
    dft_neutral_at1_e = dft_neutral_at_e[Z1];
    dft_neutral_at2_e = dft_neutral_at_e[Z2];

    n_threads = Threads.nthreads();
    simulation = Vector{SimulationSystem}();
    resize!(simulation,n_threads);
    for i in 1:n_threads
        simulation[i] = make_system_from_parsed_file(all_data[1]);
    end

    needs_casting = true;

    function cost_func(aux_X::Vector)
        aux_type = typeof(aux_X[1]);

        if needs_casting
            needs_casting = false;
            for thread_id in 1:n_threads
                cast_tot_e_coeffs_to_type!(simulation[thread_id],aux_type);
            end

            cast_tot_e_coeffs_to_type!(neutral_at1,aux_type);
            cast_tot_e_coeffs_to_type!(neutral_at2,aux_type);
        end

        set_fitted_tot_e_coeffs!(neutral_at1,Z1,Z2,aux_X);
        set_fitted_tot_e_coeffs!(neutral_at2,Z1,Z2,aux_X);

        model_neutral_at1_e = total_energy(neutral_at1);
        model_neutral_at2_e = total_energy(neutral_at2);

        ret_val = zeros(aux_type,n_threads);
        @threads for thread_id in 1:n_threads
            set_fitted_tot_e_coeffs!(simulation[thread_id],Z1,Z2,aux_X);

            for i in thread_id:n_threads:length(all_data)
                set_diatomic_system_to_parsed_file!(
                    simulation[thread_id],all_data[i]);

                dft_tot_e = all_data[i].total_energy;
                model_tot_e = total_energy(simulation[thread_id]);
                e_diff = model_tot_e - dft_tot_e;
                e_diff -= model_neutral_at1_e + model_neutral_at2_e;
                e_diff += dft_neutral_at1_e + dft_neutral_at2_e;
                ret_val[thread_id] += e_diff^2;
            end
        end

        return sum(ret_val);
    end

    aux_X = zeros(Float64,num_vars);
    sol = Optim.optimize(cost_func, aux_X[:], LBFGS(), autodiff=:forward,
        Optim.Options(show_trace=true));
    aux_X = Optim.minimizer(sol);

    simulation = make_system_from_parsed_file(all_data[1]);
    set_fitted_tot_e_coeffs!(simulation,Z1,Z2,aux_X);

    save_fitted_coeffs(simulation);
end
