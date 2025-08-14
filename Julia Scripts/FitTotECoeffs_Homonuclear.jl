using Optim, Plots;
using Printf, LinearAlgebra;
using LaTeXStrings, Latexify, Measures;
using Base.Threads;

include("FitHomonuclearCoeffs_General.jl")

which_atomic_numbers = [1,6,7,8];
for atomic_number in which_atomic_numbers
    neutral_data, cation_data, anion_data = 
        read_all_sanitized_data(atomic_number);
    all_data = vcat(neutral_data, cation_data, anion_data);

    num_1b_coeffs = 4;
    num_2b_coeffs = 4;
    num_vars = 2*num_1b_coeffs + 2*num_2b_coeffs;
    
    neutral_at = make_atom_system(atomic_number,0);
    cation_at = make_atom_system(atomic_number,1);
    anion_at = make_atom_system(atomic_number,-1);

    dft_neutral_at_e, dft_cation_at_e, dft_anion_at_e =
        get_reference_atom_total_energy();

    dft_neutral_at_e = dft_neutral_at_e[atomic_number];
    dft_cation_at_e = dft_cation_at_e[atomic_number];
    dft_anion_at_e = dft_anion_at_e[atomic_number];

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

            cast_tot_e_coeffs_to_type!(neutral_at,aux_type);
            cast_tot_e_coeffs_to_type!(cation_at,aux_type);
            cast_tot_e_coeffs_to_type!(anion_at,aux_type);
        end

        set_fitted_tot_e_coeffs!(neutral_at,atomic_number,aux_X);
        set_fitted_tot_e_coeffs!(cation_at,atomic_number,aux_X);
        set_fitted_tot_e_coeffs!(anion_at,atomic_number,aux_X);

        model_neutral_at_e = total_energy(neutral_at)
        model_cation_at_e = total_energy(cation_at);
        model_anion_at_e = total_energy(anion_at);

        ret_val = zeros(aux_type,n_threads);
        @threads for thread_id in 1:n_threads
            set_fitted_tot_e_coeffs!(
                simulation[thread_id],atomic_number,aux_X);

            for i in thread_id:n_threads:length(all_data)
                set_diatomic_system_to_parsed_file!(
                    simulation[thread_id],all_data[i]);
                charge = all_data[i].charge;

                dft_tot_e = all_data[i].total_energy;
                model_tot_e = total_energy(simulation[thread_id]);
                e_diff = model_tot_e - dft_tot_e;

                if charge == 0
                    e_diff -= model_neutral_at_e + model_neutral_at_e;
                    e_diff += dft_neutral_at_e + dft_neutral_at_e;
                    ret_val[thread_id] += e_diff^2;
                elseif charge == 1
                    e_diff -= model_neutral_at_e + model_cation_at_e;
                    e_diff += dft_neutral_at_e + dft_cation_at_e;
                    ret_val[thread_id] += e_diff^2;
                else
                    e_diff -= model_neutral_at_e + model_anion_at_e;
                    e_diff += dft_neutral_at_e + dft_anion_at_e;
                    ret_val[thread_id] += e_diff^2;
                end
            end
        end

        ret_val = sum(ret_val);

        diff1 = model_neutral_at_e - model_cation_at_e;
        diff1 -= dft_neutral_at_e - dft_cation_at_e;
        ret_val += diff1^2;

        diff2 = model_neutral_at_e - model_anion_at_e;
        diff2 -= dft_neutral_at_e - dft_anion_at_e;
        ret_val += diff2^2;

        return ret_val;
    end

    aux_X = zeros(Float64,num_vars);
    sol = Optim.optimize(cost_func, aux_X[:], LBFGS(), autodiff=:forward,
        Optim.Options(show_trace=true));
    aux_X = Optim.minimizer(sol);

    simulation = make_system_from_parsed_file(all_data[1]);
    set_fitted_tot_e_coeffs!(simulation,atomic_number,aux_X);

    save_fitted_coeffs(simulation);
end
