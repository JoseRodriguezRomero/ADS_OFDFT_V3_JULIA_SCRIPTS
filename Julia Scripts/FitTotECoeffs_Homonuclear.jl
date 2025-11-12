using Optim;
using Printf, LinearAlgebra;
using Base.Threads;

include("FitCoeffs_General.jl")

# which_atomic_numbers = [1,6,7,8];
which_atomic_numbers = [8];
for atomic_number in which_atomic_numbers
    neutral_data, cation_data, anion_data = 
        read_all_sanitized_data(atomic_number);
    all_data = vcat(neutral_data, cation_data, anion_data);

    at_neutral = make_monoatomic_system(atomic_number,0);
    at_cation = make_monoatomic_system(atomic_number,1);
    at_anion = make_monoatomic_system(atomic_number,-1);

    n_threads = Threads.nthreads();
    simulation = Vector{SimulationSystem}();
    resize!(simulation,n_threads);
    for i in 1:n_threads
        simulation[i] = make_system_from_parsed_file(all_data[1]);
    end

    dft_at_neutral_e, dft_at_cation_e, dft_at_anion_e =
        get_reference_atom_total_energy();

    dft_at_neutral_e = dft_at_neutral_e[atomic_number];
    dft_at_cation_e = dft_at_cation_e[atomic_number];
    dft_at_anion_e = dft_at_anion_e[atomic_number];

    neutral_min_index = 1;
    for index in eachindex(neutral_data)
        if neutral_data[neutral_min_index].total_energy > neutral_data[index].total_energy
            neutral_min_index = index;
        end
    end

    reference_system = neutral_data[neutral_min_index];

    needs_casting = true;
    function cost_func(aux_X::Vector)
        aux_type = typeof(aux_X[1]);

        if needs_casting == true
            needs_casting = false;

            function cast_coeffs_to_type!(simulation::SimulationSystem)
                cast_tot_e_coeffs_to_type!(simulation,aux_type);
            end

            for thread_id in 1:n_threads
                cast_coeffs_to_type!(simulation[thread_id]);
            end

            cast_coeffs_to_type!(at_neutral);
            cast_coeffs_to_type!(at_cation);
            cast_coeffs_to_type!(at_anion);
        end

        function set_fitted_coeffs!(simulation::SimulationSystem)
            set_fitted_tot_e_coeffs!(simulation,atomic_number,aux_X);
        end

        set_fitted_coeffs!(at_neutral);
        set_fitted_coeffs!(at_cation);
        set_fitted_coeffs!(at_anion);

        model_at_neutral_e = total_energy(at_neutral);
        model_at_cation_e = total_energy(at_cation);
        model_at_anion_e = total_energy(at_anion);

        set_fitted_coeffs!(simulation[1]);
        set_diatomic_system_to_parsed_file!(simulation[1],reference_system);

        dft_e0 = reference_system.total_energy;
        model_e0 = total_energy(simulation[1]);

        ret_val = zeros(aux_type,n_threads);
        @threads for thread_id in 1:n_threads
            set_fitted_coeffs!(simulation[thread_id]);

            for i in thread_id:n_threads:length(all_data)
                set_diatomic_system_to_parsed_file!(
                    simulation[thread_id],all_data[i]);

                dft_tot_e = all_data[i].total_energy;
                model_tot_e = total_energy(simulation[thread_id]);

                dft_e_diff = dft_tot_e - dft_e0;
                model_e_diff = model_tot_e - model_e0;
                ret_val[thread_id] += (dft_e_diff - model_e_diff)^2;

                if i > 1
                    if all_data[i].charge != all_data[i-1].charge
                        continue;
                    end

                    d1 = all_data[i].atomic_separation;
                    d2 = all_data[i-1].atomic_separation;
                    Δd = d2 - d1;

                    set_diatomic_system_to_parsed_file!(
                        simulation[thread_id],all_data[i-1]);

                    dft_tot_e_nxt = all_data[i-1].total_energy;
                    model_tot_e_nxt = total_energy(simulation[thread_id]);

                    dft_ΔE = dft_tot_e_nxt - dft_tot_e;
                    model_ΔE = model_tot_e_nxt - model_tot_e;

                    dft_dev = dft_ΔE/Δd;
                    model_dev = model_ΔE/Δd;

                    ret_val[thread_id] += (model_dev - dft_dev)^2;
                end
            end
        end

        ret_val = sum(ret_val) / length(all_data);

        diff1 = model_at_neutral_e - dft_at_neutral_e;
        diff2 = model_at_cation_e - dft_at_cation_e;
        diff3 = model_at_anion_e - dft_at_anion_e;
        
        ret_val += (diff1 - diff2)^2;
        ret_val += (diff1 - diff3)^2;
        ret_val += (diff2 - diff3)^2;

        return ret_val;
    end

    num_vars = 13;
    aux_X = 2.0 .* (rand(Float64,num_vars) .- 0.5);

    for i in 1:2500
        cost_func_eval = cost_func(aux_X);
        new_aux_X = 2.0 .* (rand(Float64,num_vars) .- 0.5);

        if cost_func_eval <= 0.001
            println("Initial guess found!");
            break;
        end

        if cost_func(new_aux_X) < cost_func_eval
            aux_X = new_aux_X;
            print(@sprintf "Current best %18.6E \n" cost_func(aux_X));
        end
    end

    needs_casting = true;
    sol = Optim.optimize(cost_func, aux_X[:], LBFGS(), autodiff=:forward,
        Optim.Options(show_trace=true,iterations=8000));
    aux_X = Optim.minimizer(sol);

    display(aux_X);

    simulation = make_system_from_parsed_file(all_data[1]);
    
    set_fitted_tot_e_coeffs!(simulation,atomic_number,aux_X);

    save_fitted_coeffs(simulation);
end
