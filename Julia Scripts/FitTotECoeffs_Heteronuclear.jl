using Optim, Plots;
using Printf, LinearAlgebra;
using LaTeXStrings, Latexify, Measures;
using Base.Threads;

include("FitCoeffs_General.jl")

which_atomic_numbers = [[1,6],[1,8],[6,8],[7,8]];
for atomic_numbers in which_atomic_numbers
    Z1 = atomic_numbers[1];
    Z2 = atomic_numbers[2];

    neutral_data, cation_data, anion_data = 
        read_all_sanitized_data(Z1,Z2);
    all_data = vcat(neutral_data, cation_data, anion_data);

    n_threads = Threads.nthreads();
    simulation = Vector{SimulationSystem}();
    resize!(simulation,n_threads);
    for i in 1:n_threads
        simulation[i] = make_system_from_parsed_file(all_data[1]);
    end

    num_xc_2b_coeffs = 4;
    xc_2b_order = size(simulation[1].pol_e_xc_coeffs.xc_a_2b)[2] - 1;
    num_xc_coeffs = (xc_2b_order+1)*num_xc_2b_coeffs;

    num_2b_static_coeffs = 1;
    num_static_coeffs = num_2b_static_coeffs;

    num_vars = num_xc_coeffs + num_static_coeffs;
    aux_X = zeros(Float64,num_vars);
    
    neutral_at1 = make_atom_system(Z1,0);
    neutral_at2 = make_atom_system(Z2,0);

    dft_neutral_at_e, _, _ = get_reference_atom_total_energy();
    dft_neutral_at1_e = dft_neutral_at_e[Z1];
    dft_neutral_at2_e = dft_neutral_at_e[Z2];

    polarize_to_model!(neutral_data);
    polarize_to_model!(cation_data);
    polarize_to_model!(anion_data);

    needs_casting = true;
    function cost_func(aux_X::Vector)
        aux_type = typeof(aux_X[1]);

        function set_fitted_coeffs!(simulation::SimulationSystem)
            xc_coeffs = aux_X[1:num_xc_coeffs];
            static_coeffs = aux_X[(num_xc_coeffs+1):end];
            set_fitted_tot_e_coeffs!(simulation,Z1,Z2,xc_coeffs);
            set_fitted_static_coeffs!(simulation,Z1,Z2,static_coeffs);
        end

        if needs_casting
            needs_casting = false;
            function cast_coeffs_to_type!(simulation::SimulationSystem)
                cast_tot_e_coeffs_to_type!(simulation,aux_type);
                cast_static_coeffs_to_type!(simulation,aux_type);
            end

            for thread_id in 1:n_threads
                cast_coeffs_to_type!(simulation[thread_id]);
            end

            cast_coeffs_to_type!(neutral_at1);
            cast_coeffs_to_type!(neutral_at2);
        end

        set_fitted_coeffs!(neutral_at1);
        set_fitted_coeffs!(neutral_at2);

        model_neutral_at1_e = total_energy(neutral_at1);
        model_neutral_at2_e = total_energy(neutral_at2);

        ret_val = zeros(aux_type,n_threads);
        @threads for thread_id in 1:n_threads
            set_fitted_coeffs!(simulation[thread_id])

            for i in thread_id:n_threads:length(all_data)
                set_diatomic_system_to_parsed_file!(
                    simulation[thread_id],all_data[i]);

                dft_tot_e = all_data[i].total_energy;
                model_tot_e = total_energy(simulation[thread_id]);
                e_diff = model_tot_e - dft_tot_e;
                e_diff -= model_neutral_at1_e + model_neutral_at2_e;
                e_diff += dft_neutral_at1_e + dft_neutral_at2_e;
                ret_val[thread_id] += e_diff^2;

                if (i < length(all_data))
                    if (all_data[i].charge == all_data[i+1].charge)
                        set_diatomic_system_to_parsed_file!(
                            simulation[thread_id],all_data[i+1]);

                        dft_tot_e_nxt = all_data[i+1].total_energy;
                        model_tot_e_nxt = total_energy(simulation[thread_id]);

                        Δx = all_data[i+1].atomic_separation;
                        Δx -= all_data[i].atomic_separation;

                        dft_diff = dft_tot_e - dft_tot_e_nxt;
                        model_diff = model_tot_e - model_tot_e_nxt;
                        e_diff = model_diff - dft_diff;
                        ret_val[thread_id] += (e_diff/Δx)^2;
                    end
                end
            end
        end

        return sum(ret_val);
    end

    # needs_casting = false;
    # aux_X = zeros(Float64,num_vars);
    # sol = Optim.optimize(cost_func, aux_X[:], NelderMead(),
    #     Optim.Options(show_trace=true,iterations=8000));
    # aux_X = Optim.minimizer(sol);

    needs_casting = true;
    aux_X = zeros(Float64,num_vars);
    sol = Optim.optimize(cost_func, aux_X[:], LBFGS(), autodiff=:forward,
        Optim.Options(show_trace=true,iterations=2000));
    aux_X = Optim.minimizer(sol);

    display(aux_X);

    simulation = make_system_from_parsed_file(all_data[1]);

    xc_coeffs = aux_X[1:num_xc_coeffs];
    static_coeffs = aux_X[(num_xc_coeffs+1):end];
    set_fitted_tot_e_coeffs!(simulation,Z1,Z2,xc_coeffs);
    set_fitted_static_coeffs!(simulation,Z1,Z2,static_coeffs);

    save_fitted_coeffs(simulation);
end
