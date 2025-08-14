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
    
    z1_eff = atom_eff_atomic_number(simulation[1].system.molecules[1],1);
    z2_eff = atom_eff_atomic_number(simulation[1].system.molecules[1],2);

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

                aux_x[1] = ζ1;
                aux_x[2] = ζ2;
                aux_x[3] = 0.0;
                aux_diff = (aux_m * aux_x) - aux_y;

                ret_val[thread_id] += (aux_diff[1] / z1_eff - μ)^2.0;
                ret_val[thread_id] += (aux_diff[2] / z2_eff - μ)^2.0;
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
