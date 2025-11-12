using Optim;
using Printf, LinearAlgebra;
using Base.Threads;

include("FitCoeffs_General.jl")

# which_atomic_numbers = [[1,6],[1,8],[6,8],[7,8]];
which_atomic_numbers = [[6,8]];
for atomic_numbers in which_atomic_numbers
    Z1 = atomic_numbers[1];
    Z2 = atomic_numbers[2];

    local neutral_data, cation_data, anion_data = 
        read_all_sanitized_data(Z1,Z2,true);
    all_data = vcat(neutral_data, cation_data, anion_data);

    n_threads = Threads.nthreads();
    simulation = Vector{SimulationSystem}();
    resize!(simulation,n_threads);

    simulation = Vector{SimulationSystem}();
    resize!(simulation,n_threads);
    for i in 1:n_threads
        simulation[i] = make_system_from_parsed_file(all_data[1]);
    end

    z1_eff = simulation[1].system.molecules[1].atoms[1].valence_electrons;
    z2_eff = simulation[1].system.molecules[1].atoms[2].valence_electrons;

    needs_casting = true;
    function cost_func(aux_X::Vector)
        aux_type = typeof(aux_X[1]);

        function set_fitted_coeffs!(simulation::SimulationSystem)
            set_fitted_pol_e_coeffs!(simulation,Z1,Z2,aux_X);
        end

        if needs_casting
            needs_casting = false;
            for thread_id in 1:n_threads
                cast_pol_e_coeffs_to_type!(simulation[thread_id],aux_type);
            end
        end
        
        ret_val = zeros(aux_type,n_threads);
        @threads for thread_id in 1:n_threads
            set_fitted_coeffs!(simulation[thread_id]);

            atoms_μ0 = simulation[thread_id].basis_set_settings.atoms_μ0;
            μ0_1 = atoms_μ0[Z1];
            μ0_2 = atoms_μ0[Z2];

            for i in thread_id:n_threads:(length(all_data))
                set_diatomic_system_to_parsed_file!(
                    simulation[thread_id],all_data[i]);

                ζ1 = atom_polarization_coeff(
                    simulation[thread_id].system.molecules[1],1);
                ζ2 = atom_polarization_coeff(
                    simulation[thread_id].system.molecules[1],2);
                μ = simulation[thread_id].system.chemical_potential;

                aux_m, aux_y = 
                    polarization_matrix_problem(simulation[thread_id]);
                    
                aux_x = zeros(aux_type,5);
                aux_x[1] = ζ1;
                aux_x[2] = ζ2;
                aux_x[3] = μ - μ0_1;
                aux_x[4] = μ - μ0_2;
                aux_x[5] = μ;

                diff_vec = ((aux_m \ aux_y) - aux_x);
                diff_vec[1] *= z1_eff;
                diff_vec[2] *= z2_eff;

                ret_val[thread_id] += norm(diff_vec)^2;

                if i > 1
                    if all_data[i].charge != all_data[i-1].charge
                        continue;
                    end

                    d1 = all_data[i].atomic_separation;
                    d2 = all_data[i-1].atomic_separation;
                    Δd = d2 - d1;

                    set_diatomic_system_to_parsed_file!(
                        simulation[thread_id],all_data[i-1]);

                    ζ1_nxt = atom_polarization_coeff(
                        simulation[thread_id].system.molecules[1],1);
                    ζ2_nxt = atom_polarization_coeff(
                        simulation[thread_id].system.molecules[1],2);
                    μ_nxt = simulation[thread_id].system.chemical_potential;

                    aux_m_nxt, aux_y_nxt = 
                        polarization_matrix_problem(simulation[thread_id]);

                    aux_x_nxt = zeros(aux_type,5);
                    aux_x_nxt[1] = ζ1_nxt;
                    aux_x_nxt[2] = ζ2_nxt;
                    aux_x_nxt[3] = μ_nxt - μ0_1;
                    aux_x_nxt[4] = μ_nxt - μ0_2;
                    aux_x_nxt[5] = μ_nxt;

                    diff_vec_nxt = ((aux_m_nxt \ aux_y_nxt) - aux_x_nxt);
                    diff_vec_nxt[1] *= z1_eff;
                    diff_vec_nxt[2] *= z2_eff;

                    ret_val[thread_id] += (norm(diff_vec_nxt-diff_vec) / Δd)^2;
                end
            end
        end

        ret_val = sum(ret_val) / length(all_data);

        det_cond_sum = zeros(aux_type,n_threads);
        r = 0.1:((10.0-0.1)/200.0):10.0;
        @threads for thread_id in 1:n_threads
            set_diatomic_system_to_parsed_file!(
                simulation[thread_id],all_data[1]);
                simulation[thread_id].system.molecules[1].atoms[1].coordinates .= 0.0;
                simulation[thread_id].system.molecules[1].atoms[2].coordinates .= 0.0;

            for atomic_separation in r[thread_id:n_threads:end]
                simulation[thread_id].system.molecules[1].atoms[1].coordinates[3] = 
                    atomic_separation;

                aux_m, _ = polarization_matrix_problem(simulation[thread_id]);
                det_aux_m = det(aux_m);

                det_cond_sum[thread_id] += (det_aux_m - abs(det_aux_m))^2;
            end
        end
        ret_val += sum(det_cond_sum) / length(det_cond_sum);

        return ret_val;
    end

    num_vars = 4;
    aux_X = 2.0 .* (rand(Float64,num_vars) .- 0.5);

    for i in 1:2500
        cost_func_eval = cost_func(aux_X);
        new_aux_X = 2.0 .* (rand(Float64,num_vars) .- 0.5);

        if cost_func_eval <= 0.1
            println("Initial guess found!");
            break;
        end

        if cost_func(new_aux_X) < cost_func_eval
            aux_X = new_aux_X;
            print(@sprintf "Current best %18.6E \n" cost_func(aux_X));
        end
    end

    needs_casting = true;
    sol = Optim.optimize(cost_func, aux_X, LBFGS(), autodiff=:forward,
        Optim.Options(show_trace=true,iterations=2000));
    aux_X = Optim.minimizer(sol);

    display(aux_X);

    simulation = make_system_from_parsed_file(all_data[1]);
    set_fitted_pol_e_coeffs!(simulation,Z1,Z2,aux_X);

    save_fitted_coeffs(simulation);
end
