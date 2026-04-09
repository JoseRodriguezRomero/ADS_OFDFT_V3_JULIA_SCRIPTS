using Optim
using Printf, LinearAlgebra
using Base.Threads

include("FitCoeffs_General.jl")

which_atomic_numbers = [[1,6],[1,8],[6,8],[7,8]]
for atomic_numbers in which_atomic_numbers
    Z1 = atomic_numbers[1]
    Z2 = atomic_numbers[2]

    local neutral_data, cation_data, anion_data =
        read_all_sanitized_data(Z1,Z2,true)
    all_data = vcat(neutral_data, cation_data, anion_data)

    n_threads = Threads.nthreads()
    simulation = Vector{SimulationSystem}()
    resize!(simulation,n_threads)

    simulation = Vector{SimulationSystem}()
    resize!(simulation,n_threads)
    for i in 1:n_threads
        simulation[i] = make_system_from_parsed_file(all_data[1])
    end

    z1_eff = simulation[1].system.molecules[1].atoms[1].valence_electrons
    z2_eff = simulation[1].system.molecules[1].atoms[2].valence_electrons

    needs_casting = true
    function cost_func(aux_X::Vector)
        aux_type = typeof(aux_X[1])

        function set_fitted_coeffs!(simulation::SimulationSystem)
            set_fitted_pol_e_coeffs!(simulation,Z1,Z2,aux_X)
        end

        if needs_casting
            needs_casting = false
            for thread_id in 1:n_threads
                cast_pol_e_coeffs_to_type!(simulation[thread_id],aux_type)
            end
        end
        
        ret_val = zeros(aux_type,n_threads)
        @threads for thread_id in 1:n_threads
            set_fitted_coeffs!(simulation[thread_id])

            atoms_μ0 = simulation[thread_id].basis_set_settings.atoms_μ0
            μ0_1 = atoms_μ0[Z1]
            μ0_2 = atoms_μ0[Z2]

            for i in thread_id:n_threads:(length(all_data))
                set_diatomic_system_to_parsed_file!(
                    simulation[thread_id],all_data[i])

                system = simulation[thread_id].system
                molecule = system.molecules[1]
                atom_1 = molecule.atoms[1]
                atom_2 = molecule.atoms[2]

                ζ1 = copy(atom_1.polarization_coefficient)
                ζ2 = copy(atom_2.polarization_coefficient)
                μ = copy(system.chemical_potential)

                polarize_molecules!(simulation[thread_id])
                ζ1_model = copy(atom_1.polarization_coefficient)
                ζ2_model = copy(atom_2.polarization_coefficient)
                μ_model = copy(system.chemical_potential)
                    
                diff_vec = zeros(aux_type,3)
                diff_vec[1] = ζ1 - ζ1_model
                diff_vec[2] = ζ2 - ζ2_model
                diff_vec[3] = μ - μ_model

                ret_val[thread_id] += norm(diff_vec)^2

                if i > 1
                    if all_data[i].charge != all_data[i-1].charge
                        continue
                    end

                    d1 = all_data[i].atomic_separation
                    d2 = all_data[i-1].atomic_separation
                    Δd = d2 - d1

                    set_diatomic_system_to_parsed_file!(
                        simulation[thread_id],all_data[i-1])

                    system = simulation[thread_id].system
                    molecule = system.molecules[1]
                    atom_1 = molecule.atoms[1]
                    atom_2 = molecule.atoms[2]

                    ζ1_nxt = copy(atom_1.polarization_coefficient)
                    ζ2_nxt = copy(atom_2.polarization_coefficient)
                    μ_nxt = copy(system.chemical_potential)

                    polarize_molecules!(simulation[thread_id])
                    ζ1_model_nxt = copy(atom_1.polarization_coefficient)
                    ζ2_model_nxt = copy(atom_2.polarization_coefficient)
                    μ_model_nxt = copy(system.chemical_potential)
                        
                    diff_vec_nxt = zeros(aux_type,3)
                    diff_vec_nxt[1] = ζ1_nxt - ζ1_model_nxt
                    diff_vec_nxt[2] = ζ2_nxt - ζ2_model_nxt
                    diff_vec_nxt[3] = μ_nxt - μ_model_nxt
                    
                    ret_val[thread_id] += (norm(diff_vec_nxt-diff_vec) / Δd)^2
                end
            end
        end

        ret_val = sum(ret_val) / length(all_data)

        # This makes sure that the charges predicted outside of the calibration 
        # dataset are still all positive.
        charge_cond_sum = zeros(aux_type,n_threads)
        r = 0.0:((6.0-0.0)/300.0):6.0
        @threads for thread_id in 1:n_threads
            set_diatomic_system_to_parsed_file!(
                simulation[thread_id],all_data[1])
                atom_1 = simulation[thread_id].system.molecules[1].atoms[1]
                atom_2 = simulation[thread_id].system.molecules[1].atoms[2]
                atom_1.coordinates .= 0.0
                atom_2.coordinates .= 0.0

            for atomic_separation in r[thread_id:n_threads:end]
                atom_1.coordinates[3] = atomic_separation

                aux_m, aux_y = 
                    polarization_matrix_problem(simulation[thread_id])
                aux_x = (aux_m \ aux_y)[1:2]

                charge_cond_sum[thread_id] += 
                    norm(sum(aux_x) - sum(abs.(aux_x)))^2
            end
        end
        ret_val += sum(charge_cond_sum) / length(charge_cond_sum)

        return ret_val
    end

    num_vars = 4
    aux_X = 2.0 .* (rand(Float64,num_vars) .- 0.5)

    for i in 1:2500
        cost_func_eval = cost_func(aux_X)
        new_aux_X = 2.0 .* (rand(Float64,num_vars) .- 0.5)

        if cost_func_eval <= 0.1
            println("Initial guess found!")
            break
        end

        if cost_func(new_aux_X) < cost_func_eval
            aux_X = new_aux_X
            print(@sprintf "Current best %18.6E \n" cost_func(aux_X))
        end
    end

    needs_casting = true
    sol = Optim.optimize(cost_func, aux_X, LBFGS(), autodiff=:forward,
        Optim.Options(show_trace=true,iterations=2000))
    aux_X = Optim.minimizer(sol)

    display(aux_X)

    simulation = make_system_from_parsed_file(all_data[1])
    set_fitted_pol_e_coeffs!(simulation,Z1,Z2,aux_X)

    save_fitted_coeffs(simulation)
end
