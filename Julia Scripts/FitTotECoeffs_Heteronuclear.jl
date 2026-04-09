using Optim
using Printf, LinearAlgebra
using Base.Threads

include("FitCoeffs_General.jl")

which_atomic_numbers = [[1,6],[1,8],[6,8],[7,8]]
for atomic_numbers in which_atomic_numbers
    Z1 = atomic_numbers[1]
    Z2 = atomic_numbers[2]

    local neutral_data, cation_data, anion_data = 
        read_all_sanitized_data(Z1,Z2)
    all_data = vcat(neutral_data, cation_data, anion_data)

    n_threads = Threads.nthreads()
    simulation = Vector{SimulationSystem}()
    resize!(simulation,n_threads)
    for i in 1:n_threads
        simulation[i] = make_system_from_parsed_file(all_data[1])
    end

    neutral_at1 = make_monoatomic_system(Z1,0)
    neutral_at2 = make_monoatomic_system(Z2,0)

    dft_neutral_at_e, _, _ = get_reference_atom_total_energy()
    dft_neutral_at1_e = dft_neutral_at_e[Z1]
    dft_neutral_at2_e = dft_neutral_at_e[Z2]

    polarize_to_model!(neutral_data)
    polarize_to_model!(cation_data)
    polarize_to_model!(anion_data)

    neutral_min_index = 1
    for index in eachindex(neutral_data)
        if neutral_data[neutral_min_index].total_energy > neutral_data[index].total_energy
            neutral_min_index = index
        end
    end

    reference_system = neutral_data[neutral_min_index]

    needs_casting = true
    function cost_func(aux_X::Vector)
        aux_type = typeof(aux_X[1])

        if needs_casting
            needs_casting = false
            function cast_coeffs_to_type!(simulation::SimulationSystem)
                cast_tot_e_coeffs_to_type!(simulation,aux_type)
            end

            for thread_id in 1:n_threads
                cast_coeffs_to_type!(simulation[thread_id])
            end

            cast_coeffs_to_type!(neutral_at1)
            cast_coeffs_to_type!(neutral_at2)
        end

        function set_fitted_coeffs!(simulation::SimulationSystem)
            set_fitted_tot_e_coeffs!(simulation,Z1,Z2,aux_X)
        end

        set_fitted_coeffs!(neutral_at1)
        set_fitted_coeffs!(neutral_at2)

        model_at1_neutral_e = total_energy(neutral_at1)
        model_at2_neutral_e = total_energy(neutral_at2)

        set_fitted_coeffs!(simulation[1])
        set_diatomic_system_to_parsed_file!(simulation[1],reference_system)

        dft_e0 = reference_system.total_energy
        model_e0 = total_energy(simulation[1])

        ret_val = zeros(aux_type,n_threads)
        @threads for thread_id in 1:n_threads
            set_fitted_coeffs!(simulation[thread_id])

            for i in thread_id:n_threads:length(all_data)
                set_diatomic_system_to_parsed_file!(
                    simulation[thread_id],all_data[i])

                dft_tot_e = all_data[i].total_energy
                model_tot_e = total_energy(simulation[thread_id])

                dft_e_diff = dft_tot_e - dft_e0
                model_e_diff = model_tot_e - model_e0
                ret_val[thread_id] += (model_e_diff - dft_e_diff)^2

                if i > 1
                    if all_data[i].charge != all_data[i-1].charge
                        continue
                    end

                    d1 = all_data[i].atomic_separation
                    d2 = all_data[i-1].atomic_separation
                    Δd = d2 - d1

                    set_diatomic_system_to_parsed_file!(
                        simulation[thread_id],all_data[i-1])

                    dft_tot_e_nxt = all_data[i-1].total_energy
                    model_tot_e_nxt = total_energy(simulation[thread_id])

                    dft_ΔE = dft_tot_e_nxt - dft_tot_e
                    model_ΔE = model_tot_e_nxt - model_tot_e

                    dft_dev = dft_ΔE/Δd
                    model_dev = model_ΔE/Δd

                    ret_val[thread_id] += (model_dev - dft_dev)^2
                end
            end
        end

        return sum(ret_val) / length(all_data)
    end

    num_vars = 7
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
    sol = Optim.optimize(cost_func, aux_X[:], LBFGS(), autodiff=:forward,
        Optim.Options(show_trace=true,iterations=2000))
    aux_X = Optim.minimizer(sol)

    display(aux_X)

    simulation = make_system_from_parsed_file(all_data[1])
    set_fitted_tot_e_coeffs!(simulation,Z1,Z2,aux_X)

    save_fitted_coeffs(simulation)
end
