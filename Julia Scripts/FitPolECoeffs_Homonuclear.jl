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

    at_neutral = make_atom_system(atomic_number,0);
    at_cation = make_atom_system(atomic_number,1);
    at_anion = make_atom_system(atomic_number,-1);

    neutral_atom_chem_μ, cation_atom_chem_μ, anion_atom_chem_μ = 
        get_reference_atom_chemical_potential();

    neutral_atom_chem_μ = neutral_atom_chem_μ[atomic_number];
    cation_atom_chem_μ = cation_atom_chem_μ[atomic_number];
    anion_atom_chem_μ = anion_atom_chem_μ[atomic_number];

    at_neutral.system.chemical_potential = neutral_atom_chem_μ;
    at_cation.system.chemical_potential = cation_atom_chem_μ;
    at_anion.system.chemical_potential = anion_atom_chem_μ;

    all_atoms = [at_neutral, at_cation, at_anion];

    n_threads = Threads.nthreads();
    simulation = Vector{SimulationSystem}();
    resize!(simulation,n_threads);
    for i in 1:n_threads
        simulation[i] = make_system_from_parsed_file(all_data[1]);
    end

    num_1b_coeffs = 4;
    num_2b_coeffs = 4;

    num_vars = 2*num_1b_coeffs + 2*num_2b_coeffs;
    aux_X = zeros(Float64,num_vars);

    needs_casting = true;
    function cost_func(aux_X::Vector)
        aux_type = typeof(aux_X[1]);

        if needs_casting
            needs_casting = false;
            for thread_id in 1:n_threads
                cast_pol_e_coeffs_to_type!(simulation[thread_id],aux_type);
            end

            for i in eachindex(all_atoms)
                cast_pol_e_coeffs_to_type!(all_atoms[i],aux_type);
            end
        end

        ret_val = zeros(aux_type,n_threads);
        @threads for thread_id in 1:n_threads
            # Set the trial coefficients in the simulation structure.
            set_fitted_pol_e_coeffs!(simulation[thread_id],atomic_number,aux_X);
            cast_type = simulation[thread_id].cast_types.cast_to_xc_coeff_type;

            # Calculate the error when using these trial coefficients.
            for i in thread_id:n_threads:length(all_data)
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

        ret_val = sum(ret_val);
        for atom in all_atoms
            set_fitted_pol_e_coeffs!(atom,atomic_number,aux_X);
            cast_type = atom.cast_types.cast_to_xc_coeff_type;
            aux_m, aux_y = polarization_matrix_problem(atom,cast_type);

            ζ = atom.system.molecules[1].cloud_data[1,6];

            aux_x = zeros(aux_type,2);
            aux_x[1] = ζ;
            aux_x[2] = atom.system.chemical_potential;

            aux_diff = (aux_m * aux_x) - aux_y;
            ret_val += dot(aux_diff,aux_diff);
        end

        return ret_val;
    end

    sol = Optim.optimize(cost_func, aux_X, LBFGS(), autodiff=:forward,
        Optim.Options(show_trace=true,iterations=8000));
    aux_X = Optim.minimizer(sol);

    simulation = make_system_from_parsed_file(all_data[1]);
    set_fitted_pol_e_coeffs!(simulation,atomic_number,aux_X);

    save_fitted_coeffs(simulation);
end

