using Optim, Plots;
using Printf, LinearAlgebra;
using LaTeXStrings, Latexify, Measures;
using Base.Threads;

include("FitCoeffs_General.jl")

which_atomic_numbers = [1,6,7,8];
for atomic_number in which_atomic_numbers
    neutral_data, cation_data, anion_data = 
        read_all_sanitized_data(atomic_number,true);
    all_data = vcat(neutral_data, cation_data, anion_data);

    at_neutral = make_monoatomic_system(atomic_number,0);
    at_cation = make_monoatomic_system(atomic_number,1);
    at_anion = make_monoatomic_system(atomic_number,-1);

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

    atoms_μ0 = simulation[1].basis_set_settings.atoms_μ0;
    μ0 = atoms_μ0[atomic_number];

    num_vars = 14;
    aux_X = zeros(Float64,num_vars);

    # z1_eff = simulation[1].system.molecules[1].atoms[1].valence_electrons;
    # z2_eff = simulation[1].system.molecules[1].atoms[2].valence_electrons;

    needs_casting = true;
    function cost_func(aux_X::Vector)
        aux_type = typeof(aux_X[1]);

        function set_fitted_coeffs!(simulation::SimulationSystem)
            set_fitted_pol_e_coeffs!(simulation,atomic_number,aux_X);
        end

        if needs_casting == true
            needs_casting = false;

            for thread_id in 1:n_threads
                cast_pol_e_coeffs_to_type!(simulation[thread_id],aux_type);
            end

            for i in eachindex(all_atoms)
                cast_pol_e_coeffs_to_type!(all_atoms[i],aux_type);
            end
        end

        set_fitted_coeffs!(at_neutral);
        set_fitted_coeffs!(at_cation);
        set_fitted_coeffs!(at_anion);

        model_at_neutral_e = total_energy(at_neutral)
        model_at_cation_e = total_energy(at_cation);
        model_at_anion_e = total_energy(at_anion);
    
        ret_val = zeros(aux_type,n_threads);
        @threads for thread_id in 1:n_threads
            # Set the trial coefficients in the simulation structure.
            set_fitted_coeffs!(simulation[thread_id]);

            # Calculate the error when using these trial coefficients.
            for i in thread_id:n_threads:length(all_data)
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
                aux_x[3] = μ - μ0;
                aux_x[4] = μ - μ0;
                aux_x[5] = μ;
                ret_val[thread_id] += norm((aux_m \ aux_y) - aux_x)^2;
            end
        end

        ret_val = sum(ret_val);
        for atom in all_atoms
            set_fitted_coeffs!(atom);

            aux_m, aux_y = polarization_matrix_problem(atom);

            ζ = atom_polarization_coeff(atom.system.molecules[1],1);
            μ = copy(atom.system.chemical_potential);

            aux_x = zeros(aux_type,3);
            aux_x[1] = ζ;
            aux_x[2] = μ - μ0;
            aux_x[3] = μ;
            ret_val += norm((aux_m \ aux_y) - aux_x)^2;

            if atom.system.charge == 0
                ret_val += aux_x[2]^2;
            end
        end

        return ret_val;
    end

    # needs_casting = false;
    # sol = Optim.optimize(cost_func, aux_X[:], NelderMead(),
    #     Optim.Options(show_trace=true,iterations=8000));
    # aux_X = Optim.minimizer(sol);
    
    needs_casting = true;
    sol = Optim.optimize(cost_func, aux_X[:], LBFGS(), autodiff=:forward,
        Optim.Options(show_trace=true,iterations=2000));
    aux_X = Optim.minimizer(sol);

    display(aux_X);
    
    simulation = make_monoatomic_system(atomic_number,0);
    set_fitted_pol_e_coeffs!(simulation,atomic_number,aux_X);

    save_fitted_coeffs(simulation);
end
