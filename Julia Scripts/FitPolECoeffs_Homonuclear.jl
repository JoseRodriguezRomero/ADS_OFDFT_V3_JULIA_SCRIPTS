using Optim, Plots, Evolutionary, Distributed;
using Printf, LinearAlgebra;
using LaTeXStrings, Latexify, Measures;
using Base.Threads;

include("FitHomonuclearCoeffs_General.jl")

which_atomic_numbers = [1,6,7,8];
for atomic_number in which_atomic_numbers
    neutral_data, cation_data, anion_data = readAllSanitizedData(atomic_number);
    all_data = vcat(neutral_data, cation_data, anion_data);

    num_vars = 2*num_1b_coeffs + 2*num_2b_coeffs;
    aux_X = zeros(Float64,num_vars);

    at_neutral = MakeAtom(atomic_number,0);
    at_cation = MakeAtom(atomic_number,1);
    at_anion = MakeAtom(atomic_number,-1);

    SetAtomEnergy!(at_neutral);
    SetAtomEnergy!(at_cation);
    SetAtomEnergy!(at_anion);

    all_atoms = [at_neutral, at_cation, at_anion];

    # Load the previously fitted coefficients
    global pol_all_coeffs = LoadPolCoeffs();

    all_molecs = Vector{Molecule}();
    resize!(all_molecs,length(all_data));
    for i in eachindex(all_data)
        all_molecs[i] = MakeMoleculeFromParsedFile(all_data[i]);
    end

    Z1_eff = all_molecs[1].atoms_data[1,4];
    Z2_eff = all_molecs[1].atoms_data[2,4];

    function CostFunc(aux_X::Vector)
        global pol_all_coeffs;
        aux_type = typeof(aux_X[1]);
        aux_pol_all_coeffs = Vector{Matrix}();
        resize!(aux_pol_all_coeffs,2);

        aux_pol_all_coeffs[1] = aux_type.(pol_all_coeffs[1]);
        aux_pol_all_coeffs[2] = aux_type.(pol_all_coeffs[2]);
        SetFittedCoeffs!(aux_pol_all_coeffs,atomic_number,aux_X);

        n_threads = Threads.nthreads();
        ret_val = zeros(aux_type,n_threads);
        ret_val = zeros(aux_type,n_threads);
        @threads for thread_id in 1:n_threads
            for i in thread_id:n_threads:length(all_molecs)
                _, _, aux_m, aux_y = 
                    MatPolarizeMolecules(all_molecs[i],aux_pol_all_coeffs);

                aux_x = zeros(aux_type,3);
                aux_x[1] = all_molecs[i].cloud_data[1,6];
                aux_x[2] = all_molecs[i].cloud_data[4,6];
                aux_x[3] = all_molecs[i].chem_μ;

                aux_diff = (aux_m * aux_x) - aux_y;
                ret_val[thread_id] += dot(aux_diff,aux_diff);
            end
        end

        ret_val = sum(ret_val);
        for atom in all_atoms
            _, _, aux_m, aux_y = 
                MatPolarizeMolecules(atom,aux_pol_all_coeffs);

            aux_x = zeros(aux_type,2);
            aux_x[1] = atom.cloud_data[1,6];
            aux_x[2] = atom.chem_μ;

            aux_diff = (aux_m * aux_x) - aux_y;
            ret_val += dot(aux_diff,aux_diff);
        end

        return ret_val;
    end

    sol = Optim.optimize(CostFunc, aux_X[:], LBFGS(), autodiff=:forward,
        Optim.Options(show_trace=true,iterations=8000));
    aux_X = Optim.minimizer(sol);

    # Save the currently fitted coefficients
    SetFittedCoeffs!(pol_all_coeffs,atomic_number,aux_X);

    # Save the newly fitted coefficients
    SavePolCoeffs();
end

