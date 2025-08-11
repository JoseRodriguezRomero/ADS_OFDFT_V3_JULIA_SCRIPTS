using Optim, Plots, Evolutionary, Distributed;
using Printf, LinearAlgebra;
using LaTeXStrings, Latexify, Measures;

include("FitHomonuclearCoeffs_General.jl")

which_atomic_numbers = [1,6,7,8];
for atomic_number in which_atomic_numbers
    neutral_data, cation_data, anion_data = readAllSanitizedData(atomic_number);
    all_data = vcat(neutral_data, cation_data, anion_data);

    num_vars = 2*num_1b_coeffs + 2*num_2b_coeffs;
    
    # Load the previously fitted coefficients
    global all_coeffs = LoadCoeffs();

    # Optim way
    all_molecs = Vector{Molecule}();
    resize!(all_molecs,length(all_data));
    for i in eachindex(all_data)
        all_molecs[i] = MakeMoleculeFromParsedFile(all_data[i]);
    end

    neutral_at = MakeAtom(atomic_number,0);
    cation_at = MakeAtom(atomic_number,1);
    anion_at = MakeAtom(atomic_number,-1);

    SetAtomEnergy!(neutral_at);
    SetAtomEnergy!(cation_at);
    SetAtomEnergy!(anion_at);

    function CostFunc(aux_X::Vector)
        global all_coeffs, pol_all_coeffs;
        aux_type = typeof(aux_X[1]);

        aux_all_coeffs = Vector{Matrix}();
        resize!(aux_all_coeffs,2);

        aux_all_coeffs[1] = aux_type.(all_coeffs[1]);
        aux_all_coeffs[2] = aux_type.(all_coeffs[2]);
        SetFittedCoeffs!(aux_all_coeffs,atomic_number,aux_X);

        function GetRefEnergy(mol::Molecule)
            energy = -mol.energy;
            energy += NaiveEnergyFromDensity(mol,aux_all_coeffs);
            energy += XCEnergyFromDensity(mol,aux_all_coeffs);
            return energy;
        end

        # Isolated Atom energies
        neutral_at_e = GetRefEnergy(neutral_at);
        cation_at_e = GetRefEnergy(cation_at);
        anion_at_e = GetRefEnergy(anion_at);

        # Reference Molecule energies
        neutral_mol_e = GetRefEnergy(MakeMoleculeFromParsedFile(neutral_data[end]));
        cation_mol_e = GetRefEnergy(MakeMoleculeFromParsedFile(cation_data[end]));
        anion_mol_e = GetRefEnergy(MakeMoleculeFromParsedFile(anion_data[end]));

        ret_val = 0.0;
        for i in eachindex(all_molecs)
            aux_e = GetRefEnergy(all_molecs[i]);
        
            if all_molecs[i].charge == 0
                ret_val += (aux_e - neutral_mol_e)^2.0;
                ret_val += (aux_e - neutral_at_e - neutral_at_e)^2.0;
            elseif all_molecs[i].charge == 1
                ret_val += (aux_e - cation_mol_e)^2.0;
                ret_val += (aux_e - neutral_at_e - cation_at_e)^2.0;
            else
                ret_val += (aux_e - anion_mol_e)^2.0;
                ret_val += (aux_e - neutral_at_e - anion_at_e)^2.0;
            end
        end

        ret_val += (neutral_at_e - cation_at_e)^2.0;
        ret_val += (neutral_at_e - anion_at_e)^2.0;

        return ret_val;
    end

    aux_X = zeros(Float64,num_vars);
    sol = Optim.optimize(CostFunc, aux_X[:], LBFGS(), autodiff=:forward,
        Optim.Options(show_trace=true));
    aux_X = Optim.minimizer(sol);

    # Set the currently fitted coefficients
    SetFittedCoeffs!(all_coeffs,atomic_number,aux_X);

    # Save the newly fitted coefficients
    SaveCoeffs();
end
