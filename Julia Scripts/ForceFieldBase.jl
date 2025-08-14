using Optim;
using Integrals;
using ForwardDiff;
using SpecialFunctions;
using LinearAlgebra, GenericLinearAlgebra;

mutable struct Molecule
    name::String;
    atoms_data::Matrix;
    cloud_data::Matrix;
end

mutable struct MolecularSystem
    name::String;
    molecules::Vector{Molecule};
    charge::Int;
    energy::Number;
    chemical_potential::Number;
end

mutable struct EmpiricalCoefficients
    max_atomic_number::Int;

    xc_a_1b::Matrix{Number};
    xc_c_1b::Matrix{Number};
    xc_e_1b::Matrix{Number};
    xc_f_1b::Matrix{Number};

    xc_a_2b::Matrix{Number};
    xc_b_2b::Matrix{Number};
    xc_c_2b::Matrix{Number};
    xc_d_2b::Matrix{Number};
end

struct EmpiricalCoefficientMappings
    max_atomic_number::Int;
    coeff_1b_combinations::Int;
    coeff_2b_combinations::Int;
    coeff_1b_map::Dict{Int,Int};
    inv_coeff_1b_map::Dict{Int,Int};
    coeff_2b_map::Dict{Tuple{Int,Int},Int};
    inv_coeff_2b_map::Dict{Int,Tuple{Int,Int}};
end

struct BasisSetSettings
    clouds_per_atom::Int;
    max_atomic_number::Int;
    element_cloud_coeffs::Matrix{Float64};
    thomas_fermi_kes::Vector{Float64};
    von_weizsacker_kes::Vector{Float64};
end

struct SimulationCastTypes
    no_cast::Int;
    cast_to_xc_coeff_type::Int;
    cast_to_pol_coeff_type::Int;
end

mutable struct SimulationSystem
    system::MolecularSystem;
    pol_e_xc_coeffs::EmpiricalCoefficients;
    tot_e_xc_coeffs::EmpiricalCoefficients;
    coeff_mappings::EmpiricalCoefficientMappings;
    basis_set_settings::BasisSetSettings;
    cast_types::SimulationCastTypes;
end

function Base.copy(mol::Molecule)
    return Molecule(
        mol.name,
        copy(mol.atoms_data),
        copy(mol.cloud_data)
    )
end

function Molecule(type = Float64)
    name = "Kryptonite";
    atoms_data = zeros(type,0,6);
    cloud_data = zeros(type,0,6);
    return Molecule(name,atoms_data,cloud_data);
end

function MolecularSystem(type = Float64)
    name = "Kryptonite++";
    molecules = [Molecule(type)];
    charge = 0;
    energy = 0;
    chemical_potential = 0;

    return MolecularSystem(name, molecules, charge, 
        energy, chemical_potential);
end

function read_molecule(molec_data::String, type = Float64)
    file_id = open(molec_data,"r");
    molecule = Molecule(type);

    readline(file_id);

    while true
        line = readline(file_id);
        line_splitted = split(line);

        if length(line_splitted) < 4
            break;
        end

        x = type(parse(Float64,line_splitted[1]));
        y = type(parse(Float64,line_splitted[2]));
        z = type(parse(Float64,line_splitted[3]));
        q = type(parse(Float64,line_splitted[4]));
        Z = type(parse(Float64,line_splitted[5]));
        
        new_atom = zeros(type,1,6);
        new_atom[1] = x;
        new_atom[2] = y;
        new_atom[3] = z;
        new_atom[4] = q;            # Effective Atomic Number 
        new_atom[5] = Z;            # Atomic Number
        new_atom[6] = type(1.0);    # Polarization coefficient

        molecule.atoms_data = vcat(molecule.atoms_data,new_atom);
    end

    readline(file_id);
    while true
        line = readline(file_id);
        line_splitted = split(line);

        if length(line_splitted) < 5
            break;
        end

        x = type(parse(Float64,line_splitted[1]));
        y = type(parse(Float64,line_splitted[2]));
        z = type(parse(Float64,line_splitted[3]));
        c = type(parse(Float64,line_splitted[4]));
        λ = type(parse(Float64,line_splitted[5]));

        new_cloud = zeros(type,1,6)
        new_cloud[1] = x;
        new_cloud[2] = y;
        new_cloud[3] = z;
        new_cloud[4] = c;         # Magnitude
        new_cloud[5] = λ;         # Decay coefficient
        new_cloud[6] = type(1.0); # Polarization coefficient

        molecule.cloud_data = vcat(molecule.cloud_data,new_cloud);
    end

    close(file_id);
    return molecule;
end

function number_of_atoms(molecule::Molecule)
    # Returns the number of atoms in the molecule.
    return length(molecule.atoms_data[:,1]);
end

function number_of_clouds(molecule::Molecule)
    # Returns the number of clouds in the molecule.
    return length(molecule.cloud_data[:,1]);
end

function clouds_per_atom(molecule::Molecule)
    # Returns the ratio between clouds and atoms in the molecule.
    return round(Int,number_of_clouds(molecule) / number_of_atoms(molecule));
end

function atom_atomic_number(molecule::Molecule, i::Int)
    # Returns the atomic number of the iᵗʰ atom.
    return round(Int, molecule.atoms_data[i,5]);
end

function atom_eff_atomic_number(molecule::Molecule, i::Int)
    # Returns the effective atomic number of the iᵗʰ atom.
    return round(Int, ForwardDiff.value(molecule.atoms_data[i,4]));
end

function atom_position(molecule::Molecule, i::Int)
    # Returns the position vector of the iᵗʰ atom.
    return molecule.atoms_data[i,1:3];
end

function atom_polarization_coeff(molecule::Molecule, i::Int)
    # Returns the polarization coefficient of the iᵗʰ atom.
    return molecule.atoms_data[i,6];

    # _clouds_per_atom_ = clouds_per_atom(molecule);
    # return molecule.cloud_data[i*_clouds_per_atom_,6];
end

function cloud_amplitude_coeff(molecule::Molecule, i::Int)
    # Returns the amplitude coefficient of the iᵗʰ cloud.
    return molecule.cloud_data[i,4];
end

function cloud_decay_coeff(molecule::Molecule, i::Int)
    # Returns the decay coefficient of the iᵗʰ cloud.
    return molecule.cloud_data[i,5];
end

function cloud_polarization_coeff(molecule::Molecule, i::Int)
    # Returns the polarization coefficient of the iᵗʰ cloud.
    return molecule.cloud_data[i,6];
end

function cloud_position(molecule::Molecule, i::Int)
    # Returns the position vector of the iᵗʰ cloud.
    return molecule.cloud_data[i,1:3];
end

function cloud_atomic_number(molecule::Molecule, i::Int)
    # Returns the atomic number of the atom that the iᵗʰ cloud belongs to.
    atom_index = ceil(Int,i/clouds_per_atom(molecule));
    return atom_atomic_number(molecule,atom_index);
end

function set_atom_pol_coeff!(molecule::Molecule, ζ::Number, i::Int)
    # Sets the atom iᵗʰ atom's polarization coefficient and its corresponding 
    # electron cloud.
    molecule.atoms_data[i,6] = ζ;

    i0 = (i-1)*clouds_per_atom(molecule);
    for i in 1:clouds_per_atom(molecule)
        molecule.cloud_data[i0+i,6] = ζ;
    end
end

function set_atom_partial_charge!(molecule::Molecule, ρ::Number, i::Int)
    # Sets the atom iᵗʰ atom's polarization coefficient and its corresponding 
    # electron cloud so that they match the desired partial charge.
    z_eff = atom_eff_atomic_number(molecule,i);
    ζ = (z_eff - ρ) / z_eff;
    set_atom_pol_coeff!(molecule,ζ,i);
end

function set_mol_pol_coeffs!(molecule::Molecule, ζ::Vector)
    # Sets the polarization coefficients of the atoms of the molecules. If the
    # coefficient vector's size is greater than the number of atoms the 
    # program will crash.
    for i in eachindex(ζ)
        set_atom_pol_coeff!(molecule,ζ[i],i);
    end
end

function set_mol_partial_charges!(molecule::Molecule, ρ::Vector)
    # Sets the polarization coefficients of the atoms of the molecules. If the
    # coefficient vector's size is greater than the number of atoms the 
    # program will crash.
    for i in eachindex(ζ)
        set_atom_partial_charge!(molecule,ρ[i],i);
    end
end

function set_atom_position!(molecule::Molecule, r::Vector, i::Int)
    # Sets the position of the iᵗʰ atom and its corresponding electron cloud.
    molecule.atoms_data[i,1:3] = r;

    i0 = (i-1)*clouds_per_atom(molecule);
    for i in 1:clouds_per_atom(molecule)
        molecule.cloud_data[i0+i,1:3] = r;
    end
end

function naive_coulomb_integral(λ::Real, d::Real)
    if (abs(d) < 1.0E-6)
        return 2*sqrt(λ/π);
    end

    erf_arg = d*sqrt(λ);
    if erf_arg > 20.0
        return 1.0/d;
    end

    return erf(erf_arg)/d;
end 

function xc_sph(λ::Real, d::Real)
    if (abs(d) < 1.0E-6)
        return -((4*λ^(3/2)*exp(-λ))/sqrt(π));
    end

    erf_arg = d*sqrt(λ);
    if erf_arg > 20.0
        return 0.0;
    end

    return (1/d)*sqrt(λ/π)*((1-exp(4*λ*d))/(exp(λ*(d+1)^2)));
end

function xc_cyl(λ::Real, d::Real)
    erf_arg = d*sqrt(λ);
    if erf_arg > 20.0
        return 0.0;
    end

    return (1/d)*sqrt(λ/π)*((((exp(2*λ*d)-1)^2)/(exp(λ*(d+1)^2))) + 
        ((2*(exp(-λ)-1))/(exp(λ*d^2))));
end

function ee_energy(mol1::Molecule, mol2::Molecule, i::Int, j::Int)
    # Returns the naive and xc electron-electron energy functionals between  
    # the iᵗʰ cloud of mol1 and the jᵗʰ cloud of mol2.
    c1 = cloud_amplitude_coeff(mol1,i);
    λ1 = cloud_decay_coeff(mol1,i)
    r1 = cloud_position(mol1,i);

    c2 = cloud_amplitude_coeff(mol2,j);
    λ2 = cloud_decay_coeff(mol2,j)
    r2 = cloud_position(mol2,j);

    λ = (λ1*λ2)/(λ1+λ2);
    d = norm(r1-r2);
    c = c1*c2;

    ee_naive = c*naive_coulomb_integral(λ,d);
    ee_xc_sph = c*xc_sph(λ,d);
    ee_xc_cyl = c*xc_cyl(λ,d);

    if (d < 1.0E-6)
        ee_xc_cyl = Inf;
    end

    return ee_naive, ee_xc_sph, ee_xc_cyl;
end

function en_energy(mol1::Molecule, mol2::Molecule, i::Int, j::Int)
    # Returns the naive and xc electron-nuclei energy functionals between  
    # the iᵗʰ cloud of mol1 and the jᵗʰ atom of mol2.
    c1 = cloud_amplitude_coeff(mol1,i);
    λ1 = cloud_decay_coeff(mol1,i)
    r1 = cloud_position(mol1,i);

    z2_eff = atom_eff_atomic_number(mol2,j);
    r2 = atom_position(mol2,j);

    λ = λ1;
    d = norm(r1-r2);
    c = c1*z2_eff;

    en_naive = c*naive_coulomb_integral(λ,d);
    en_xc_sph = c*xc_sph(λ,d);
    en_xc_cyl = c*xc_cyl(λ,d);

    if (d < 1.0E-6)
        en_xc_cyl = Inf;
    end

    return en_naive, en_xc_sph, en_xc_cyl;
end

function move_molecule!(molecule::Molecule, ΔR::Vector{Number})
    # Moves the whole molecule by the vector ΔR.
    molecule.atoms_data[:,1:3] .+= ΔR;
    molecule.cloud_data[:,1:3] .+= ΔR;
end

function rotate_molecule!(molecule::Molecule, Δθ::Matrix{Number})
    # Rotates the whole molecule by the rotation matrix Δθ.
    molecule.atoms_data[:,1:3] = molecule.atoms_data[:,1:3]*Δθ;
    molecule.cloud_data[:,1:3] = molecule.cloud_data[:,1:3]*Δθ;
end

function polarization_matrix_problem(simulation::SimulationSystem, 
        cast_type::Int = simulation.cast_types.no_cast)
    # Calculates the matrix problem that needs to be solved for the 
    # polarization lagrangian to be stationary.
    xc_a_1b = simulation.pol_e_xc_coeffs.xc_a_1b;
    xc_c_1b = simulation.pol_e_xc_coeffs.xc_c_1b;
    xc_e_1b = simulation.pol_e_xc_coeffs.xc_e_1b;
    xc_f_1b = simulation.pol_e_xc_coeffs.xc_f_1b;

    xc_a_2b = simulation.pol_e_xc_coeffs.xc_a_2b;
    xc_b_2b = simulation.pol_e_xc_coeffs.xc_b_2b;
    xc_c_2b = simulation.pol_e_xc_coeffs.xc_c_2b;
    xc_d_2b = simulation.pol_e_xc_coeffs.xc_d_2b;

    coeff_1b_map = simulation.coeff_mappings.coeff_1b_map;
    coeff_2b_map = simulation.coeff_mappings.coeff_2b_map;

    unpol_tf_ke = simulation.basis_set_settings.thomas_fermi_kes;
    unpol_vw_ke = simulation.basis_set_settings.von_weizsacker_kes;

    molecules = simulation.system.molecules;
    charge = simulation.system.charge;

    aux_type = Float64;

    if cast_type == simulation.cast_types.cast_to_xc_coeff_type
        aux_type = typeof(xc_a_1b[1]);
    elseif cast_type == simulation.cast_types.cast_to_pol_coeff_type
        aux_type = typeof(simulation.system.chemical_potential);
    end

    function comb_coeffs_1b(coeffs::AbstractMatrix, Z::Int, ζ::Number)
        # Calculates the XC coefficient based on the polarization coefficients
        # (One body cases).
        b = coeffs[coeff_1b_map[Z],1];
        m = coeffs[coeff_1b_map[Z],2];

        return 2*m*ζ + b;
    end

    function comb_coeffs_2b(coeffs::AbstractMatrix, Z1::Int, Z2::Int,
        ζ1::Number, ζ2::Number)
        # Calculates the XC coefficient based on the polarization coefficients.
        # (One body cases)
        b = coeffs[coeff_2b_map[(Z1,Z2)],1];
        m1 = coeffs[coeff_2b_map[(Z1,Z2)],2];
        m2 = coeffs[coeff_2b_map[(Z1,Z2)],3];

        if (Z1 < Z2)
            return 2*m1*ζ1 + m2*ζ2 + b;
        else
            return m1*ζ2 + 2*m2*ζ1 + b;
        end
    end

    tot_num_atoms = 0;
    tot_num_clouds = 0;
    atom_ind_base = Vector{Int}();
    for molecule in molecules
        if isempty(atom_ind_base)
            push!(atom_ind_base,1);
            tot_num_atoms += number_of_atoms(molecule);
            tot_num_clouds += number_of_clouds(molecule);
        else
            aux_num_atoms = number_of_atoms(molecule);
            aux_num_clouds = number_of_clouds(molecule);
            
            push!(atom_ind_base, tot_num_atoms + 1);
            tot_num_atoms += aux_num_atoms;
            tot_num_clouds += aux_num_clouds;
        end
    end

    num_vars = tot_num_atoms + 1;
    aux_M = zeros(aux_type, num_vars, num_vars);
    aux_Y = zeros(aux_type, num_vars, 1);
    aux_Y[end] -= charge;

    clouds_per_atom = simulation.basis_set_settings.clouds_per_atom;

    for ii in eachindex(molecules)
        molecule1 = molecules[ii];
        num_atoms1 = number_of_atoms(molecule1);
        num_clouds1 = number_of_clouds(molecule1);

        # CHEMICAL POTENTIAL
        # Lagrange multiplier to keep the number of electrons constant in the 
        # new fictitious electron density.
        ind1 = atom_ind_base[ii];
        ind2 = ind1 + num_atoms1 - 1;

        aux_M[ind1:ind2,end] = -molecule1.atoms_data[:,4];
        aux_M[end,ind1:ind2] = molecule1.atoms_data[:,4];
        aux_Y[end] += round(Int, sum(molecule1.atoms_data[:,4]));

        # Kinetic Energy Contributions
        for i in 1:num_atoms1
            z1 = atom_atomic_number(molecule1,i);
            ζ1 = atom_polarization_coeff(molecule1,i);

            ii0 = atom_ind_base[ii] + i - 1;

            tf_fit_coeff = aux_type(0.0);
            tf_fit_coeff += (5.0/3.0)*xc_e_1b[z1,1]*(ζ1^(2.0/3.0));
            tf_fit_coeff += (8.0/3.0)*xc_e_1b[z1,2]*(ζ1^(5.0/3.0));

            vw_fit_coeff = aux_type(0.0);
            vw_fit_coeff += 1.0*xc_f_1b[z1,1]*(ζ1^(0.0));
            vw_fit_coeff += 2.0*xc_f_1b[z1,2]*(ζ1^(1.0));

            aux_M[ii0,ii0] += tf_fit_coeff*unpol_tf_ke[z1];
            aux_M[ii0,ii0] += vw_fit_coeff*unpol_vw_ke[z1];
        end

        # Fill the rest of the matrix
        for jj in eachindex(molecules)
            molecule2 = molecules[jj];
            num_atoms2 = number_of_atoms(molecule2);
            num_clouds2 = number_of_clouds(molecule2);

            # cloud-nuclei
            for i in 1:num_clouds1
                for j in 1:num_atoms2
                    Z1 = cloud_atomic_number(molecule1,i);
                    Z2 = atom_atomic_number(molecule2,j);

                    ζ1 = cloud_polarization_coeff(molecule1,i);
                    ζ2 = atom_polarization_coeff(molecule2,j);

                    en_naive, en_xc_sph, en_ex_cyl = 
                        en_energy(molecule1,molecule2,i,j);
                    
                    ii0 = atom_ind_base[ii] + ceil(Int,i/clouds_per_atom) - 1;

                    # Naive contribution
                    aux_Y[ii0] += en_naive;
                    
                    # XC contributions
                    if isinf(en_ex_cyl)
                        xc_coeff_1 = comb_coeffs_1b(xc_c_1b,Z1,ζ1);

                        aux_Y[ii0] -= xc_coeff_1*en_xc_sph;
                    else
                        xc_coeff_1 = comb_coeffs_2b(xc_c_2b,Z1,Z2,ζ1,ζ2);
                        xc_coeff_2 = comb_coeffs_2b(xc_d_2b,Z1,Z2,ζ1,ζ2);
                        
                        aux_Y[ii0] -= xc_coeff_1*en_xc_sph;
                        aux_Y[ii0] -= xc_coeff_2*en_ex_cyl;
                    end
                end
            end

            # cloud-cloud
            for i in 1:num_clouds1
                for j in 1:num_clouds2
                    Z1 = cloud_atomic_number(molecule1,i);
                    Z2 = cloud_atomic_number(molecule2,j);

                    ζ1 = cloud_polarization_coeff(molecule1,i);
                    ζ2 = cloud_polarization_coeff(molecule2,j);

                    ee_naive, ee_xc_sph, ee_xc_cyl = 
                        ee_energy(molecule1,molecule2,i,j);

                    ii0 = atom_ind_base[ii] + ceil(Int,i/clouds_per_atom) - 1;
                    jj0 = atom_ind_base[jj] + ceil(Int,j/clouds_per_atom) - 1;

                    # Naive contribution
                    aux_M[ii0,jj0] += ee_naive;

                    # XC contributions
                    if isinf(ee_xc_cyl)
                        xc_coeff_1 = comb_coeffs_1b(xc_a_1b,Z1,ζ1);
                        
                        aux_M[ii0,jj0] += xc_coeff_1*ee_xc_sph;
                    else
                        xc_coeff_1 = comb_coeffs_2b(xc_a_2b,Z1,Z2,ζ1,ζ2);
                        xc_coeff_2 = comb_coeffs_2b(xc_b_2b,Z1,Z2,ζ1,ζ2);

                        aux_M[ii0,jj0] += xc_coeff_1*ee_xc_sph;
                        aux_M[ii0,jj0] += xc_coeff_2*ee_xc_cyl;
                    end
                end
            end
        end
    end

    return aux_M, aux_Y;
end

function scf_min_func(simulation::SimulationSystem)
    molecules = simulation.system.molecules;
    clouds_per_atom = simulation.basis_set_settings.clouds_per_atom;

    tot_num_atoms = 0;
    tot_num_clouds = 0;
    atom_ind_base = [];
    for molecule in molecules
        if isempty(atom_ind_base)
            push!(atom_ind_base,1);
            tot_num_atoms += number_of_atoms(molecule);
            tot_num_clouds += number_of_clouds(molecule);
        else
            aux_num_atoms = number_of_atoms(molecule);
            aux_num_clouds = number_of_clouds(molecule);
            
            push!(atom_ind_base, tot_num_atoms + 1);
            tot_num_atoms += aux_num_atoms;
            tot_num_clouds += aux_num_clouds;
        end
    end

    # Make the initial guess of the polarization coefficients and the 
    # chemical potential.
    num_vars = tot_num_atoms + 1;
    x0 = ones(Float64,num_vars);
    x0[end] = 0.0;

    # x0 = zeros(Float64,0);
    # for molecule in molecules
    #     sqrt_ζ = sqrt.(molecule.cloud_data[1:clouds_per_atom:end,6]);
    #     x0 = vcat(x0, sqrt_ζ);
    # end
    # x0 = vcat(x0,simulation.system.chemical_potential);

    needs_casting = true;
    old_system = deepcopy(simulation.system);

    function min_func(aux_X::AbstractVector)
        # Error function that if the polarization coefficients and chemical 
        # potential are exact will always be equal to zero. This function 
        # always returns non-negative values.

        aux_type = typeof(aux_X[1]);
        if (needs_casting == true)
            needs_casting = false;

            for i in eachindex(simulation.system.molecules)
                simulation.system.molecules[i].cloud_data =
                    aux_type.(simulation.system.molecules[i].cloud_data);

                simulation.system.molecules[i].atoms_data =
                    aux_type.(simulation.system.molecules[i].atoms_data);
            end

            simulation.system.chemical_potential = 
                aux_type(simulation.system.chemical_potential);
        end
        
        # This ensures that the polarization coefficients are limited 
        # to non-negative numbers.
        aux_X[1:(end-1)] .^= 2.0;

        # Copy the trial polarization coefficients to the deepcopied 
        # vector of molecules.
        for i in eachindex(simulation.system.molecules)
            i0 = atom_ind_base[i];
            num_atoms = number_of_atoms(simulation.system.molecules[i]);

            ζ = aux_X[i0:(i0+num_atoms-1)];
            set_mol_pol_coeffs!(simulation.system.molecules[i],ζ);
        end
        simulation.system.chemical_potential = aux_X[end];

        # Get the system of equations associated with making the  
        # polarization Lagrangian stationary.
        aux_M, aux_Y = polarization_matrix_problem(simulation,
            simulation.cast_types.cast_to_pol_coeff_type);
        ret_vec = aux_M * aux_X - aux_Y;

        return dot(ret_vec,ret_vec);
    end

    cost_trace = zeros(Float64,0);

    # Callback function
    function trace_callback(state)
        push!(cost_trace, state.value);
        return false;
    end

    sol = Optim.optimize(min_func, x0, LBFGS(), autodiff=:forward,
        Optim.Options(show_trace=false, callback=trace_callback,
        iterations=8000));
    sol = Optim.minimizer(sol);
    sol[1:(end-1)] .^= 2.0;

    simulation.system = old_system;
    return cost_trace, sol;
end

function polarize_molecules!(simulation::SimulationSystem)
    # Calculates and sets the polarization coefficients of the atoms.
    _, minimizer = scf_min_func(simulation);

    tot_num_atoms = 0;
    tot_num_clouds = 0;
    atom_ind_base = [];
    for molecule in simulation.system.molecules
        if isempty(atom_ind_base)
            push!(atom_ind_base,1);
            tot_num_atoms += number_of_atoms(molecule);
            tot_num_clouds += number_of_clouds(molecule);
        else
            aux_num_atoms = number_of_atoms(molecule);
            aux_num_clouds = number_of_clouds(molecule);
            
            push!(atom_ind_base, tot_num_atoms + 1);
            tot_num_atoms += aux_num_atoms;
            tot_num_clouds += aux_num_clouds;
        end
    end

    clouds_per_atom = simulation.basis_set_settings.clouds_per_atom;
    simulation.system.chemical_potential = minimizer[end];

    for i in eachindex(simulation.system.molecules)
        i0 = atom_ind_base[i];
        num_atoms = number_of_atoms(simulation.system.molecules[i]);

        for ii in 1:num_atoms
            ζ = minimizer[i0+ii-1];
            simulation.system.molecules[i].atoms_data[ii,6] = ζ;

            for jj in 1:clouds_per_atom
                jj0 = (ii-1)*clouds_per_atom + jj;
                simulation.system.molecules[i].cloud_data[jj0,6] = ζ;
            end
        end
    end
end

function system_energies(simulation::SimulationSystem)
    # Returns the kinetic, naive and xc energy contributions of our 
    # simplified model.
    molecules = simulation.system.molecules;
    num_molecules = length(molecules);

    xc_a_1b = simulation.tot_e_xc_coeffs.xc_a_1b;
    xc_c_1b = simulation.tot_e_xc_coeffs.xc_c_1b;
    xc_e_1b = simulation.tot_e_xc_coeffs.xc_e_1b;
    xc_f_1b = simulation.tot_e_xc_coeffs.xc_f_1b;

    xc_a_2b = simulation.tot_e_xc_coeffs.xc_a_2b;
    xc_b_2b = simulation.tot_e_xc_coeffs.xc_b_2b;
    xc_c_2b = simulation.tot_e_xc_coeffs.xc_c_2b;
    xc_d_2b = simulation.tot_e_xc_coeffs.xc_d_2b;

    unpol_tf_ke = simulation.basis_set_settings.thomas_fermi_kes;
    unpol_vw_ke = simulation.basis_set_settings.von_weizsacker_kes;

    coeff_1b_map = simulation.coeff_mappings.coeff_1b_map;
    coeff_2b_map = simulation.coeff_mappings.coeff_2b_map;

    clouds_per_atom = simulation.basis_set_settings.clouds_per_atom;

    function comb_coeffs_1b(coeffs::AbstractMatrix, Z::Int, ζ::Number)
        # Calculates the XC coefficient based on the polarization coefficients.
        # (One body cases)
        b = coeffs[coeff_1b_map[Z],1];
        m = coeffs[coeff_1b_map[Z],2];

        return m*ζ + b;
    end

    function comb_coeffs_2b(coeffs::AbstractMatrix, Z1::Int, Z2::Int,
        ζ1::Number, ζ2::Number)
        # Calculates the XC coefficient based on the polarization coefficients.
        # (One body cases)
        b = coeffs[coeff_2b_map[(Z1,Z2)],1];
        m1 = coeffs[coeff_2b_map[(Z1,Z2)],2];
        m2 = coeffs[coeff_2b_map[(Z1,Z2)],3];

        if (Z1 < Z2)
            return m1*ζ1 + m2*ζ2 + b;
        else
            return m1*ζ2 + m2*ζ1 + b;
        end
    end
    
    aux_type = typeof(xc_a_1b[1]);
    xc_energy = aux_type(0.0);
    naive_energy = aux_type(0.0);
    kinetic_energy = aux_type(0.0);

    for ii in 1:num_molecules
        molecule1 = molecules[ii];
        num_atoms1 = number_of_atoms(molecule1);
        num_clouds1 = number_of_clouds(molecule1);

        # Kinetic Energy
        for i in 1:num_atoms1
            Z1 = atom_atomic_number(molecule1,i);
            ζ1 = atom_polarization_coeff(molecule1,i);

            tf_fit_coeff = 0.0;
            tf_fit_coeff += xc_e_1b[Z1,1]*(ζ1^(5.0/3.0));
            tf_fit_coeff += xc_e_1b[Z1,2]*(ζ1^(8.0/3.0));

            vw_fit_coeff = 0.0;
            vw_fit_coeff += xc_f_1b[Z1,1]*(ζ1^(1.0));
            vw_fit_coeff += xc_f_1b[Z1,2]*(ζ1^(2.0));

            kinetic_energy += tf_fit_coeff*unpol_tf_ke[Z1];
            kinetic_energy += vw_fit_coeff*unpol_vw_ke[Z1];
        end

        # Intramolecular interactions
        # nuclei-nuclei
        for i in 1:num_atoms1
            for j in (i+1):num_atoms1
                z1_eff = atom_eff_atomic_number(molecule1,i);
                z2_eff = atom_eff_atomic_number(molecule1,j);

                r1 = atom_position(molecule1,i);
                r2 = atom_position(molecule1,j);
                aux_dist = norm(r1-r2);

                naive_energy += z1_eff*z2_eff*(1.0/aux_dist);
            end
        end

        # nuclei-cloud 
        for i in 1:num_atoms1
            for j in 1:num_clouds1
                Z1 = atom_atomic_number(molecule1,i);
                Z2 = cloud_atomic_number(molecule1,j);

                ζ1 = atom_polarization_coeff(molecule1,i);
                ζ2 = cloud_polarization_coeff(molecule1,j);

                en_naive, en_xc_sph, en_xc_cyl = 
                    en_energy(molecule1,molecule1,j,i);

                # XC contributions
                if isinf(en_xc_cyl)
                    # polarization
                    en_naive *= ζ2;
                    en_xc_sph *= ζ2;
                    xc_coeff_1 = comb_coeffs_1b(xc_c_1b,Z1,ζ1);

                    xc_energy += xc_coeff_1*en_xc_sph;
                else
                    # polarization
                    en_naive *= ζ2;
                    en_xc_sph *= ζ2;
                    en_xc_cyl *= ζ2;
                    xc_coeff_1 = comb_coeffs_2b(xc_c_2b,Z1,Z2,ζ1,ζ2);
                    xc_coeff_2 = comb_coeffs_2b(xc_d_2b,Z1,Z2,ζ1,ζ2);
                    
                    xc_energy += xc_coeff_1*en_xc_sph;
                    xc_energy += xc_coeff_2*en_xc_cyl;
                end

                naive_energy -= en_naive;
            end
        end

        # cloud-cloud
        for i in 1:num_clouds1
            for j in i:num_clouds1
                Z1 = cloud_atomic_number(molecule1,i);
                Z2 = cloud_atomic_number(molecule1,j);

                ζ1 = cloud_polarization_coeff(molecule1,i);
                ζ2 = cloud_polarization_coeff(molecule1,j);

                ee_naive, ee_xc_sph, ee_xc_cyl = 
                    ee_energy(molecule1,molecule1,i,j);

                ii0 = ceil(Int,i/clouds_per_atom);
                jj0 = ceil(Int,j/clouds_per_atom);

                if (ii0 == jj0) && (i == j)
                    ee_naive *= 0.5;
                end                

                # XC contributions
                if isinf(ee_xc_cyl)
                    # polarization
                    ee_naive *= ζ1*ζ2;
                    ee_xc_sph *= ζ1*ζ2;
                    xc_coeff_1 = comb_coeffs_1b(xc_a_1b,Z1,ζ1);

                    xc_energy += xc_coeff_1*ee_xc_sph;
                else
                    # polarization
                    ee_naive *= ζ1*ζ2;
                    ee_xc_sph *= ζ1*ζ2;
                    ee_xc_cyl *= ζ1*ζ2;
                    xc_coeff_1 = comb_coeffs_2b(xc_a_2b,Z1,Z2,ζ1,ζ2);
                    xc_coeff_2 = comb_coeffs_2b(xc_b_2b,Z1,Z2,ζ1,ζ2);
                    
                    xc_energy += xc_coeff_1*ee_xc_sph;
                    xc_energy += xc_coeff_2*ee_xc_cyl;
                end

                naive_energy += ee_naive;
            end
        end

        # Intermolecular interactions
        for jj in (ii+1):num_molecules
            molecule2 = molecules[jj];
            num_atoms2 = number_of_atoms(molecule2);
            num_clouds2 = number_of_clouds(molecule2);

            # nuclei-nuclei
            for i in 1:num_atoms1
                for j in 1:num_atoms2
                    z1_eff = atom_eff_atomic_number(molecule1,i);
                    z2_eff = atom_eff_atomic_number(molecule2,j);

                    r1 = atom_position(molecule1,i);
                    r2 = atom_position(molecule2,j);
                    aux_dist = norm(r1-r2);

                    naive_energy += z1_eff*z2_eff*(1.0/aux_dist);
                end
            end

            # nuclei-cloud 
            for i in 1:num_atoms1
                for j in 1:num_clouds2
                    Z1 = atom_atomic_number(molecule1,i);
                    Z2 = cloud_atomic_number(molecule2,j);

                    ζ1 = atom_polarization_coeff(molecule1,i);
                    ζ2 = cloud_polarization_coeff(molecule2,j);

                    en_naive, en_xc_sph, en_xc_cyl = 
                        en_energy(molecule2,molecule1,j,i);

                    # polarization
                    en_naive *= ζ2;
                    en_xc_sph *= ζ2;
                    en_xc_cyl *= ζ2;

                    naive_energy -= en_naive;

                    # XC contributions
                    xc_coeff_1 = comb_coeffs_2b(xc_c_2b,Z1,Z2,ζ1,ζ2);
                    xc_coeff_2 = comb_coeffs_2b(xc_d_2b,Z1,Z2,ζ1,ζ2);
                    
                    xc_energy += xc_coeff_1*en_xc_sph;
                    xc_energy += xc_coeff_2*en_xc_cyl;
                end
            end

            for i in 1:num_clouds1
                for j in 1:num_atoms2
                    Z1 = cloud_atomic_number(molecule1,i);
                    Z2 = atom_atomic_number(molecule2,j);

                    ζ1 = cloud_polarization_coeff(molecule1,i);
                    ζ2 = atom_polarization_coeff(molecule2,j);

                    en_naive, en_xc_sph, en_xc_cyl = 
                        en_energy(molecule1,molecule2,i,j);

                    # polarization
                    en_naive *= ζ1;
                    en_xc_sph *= ζ1;
                    en_xc_cyl *= ζ1;

                    naive_energy -= en_naive;

                    # XC contributions
                    xc_coeff_1 = comb_coeffs_2b(xc_c_2b,Z1,Z2,ζ1,ζ2);
                    xc_coeff_2 = comb_coeffs_2b(xc_d_2b,Z1,Z2,ζ1,ζ2);
                    
                    xc_energy += xc_coeff_1*en_xc_sph;
                    xc_energy += xc_coeff_2*en_xc_cyl;
                end
            end

            # cloud-cloud
            for i in 1:num_clouds1
                for j in 1:num_clouds2
                    Z1 = cloud_atomic_number(molecule1,i);
                    Z2 = cloud_atomic_number(molecule2,j);

                    ζ1 = cloud_polarization_coeff(molecule1,i);
                    ζ2 = cloud_polarization_coeff(molecule2,j);

                    ee_naive, ee_xc_sph, ee_xc_cyl = 
                        ee_energy(molecule1,molecule2,i,j);

                    # polarization
                    ee_naive *= ζ1*ζ2;
                    ee_xc_sph *= ζ1*ζ2;
                    ee_xc_cyl *= ζ1*ζ2;

                    naive_energy += ee_naive;

                    # XC contributions
                    xc_coeff_1 = comb_coeffs_2b(xc_a_2b,Z1,Z2,ζ1,ζ2);
                    xc_coeff_2 = comb_coeffs_2b(xc_b_2b,Z1,Z2,ζ1,ζ2);

                    xc_energy += xc_coeff_1*ee_xc_sph;
                    xc_energy += xc_coeff_2*ee_xc_cyl;
                end
            end
        end
    end

    return xc_energy, naive_energy, kinetic_energy;
end

function total_energy(simulation::SimulationSystem)
    xc_energy, naive_energy, kinetic_energy = system_energies(simulation);
    return xc_energy + naive_energy + kinetic_energy;
end
