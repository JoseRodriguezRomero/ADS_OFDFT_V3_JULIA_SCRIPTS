using Optim
using Integrals
using ForwardDiff
using SpecialFunctions
using LinearAlgebra, GenericLinearAlgebra

include("ForceFieldBase.jl")
include("ForceFieldBaseBasisSet.jl")
include("ForceFieldBaseFitCoeffs.jl")

function naive_coulomb_integral(λ::Real, d::Real)
    if (abs(d) < atoms_dist_cutoff())
        return 2*sqrt(λ/π)
    end

    if isinf(λ)
        return 1.0/d
    end

    erf_arg = d*sqrt(λ)
    if erf_arg > erf_arg_cutoff()
        return 1.0/d
    end

    return erf(erf_arg)/d
end 

function xc_sph(λ::Real, d::Real)
    ret_val = 2*sqrt(λ/π)*(exp(-2*λ*d-λ)-1)/exp(λ*d^2)

    if isnan(ret_val)
        return 0.0
    end

    return ret_val
end

function xc_cyl(λ::Real, d::Real)
    ret_val = (4*exp(-(d^2*λ)-2*d*λ-λ)*(exp(2*d*λ+λ)*d-d-1)*λ^(3/2))/sqrt(π)

    if isnan(ret_val)
        return 0.0
    end

    return ret_val
end

function morse_u(depth::Real, stiffness_parameter::Real, 
    equilibrium_distance::Real, atomic_separation::Real)
    A = depth
    B = stiffness_parameter
    C = equilibrium_distance

    r = atomic_separation

    return A * ((1 - exp(-B * (r - C)))^2 - 1)
end

function smear_comb(λ1::Number, λ2::Number)
    if isinf(λ1) && isinf(λ2)
        return ∞
    end

    if isinf(λ1)
        return λ2
    end

    if isinf(λ2)
        return λ1
    end

    return (λ1 * λ2) / (λ1 + λ2)
end

"""
    unpol_ee_energy(atom_1::Atom, atom_2::Atom) → (e_naive::Real, e_xc_sph::Real, e_xc_cyl::Real)

Calculates the unpolarized electron-electron interactions between the valence 
electron shells of atoms atom_1 and atom_2.
"""
function unpol_ee_energy(atom_1::Atom, atom_2::Atom)
    atom_1_cloud = atom_1.electron_cloud_shells[end]
    atom_1_basis_size = length(atom_1_cloud.basis_function_amplitude)

    atom_2_cloud = atom_2.electron_cloud_shells[end]
    atom_2_basis_size = length(atom_2_cloud.basis_function_amplitude)

    d = norm(atom_1.coordinates - atom_2.coordinates)

    e_naive = 0.0
    e_xc_sph = 0.0
    e_xc_cyl = 0.0

    for i in 1:atom_1_basis_size
        for j in 1:atom_2_basis_size
            c1 = atom_1_cloud.basis_function_amplitude[i]
            λ1 = atom_1_cloud.basis_function_decay[i]

            c2 = atom_2_cloud.basis_function_amplitude[j]
            λ2 = atom_2_cloud.basis_function_decay[j]

            c = c1 * c2
            λ = smear_comb(λ1,λ2)

            e_xc_sph += c*xc_sph(λ,d)
            e_xc_cyl += c*xc_cyl(λ,d)
            e_naive += c*naive_coulomb_integral(λ,d)
        end
    end

    return e_naive, e_xc_sph, e_xc_cyl
end

"""
    ee_energy(atom_1::Atom, atom_2::Atom) → (e_naive::Real, e_xc_sph::Real, e_xc_cyl::Real)

Calculates the polarized electron-electron interactions between the electron 
shells of atoms atom_1 and atom_2.
"""
function ee_energy(atom_1::Atom, atom_2::Atom)
    ζ1 = atom_1.polarization_coefficient
    ζ2 = atom_2.polarization_coefficient

    e_naive, e_xc_sph, e_xc_cyl = unpol_ee_energy(atom_1,atom_2)

    e_naive = ζ1*ζ2*e_naive
    e_xc_sph = ζ1*ζ2*e_xc_sph
    e_xc_cyl = ζ1*ζ2*e_xc_cyl

    return e_naive, e_xc_sph, e_xc_cyl
end

"""
    unpol_en_energy(atom_1::Atom, atom_2::Atom) → (e_naive::Real, e_xc_sph::Real, e_xc_cyl::Real)

Calculates the unpolarized interactions between the valence electron shell of 
atom atom_1 and the nuclei of atom atom_2. The charge of the nuclei of atom_2 is 
replaced with its number of valence electrons.
"""
function unpol_en_energy(atom_1::Atom, atom_2::Atom)
    atom_1_cloud = atom_1.electron_cloud_shells[end]
    atom_1_basis_size = length(atom_1_cloud.basis_function_amplitude)

    d = norm(atom_1.coordinates - atom_2.coordinates)

    e_naive = 0
    e_xc_sph = 0
    e_xc_cyl = 0

    for i in 1:atom_1_basis_size
        c1 = atom_1_cloud.basis_function_amplitude[i]
        λ1 = atom_1_cloud.basis_function_decay[i]

        q2 = atom_2.valence_electrons

        c = c1 * q2
        λ = λ1

        e_xc_sph += c*xc_sph(λ,d)
        e_xc_cyl += c*xc_cyl(λ,d)
        e_naive -= c*naive_coulomb_integral(λ,d)
    end

    return e_naive, e_xc_sph, e_xc_cyl
end

"""
    unpol_ne_energy(atom_1::Atom, atom_2::Atom) → (e_naive::Real, e_xc_sph::Real, e_xc_cyl::Real)

Calculates the unpolarized interactions between the valence electron shell of 
atom atom_2 and the nuclei of atom atom_1. The charge of the nuclei of atom_1 is 
replaced with its number of valence electrons.
"""
function unpol_ne_energy(atom_1::Atom, atom_2::Atom)
    e_naive, e_xc_sph, e_xc_cyl = unpol_en_energy(atom_2,atom_1)
    return e_naive, e_xc_sph, e_xc_cyl
end

"""
    en_energy(atom_1::Atom, atom_2::Atom) → (e_naive::Real, e_xc_sph::Real, e_xc_cyl::Real)

Calculates the polarized interactions between the electron shells of atom atom_1 
and the nuclei of atom atom_2.
"""
function en_energy(atom_1::Atom, atom_2::Atom)
    ζ1 = atom_1.polarization_coefficient

    e_naive, e_xc_sph, e_xc_cyl = unpol_en_energy(atom_1,atom_2)

    e_naive = ζ1*e_naive
    e_xc_sph = ζ1*e_xc_sph
    e_xc_cyl = ζ1*e_xc_cyl

    return e_naive, e_xc_sph, e_xc_cyl
end

"""
    ne_energy(atom_1::Atom, atom_2::Atom) → (e_naive::Real, e_xc_sph::Real, e_xc_cyl::Real)

Calculates the polarized interactions between the nuclei of atom atom_1 and the 
electron shells of atom atom_2.
"""
function ne_energy(atom_1::Atom, atom_2::Atom)
    e_naive, e_xc_sph, e_xc_cyl = en_energy(atom_2,atom_1)
    return e_naive, e_xc_sph, e_xc_cyl
end

"""
    nn_energy(atom_1::Atom, atom_2::Atom) → e_naive::Real

Calculates the Couloumb interaction between the nuclei of atoms atom_1 and 
atom_2. The atomic numbers are replaced with the number of valence electrons.
"""
function nn_energy(atom_1::Atom, atom_2::Atom)
    q1 = atom_1.valence_electrons
    q2 = atom_2.valence_electrons

    d = norm(atom_1.coordinates - atom_2.coordinates)
    e_naive = (q1*q2)/d

    return e_naive
end

"""
    rotate_molecule!(molecule::Molecule, Δθ::Matrix{Number})

Rotates the whole molecule by the rotation matrix Δθ.
"""
function rotate_molecule!(molecule::Molecule, Δθ::Matrix{Number})
    for i in eachindex(molecule.atoms)
        molecule.atoms[i].position = Δθ*(molecule.atoms[i].position)
    end

    return
end

"""
    polarization_matrix_problem(simulation::SimulationSystem) → (aux_M::Matrix{Real}, aux_Y::Matrix{Real})

Calculates the matrix problem that needs to be solved for the polarization 
Lagrangian to be stationary.
"""
function polarization_matrix_problem(simulation::SimulationSystem)
    molecules = simulation.system.molecules
    charge = simulation.system.charge

    coeffs = simulation.pol_e_coeffs

    xc_a_1b = coeffs.xc_coeffs.xc_a_1b
    xc_b_1b = coeffs.xc_coeffs.xc_b_1b
    xc_c_1b = coeffs.xc_coeffs.xc_c_1b
    xc_d_1b = coeffs.xc_coeffs.xc_d_1b
    ke_e_1b = coeffs.ke_coeffs.ke_e_1b
    ke_f_1b = coeffs.ke_coeffs.ke_f_1b
    
    xc_a_2b = coeffs.xc_coeffs.xc_a_2b
    xc_b_2b = coeffs.xc_coeffs.xc_b_2b
    xc_c_2b = coeffs.xc_coeffs.xc_c_2b
    xc_d_2b = coeffs.xc_coeffs.xc_d_2b

    aux_type = typeof(xc_a_1b[1])
    atoms_μ0 = simulation.basis_set_settings.atoms_μ0

    tot_num_atoms = 0
    atom_ind_base = Vector{Int}()
    for molecule in molecules
        if isempty(atom_ind_base)
            push!(atom_ind_base,1)
            tot_num_atoms += number_of_atoms(molecule)
        else            
            push!(atom_ind_base, tot_num_atoms + 1)
            tot_num_atoms += number_of_atoms(molecule)
        end
    end

    num_vars = 2*tot_num_atoms + 1
    aux_M = zeros(aux_type, num_vars, num_vars)
    aux_Y = zeros(aux_type, num_vars, 1)
    aux_Y[end] -= charge

    function compute_chemical_potential(atom::Atom, atom_index::Int)
        aux_M[atom_index,end] -= atom.valence_electrons
        aux_M[end,atom_index] += atom.valence_electrons
        aux_Y[end] += atom.valence_electrons

        return
    end

    function compute_chemical_potential_shift(atom::Atom, atom_index::Int)
        atom_shift_index = tot_num_atoms + atom_index
        μ0 = atoms_μ0[atom.atomic_number]

        aux_M[atom_shift_index,end] += atom.valence_electrons
        aux_M[atom_shift_index,atom_shift_index] -= atom.valence_electrons
        aux_M[atom_index,atom_shift_index] -= atom.valence_electrons
        aux_Y[atom_shift_index] = -μ0

        return
    end

    function compute_kinetic_energies(atom::Atom, atom_index::Int)
        z = atom.atomic_number
        
        valence_shell = atom.electron_cloud_shells[end]
        unpol_tf_ke = valence_shell.thomas_fermi_ke
        unpol_vw_ke = valence_shell.von_weizsacker_ke

        tf_fit_coeff_1 = aux_type(0.0)
        tf_fit_coeff_2 = aux_type(0.0)
        tf_fit_coeff_1 += (10.0/9.0)*ke_e_1b[z]*unpol_tf_ke
        tf_fit_coeff_2 += (5.0/9.0)*ke_e_1b[z]*unpol_tf_ke

        aux_M[atom_index,atom_index] += tf_fit_coeff_1
        aux_Y[atom_index] -= tf_fit_coeff_2

        vw_fit_coeff = aux_type(0.0)
        vw_fit_coeff += ke_f_1b[z]*unpol_vw_ke
        aux_Y[atom_index] -= vw_fit_coeff

        return
    end

    function compute_en_energies(atom_1::Atom, atom_2::Atom,
        atom_index_1::Int, atom_index_2::Int)
        z1 = atom_1.atomic_number
        z2 = atom_2.atomic_number
        
        e_naive, e_xc_sph, e_xc_cyl = unpol_en_energy(atom_1,atom_2)

        energy_sum_12 = e_naive
        energy_sum_12 += e_xc_sph*xc_c_2b[(z1,z2)]
        energy_sum_12 += e_xc_cyl*xc_d_2b[(z1,z2)]

        e_naive, e_xc_sph, e_xc_cyl = unpol_en_energy(atom_2,atom_1)
            
        energy_sum_21 = e_naive
        energy_sum_21 += e_xc_sph*xc_c_2b[(z2,z1)]
        energy_sum_21 += e_xc_cyl*xc_d_2b[(z2,z1)]

        aux_Y[atom_index_1] -= energy_sum_12
        aux_Y[atom_index_2] -= energy_sum_21

        return
    end

    function compute_en_energies(atom::Atom, atom_index::Int)
        z = atom.atomic_number
        e_naive, e_xc_sph, e_xc_cyl = unpol_en_energy(atom,atom)

        energy_sum = e_naive
        energy_sum += e_xc_sph*xc_c_1b[z]
        energy_sum += e_xc_cyl*xc_d_1b[z]

        aux_Y[atom_index] -= energy_sum

        return
    end

    function compute_ee_energies(atom_1::Atom, atom_2::Atom,
        atom_index_1::Int, atom_index_2::Int)
        z1 = atom_1.atomic_number
        z2 = atom_2.atomic_number

        e_naive, e_xc_sph, e_xc_cyl = unpol_ee_energy(atom_1,atom_2)

        energy_sum = e_naive
        energy_sum += e_xc_sph*xc_a_2b[(z1,z2)]
        energy_sum += e_xc_cyl*xc_b_2b[(z1,z2)]

        aux_M[atom_index_1,atom_index_2] += energy_sum
        aux_M[atom_index_2,atom_index_1] += energy_sum

        return
    end

    function compute_ee_energies(atom::Atom, atom_index::Int)
        z = atom.atomic_number
        e_naive, e_xc_sph, e_xc_cyl = unpol_ee_energy(atom,atom)

        energy_sum = e_naive
        energy_sum += e_xc_sph*xc_a_1b[z]
        energy_sum += e_xc_cyl*xc_b_1b[z]

        aux_M[atom_index,atom_index] += energy_sum
        
        return
    end

    function compute_one_body_terms()
        for molecule_index in eachindex(molecules)
            molecule = molecules[molecule_index]

            for i in eachindex(molecule.atoms)
                atom = molecule.atoms[i]
                atom_index = atom_ind_base[molecule_index] + i - 1

                compute_en_energies(atom,atom_index)
                compute_ee_energies(atom,atom_index)
                compute_kinetic_energies(atom,atom_index)
                compute_chemical_potential(atom,atom_index)
                compute_chemical_potential_shift(atom,atom_index)
            end
        end

        return
    end

    function compute_two_body_intramolecular_terms()
        for molecule_index in eachindex(molecules)
            molecule = molecules[molecule_index]

            for i in eachindex(molecule.atoms)
                for j in (i+1):length(molecule.atoms)
                    atom_1 = molecule.atoms[i]
                    atom_2 = molecule.atoms[j]

                    atom_index_1 = atom_ind_base[molecule_index] + i - 1
                    atom_index_2 = atom_ind_base[molecule_index] + j - 1

                    compute_en_energies(atom_1,atom_2,atom_index_1,atom_index_2)
                    compute_ee_energies(atom_1,atom_2,atom_index_1,atom_index_2)
                end
            end
        end

        return
    end

    function compute_two_body_intermolecular_terms()
        for molecule_index_1 in eachindex(molecules)
            for molecule_index_2 in (molecule_index_1+1):length(molecules)
                molecule_1 = molecules[molecule_index_1]
                molecule_2 = molecules[molecule_index_2]

                for i in eachindex(molecule_1)
                    for j in eachindex(molecule_2)
                        atom_1 = molecule_1.atoms[i]
                        atom_2 = molecule_2.atoms[j]

                        atom_index_1 = atom_ind_base[molecule_index_1] + i - 1
                        atom_index_2 = atom_ind_base[molecule_index_2] + j - 1

                        compute_en_energies(atom_1,atom_2,atom_index_1,atom_index_2)
                        compute_ee_energies(atom_1,atom_2,atom_index_1,atom_index_2)
                    end
                end
            end
        end

        return
    end

    function compute_two_body_terms()
        compute_two_body_intramolecular_terms()
        compute_two_body_intermolecular_terms()

        return
    end
    
    compute_one_body_terms()
    compute_two_body_terms()

    return aux_M, aux_Y
end

"""
    polarize_molecules!(simulation::SimulationSystem)

Calculates and sets the polarization coefficients of the atoms.
"""
function polarize_molecules!(simulation::SimulationSystem)
    aux_m, aux_y = polarization_matrix_problem(simulation)
    minimizer = aux_m \ aux_y

    molecules = simulation.system.molecules

    tot_num_atoms = 0
    atom_ind_base = []
    for molecule in simulation.system.molecules
        if isempty(atom_ind_base)
            push!(atom_ind_base,1)
            tot_num_atoms += number_of_atoms(molecule)
        else            
            push!(atom_ind_base, tot_num_atoms + 1)
            tot_num_atoms += number_of_atoms(molecule)
        end
    end

    # Set the polarization coefficients
    for molecule_index in eachindex(molecules)
        molecule = molecules[molecule_index]

        for i in eachindex(molecule.atoms)
            atom_index = atom_ind_base[molecule_index] + i - 1

            ζ = minimizer[atom_index]
            molecules[molecule_index].atoms[i].polarization_coefficient = ζ
        end
    end

    # Set the chemical potential
    simulation.system.chemical_potential = minimizer[end]
    return
end

"""
    system_energies(simulation::SimulationSystem) → (naive_energy::Real, kinetic_energy::Real, xc_energy::Real, non_polarizable_energy::Real)

Returns the kinetic, naive and xc energy contributions of our simplified model.
"""
function system_energies(simulation::SimulationSystem)
    molecules = simulation.system.molecules

    coeffs = simulation.tot_e_coeffs

    xc_a_1b = coeffs.xc_coeffs.xc_a_1b
    xc_b_1b = coeffs.xc_coeffs.xc_b_1b
    xc_c_1b = coeffs.xc_coeffs.xc_c_1b
    xc_d_1b = coeffs.xc_coeffs.xc_d_1b
    ke_e_1b = coeffs.ke_coeffs.ke_e_1b
    ke_f_1b = coeffs.ke_coeffs.ke_f_1b

    xc_a_2b = coeffs.xc_coeffs.xc_a_2b
    xc_b_2b = coeffs.xc_coeffs.xc_b_2b
    xc_c_2b = coeffs.xc_coeffs.xc_c_2b
    xc_d_2b = coeffs.xc_coeffs.xc_d_2b

    morse_depth = coeffs.non_polarizable_coeffs.depth
    morse_stiffness_parameter = coeffs.non_polarizable_coeffs.stiffness_parameter
    morse_equilibrium_distance = coeffs.non_polarizable_coeffs.equilibrium_distance
    
    aux_type = typeof(xc_a_1b[1])
    xc_energy = aux_type(0.0)
    naive_energy = aux_type(0.0)
    kinetic_energy = aux_type(0.0)
    non_polarizable_energy = aux_type(0.0)

    function compute_nn_energy(atom_1::Atom, atom_2::Atom)
        naive_energy += nn_energy(atom_1,atom_2)

        return
    end

    function compute_en_energies(atom_1::Atom, atom_2::Atom)
        z1 = atom_1.atomic_number
        z2 = atom_2.atomic_number

        en_naive, en_xc_sph, en_xc_cyl = en_energy(atom_1,atom_2)
        naive_energy += en_naive
        xc_energy += en_xc_sph*xc_c_2b[(z1,z2)]
        xc_energy += en_xc_cyl*xc_d_2b[(z1,z2)]

        en_naive, en_xc_sph, en_xc_cyl = en_energy(atom_2,atom_1)
        naive_energy += en_naive
        xc_energy += en_xc_sph*xc_c_2b[(z2,z1)]
        xc_energy += en_xc_cyl*xc_d_2b[(z2,z1)]

        return
    end

    function compute_en_energies(atom_1::Atom)
        z1 = atom_1.atomic_number

        en_naive, en_xc_sph, en_xc_cyl = en_energy(atom_1,atom_1)

        naive_energy += en_naive
        xc_energy += en_xc_sph*xc_c_1b[z1]
        xc_energy += en_xc_cyl*xc_d_1b[z1]
    end

    function compute_ee_energies(atom_1::Atom, atom_2::Atom)
        z1 = atom_1.atomic_number
        z2 = atom_2.atomic_number
        
        ee_naive, ee_xc_sph, ee_xc_cyl = ee_energy(atom_1,atom_2)
        
        naive_energy += ee_naive
        xc_energy += ee_xc_sph*xc_a_2b[(z1,z2)]
        xc_energy += ee_xc_cyl*xc_b_2b[(z1,z2)]

        return
    end

    function compute_ee_energies(atom_1::Atom)
        z1 = atom_1.atomic_number
        
        ee_naive, ee_xc_sph, ee_xc_cyl = ee_energy(atom_1,atom_1)

        naive_energy += ee_naive / 2.0
        xc_energy += ee_xc_sph*xc_a_1b[z1]
        xc_energy += ee_xc_cyl*xc_b_1b[z1]

        return
    end

    function compute_kinetic_energies(atom_1::Atom)
        z1 = atom_1.atomic_number
        ζ1 = atom_1.polarization_coefficient
        e_shell = atom_1.electron_cloud_shells[end]

        unpol_tf_ke = e_shell.thomas_fermi_ke
        unpol_vw_ke = e_shell.von_weizsacker_ke

        tf_fit_coeff = aux_type(0.0)

        tf_fit_coeff += (2.0/3.0)*ζ1*ζ1*ke_e_1b[z1]
        tf_fit_coeff += (1.0/3.0)*ζ1*ke_e_1b[z1]

        vw_fit_coeff = aux_type(0.0)
        vw_fit_coeff += ζ1*ke_f_1b[z1]

        kinetic_energy += tf_fit_coeff*unpol_tf_ke
        kinetic_energy += vw_fit_coeff*unpol_vw_ke

        return
    end

    function compute_morse_potential(atom_1::Atom, atom_2::Atom)
        z1 = atom_1.atomic_number
        z2 = atom_2.atomic_number

        d = norm(atom_1.coordinates - atom_2.coordinates)

        A = morse_depth[(z1,z2)]
        B = morse_stiffness_parameter[(z1,z2)]
        C = morse_equilibrium_distance[(z1,z2)]
        non_polarizable_energy += morse_u(A,B,C,d)

        return
    end

    function compute_one_body_energies()
        for molecule in molecules
            for atom in molecule.atoms
                compute_kinetic_energies(atom)
                compute_en_energies(atom)
                compute_ee_energies(atom)
            end
        end

        return
    end

    function compute_two_body_intramolecular_energies()
        for molecule in molecules
            for atom_index_1 in eachindex(molecule.atoms)
                for atom_index_2 in (atom_index_1+1):length(molecule.atoms)
                    atom_1 = molecule.atoms[atom_index_1]
                    atom_2 = molecule.atoms[atom_index_2]

                    compute_nn_energy(atom_1,atom_2)
                    compute_en_energies(atom_1,atom_2)
                    compute_ee_energies(atom_1,atom_2)
                    compute_morse_potential(atom_1,atom_2)
                end
            end
        end

        return
    end

    function compute_two_body_intermolecular_energies()
        for molecule_index_1 in eachindex(molecules)
            for molecule_index_2 in (molecule_index_1+1):length(molecules)
                molecule_1 = molecules[molecule_index_1]
                molecule_2 = molecules[molecule_index_2]

                for atom_1 in molecule_1.atoms
                    for atom_2 in molecule_2.atoms
                        compute_nn_energy(atom_1,atom_2)
                        compute_en_energies(atom_1,atom_2)
                        compute_ee_energies(atom_1,atom_2)
                        compute_morse_potentials(atom_1,atom_2)
                    end
                end
            end
        end

        return
    end

    function compute_two_body_energies()
        compute_two_body_intramolecular_energies()
        compute_two_body_intermolecular_energies()

        return
    end

    compute_one_body_energies()
    compute_two_body_energies()

    return naive_energy, kinetic_energy, xc_energy, non_polarizable_energy
end

"""
    total_energy(simulation::SimulationSystem) → total_energy::Real

Returns the sum of all three energy contributions.
"""
function total_energy(simulation::SimulationSystem)
    naive_energy, kinetic_energy, xc_energy, non_polarizable_energy = system_energies(simulation)
    return naive_energy + kinetic_energy + xc_energy + non_polarizable_energy
end

function save_fitted_coeffs(simulation::SimulationSystem)
    save_1b_pol_e_coeffs(simulation)
    save_2b_pol_e_coeffs(simulation)

    save_1b_tot_e_coeffs(simulation)
    save_2b_tot_e_coeffs(simulation)

    return
end

function initialize_simulation_environment()
    basis_set_settings = load_basis_set()
    
    xc_a_1b = Dict{Int,Float64}()
    xc_b_1b = Dict{Int,Float64}()
    xc_c_1b = Dict{Int,Float64}()
    xc_d_1b = Dict{Int,Float64}()
    ke_e_1b = Dict{Int,Float64}()
    ke_f_1b = Dict{Int,Float64}()
    
    xc_a_2b = Dict{Tuple{Int,Int},Number}()
    xc_b_2b = Dict{Tuple{Int,Int},Number}()
    xc_c_2b = Dict{Tuple{Int,Int},Number}()
    xc_d_2b = Dict{Tuple{Int,Int},Number}()

    morse_depth = Dict{Tuple{Int,Int},Number}()
    morse_stiffness_parameter = Dict{Tuple{Int,Int},Number}()
    morse_equilibrium_distance = Dict{Tuple{Int,Int},Number}()

    xc_coeffs = EmpiricalXCCoefficients(xc_a_1b,xc_b_1b,xc_c_1b,xc_d_1b,xc_a_2b,xc_b_2b,xc_c_2b,xc_d_2b)
    ke_coeffs = EmpiricalKECoefficients(ke_e_1b,ke_f_1b)
    non_polarizable_coeffs = EmpiricalMorseCoefficients(morse_depth,morse_stiffness_parameter,morse_equilibrium_distance)

    tot_e_coeffs = deepcopy(TotalEnergyCoefficients(xc_coeffs,ke_coeffs,non_polarizable_coeffs))
    pol_e_coeffs = deepcopy(PolarizationEnergyCoefficients(xc_coeffs,ke_coeffs))

    system = MolecularSystem()
    simulation = SimulationSystem(system,tot_e_coeffs,pol_e_coeffs,basis_set_settings)

    load_fitted_coeffs!(simulation)

    return simulation
end

"""
    make_monoatomic_system(Z::Int, charge::Number) → simulation::SimulationSystem
    make_monoatomic_system(element::String,charge::Int) → simulation::SimulationSystem

Makes a molecule object with a single atom whose atomic number is Z.
"""
function make_monoatomic_system(Z::Int, charge::Number)
    simulation = initialize_simulation_environment()
    simulation.system.charge = charge
    reference_atoms = simulation.basis_set_settings.reference_atoms

    molecule = Molecule()
    molecule.atoms = [copy(reference_atoms[Z])]

    z_eff = molecule.atoms[1].valence_electrons
    molecule.atoms[1].polarization_coefficient = (z_eff - charge) / z_eff

    simulation.system.molecules = [molecule]
    return simulation
end

function make_monoatomic_system(element::String,charge::Int)
    Z = get_atomic_number(element)
    return make_monoatomic_system(Z,charge)
end

"""
    make_diatomic_molecule(reference_atoms::Dict{Int,Atom}, Z1::Int, Z2::Int, d::Number, charge::Int) → molecule::Molecule

Makes a diatomic molecule with a diatomic separation d in Bohr with the 
specified charge.
"""
function make_diatomic_molecule(reference_atoms::Dict{Int,Atom}, Z1::Int, Z2::Int, d::Number, charge::Int)
    atom_1 = copy(reference_atoms[Z1])
    atom_2 = copy(reference_atoms[Z2])

    z1_eff = atom_1.valence_electrons
    z2_eff = atom_2.valence_electrons

    charge1 = charge*z1_eff/(z1_eff+z2_eff)
    charge2 = charge*z2_eff/(z1_eff+z2_eff)

    atom_1.polarization_coefficient = (z1_eff - charge1) / z1_eff
    atom_2.polarization_coefficient = (z2_eff - charge2) / z2_eff

    atom_1.coordinates = [0.0, 0.0, +d/2.0]
    atom_2.coordinates = [0.0, 0.0, -d/2.0]

    molecule = Molecule()
    molecule.atoms = [atom_1,atom_2]

    return molecule
end

"""
    make_triatomic_molecule(reference_atoms::Dict{Int,Atom}, Z1::Int, Z2::Int, Z3::Int, distance_12::Number, distance_13::Number, angle_123::Number) → molecule::Molecule

Makes a triatomic molecule with a diatomic separation distance_12 in Bohr 
between atoms 1 and 2, distance_13 in Bohr between atoms 2 and 3, all of them 
describing an angle angle_213 in degrees.
"""
function make_triatomic_molecule(reference_atoms::Dict{Int,Atom}, Z1::Int, Z2::Int, Z3::Int, distance_12::Number, distance_13::Number, angle_123::Number)
    atom_1 = copy(reference_atoms[Z1])
    atom_2 = copy(reference_atoms[Z2])
    atom_3 = copy(reference_atoms[Z3])

    x1 = 0.0
    y1 = 0.0
    z1 = 0.0

    x2 = distance_12
    y2 = 0.0
    z2 = 0.0

    x3 = distance_13 * cos(angle_123 * (π / 180.0))
    y3 = distance_13 * sin(angle_123 * (π / 180.0))
    z3 = 0.0

    atom_1.coordinates = [x1, y1, z1]
    atom_2.coordinates = [x2, y2, z2]
    atom_3.coordinates = [x3, y3, z3]

    molecule = Molecule()
    molecule.atoms = [atom_1,atom_2,atom_3]

    return molecule
end

"""
    make_diatomic_system(Z1::Int, Z2::Int, d::Number, charge::Int) → simulation::SimulationSystem
    make_diatomic_system(element1::String, element2::String,d::Number,charge::Int) → simulation::SimulationSystem
    make_diatomic_system(Z::Int,d::Number,charge::Int) → simulation::SimulationSystem
    make_diatomic_system(element::String,d::Number,charge::Int) → simulation::SimulationSystem

Makes a simulation structure with a diatomic molecule with a diatomic separation 
d in Bohr with the specified charge.
"""
function make_diatomic_system(Z1::Int, Z2::Int, d::Number, charge::Int)
    simulation = initialize_simulation_environment()
    simulation.system.charge = charge
    reference_atoms = simulation.basis_set_settings.reference_atoms

    molecule = make_diatomic_molecule(reference_atoms,Z1,Z2,d,charge)

    simulation.system.molecules = [molecule]
    return simulation
end

function make_diatomic_system(
    element1::String, element2::String,d::Number,charge::Int)
    Z1 = get_atomic_number(element1)
    Z2 = get_atomic_number(element2)
    return make_diatomic_system(Z1,Z2,d,charge)
end

"""
    make_diatomic_system(Z::Int,d::Number,charge::Int) → simulation::SimulationSystem

Makes a homonuclear diatomic molecule with a diatomic separation d in Bohr with 
the specified charge.
"""
function make_diatomic_system(Z::Int,d::Number,charge::Int)
    return make_diatomic_system(Z,Z,d,charge)
end

"""
    make_diatomic_system(element::String,d::Number,charge::Int) → simulation::SimulationSystem

Makes a homonuclear diatomic molecule with a diatomic separation d in Bohr with 
the specified charge.
"""
function make_diatomic_system(element::String,d::Number,charge::Int)
    return make_diatomic_system(element,element,d,charge)
end

"""
    full_model_reset()

Sets all the empirical coefficients to zero and recalculates the Thomas-Fermi 
and von Weizacker kinetic energies for the atoms allowed in the basis-set.
"""
function full_model_reset()
    reset_fitted_coeffs()
    load_basis_set(true)

    return
end
