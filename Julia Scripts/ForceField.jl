using Integrals;
using ForwardDiff;
using SpecialFunctions;
using LinearAlgebra, GenericLinearAlgebra;

include("ForceFieldBase.jl")
include("ForceFieldBaseXCCoeffs.jl");
include("ForceFieldBaseBasisSet.jl");

function naive_coulomb_integral(λ::Real, d::Real)
    if (abs(d) < atoms_dist_cutoff())
        return 2*sqrt(λ/π);
    end

    erf_arg = d*sqrt(λ);
    if erf_arg > erf_arg_cutoff()
        return 1.0/d;
    end

    return erf(erf_arg)/d;
end 

function xc_sph(λ::Real, d::Real)
    if (abs(d) < atoms_dist_cutoff())
        return 0.0;
    end

    ret_val = sqrt(λ/π)*(1-exp(4*λ*d))/exp(λ*(d+1)^2);

    if isnan(ret_val)
        return 0.0;
    end

    return ret_val;
end

function xc_cyl(λ::Real, d::Real)
    if (abs(d) < atoms_dist_cutoff())
        return (2*exp(-λ)*sqrt(λ))/sqrt(π)-(2*sqrt(λ))/sqrt(π);
    end

    ret_val = sqrt(λ/π)*((exp(2*λ*d)-exp(λ))^2-exp(2*λ)+1)/exp(λ*(d+1)^2);

    if isnan(ret_val)
        return 0.0;
    end

    return ret_val;
end

function unpol_ee_energy(at1::Atom, at2::Atom, at1_shell::Int, at2_shell::Int)
    # Calculates the unpolarized electron-electron interactions between the 
    # specified electron shells of atoms at1 and at2.
    at1_cloud = at1.electron_cloud_shells[at1_shell];
    at1_basis_size = length(at1_cloud.basis_function_amplitude);

    at2_cloud = at2.electron_cloud_shells[at2_shell];
    at2_basis_size = length(at2_cloud.basis_function_amplitude);

    d = norm(at1.coordinates - at2.coordinates);

    e_naive = 0;
    e_xc_sph = 0;
    e_xc_cyl = 0;

    for i in 1:at1_basis_size
        for j in 1:at2_basis_size
            c1 = at1_cloud.basis_function_amplitude[i];
            λ1 = at1_cloud.basis_function_decay[i];

            c2 = at2_cloud.basis_function_amplitude[j];
            λ2 = at2_cloud.basis_function_decay[j];

            c = c1 * c2;
            λ = (λ1 * λ2) / (λ1 + λ2);

            e_xc_sph += c*xc_sph(λ,d);
            e_xc_cyl += c*xc_cyl(λ,d);
            e_naive += c*naive_coulomb_integral(λ,d);
        end
    end

    return e_naive, e_xc_sph, e_xc_cyl;
end

function unpol_ee_energy_vv(at1::Atom, at2::Atom)
    # Calculates the unpolarized electron-electron interactions between the 
    # valence electron shells of atoms at1 and at2.
    at1_shell = length(at1.electron_cloud_shells);
    at2_shell = length(at2.electron_cloud_shells);

    e_naive, e_xc_sph, e_xc_cyl = 
        unpol_ee_energy(at1,at2,at1_shell,at2_shell);
    
    return e_naive, e_xc_sph, e_xc_cyl;
end

function unpol_ee_energy_vc(at1::Atom, at2::Atom)
    # Calculates the unpolarized electron-electron interactions between the 
    # valence electron shell of atom at1 and the core electron shells of atom 
    # at2.
    at1_shell = length(at1.electron_cloud_shells);
    at2_shells = 1:(length(at2.electron_cloud_shells)-1);

    e_naive = 0;
    e_xc_sph = 0;
    e_xc_cyl = 0;

    for at2_shell in at2_shells
        i_e_naive, i_e_xc_sph, i_e_xc_cyl = 
            unpol_ee_energy(at1,at2,at1_shell,at2_shell);
        
        e_naive += i_e_naive;
        e_xc_sph += i_e_xc_sph;
        e_xc_cyl += i_e_xc_cyl;
    end

    return e_naive, e_xc_sph, e_xc_cyl;
end

function unpol_ee_energy_cv(at1::Atom, at2::Atom)
    # Calculates the unpolarized electron-electron interactions between the 
    # core electron shells of atom at1 and the valence electron shell of atom 
    # at2.
    e_naive, e_xc_sph, e_xc_cyl = unpol_ee_energy_vc(at2,at1);
    return e_naive, e_xc_sph, e_xc_cyl;
end

function unpol_ee_energy_cc(at1::Atom, at2::Atom)
    # Calculates the unpolarized electron-electron interactions between the 
    # core electron shells of atoms at1 and at2.
    at1_shells = 1:(length(at1.electron_cloud_shells)-1);
    at2_shells = 1:(length(at2.electron_cloud_shells)-1);

    e_naive = 0;
    e_xc_sph = 0;
    e_xc_cyl = 0;

    for at1_shell in at1_shells
        for at2_shell in at2_shells
            ij_e_naive, ij_e_xc_sph, ij_e_xc_cyl = 
                unpol_ee_energy(at1,at2,at1_shell,at2_shell);
            
            e_naive += ij_e_naive;
            e_xc_sph += ij_e_xc_sph;
            e_xc_cyl += ij_e_xc_cyl;
        end
    end

    return e_naive, e_xc_sph, e_xc_cyl;
end

function ee_energy(at1::Atom, at2::Atom)
    # Calculates the polarized electron-electron interactions between the 
    # electron shells of atoms at1 and at2.
    ζ1 = at1.polarization_coefficient;
    ζ2 = at2.polarization_coefficient;

    e_naive_vv, e_xc_sph_vv, e_xc_cyl_vv = unpol_ee_energy_vv(at1,at2);
    e_naive_vc, e_xc_sph_vc, e_xc_cyl_vc = unpol_ee_energy_vc(at1,at2);
    e_naive_cv, e_xc_sph_cv, e_xc_cyl_cv = unpol_ee_energy_cv(at1,at2);
    e_naive_cc, e_xc_sph_cc, e_xc_cyl_cc = unpol_ee_energy_cc(at1,at2);

    e_naive = 
        ζ1*ζ2*e_naive_vv + ζ1*e_naive_vc + ζ2*e_naive_cv + e_naive_cc;
    e_xc_sph = 
        ζ1*ζ2*e_xc_sph_vv + ζ1*e_xc_sph_vc + ζ2*e_xc_sph_cv + e_xc_sph_cc;
    e_xc_cyl = 
        ζ1*ζ2*e_xc_cyl_vv + ζ1*e_xc_cyl_vc + ζ2*e_xc_cyl_cv + e_xc_cyl_cc;

    return e_naive, e_xc_sph, e_xc_cyl;
end

function unpol_en_energy(at1::Atom, at2::Atom, at1_shell::Int)
    # Calculates the unpolarized interactions between the specified electron 
    # shell of atom at1 and the nuclei of atom at2.
    at1_cloud = at1.electron_cloud_shells[at1_shell];
    at1_basis_size = length(at1_cloud.basis_function_amplitude);

    d = norm(at1.coordinates - at2.coordinates);

    e_naive = 0;
    e_xc_sph = 0;
    e_xc_cyl = 0;

    for i in 1:at1_basis_size
        c1 = at1_cloud.basis_function_amplitude[i];
        λ1 = at1_cloud.basis_function_decay[i];

        q2 = at2.atomic_number;

        c = c1 * q2;
        λ = λ1;

        e_xc_sph += c*xc_sph(λ,d);
        e_xc_cyl += c*xc_cyl(λ,d);
        e_naive -= c*naive_coulomb_integral(λ,d);
    end

    return e_naive, e_xc_sph, e_xc_cyl;
end

function unpol_en_energy_v(at1::Atom, at2::Atom)
    # Calculates the unpolarized interactions between the valence electron  
    # shell of atom at1 and the nuclei of atom at2.
    at1_shell = length(at1.electron_cloud_shells);
    e_naive, e_xc_sph, e_xc_cyl = unpol_en_energy(at1,at2,at1_shell);
    return e_naive, e_xc_sph, e_xc_cyl;
end

function unpol_en_energy_c(at1::Atom, at2::Atom)
    # Calculates the unpolarized interactions between the core electron shells
    # of atom at1 and the nuclei of atom at2.
    at1_shells = 1:(length(at1.electron_cloud_shells)-1);

    e_naive = 0;
    e_xc_sph = 0;
    e_xc_cyl = 0;

    for at1_shell in at1_shells
        i_e_naive, i_e_xc_sph, i_e_xc_cyl = 
            unpol_en_energy(at1,at2,at1_shell);
        
        e_naive += i_e_naive;
        e_xc_sph += i_e_xc_sph;
        e_xc_cyl += i_e_xc_cyl;
    end

    return e_naive, e_xc_sph, e_xc_cyl;
end

function en_energy(at1::Atom, at2::Atom)
    # Calculates the polarized interactions between the electron shells of atom 
    # at1 and the nuclei of atom at2.
    ζ1 = at1.polarization_coefficient;

    e_naive_v, e_xc_sph_v, e_xc_cyl_v = unpol_en_energy_v(at1,at2);
    e_naive_c, e_xc_sph_c, e_xc_cyl_c = unpol_en_energy_c(at1,at2);

    e_naive = ζ1*e_naive_v + e_naive_c;
    e_xc_sph = ζ1*e_xc_sph_v + e_xc_sph_c;
    e_xc_cyl = ζ1*e_xc_cyl_v + e_xc_cyl_c;

    return e_naive, e_xc_sph, e_xc_cyl;
end

function ne_energy(at1::Atom, at2::Atom)
    # Calculates the polarized interactions between the nuclei of atom at1 and 
    # the electron shells of atom at2.
    e_naive, e_xc_sph, e_xc_cyl = en_energy(at2,at1);
    return e_naive, e_xc_sph, e_xc_cyl;
end

function nn_energy(at1::Atom, at2::Atom)
    # Calculates the Couloumb interaction between the nuclei of atoms at1 and 
    # at2.
    q1 = at1.atomic_number;
    q2 = at2.atomic_number;

    d = norm(at1.coordinates - at2.coordinates);
    e_naive = (q1*q2)/d;

    return e_naive;
end

function rotate_molecule!(molecule::Molecule, Δθ::Matrix{Number})
    # Rotates the whole molecule by the rotation matrix Δθ.
    for i in eachindex(molecule.atoms)
        molecule.atoms[i].position = Δθ*(molecule.atoms[i].position);
    end
end

function polarization_matrix_problem(simulation::SimulationSystem)
    # Calculates the matrix problem that needs to be solved for the 
    # polarization lagrangian to be stationary.
    xc_coeffs = simulation.pol_e_xc_coeffs;

    molecules = simulation.system.molecules;
    charge = simulation.system.charge;

    aux_type = typeof(xc_coeffs.xc_a_1b[1].m);
    atoms_μ0 = simulation.basis_set_settings.atoms_μ0;

    tot_num_atoms = 0;
    atom_ind_base = Vector{Int}();
    for molecule in molecules
        if isempty(atom_ind_base)
            push!(atom_ind_base,1);
            tot_num_atoms += number_of_atoms(molecule);
        else            
            push!(atom_ind_base, tot_num_atoms + 1);
            tot_num_atoms += number_of_atoms(molecule);
        end
    end

    num_vars = 2*tot_num_atoms + 1;
    aux_M = zeros(aux_type, num_vars, num_vars);
    aux_Y = zeros(aux_type, num_vars, 1);
    aux_Y[end] -= charge;

    for ii in eachindex(molecules)
        molecule1 = molecules[ii];
        num_atoms1 = number_of_atoms(molecule1);

        # Chemical Potential
        # Lagrange multiplier to keep the number of electrons constant in the 
        # new fictitious electron density.
        for i in eachindex(molecule1.atoms)
            i_ind = atom_ind_base[ii] + i - 1;
            aux_M[i_ind,end] -= molecule1.atoms[i].valence_electrons;
            aux_M[end,i_ind] += molecule1.atoms[i].valence_electrons;
            aux_Y[end] += molecule1.atoms[i].valence_electrons;
        end

        # Atomic Chemical Potential Shifts
        for i in 1:num_atoms1
            z1 = molecule1.atoms[i].atomic_number;
            μ0 = atoms_μ0[z1];

            Δμ_atom_index = tot_num_atoms+i;
            aux_M[Δμ_atom_index,end] = -1.0;
            aux_M[Δμ_atom_index,Δμ_atom_index] = 1.0;
            aux_M[i,Δμ_atom_index] = -molecule1.atoms[i].valence_electrons;
            aux_Y[Δμ_atom_index] = -μ0;
        end

        # Kinetic Energy Contributions
        for i in 1:num_atoms1
            z1 = molecule1.atoms[i].atomic_number;
            ii0 = atom_ind_base[ii] + i - 1;

            valence_shell = molecule1.atoms[i].electron_cloud_shells[end];
            unpol_tf_ke = valence_shell.thomas_fermi_ke;
            unpol_vw_ke = valence_shell.von_weizsacker_ke;

            tf_fit_coeff_1 = aux_type(0.0);
            tf_fit_coeff_2 = aux_type(0.0);
            tf_fit_coeff_1 += (10.0/9.0)*xc_coeffs.xc_e_1b[z1];
            tf_fit_coeff_2 += (5.0/9.0)*xc_coeffs.xc_e_1b[z1];
            aux_M[ii0,ii0] += tf_fit_coeff_1*unpol_tf_ke;
            aux_Y[ii0] -= tf_fit_coeff_2*unpol_tf_ke;

            vw_fit_coeff = aux_type(0.0);
            vw_fit_coeff += xc_coeffs.xc_f_1b[z1];
            aux_Y[ii0] -= vw_fit_coeff*unpol_vw_ke;
        end

        # Fill the rest of the matrix
        for jj in eachindex(molecules)
            molecule2 = molecules[jj];
            num_atoms2 = number_of_atoms(molecule2);

            for i in 1:num_atoms1
                for j in 1:num_atoms2
                    at1 = molecule1.atoms[i];
                    at2 = molecule2.atoms[j];

                    d = norm(at1.coordinates - at2.coordinates);

                    z1 = at1.atomic_number;
                    z2 = at2.atomic_number;

                    ii0 = atom_ind_base[ii] + i - 1;
                    jj0 = atom_ind_base[jj] + j - 1;

                    # cloud-nuclei (valence-nuclei)
                    en_naive, en_xc_sph, en_xc_cyl = 
                        unpol_en_energy_v(at1,at2);

                    if (i == j) && (ii == jj)
                        aux_Y[ii0] -= en_naive;

                        sph_coeff = get_xc_coeff(d,xc_coeffs.xc_c_1b[z1]);
                        cyl_coeff = get_xc_coeff(d,xc_coeffs.xc_d_1b[z1]);
                        aux_Y[ii0] -= en_xc_sph*sph_coeff;
                        aux_Y[ii0] -= en_xc_cyl*cyl_coeff;
                    else
                        aux_Y[ii0] -= en_naive;

                        sph_coeff = get_xc_coeff(d,xc_coeffs.xc_c_2b[(z1,z2)]);
                        cyl_coeff = get_xc_coeff(d,xc_coeffs.xc_d_2b[(z1,z2)]);
                        aux_Y[ii0] -= en_xc_sph*sph_coeff;
                        aux_Y[ii0] -= en_xc_cyl*cyl_coeff;
                    end

                    # cloud-cloud (valence-core)
                    ee_naive, ee_xc_sph, ee_xc_cyl = 
                        unpol_ee_energy_vc(at1,at2);

                    if (i == j) && (ii == jj)
                        aux_Y[ii0] -= ee_naive;
                        
                        sph_coeff = get_xc_coeff(d,xc_coeffs.xc_a_1b[z1]);
                        cyl_coeff = get_xc_coeff(d,xc_coeffs.xc_b_1b[z1]);
                        aux_Y[ii0] -= ee_xc_sph*sph_coeff;
                        aux_Y[ii0] -= ee_xc_cyl*cyl_coeff;
                    else
                        aux_Y[ii0] -= ee_naive;

                        sph_coeff = get_xc_coeff(d,xc_coeffs.xc_a_2b[(z1,z2)]);
                        cyl_coeff = get_xc_coeff(d,xc_coeffs.xc_b_2b[(z1,z2)]);
                        aux_Y[ii0] -= ee_xc_sph*sph_coeff;
                        aux_Y[ii0] -= ee_xc_cyl*cyl_coeff;
                    end

                    # cloud-cloud (valence-valence)
                    ee_naive, ee_xc_sph, ee_xc_cyl = 
                        unpol_ee_energy_vv(at1,at2);

                    if (i == j) && (ii == jj)
                        aux_M[ii0,jj0] += ee_naive;

                        sph_coeff = get_xc_coeff(d,xc_coeffs.xc_a_1b[z1]);
                        cyl_coeff = get_xc_coeff(d,xc_coeffs.xc_b_1b[z1]);
                        aux_M[ii0,jj0] += ee_xc_sph*sph_coeff;
                        aux_M[ii0,jj0] += ee_xc_cyl*cyl_coeff;
                    else
                        aux_M[ii0,jj0] += ee_naive;

                        sph_coeff = get_xc_coeff(d,xc_coeffs.xc_a_2b[(z1,z2)]);
                        cyl_coeff = get_xc_coeff(d,xc_coeffs.xc_b_2b[(z1,z2)]);
                        aux_M[ii0,jj0] += ee_xc_sph*sph_coeff;
                        aux_M[ii0,jj0] += ee_xc_cyl*cyl_coeff;
                    end
                end
            end
        end
    end

    return aux_M, aux_Y;
end

function polarize_molecules!(simulation::SimulationSystem)
    # Calculates and sets the polarization coefficients of the atoms.
    aux_m, aux_y = polarization_matrix_problem(simulation);
    minimizer = aux_m \ aux_y;

    # display(hcat(aux_m,aux_y));

    tot_num_atoms = 0;
    atom_ind_base = [];
    for molecule in simulation.system.molecules
        if isempty(atom_ind_base)
            push!(atom_ind_base,1);
            tot_num_atoms += number_of_atoms(molecule);
        else            
            push!(atom_ind_base, tot_num_atoms + 1);
            tot_num_atoms += number_of_atoms(molecule);
        end
    end

    simulation.system.chemical_potential = minimizer[end];

    for i in eachindex(simulation.system.molecules)
        i0 = atom_ind_base[i];
        atoms = simulation.system.molecules[i].atoms;

        for ii in eachindex(atoms)
            ζ = minimizer[i0+ii-1];
            atoms[ii].polarization_coefficient = ζ;
        end

        simulation.system.molecules[i].atoms = atoms;
    end
end

function system_energies(simulation::SimulationSystem)
    # Returns the kinetic, naive and xc energy contributions of our 
    # simplified model.
    xc_coeffs = simulation.tot_e_xc_coeffs;
    molecules = simulation.system.molecules;
    num_molecules = length(molecules);
    
    aux_type = typeof(xc_coeffs.xc_a_1b[1].m);
    xc_energy = aux_type(0.0);
    naive_energy = aux_type(0.0);
    kinetic_energy = aux_type(0.0);

    function compute_nn_energy(at1::Atom, at2::Atom)
        d = norm(at1.coordinates - at2.coordinates);

        if d > atoms_dist_cutoff()
            naive_energy += nn_energy(at1,at2);
        end

        return;
    end

    function compute_en_naive_xc_energies(at1::Atom, at2::Atom)
        z1 = at1.atomic_number;
        z2 = at2.atomic_number;

        d = norm(at1.coordinates - at2.coordinates);
        en_naive, en_xc_sph, en_xc_cyl = en_energy(at1,at2);

        if d < atoms_dist_cutoff()
            naive_energy += en_naive;

            sph_coeff = get_xc_coeff(d,xc_coeffs.xc_c_1b[z1]);
            cyl_coeff = get_xc_coeff(d,xc_coeffs.xc_d_1b[z1]);
            xc_energy += en_xc_sph*sph_coeff;
            xc_energy += en_xc_cyl*cyl_coeff;
        else
            naive_energy += en_naive;

            sph_coeff = get_xc_coeff(d,xc_coeffs.xc_c_2b[(z1,z2)]);
            cyl_coeff = get_xc_coeff(d,xc_coeffs.xc_d_2b[(z1,z2)]);
            xc_energy += en_xc_sph*sph_coeff;
            xc_energy += en_xc_cyl*cyl_coeff;
        end

        return;
    end

    function compute_ne_naive_xc_energies(at1::Atom, at2::Atom)
        d = norm(at1.coordinates - at2.coordinates);

        if d > atoms_dist_cutoff()
            compute_en_naive_xc_energies(at2,at1);
        end

        return;
    end

    function compute_ee_naive_xc_energies(at1::Atom, at2::Atom)
        z1 = at1.atomic_number;
        z2 = at2.atomic_number;
        
        d = norm(at1.coordinates - at2.coordinates);

        # valence-valence
        ee_naive, ee_xc_sph, ee_xc_cyl = ee_energy(at1,at2);

        if d < atoms_dist_cutoff()
            naive_energy += ee_naive;

            sph_coeff = get_xc_coeff(d,xc_coeffs.xc_a_1b[z1]);
            cyl_coeff = get_xc_coeff(d,xc_coeffs.xc_b_1b[z1]);
            xc_energy += ee_xc_sph*sph_coeff;
            xc_energy += ee_xc_cyl*cyl_coeff;
        else
            naive_energy += ee_naive;

            sph_coeff = get_xc_coeff(d,xc_coeffs.xc_a_2b[(z2,z2)]);
            cyl_coeff = get_xc_coeff(d,xc_coeffs.xc_b_2b[(z2,z2)]);
            xc_energy += ee_xc_sph*sph_coeff;
            xc_energy += ee_xc_cyl*cyl_coeff;
        end

        return;
    end

    for ii in 1:num_molecules
        molecule1 = molecules[ii];
        num_atoms1 = number_of_atoms(molecule1);

        # Kinetic Energy
        for i in 1:num_atoms1
            at1 = molecule1.atoms[i];
            num_shells = length(at1.electron_cloud_shells);
            for j in 1:num_shells
                z1 = atom_atomic_number(molecule1,i);
                ζ1 = atom_polarization_coeff(molecule1,i);
                e_shell = at1.electron_cloud_shells[j];

                unpol_tf_ke = e_shell.thomas_fermi_ke;
                unpol_vw_ke = e_shell.von_weizsacker_ke;

                tf_fit_coeff = aux_type(0.0);

                if i < num_shells
                    tf_fit_coeff += (2.0/3.0)*xc_coeffs.xc_e_1b[z1];
                    tf_fit_coeff += (1.0/3.0)*xc_coeffs.xc_e_1b[z1];
                else
                    tf_fit_coeff += (2.0/3.0)*ζ1*ζ1*xc_coeffs.xc_e_1b[z1];
                    tf_fit_coeff += (1.0/3.0)*ζ1*xc_coeffs.xc_e_1b[z1];
                end

                vw_fit_coeff = aux_type(0.0);

                if i < num_shells
                    vw_fit_coeff += xc_coeffs.xc_f_1b[z1];
                else
                    vw_fit_coeff += ζ1*xc_coeffs.xc_f_1b[z1];
                end

                kinetic_energy += tf_fit_coeff*unpol_tf_ke;
                kinetic_energy += vw_fit_coeff*unpol_vw_ke;
            end
        end

        # Intramolecular interactions
        for i in 1:num_atoms1
            for j in i:num_atoms1
                at1 = molecule1.atoms[i];
                at2 = molecule1.atoms[j];

                compute_nn_energy(at1,at2);
                compute_en_naive_xc_energies(at1,at2);
                compute_ne_naive_xc_energies(at1,at2);
                compute_ee_naive_xc_energies(at1,at2);
            end
        end

        # Intermolecular interactions
        for jj in (ii+1):num_molecules
            molecule2 = molecules[jj];
            num_atoms2 = number_of_atoms(molecule2);

            for i in 1:num_atoms1
                for j in 1:num_atoms2
                    at1 = molecule1.atoms[i];
                    at2 = molecule2.atoms[j];

                    compute_nn_energy(at1,at2);
                    compute_en_naive_xc_energies(at1,at2);
                    compute_ne_naive_xc_energies(at1,at2);
                    compute_ee_naive_xc_energies(at1,at2);
                end
            end            
        end
    end

    return naive_energy, kinetic_energy, xc_energy;
end

function total_energy(simulation::SimulationSystem)
    # Returns the sum of all three energy contributions.
    naive_energy, kinetic_energy, xc_energy = system_energies(simulation);
    return naive_energy + kinetic_energy + xc_energy;
end

function save_fitted_coeffs(simulation::SimulationSystem)
    save_1b_pol_e_coeffs(simulation);
    save_2b_pol_e_coeffs(simulation);

    save_1b_tot_e_coeffs(simulation);
    save_2b_tot_e_coeffs(simulation);

    return;
end

function initialize_simulation_environment()
    basis_set_settings = load_basis_set();
    max_atomic_number = basis_set_settings.max_atomic_number;
    
    xc_a_1b = Dict{Int,XCCoefficient}();
    xc_b_1b = Dict{Int,XCCoefficient}();
    xc_c_1b = Dict{Int,XCCoefficient}();
    xc_d_1b = Dict{Int,XCCoefficient}();
    xc_e_1b = Dict{Int,Float64}();
    xc_f_1b = Dict{Int,Float64}();
    
    xc_a_2b = Dict{Tuple{Int,Int},XCCoefficient}();
    xc_b_2b = Dict{Tuple{Int,Int},XCCoefficient}();
    xc_c_2b = Dict{Tuple{Int,Int},XCCoefficient}();
    xc_d_2b = Dict{Tuple{Int,Int},XCCoefficient}();

    pol_e_xc_coeffs = EmpiricalXCCoefficients(max_atomic_number,
        xc_a_1b,xc_b_1b,xc_c_1b,xc_d_1b,xc_e_1b,xc_f_1b,
        xc_a_2b,xc_b_2b,xc_c_2b,xc_d_2b);
    tot_e_xc_coeffs = EmpiricalXCCoefficients(max_atomic_number,
        xc_a_1b,xc_b_1b,xc_c_1b,xc_d_1b,xc_e_1b,xc_f_1b,
        xc_a_2b,xc_b_2b,xc_c_2b,xc_d_2b);

    system = MolecularSystem();
    simulation = SimulationSystem(system,pol_e_xc_coeffs,
        tot_e_xc_coeffs,basis_set_settings);

    load_fitted_coeffs!(simulation);

    return simulation;
end

function make_monoatomic_system(Z::Int,charge::Int)
    # Makes a molecule object with a single atom whose atomic number is Z.
    simulation = initialize_simulation_environment();
    simulation.system.charge = charge;
    reference_atoms = simulation.basis_set_settings.reference_atoms;

    molecule = Molecule();
    molecule.atoms = [copy(reference_atoms[Z])];

    z_eff = molecule.atoms[1].valence_electrons;
    molecule.atoms[1].polarization_coefficient = (z_eff - charge) / z_eff;

    simulation.system.molecules = [molecule];
    return simulation;
end

function make_monoatomic_system(element::String,charge::Int)
    Z = get_atomic_number(element);
    return make_monoatomic_system(Z,charge);
end

function make_diatomic_molecule(
    reference_atoms::Dict{Int,Atom}, Z1::Int,Z2::Int,d::Number,charge::Int)
    # Makes a diatomic molecule with a diatomic separation d in Bohr with the 
    # specified charge.
    at1 = copy(reference_atoms[Z1]);
    at2 = copy(reference_atoms[Z2]);

    z1_eff = at1.valence_electrons;
    z2_eff = at2.valence_electrons;

    charge1 = charge*z1_eff/(z1_eff+z2_eff);
    charge2 = charge*z2_eff/(z1_eff+z2_eff);

    at1.polarization_coefficient = (z1_eff - charge1) / z1_eff;
    at2.polarization_coefficient = (z2_eff - charge2) / z2_eff;

    at1.coordinates = [0.0, 0.0, +d/2.0];
    at2.coordinates = [0.0, 0.0, -d/2.0];

    molecule = Molecule();
    molecule.atoms = [at1,at2];

    return molecule;
end

function make_diatomic_system(Z1::Int,Z2::Int,d::Number,charge::Int)
    # Makes a simulation structure with a diatomic molecule with a diatomic 
    # separation d in Bohr with the specified charge.
    simulation = initialize_simulation_environment();
    simulation.system.charge = charge;
    reference_atoms = simulation.basis_set_settings.reference_atoms;

    molecule = make_diatomic_molecule(reference_atoms,Z1,Z2,d,charge);

    simulation.system.molecules = [molecule];
    return simulation;
end

function make_diatomic_system(
    element1::String, element2::String,d::Number,charge::Int)
    # Makes a diatomic molecule with a diatomic separation d in Bohr with the 
    # specified charge.
    Z1 = get_atomic_number(element1);
    Z2 = get_atomic_number(element2);
    return make_diatomic_system(Z1,Z2,d,charge);
end

function make_diatomic_system(Z::Int,d::Number,charge::Int)
    # Makes a homonuclear diatomic molecule with a diatomic separation d in 
    # Bohr with the specified charge.
    return make_diatomic_system(Z,Z,d,charge);
end

function make_diatomic_system(element::String,d::Number,charge::Int)
    # Makes a homonuclear diatomic molecule with a diatomic separation d in 
    # Bohr with the specified charge.
    return make_diatomic_system(element,element,d,charge);
end

function full_model_reset()
    # Sets all the empirical coefficients to zero and recalculates the 
    # Thomas-Fermi and von Weizacker kinetic energies for the atoms allowed 
    # in the basis-set.
    reset_xc_coeffs();
    load_basis_set(true);

    return;
end
