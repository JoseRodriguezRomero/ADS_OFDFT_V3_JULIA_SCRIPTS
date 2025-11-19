using Interpolations;
using Base.Threads;

include("ForceField.jl")

mutable struct ParsedFile
    atomic_separation::Float64;
    total_energy::Float64;
    partial_charge_1::Float64;
    partial_charge_2::Float64;
    atomic_number_1::Int;
    atomic_number_2::Int;
    HOMO_energy::Float64;
    LUMO_energy::Float64;
    chemical_potential::Float64;
    charge::Int;
end

function read_scanned_data(file_name::String)
    file_data = Vector{ParsedFile}();
    lines = readlines(file_name);

    atom_elem_label_1 = split(lines[2])[3];
    atom_elem_label_2 = split(lines[2])[5];

    atom_elem_label_1 = String(atom_elem_label_1[2:end]);
    atom_elem_label_2 = String(atom_elem_label_2[2:end]);

    atomic_number_1 = get_atomic_number(atom_elem_label_1);
    atomic_number_2 = get_atomic_number(atom_elem_label_2);

    for i in eachindex(lines)
        if i < 3
            continue;
        end

        line_splitted = split(lines[i]);

        atomic_separation = parse(Float64,line_splitted[1]);
        total_energy = parse(Float64,line_splitted[2]);
        partial_charge_1 = parse(Float64,line_splitted[3]);
        partial_charge_2 = parse(Float64,line_splitted[4]);
        HOMO_energy = parse(Float64,line_splitted[5]);
        LUMO_energy = parse(Float64,line_splitted[6]);
        chemical_potential = 0.0;
        charge = round(Int, partial_charge_1 + partial_charge_2);

        angstrom_to_bohr = 1.88973;
        atomic_separation *= angstrom_to_bohr;

        parsed_file = ParsedFile(atomic_separation, total_energy, 
            partial_charge_1, partial_charge_2, atomic_number_1,
            atomic_number_2, HOMO_energy, LUMO_energy, chemical_potential, 
            charge);
        
        push!(file_data, parsed_file);
    end

    return file_data;
end

function read_all_scanned_data(Z1::Int, Z2::Int)
    neutral_file_name = "../Reference Data/diatomic_data/Neutral/ParsedData/";
    cation_file_name = "../Reference Data/diatomic_data/Cation/ParsedData/";
    anion_file_name = "../Reference Data/diatomic_data/Anion/ParsedData/";

    z1_elem_symbol = get_element_symbol(Z1);
    z2_elem_symbol = get_element_symbol(Z2);
    parsed_file_name = "scan_"*z1_elem_symbol*z2_elem_symbol*".txt";

    neutral_file_name *= parsed_file_name;
    cation_file_name *= parsed_file_name;
    anion_file_name *= parsed_file_name;

    neutral_data = read_scanned_data(neutral_file_name);
    cation_data = read_scanned_data(cation_file_name);
    anion_data = read_scanned_data(anion_file_name);

    return neutral_data, cation_data, anion_data;
end

function read_all_scanned_data(Z::Int)
    return read_all_scanned_data(Z,Z);
end

function get_homo_energies(data::Vector{ParsedFile})
    r = zeros(Float64,length(data));
    e = zeros(Float64,length(data));

    for i in eachindex(data)
        r[i] = data[i].atomic_separation;
        e[i] = data[i].HOMO_energy;
    end

    return r, e;
end

function get_lumo_energies(data::Vector{ParsedFile})
    r = zeros(Float64,length(data));
    e = zeros(Float64,length(data));

    for i in eachindex(data)
        r[i] = data[i].atomic_separation;
        e[i] = data[i].LUMO_energy;
    end

    return r, e;
end

function get_total_energies(data::Vector{ParsedFile})
    r = zeros(Float64,length(data));
    e = zeros(Float64,length(data));

    for i in eachindex(data)
        r[i] = data[i].atomic_separation;
        e[i] = data[i].total_energy;
    end

    return r, e;
end

function make_system_from_parsed_file(parsed_file::ParsedFile)
    # Creates a simulation for a molecule using the data in parsed_file.
    Z1 = Int(parsed_file.atomic_number_1);
    Z2 = Int(parsed_file.atomic_number_2);
    d = parsed_file.atomic_separation;
    charge = parsed_file.charge;
    chemical_potential = parsed_file.chemical_potential;

    simulation = make_diatomic_system(Z1,Z2,d,charge);

    z1_eff = simulation.system.molecules[1].atoms[1].valence_electrons;
    z2_eff = simulation.system.molecules[1].atoms[2].valence_electrons;

    ζ1 = (z1_eff - parsed_file.partial_charge_1) / z1_eff;
    ζ2 = (z2_eff - parsed_file.partial_charge_2) / z2_eff;
    
    simulation.system.chemical_potential = chemical_potential;
    simulation.system.molecules[1].atoms[1].polarization_coefficient = ζ1;
    simulation.system.molecules[1].atoms[2].polarization_coefficient = ζ2;

    return simulation;
end

function set_diatomic_system_to_parsed_file!(
    simulation::SimulationSystem, parsed_file::ParsedFile)
    # Sets a simulation for a molecule using the data in parsed_file.
    Z1 = parsed_file.atomic_number_1;
    Z2 = parsed_file.atomic_number_2;
    q1 = parsed_file.partial_charge_1;
    q2 = parsed_file.partial_charge_2;
    d = parsed_file.atomic_separation;
    charge = parsed_file.charge;
    reference_atoms = simulation.basis_set_settings.reference_atoms;
    molecule = make_diatomic_molecule(reference_atoms,Z1,Z2,d,charge);

    z1_eff = molecule.atoms[1].valence_electrons;
    z2_eff = molecule.atoms[2].valence_electrons;

    molecule.atoms[1].polarization_coefficient = (z1_eff-q1)/z1_eff;
    molecule.atoms[2].polarization_coefficient = (z2_eff-q2)/z2_eff;

    simulation.system.molecules = [molecule];

    simulation.system.charge = parsed_file.charge;
    simulation.system.energy = parsed_file.total_energy;
    simulation.system.chemical_potential = parsed_file.chemical_potential;

    return;
end

function set_fitted_coeffs!(coeffs::PolarizationEnergyCoefficients, 
    Z::Int, fitted_coeffs::Vector)
    # Store the 1B fitted coefficients
    xc_a_1b_ind = 1;
    xc_b_1b_ind = 2;
    xc_c_1b_ind = 3;
    xc_d_1b_ind = 4;
    ke_e_1b_ind = 5;
    ke_f_1b_ind = 6;

    coeffs.xc_coeffs.xc_a_1b[Z] = fitted_coeffs[xc_a_1b_ind];
    coeffs.xc_coeffs.xc_b_1b[Z] = fitted_coeffs[xc_b_1b_ind];
    coeffs.xc_coeffs.xc_c_1b[Z] = fitted_coeffs[xc_c_1b_ind];
    coeffs.xc_coeffs.xc_d_1b[Z] = fitted_coeffs[xc_d_1b_ind];
    coeffs.ke_coeffs.ke_e_1b[Z] = fitted_coeffs[ke_e_1b_ind];
    coeffs.ke_coeffs.ke_f_1b[Z] = fitted_coeffs[ke_f_1b_ind];

    # Store the 2B fitted coefficients
    xc_a_2b_ind         =  7;
    xc_b_2b_ind         =  8;
    xc_c_2b_ind         =  9;
    xc_d_2b_ind         = 10;

    coeffs.xc_coeffs.xc_a_2b[(Z,Z)] = fitted_coeffs[xc_a_2b_ind];
    coeffs.xc_coeffs.xc_b_2b[(Z,Z)] = fitted_coeffs[xc_b_2b_ind];
    coeffs.xc_coeffs.xc_c_2b[(Z,Z)] = fitted_coeffs[xc_c_2b_ind];
    coeffs.xc_coeffs.xc_d_2b[(Z,Z)] = fitted_coeffs[xc_d_2b_ind];

    return;
end

function set_fitted_coeffs!(coeffs::TotalEnergyCoefficients, 
    Z::Int, fitted_coeffs::Vector)
    # Store the 1B fitted coefficients
    xc_a_1b_ind = 1;
    xc_b_1b_ind = 2;
    xc_c_1b_ind = 3;
    xc_d_1b_ind = 4;
    ke_e_1b_ind = 5;
    ke_f_1b_ind = 6;

    coeffs.xc_coeffs.xc_a_1b[Z] = fitted_coeffs[xc_a_1b_ind];
    coeffs.xc_coeffs.xc_b_1b[Z] = fitted_coeffs[xc_b_1b_ind];
    coeffs.xc_coeffs.xc_c_1b[Z] = fitted_coeffs[xc_c_1b_ind];
    coeffs.xc_coeffs.xc_d_1b[Z] = fitted_coeffs[xc_d_1b_ind];
    coeffs.ke_coeffs.ke_e_1b[Z] = fitted_coeffs[ke_e_1b_ind];
    coeffs.ke_coeffs.ke_f_1b[Z] = fitted_coeffs[ke_f_1b_ind];

    # Store the 2B fitted coefficients
    xc_a_2b_ind                            =  7;
    xc_b_2b_ind                            =  8;
    xc_c_2b_ind                            =  9;
    xc_d_2b_ind                            = 10;
    non_pol_morse_depth_ind                = 11;
    non_pol_morse_stiffness_parameter_ind  = 12;
    non_pol_morse_equilibrium_distance_ind = 13;

    coeffs.xc_coeffs.xc_a_2b[(Z,Z)] = fitted_coeffs[xc_a_2b_ind];
    coeffs.xc_coeffs.xc_b_2b[(Z,Z)] = fitted_coeffs[xc_b_2b_ind];
    coeffs.xc_coeffs.xc_c_2b[(Z,Z)] = fitted_coeffs[xc_c_2b_ind];
    coeffs.xc_coeffs.xc_d_2b[(Z,Z)] = fitted_coeffs[xc_d_2b_ind];

    coeffs.non_polarizable_coeffs.depth[(Z,Z)] = 
        fitted_coeffs[non_pol_morse_depth_ind]^2;
    coeffs.non_polarizable_coeffs.stiffness_parameter[(Z,Z)] = 
        fitted_coeffs[non_pol_morse_stiffness_parameter_ind]^2;
    coeffs.non_polarizable_coeffs.equilibrium_distance[(Z,Z)] = 
        fitted_coeffs[non_pol_morse_equilibrium_distance_ind]^2;

    return;
end

function set_fitted_tot_e_coeffs!(simulation::SimulationSystem,
    Z::Int, fitted_coeffs::Vector)
    set_fitted_coeffs!(simulation.tot_e_coeffs,Z,fitted_coeffs);
    return;
end

function set_fitted_pol_e_coeffs!(simulation::SimulationSystem,
    Z::Int, fitted_coeffs::Vector)
    set_fitted_coeffs!(simulation.pol_e_coeffs,Z,fitted_coeffs);
    return;
end

function set_fitted_coeffs!(coeffs::PolarizationEnergyCoefficients,
    Z1::Int, Z2::Int, fitted_coeffs::Vector)
    # Store the 2B fitted coefficients
    xc_a_2b_ind         =  1;
    xc_b_2b_ind         =  2;
    xc_c_2b_ind         =  3;
    xc_d_2b_ind         =  4;

    coeffs.xc_coeffs.xc_a_2b[(Z1,Z2)] = fitted_coeffs[xc_a_2b_ind];
    coeffs.xc_coeffs.xc_b_2b[(Z1,Z2)] = fitted_coeffs[xc_b_2b_ind];
    coeffs.xc_coeffs.xc_c_2b[(Z1,Z2)] = fitted_coeffs[xc_c_2b_ind];
    coeffs.xc_coeffs.xc_d_2b[(Z1,Z2)] = fitted_coeffs[xc_d_2b_ind];

    return;
end

function set_fitted_coeffs!(coeffs::TotalEnergyCoefficients,
    Z1::Int, Z2::Int, fitted_coeffs::Vector)
    # Store the 2B fitted coefficients
    xc_a_2b_ind                            = 1;
    xc_b_2b_ind                            = 2;
    xc_c_2b_ind                            = 3;
    xc_d_2b_ind                            = 4;
    non_pol_morse_depth_ind                = 5;
    non_pol_morse_stiffness_parameter_ind  = 6;
    non_pol_morse_equilibrium_distance_ind = 7;

    coeffs.xc_coeffs.xc_a_2b[(Z1,Z2)] = fitted_coeffs[xc_a_2b_ind];
    coeffs.xc_coeffs.xc_b_2b[(Z1,Z2)] = fitted_coeffs[xc_b_2b_ind];
    coeffs.xc_coeffs.xc_c_2b[(Z1,Z2)] = fitted_coeffs[xc_c_2b_ind];
    coeffs.xc_coeffs.xc_d_2b[(Z1,Z2)] = fitted_coeffs[xc_d_2b_ind];

    coeffs.non_polarizable_coeffs.depth[(Z1,Z2)] = 
        fitted_coeffs[non_pol_morse_depth_ind]^2;
    coeffs.non_polarizable_coeffs.stiffness_parameter[(Z1,Z2)] = 
        fitted_coeffs[non_pol_morse_stiffness_parameter_ind]^2;
    coeffs.non_polarizable_coeffs.equilibrium_distance[(Z1,Z2)] = 
        fitted_coeffs[non_pol_morse_equilibrium_distance_ind]^2;

    return;
end

function set_fitted_pol_e_coeffs!(simulation::SimulationSystem,
    Z1::Int, Z2::Int, fitted_coeffs::Vector)
    set_fitted_coeffs!(simulation.pol_e_coeffs,Z1,Z2,fitted_coeffs);
    set_fitted_coeffs!(simulation.pol_e_coeffs,Z2,Z1,fitted_coeffs);
    return;
end

function set_fitted_tot_e_coeffs!(simulation::SimulationSystem,
    Z1::Int, Z2::Int, fitted_coeffs::Vector)
    set_fitted_coeffs!(simulation.tot_e_coeffs,Z1,Z2,fitted_coeffs);
    set_fitted_coeffs!(simulation.tot_e_coeffs,Z2,Z1,fitted_coeffs);
    return;
end

function get_reference_atom_total_energy()
    base_dir = "../Reference Data/atoms_data/";
    neutral_file_ID = base_dir*"Neutral/elem_total_energy.txt";
    cation_file_ID = base_dir*"Cation/elem_total_energy.txt";
    anion_file_ID = base_dir*"Anion/elem_total_energy.txt";

    function ReadFile(file_ID::String)
        lines = readlines(file_ID);
        energies = zeros(Float64,length(lines));

        for i in eachindex(lines)
            energies[i] = parse(Float64,lines[i]);
        end

        return energies;
    end

    neutral_energies = ReadFile(neutral_file_ID);
    cation_energies = ReadFile(cation_file_ID);
    anion_energies = ReadFile(anion_file_ID);

    return neutral_energies, cation_energies, anion_energies;
end

function get_reference_atom_chemical_potential()
    neutral_energies, cation_energies, anion_energies = 
        get_reference_atom_total_energy();

    m = zeros(Float64,3,3);
    y = zeros(Float64,3);

    neutral_μ = zeros(Float64,length(neutral_energies));
    cation_μ = zeros(Float64,length(cation_energies));
    anion_μ = zeros(Float64,length(anion_energies));

    for i in eachindex(neutral_energies)
        y[1] = neutral_energies[i];
        m[1,1] = 0.0^2;
        m[1,2] = 0.0^1;
        m[1,3] = 1.0;

        y[2] = cation_energies[i];
        m[2,1] = 1.0^2;
        m[2,2] = 1.0^1;
        m[2,3] = 1.0;

        y[3] = anion_energies[i];
        m[3,1] = (-1.0)^2;
        m[3,2] = (-1.0)^1;
        m[3,3] = 1.0;

        m_coeffs = (m \ y);

        function fitted_chemical_potential(charge::Number)
            return -(2*m_coeffs[1]*charge + m_coeffs[2]);
        end

        neutral_μ[i] = fitted_chemical_potential(0);
        cation_μ[i] = fitted_chemical_potential(1);
        anion_μ[i] = fitted_chemical_potential(-1);
    end

    return neutral_μ, cation_μ, anion_μ;
end

function sanitize_data(scanned_data::Vector{ParsedFile})
    mid_ind = round(Int32, length(scanned_data)/2.0);

    i0 = 1;
    i1 = length(scanned_data);

    diff_threshold = 0.05;

    for i in reverse((i0+1):mid_ind)
        d0 = scanned_data[i+1].atomic_separation;
        u0 = scanned_data[i+1].total_energy;

        d1 = scanned_data[i].atomic_separation;
        u1 = scanned_data[i].total_energy;

        d2 = scanned_data[i-1].atomic_separation;
        u2 = scanned_data[i-1].total_energy;

        m = (u1-u0)/(d1-d0);
        b = -((d0*u1-d1*u0)/(d1-d0));

        U2 = m*d2 + b;

        if abs(U2-u2) > diff_threshold
            i0 = i;
            break;
        end
    end

    for i in (mid_ind:(i1-1))
        d0 = scanned_data[i+1].atomic_separation;
        u0 = scanned_data[i+1].total_energy;

        d1 = scanned_data[i].atomic_separation;
        u1 = scanned_data[i].total_energy;

        d2 = scanned_data[i-1].atomic_separation;
        u2 = scanned_data[i-1].total_energy;

        m = (u1-u0)/(d1-d0);
        b = -((d0*u1-d1*u0)/(d1-d0));

        U2 = m*d2 + b;

        if abs(U2-u2) > diff_threshold
            i1 = i;
            break;
        end
    end

    return scanned_data[i0:i1];
end

function read_all_sanitized_data(Z1::Int, Z2::Int, clip_data::Bool = false)
    # If clip_data is true, then, the data is interpoated in a uniform grid
    # in the largest common interatomic range in the neutral, cation and 
    # anion energy scans. 
    neutral_data, cation_data, anion_data = read_all_scanned_data(Z1,Z2);

    neutral_data = sanitize_data(neutral_data);
    cation_data = sanitize_data(cation_data);
    anion_data = sanitize_data(anion_data);

    if clip_data == false
        return neutral_data, cation_data, anion_data;
    end

    r0 = [
        neutral_data[1].atomic_separation,
        cation_data[1].atomic_separation,
        anion_data[1].atomic_separation
    ];

    r1 = [
        neutral_data[end].atomic_separation,
        cation_data[end].atomic_separation,
        anion_data[end].atomic_separation
    ];

    r0 = maximum(r0);
    r1 = minimum(r1);
    num_samples = 75;

    function interpolate_between(data::Vector{ParsedFile}, 
        r0::Number,r1::Number,num_samples::Int)
        # Auxiliary function to create interpolated parsed files structs, so 
        # that neutral, cation and anion classes share the same interatomic 
        # separation gridspace.
        atomic_number_1 = data[1].atomic_number_1;
        atomic_number_2 = data[1].atomic_number_2;
        charge = data[1].charge;

        v_atomic_separation = zeros(Float64,0);
        v_total_energy = zeros(Float64,0);
        v_partial_charge_1 = zeros(Float64,0);
        v_partial_charge_2 = zeros(Float64,0);
        v_HOMO_energy = zeros(Float64,0);
        v_LUMO_energy = zeros(Float64,0);

        for i in eachindex(data)
            append!(v_atomic_separation,data[i].atomic_separation);
            append!(v_total_energy,data[i].total_energy);
            append!(v_partial_charge_1,data[i].partial_charge_1);
            append!(v_partial_charge_2,data[i].partial_charge_2);
            append!(v_HOMO_energy,data[i].HOMO_energy);
            append!(v_LUMO_energy,data[i].LUMO_energy);
        end

        itp_atomic_separation = interpolate(v_atomic_separation, 
            BSpline(Cubic(Line())), OnGrid());
        itp_total_energy = interpolate(v_total_energy, 
            BSpline(Cubic(Line())), OnGrid());
        itp_partial_charge_1 = interpolate(v_partial_charge_1, 
            BSpline(Cubic(Line())), OnGrid());
        itp_partial_charge_2 = interpolate(v_partial_charge_2, 
            BSpline(Cubic(Line())), OnGrid());
        itp_HOMO_energy = interpolate(v_HOMO_energy, 
            BSpline(Cubic(Line())), OnGrid());
        itp_LUMO_energy = interpolate(v_LUMO_energy, 
            BSpline(Cubic(Line())), OnGrid());

        R0 = data[1].atomic_separation;
        R1 = data[end].atomic_separation;

        i0 = 1 + ((r0 - R0) / (R1 - R0)) * (length(data) - 1);
        i1 = 1 + ((r1 - R0) / (R1 - R0)) * (length(data) - 1);
        i_itp = collect(i0:((i1-i0)/(num_samples-1)):i1);
        itp_data = Vector{ParsedFile}();
        resize!(itp_data,length(i_itp));

        for i in eachindex(i_itp)
            curr_i = i_itp[i];
            atomic_separation = itp_atomic_separation(curr_i);
            total_energy = itp_total_energy(curr_i);
            partial_charge_1 = itp_partial_charge_1(curr_i);
            partial_charge_2 = itp_partial_charge_2(curr_i);
            HOMO_energy = itp_HOMO_energy(curr_i);
            LUMO_energy = itp_LUMO_energy(curr_i);
            chemical_potential = 0.0;
            new_file = ParsedFile(atomic_separation,total_energy,
            partial_charge_1,partial_charge_2,atomic_number_1, atomic_number_2,
            HOMO_energy, LUMO_energy,chemical_potential,charge);
            itp_data[i] = new_file;
        end

        return itp_data;
    end

    function set_chemical_potentials!(neutral_data::Vector{ParsedFile},
        cation_data::Vector{ParsedFile}, anion_data::Vector{ParsedFile})
        # Sets the chemical potential of the reference files by fitting a 
        # quadratic with respect the total energy as a function of the 
        # system charge.
        for i in eachindex(neutral_data)
            neutral_d = neutral_data[i].atomic_separation;
            cation_d = cation_data[i].atomic_separation;
            anion_d = anion_data[i].atomic_separation;

            error_distance_threshold = 1.0E-4;
            cond1 = abs(neutral_d - cation_d) <= error_distance_threshold;
            cond2 = abs(neutral_d - anion_d) <= error_distance_threshold;

            if (cond1 == false) || (cond2 == false)
                display([neutral_d,cation_d,anion_d])
                println("Incorrect data!");
                return;
            end

            neutral_tot_e = neutral_data[i].total_energy;
            cation_tot_e = cation_data[i].total_energy;
            anion_tot_e = anion_data[i].total_energy;

            m = zeros(Float64,3,3);
            y = zeros(Float64,3);

            y[1] = neutral_tot_e;
            m[1,1] = 0.0^2;
            m[1,2] = 0.0^1;
            m[1,3] = 1.0;

            y[2] = cation_tot_e;
            m[2,1] = 1.0^2;
            m[2,2] = 1.0^1;
            m[2,3] = 1.0;

            y[3] = anion_tot_e;
            m[3,1] = (-1.0)^2;
            m[3,2] = (-1.0)^1;
            m[3,3] = 1.0;

            m_coeffs = -(m \ y);

            function fitted_chemical_potential(charge::Number)
                return 2*charge*m_coeffs[1] + m_coeffs[2];
            end

            function set_to_fitted_chemical_potential(file::ParsedFile)
                atomic_separation = file.atomic_separation;
                total_energy = file.total_energy;
                partial_charge_1 = file.partial_charge_1;
                partial_charge_2 = file.partial_charge_2;
                atomic_number_1 = file.atomic_number_1;
                atomic_number_2 = file.atomic_number_2;
                HOMO_energy = file.HOMO_energy;
                LUMO_energy = file.LUMO_energy;
                chemical_potential = fitted_chemical_potential(file.charge);
                charge = file.charge;

                return ParsedFile(atomic_separation,total_energy,
                    partial_charge_1,partial_charge_2,atomic_number_1,
                    atomic_number_2,HOMO_energy,LUMO_energy,chemical_potential,
                    charge);;
            end

            neutral_data[i] = set_to_fitted_chemical_potential(neutral_data[i]);
            cation_data[i] = set_to_fitted_chemical_potential(cation_data[i]);
            anion_data[i] = set_to_fitted_chemical_potential(anion_data[i]);
        end
    end

    neutral_data = interpolate_between(neutral_data,r0,r1,num_samples);
    cation_data = interpolate_between(cation_data,r0,r1,num_samples);
    anion_data = interpolate_between(anion_data,r0,r1,num_samples);
    set_chemical_potentials!(neutral_data,cation_data,anion_data);

    return neutral_data, cation_data, anion_data;
end

function read_all_sanitized_data(Z::Int, clip_data::Bool = false)
    return read_all_sanitized_data(Z,Z,clip_data);
end

function plot_reference_chemical_potential(Z1::Int, Z2::Int)
    neutral_data, cation_data, anion_data = read_all_sanitized_data(Z1,Z2,true);

    function get_data(data::Vector{ParsedFile})
        r = zeros(Float64,0);
        μ = zeros(Float64,0);

        for file_data in data
            r_i = file_data.atomic_separation;
            μ_i = file_data.chemical_potential;

            r = vcat(r,r_i);
            μ = vcat(μ,μ_i);
        end

        return r, μ;
    end

    neutral_r, neutral_μ = get_data(neutral_data);
    cation_r, cation_μ = get_data(cation_data);
    anion_r, anion_μ = get_data(anion_data);

    p = plot(neutral_r,neutral_μ,label="neutral",marker=true);
    plot!(cation_r,cation_μ,label="cation",marker=true);
    plot!(anion_r,anion_μ,label="anion",marker=true);
    plot!(xlabel="Atomic Separation [Borh]");
    plot!(ylabel="Chemical Potential [Hertree]");

    return p;
end

function polarize_to_model!(data::Vector{ParsedFile})
    n_threads = Threads.nthreads();
    @threads for thread_id in 1:n_threads
        simulation = make_system_from_parsed_file(data[1]);
        for i in thread_id:n_threads:length(data)
            set_diatomic_system_to_parsed_file!(simulation,data[i]);
            polarize_molecules!(simulation);

            mol = simulation.system.molecules[1];
            z1_eff = mol.atoms[1].valence_electrons;
            z2_eff = mol.atoms[2].valence_electrons;
            ζ1 = mol.atoms[1].polarization_coefficient;
            ζ2 = mol.atoms[2].polarization_coefficient;
            μ = simulation.system.chemical_potential;

            data[i].partial_charge_1 = z1_eff * (1.0 - ζ1);
            data[i].partial_charge_2 = z2_eff * (1.0 - ζ2);
            data[i].chemical_potential = μ;
        end
    end

    return;
end

function plot_reference_total_energy(Z1::Int, Z2::Int)
    neutral_data, cation_data, anion_data = read_all_sanitized_data(Z1,Z2,true);

    function get_data(data::Vector{ParsedFile})
        r = zeros(Float64,0);
        e = zeros(Float64,0);

        for file_data in data
            r_i = file_data.atomic_separation;
            e_i = file_data.total_energy / 2.0;

            r = vcat(r,r_i);
            e = vcat(e,e_i);
        end

        return r, e;
    end

    neutral_r, neutral_e = get_data(neutral_data);
    cation_r, cation_e = get_data(cation_data);
    anion_r, anion_e = get_data(anion_data);

    p = plot(neutral_r,neutral_e,label="neutral",marker=true);
    plot!(cation_r,cation_e,label="cation",marker=true);
    plot!(anion_r,anion_e,label="anion",marker=true);
    plot!(xlabel="Atomic Separation [Borh]");
    plot!(ylabel="Total Energy [Hertree]");

    return p;
end

function convert_to_type!(coeffs::Dict{Int,Number}, which_type)
    for key in keys(coeffs)
        coeffs[key] = convert(which_type,coeffs[key]);
    end

    return;
end

function convert_to_type!(coeffs::Dict{Tuple{Int,Int},Number}, which_type)
    for key in keys(coeffs)
        coeffs[key] = convert(which_type,coeffs[key]);
    end

    return;
end

function cast_coeffs_to_type!(
    coeffs::PolarizationEnergyCoefficients, which_type)
    convert_to_type!(coeffs.xc_coeffs.xc_a_1b,which_type);
    convert_to_type!(coeffs.xc_coeffs.xc_b_1b,which_type);
    convert_to_type!(coeffs.xc_coeffs.xc_c_1b,which_type);
    convert_to_type!(coeffs.xc_coeffs.xc_d_1b,which_type);
    convert_to_type!(coeffs.ke_coeffs.ke_e_1b,which_type);
    convert_to_type!(coeffs.ke_coeffs.ke_f_1b,which_type);

    convert_to_type!(coeffs.xc_coeffs.xc_a_2b,which_type);
    convert_to_type!(coeffs.xc_coeffs.xc_b_2b,which_type);
    convert_to_type!(coeffs.xc_coeffs.xc_c_2b,which_type);
    convert_to_type!(coeffs.xc_coeffs.xc_d_2b,which_type);

    return;
end

function cast_coeffs_to_type!(
    coeffs::TotalEnergyCoefficients, which_type)
    convert_to_type!(coeffs.xc_coeffs.xc_a_1b,which_type);
    convert_to_type!(coeffs.xc_coeffs.xc_b_1b,which_type);
    convert_to_type!(coeffs.xc_coeffs.xc_c_1b,which_type);
    convert_to_type!(coeffs.xc_coeffs.xc_d_1b,which_type);
    convert_to_type!(coeffs.ke_coeffs.ke_e_1b,which_type);
    convert_to_type!(coeffs.ke_coeffs.ke_f_1b,which_type);

    convert_to_type!(coeffs.xc_coeffs.xc_a_2b,which_type);
    convert_to_type!(coeffs.xc_coeffs.xc_b_2b,which_type);
    convert_to_type!(coeffs.xc_coeffs.xc_c_2b,which_type);
    convert_to_type!(coeffs.xc_coeffs.xc_d_2b,which_type);

    convert_to_type!(
        coeffs.non_polarizable_coeffs.depth,which_type);
    convert_to_type!(
        coeffs.non_polarizable_coeffs.stiffness_parameter,which_type);
    convert_to_type!(
        coeffs.non_polarizable_coeffs.equilibrium_distance,which_type);
    return;
end

function cast_tot_e_coeffs_to_type!(simulation::SimulationSystem, which_type)
    cast_coeffs_to_type!(simulation.tot_e_coeffs,which_type);
    return;
end

function cast_pol_e_coeffs_to_type!(simulation::SimulationSystem, which_type)
    cast_coeffs_to_type!(simulation.pol_e_coeffs,which_type);
    return;
end
