include("ForceFieldBaseECP.jl")

struct ParsedFile
    atomic_separation::Float64;
    total_energy::Float64;
    partial_charge_1::Float64;
    partial_charge_2::Float64;
    atomic_number_1::Int;
    atomic_number_2::Int;
    HOMO_energy::Float64;
    LUMO_energy::Float64;
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
        charge = round(Int, partial_charge_1 + partial_charge_2);

        angstrom_to_bohr = 1.88973;
        atomic_separation *= angstrom_to_bohr;

        parsed_file = ParsedFile(atomic_separation, total_energy, 
            partial_charge_1, partial_charge_2, atomic_number_1,
            atomic_number_2, HOMO_energy, LUMO_energy, charge);
        
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

    system = make_diatomic_system(Z1,Z2,d,charge);
    return system;
end

function set_diatomic_system_to_parsed_file!(simulation::SimulationSystem,
    parsed_file::ParsedFile)
    # Sets a simulation for a molecule using the data in parsed_file.
    Z1 = parsed_file.atomic_number_1;
    Z2 = parsed_file.atomic_number_2;
    d = parsed_file.atomic_separation;
    charge = parsed_file.charge;

    element_cloud_coeffs = simulation.basis_set_settings.element_cloud_coeffs;
    molecule = Molecule();

    molecule.atoms_data = zeros(Float64,2,6);
    molecule.cloud_data = zeros(Float64,6,6);

    z1_eff = round(sum(element_cloud_coeffs[Z1,1:2:end]));
    z2_eff = round(sum(element_cloud_coeffs[Z2,1:2:end]));

    charge1 = charge*z1_eff/(z1_eff+z2_eff);
    charge2 = charge*z2_eff/(z1_eff+z2_eff);

    molecule.atoms_data[1,5] = Z1;
    molecule.atoms_data[1,4] = z1_eff;
    molecule.atoms_data[1,6] = (z1_eff - charge1) / z1_eff;

    molecule.atoms_data[2,5] = Z2;
    molecule.atoms_data[2,4] = z2_eff;
    molecule.atoms_data[2,6] = (z2_eff - charge2) / z2_eff;

    molecule.cloud_data[1:3,4] .= element_cloud_coeffs[Z1,1:2:end];
    molecule.cloud_data[1:3,5] .= element_cloud_coeffs[Z1, 2:2:end];
    molecule.cloud_data[1:3,6] .= (z1_eff - charge1) / z1_eff;

    molecule.cloud_data[4:6,4] .= element_cloud_coeffs[Z2,1:2:end];
    molecule.cloud_data[4:6,5] .= element_cloud_coeffs[Z2, 2:2:end];
    molecule.cloud_data[4:6,6] .= (z2_eff - charge2) / z2_eff;

    molecule.atoms_data[1,3] = +d/2.0;
    molecule.atoms_data[2,3] = -d/2.0;

    molecule.cloud_data[1:3,3] .= +d/2.0;
    molecule.cloud_data[4:6,3] .= -d/2.0;

    chemical_potential = 
        (parsed_file.HOMO_energy + parsed_file.LUMO_energy)/2.0;

    simulation.system.molecules = [molecule];
    simulation.system.chemical_potential = chemical_potential;
end

function set_fitted_coeffs!(coeffs::EmpiricalCoefficients,
    mappings::EmpiricalCoefficientMappings, Z::Int, fitted_coeffs::Vector)
    num_1b_coeffs = 4;
    num_2b_coeffs = 4;

    # Separate aux_X between 1B and 2B coefficient components.
    i0_1B = 1;
    i1_1B = i0_1B + 2*num_1b_coeffs - 1;

    i0_2B = i1_1B + 1;
    i1_2B = i0_2B + 2*num_2b_coeffs - 1;

    aux_X_1B = fitted_coeffs[i0_1B:i1_1B];
    aux_X_2B = fitted_coeffs[i0_2B:i1_2B];

    # Store the 1B fitted coefficients
    xc_a_1b_b_index = 1;
    xc_a_1b_m_index = 2;
    xc_c_1b_b_index = 3;
    xc_c_1b_m_index = 4;
    xc_e_1b_b_index = 5;
    xc_e_1b_m_index = 6;
    xc_f_1b_b_index = 7;
    xc_f_1b_m_index = 8;

    coeff_1b_map = mappings.coeff_1b_map[Z];
    coeffs.xc_a_1b[coeff_1b_map,1] = aux_X_1B[xc_a_1b_b_index];
    coeffs.xc_a_1b[coeff_1b_map,2] = aux_X_1B[xc_a_1b_m_index];
    coeffs.xc_c_1b[coeff_1b_map,1] = aux_X_1B[xc_c_1b_b_index];
    coeffs.xc_c_1b[coeff_1b_map,2] = aux_X_1B[xc_c_1b_m_index];
    coeffs.xc_e_1b[coeff_1b_map,1] = aux_X_1B[xc_e_1b_b_index];
    coeffs.xc_e_1b[coeff_1b_map,2] = aux_X_1B[xc_e_1b_m_index];
    coeffs.xc_f_1b[coeff_1b_map,1] = aux_X_1B[xc_f_1b_b_index];
    coeffs.xc_f_1b[coeff_1b_map,2] = aux_X_1B[xc_f_1b_m_index];

    # Store the 2B fitted coefficients
    xc_a_2b_b_index = 1;
    xc_a_2b_m_index = 2;
    xc_b_2b_b_index = 3;
    xc_b_2b_m_index = 4;
    xc_c_2b_b_index = 5;
    xc_c_2b_m_index = 6;
    xc_d_2b_b_index = 7;
    xc_d_2b_m_index = 8;

    coeff_2b_map = mappings.coeff_2b_map[(Z,Z)];
    coeffs.xc_a_2b[coeff_2b_map,1] = aux_X_2B[xc_a_2b_b_index];
    coeffs.xc_a_2b[coeff_2b_map,2] = aux_X_2B[xc_a_2b_m_index];
    coeffs.xc_a_2b[coeff_2b_map,3] = aux_X_2B[xc_a_2b_m_index];
    coeffs.xc_b_2b[coeff_2b_map,1] = aux_X_2B[xc_b_2b_b_index];
    coeffs.xc_b_2b[coeff_2b_map,2] = aux_X_2B[xc_b_2b_m_index];
    coeffs.xc_b_2b[coeff_2b_map,3] = aux_X_2B[xc_b_2b_m_index];
    coeffs.xc_c_2b[coeff_2b_map,1] = aux_X_2B[xc_c_2b_b_index];
    coeffs.xc_c_2b[coeff_2b_map,2] = aux_X_2B[xc_c_2b_m_index];
    coeffs.xc_c_2b[coeff_2b_map,3] = aux_X_2B[xc_c_2b_m_index];
    coeffs.xc_d_2b[coeff_2b_map,1] = aux_X_2B[xc_d_2b_b_index];
    coeffs.xc_d_2b[coeff_2b_map,2] = aux_X_2B[xc_d_2b_m_index];
    coeffs.xc_d_2b[coeff_2b_map,3] = aux_X_2B[xc_d_2b_m_index];
end

function set_fitted_pol_e_coeffs!(simulation::SimulationSystem,
    Z::Int, fitted_coeffs::Vector)
    set_fitted_coeffs!(simulation.pol_e_xc_coeffs,
        simulation.coeff_mappings,Z,fitted_coeffs);
end

function set_fitted_tot_e_coeffs!(simulation::SimulationSystem,
    Z::Int, fitted_coeffs::Vector)
    set_fitted_coeffs!(simulation.tot_e_xc_coeffs,
        simulation.coeff_mappings,Z,fitted_coeffs);
end

function set_fitted_coeffs!(coeffs::EmpiricalCoefficients,
    mappings::EmpiricalCoefficientMappings, Z1::Int, Z2::Int, 
    fitted_coeffs::Vector)
    # Store the 2B fitted coefficients
    xc_a_2b_b_index  = 1;
    xc_a_2b_m1_index = 2;
    xc_a_2b_m2_index = 3;
    xc_b_2b_b_index  = 4;
    xc_b_2b_m1_index = 5;
    xc_b_2b_m2_index = 6;
    xc_c_2b_b_index  = 7;
    xc_c_2b_m1_index = 8;
    xc_c_2b_m2_index = 9;
    xc_d_2b_b_index  = 10;
    xc_d_2b_m1_index = 11;
    xc_d_2b_m2_index = 12;

    coeff_2b_map = mappings.coeff_2b_map[(Z1,Z2)];
    coeffs.xc_a_2b[coeff_2b_map,1] = fitted_coeffs[xc_a_2b_b_index];
    coeffs.xc_a_2b[coeff_2b_map,2] = fitted_coeffs[xc_a_2b_m1_index];
    coeffs.xc_a_2b[coeff_2b_map,3] = fitted_coeffs[xc_a_2b_m2_index];
    coeffs.xc_b_2b[coeff_2b_map,1] = fitted_coeffs[xc_b_2b_b_index];
    coeffs.xc_b_2b[coeff_2b_map,2] = fitted_coeffs[xc_b_2b_m1_index];
    coeffs.xc_b_2b[coeff_2b_map,3] = fitted_coeffs[xc_b_2b_m2_index];
    coeffs.xc_c_2b[coeff_2b_map,1] = fitted_coeffs[xc_c_2b_b_index];
    coeffs.xc_c_2b[coeff_2b_map,2] = fitted_coeffs[xc_c_2b_m1_index];
    coeffs.xc_c_2b[coeff_2b_map,3] = fitted_coeffs[xc_c_2b_m2_index];
    coeffs.xc_d_2b[coeff_2b_map,1] = fitted_coeffs[xc_d_2b_b_index];
    coeffs.xc_d_2b[coeff_2b_map,2] = fitted_coeffs[xc_d_2b_m1_index];
    coeffs.xc_d_2b[coeff_2b_map,3] = fitted_coeffs[xc_d_2b_m2_index];
end

function set_fitted_pol_e_coeffs!(simulation::SimulationSystem,
    Z1::Int, Z2::Int, fitted_coeffs::Vector)
    set_fitted_coeffs!(simulation.pol_e_xc_coeffs,
        simulation.coeff_mappings,Z1,Z2,fitted_coeffs);
end

function set_fitted_tot_e_coeffs!(simulation::SimulationSystem,
    Z1::Int, Z2::Int, fitted_coeffs::Vector)
    set_fitted_coeffs!(simulation.tot_e_xc_coeffs,
        simulation.coeff_mappings,Z1,Z2,fitted_coeffs);
end

function get_reference_atom_chemical_potential()
    base_dir = "../Reference Data/atoms_data/";
    neutral_file_ID = base_dir*"Neutral/elem_orb_energies.txt";
    cation_file_ID = base_dir*"Cation/elem_orb_energies.txt";
    anion_file_ID = base_dir*"Anion/elem_orb_energies.txt";

    function ReadFile(file_ID::String)
        lines = readlines(file_ID);
        energies = zeros(Float64,length(lines));

        for i in eachindex(lines)
            line_splitted = split(lines[i]);
            homo = parse(Float64,line_splitted[1]);
            lumo = parse(Float64,line_splitted[2]);

            if isnan(homo)
                energies[i] = lumo;
            else
                energies[i] = (homo+lumo)/2.0;
            end
        end

        return energies;
    end

    neutral_chem_μ = ReadFile(neutral_file_ID);
    cation_chem_μ = ReadFile(cation_file_ID);
    anion_chem_μ = ReadFile(anion_file_ID);

    return neutral_chem_μ, cation_chem_μ, anion_chem_μ;
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

function read_all_sanitized_data(Z1::Int, Z2::Int)
    neutral_data, cation_data, anion_data = read_all_scanned_data(Z1,Z2);

    neutral_data = sanitize_data(neutral_data);
    cation_data = sanitize_data(cation_data);
    anion_data = sanitize_data(anion_data);

    return neutral_data, cation_data, anion_data;
end

function read_all_sanitized_data(Z::Int)
    return read_all_sanitized_data(Z,Z);
end

function cast_coeffs_to_type!(coeffs::EmpiricalCoefficients, which_type)
    coeffs.xc_a_1b = which_type.(coeffs.xc_a_1b);
    coeffs.xc_c_1b = which_type.(coeffs.xc_c_1b);
    coeffs.xc_e_1b = which_type.(coeffs.xc_e_1b);
    coeffs.xc_f_1b = which_type.(coeffs.xc_f_1b);

    coeffs.xc_a_2b = which_type.(coeffs.xc_a_2b);
    coeffs.xc_b_2b = which_type.(coeffs.xc_b_2b);
    coeffs.xc_c_2b = which_type.(coeffs.xc_c_2b);
    coeffs.xc_d_2b = which_type.(coeffs.xc_d_2b);
end

function cast_pol_e_coeffs_to_type!(simulation::SimulationSystem, which_type)
    cast_coeffs_to_type!(simulation.pol_e_xc_coeffs,which_type);
end

function cast_tot_e_coeffs_to_type!(simulation::SimulationSystem, which_type)
    cast_coeffs_to_type!(simulation.tot_e_xc_coeffs,which_type);
end
