include("ForceFieldBase.jl")
include("ReadFittedCoeffs.jl")

struct ParsedFile
    atomic_separation::Float64;
    total_energy::Float64;
    partial_charge_1::Float64;
    partial_charge_2::Float64;
    atomic_number_1::Int32;
    atomic_number_2::Int32;
    HOMO_energy::Float64;
    LUMO_energy::Float64;
    charge::Int64;
end

function getElementSymbol(atomic_number::Int)
    symbols = [
        "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne",
        "Na", "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar"
    ]

    if 1 <= atomic_number <= 18
        return symbols[atomic_number]
    else
        return "Invalid atomic number. Please enter a number from 1 to 18."
    end
end

function getAtomicNumber(symbol::String)
    symbols = Dict(
        "H" => 1,   "He" => 2,  "Li" => 3,  "Be" => 4,  "B" => 5,   "C" => 6,
        "N" => 7,   "O" => 8,   "F" => 9,   "Ne" => 10, "Na" => 11, "Mg" => 12,
        "Al" => 13, "Si" => 14, "P" => 15,  "S" => 16,  "Cl" => 17, "Ar" => 18
    )

    symbol = strip(symbol);
    if haskey(symbols, symbol)
        return symbols[symbol]
    else
        return "Invalid element symbol. Please enter a symbol from H to Ar."
    end
end

function readScannedData(file_name::String)
    file_data = Vector{ParsedFile}();
    lines = readlines(file_name);

    atom_elem_label_1 = split(lines[2])[3];
    atom_elem_label_2 = split(lines[2])[5];

    atom_elem_label_1 = String(atom_elem_label_1[2:end]);
    atom_elem_label_2 = String(atom_elem_label_2[2:end]);

    atomic_number_1 = getAtomicNumber(atom_elem_label_1);
    atomic_number_2 = getAtomicNumber(atom_elem_label_2);

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

        parsed_file = ParsedFile(atomic_separation, total_energy, 
            partial_charge_1, partial_charge_2, atomic_number_1,
            atomic_number_2, HOMO_energy, LUMO_energy, charge);
        
        push!(file_data, parsed_file);
    end

    return file_data;
end

function readAllScannedData(Z1::Int, Z2::Int)
    neutral_file_name = "../Reference Data/diatomic_data/Neutral/ParsedData/";
    cation_file_name = "../Reference Data/diatomic_data/Cation/ParsedData/";
    anion_file_name = "../Reference Data/diatomic_data/Anion/ParsedData/";

    z1_elem_symbol = getElementSymbol(Z1);
    z2_elem_symbol = getElementSymbol(Z2);
    parsed_file_name = "scan_"*z1_elem_symbol*z2_elem_symbol*".txt";

    neutral_file_name *= parsed_file_name;
    cation_file_name *= parsed_file_name;
    anion_file_name *= parsed_file_name;

    neutral_data = readScannedData(neutral_file_name);
    cation_data = readScannedData(cation_file_name);
    anion_data = readScannedData(anion_file_name);

    return neutral_data, cation_data, anion_data;
end

function readAllScannedData(Z::Int)
    return readAllScannedData(Z,Z);
end

function GetHOMOEnergies(data::Vector{ParsedFile})
    r = zeros(Float64,length(data));
    e = zeros(Float64,length(data));

    for i in eachindex(data)
        r[i] = data[i].atomic_separation;
        e[i] = data[i].HOMO_energy;
    end

    return r, e;
end

function GetLUMOEnergies(data::Vector{ParsedFile})
    r = zeros(Float64,length(data));
    e = zeros(Float64,length(data));

    for i in eachindex(data)
        r[i] = data[i].atomic_separation;
        e[i] = data[i].LUMO_energy;
    end

    return r, e;
end

function GetTotalEnergies(data::Vector{ParsedFile})
    r = zeros(Float64,length(data));
    e = zeros(Float64,length(data));

    for i in eachindex(data)
        r[i] = data[i].atomic_separation;
        e[i] = data[i].total_energy;
    end

    return r, e;
end

function MakeMoleculeFromParsedFile(parsed_file::ParsedFile)
    # Creates a molecule structure using the data in parsed_file.
    atom1 = MakeAtom(parsed_file.atomic_number_1,0);
    atom2 = MakeAtom(parsed_file.atomic_number_2,0);

    angstrom_to_bohr = 1.88973;
    atom2.atoms_data[3] = angstrom_to_bohr * parsed_file.atomic_separation;
    atom2.cloud_data[:,3] .= angstrom_to_bohr * parsed_file.atomic_separation;

    mol =  Molecule();
    mol.charge = parsed_file.charge;
    mol.energy = parsed_file.total_energy;
    mol.chem_μ = (parsed_file.HOMO_energy + parsed_file.LUMO_energy)/2.0;

    mol.cloud_data = vcat(atom1.cloud_data, atom2.cloud_data);
    mol.atoms_data = vcat(atom1.atoms_data, atom2.atoms_data);

    Zeff_1 = atom1.atoms_data[1,4];
    ζ1 = (Zeff_1 - parsed_file.partial_charge_1) / Zeff_1;

    Zeff_2 = atom2.atoms_data[1,4];
    ζ2 = (Zeff_2 - parsed_file.partial_charge_2) / Zeff_2;

    mol.cloud_data[1:3,6] .= ζ1;
    mol.cloud_data[4:6,6] .= ζ2;

    return mol;
end

function SetFittedCoeffs!(all_coeffs::Vector, Z::Integer, fitted_coeffs::Vector)
    # Separate aux_X between 1B and 2B coefficient components.
    i0_1B = 1;
    i1_1B = i0_1B + 2*num_1b_coeffs - 1;

    i0_2B = i1_1B + 1;
    i1_2B = i0_2B + 2*num_2b_coeffs - 1;

    aux_X_1B = fitted_coeffs[i0_1B:i1_1B];
    aux_X_2B = fitted_coeffs[i0_2B:i1_2B];

    # Store the 1B fitted coefficients
    for i in 1:num_1b_coeffs
        aux_ind = collect(1:length(all_coeffs[1][:,1]));
        aux_ind = aux_ind[i:num_1b_coeffs:end];
        aux_ind = aux_ind[Z];

        all_coeffs[1][aux_ind,1] = aux_X_1B[(i-1)*2+1];
        all_coeffs[1][aux_ind,2] = aux_X_1B[(i-1)*2+2];
    end

    # Store the 2B fitted coefficients
    for i in 1:num_2b_coeffs
        aux_ind = collect(1:length(all_coeffs[2][:,1]));
        aux_ind = aux_ind[i:num_2b_coeffs:end];
        map_key = (Z,Z);
        aux_ind = aux_ind[coeff_map_2b[(map_key)]];

        all_coeffs[2][aux_ind,1] = aux_X_2B[(i-1)*2+1];
        all_coeffs[2][aux_ind,2] = aux_X_2B[(i-1)*2+2];
        all_coeffs[2][aux_ind,3] = aux_X_2B[(i-1)*2+2];
    end
end

function SetFittedCoeffs!(all_coeffs::Vector, Z1::Integer, Z2::Integer, 
    fitted_coeffs::Vector)
    # Store the 2B fitted coefficients
    for i in 1:num_2b_coeffs
        aux_ind = collect(1:length(all_coeffs[2][:,1]));
        aux_ind = aux_ind[i:num_2b_coeffs:end];
        map_key = (Z1,Z2);
        aux_ind = aux_ind[coeff_map_2b[(map_key)]];

        all_coeffs[2][aux_ind,1] = fitted_coeffs[(i-1)*3+1];
        all_coeffs[2][aux_ind,2] = fitted_coeffs[(i-1)*3+2];
        all_coeffs[2][aux_ind,3] = fitted_coeffs[(i-1)*3+3];
    end
end

function GetAtomChemμEnergies()
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

function GetAtomEnergies()
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

function SanitizeData(scanned_data::Vector{ParsedFile})
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

function readAllSanitizedData(Z1::Integer, Z2::Integer)
    neutral_data, cation_data, anion_data = readAllScannedData(Z1,Z2);

    neutral_data = SanitizeData(neutral_data);
    cation_data = SanitizeData(cation_data);
    anion_data = SanitizeData(anion_data);

    return neutral_data, cation_data, anion_data;
end

function readAllSanitizedData(Z::Integer)
    return readAllSanitizedData(Z,Z);
end

neutral_energies, cation_energies, anion_energies = GetAtomEnergies();
neutral_chem_μ, cation_chem_μ, anion_chem_μ = GetAtomChemμEnergies();

function SetAtomEnergy!(atom::Molecule)
    # Sets the energy of the atom to its B3LYP/aug-cc-pVTZ value.
    global neutral_energies, cation_energies, anion_energies;
    global neutral_chem_μ, cation_chem_μ, anion_chem_μ;

    atomic_num = round(Int,atom.atoms_data[1,5]);
    if atom.charge == 0
        atom.energy = neutral_energies[atomic_num];
        atom.chem_μ = neutral_chem_μ[atomic_num];
    elseif atom.charge == 1
        atom.energy = cation_energies[atomic_num];
        atom.chem_μ = cation_chem_μ[atomic_num];
    else
        atom.energy = anion_energies[atomic_num];
        atom.chem_μ = anion_chem_μ[atomic_num];
    end
end
