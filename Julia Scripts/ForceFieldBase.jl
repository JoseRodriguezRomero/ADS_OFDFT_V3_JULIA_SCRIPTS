using Printf;

mutable struct ElectronCloud
    basis_function_amplitude::Vector{Float64};
    basis_function_decay::Vector{Float64};
    thomas_fermi_ke::Float64;
    von_weizsacker_ke::Float64;
end

mutable struct Atom
    atomic_number::Int;
    coordinates::Vector;
    core_electrons::Int;
    valence_electrons::Int;
    polarization_coefficient::Number;
    electron_cloud_shells::Vector{ElectronCloud};
end

mutable struct Molecule
    name::String;
    atoms::Vector{Atom};
end

mutable struct ExponentialCoefficients
    amplitude::Number;
    decay::Number;
end

mutable struct XCCoeff2B
    m::Number;
    b::Number;
end

mutable struct EmpiricalXCCoefficients
    xc_a_1b::Dict{Int,Number};
    xc_b_1b::Dict{Int,Number};
    xc_c_1b::Dict{Int,Number};
    xc_d_1b::Dict{Int,Number};

    xc_a_2b::Dict{Tuple{Int,Int},XCCoeff2B};
    xc_b_2b::Dict{Tuple{Int,Int},XCCoeff2B};
    xc_c_2b::Dict{Tuple{Int,Int},XCCoeff2B};
    xc_d_2b::Dict{Tuple{Int,Int},XCCoeff2B};
end

mutable struct EmpiricalKECoefficients
    ke_e_1b::Dict{Int,Number};
    ke_f_1b::Dict{Int,Number};
end

mutable struct TotalEnergyEmpiricalStaticCoefficients
    core_core::Dict{Tuple{Int,Int},ExponentialCoefficients};
    core_valence::Dict{Tuple{Int,Int},ExponentialCoefficients};
end

mutable struct PolarizationEnergyEmpiricalStaticCoefficients
    core_valence::Dict{Tuple{Int,Int},ExponentialCoefficients};
end

mutable struct BasisSetSettings
    max_atomic_number::Int;
    reference_atoms::Dict{Int,Atom};
    atoms_μ0::Vector{Number};
end

mutable struct MolecularSystem
    name::String;
    molecules::Vector{Molecule};
    charge::Number;
    energy::Number;
    chemical_potential::Number;
end

mutable struct TotalEnergyCoefficients
    max_atomic_number::Int;
    tot_e_xc_coeffs::EmpiricalXCCoefficients;
    tot_e_ke_coeffs::EmpiricalKECoefficients;
    tot_e_static_coeffs::TotalEnergyEmpiricalStaticCoefficients;
end

mutable struct PolarizationEnergyCoefficients
    max_atomic_number::Int;
    pol_e_xc_coeffs::EmpiricalXCCoefficients;
    pol_e_ke_coeffs::EmpiricalKECoefficients;
    pol_e_static_coeffs::PolarizationEnergyEmpiricalStaticCoefficients;
end

mutable struct SimulationSystem
    system::MolecularSystem;
    tot_e_coeffs::TotalEnergyCoefficients;
    pol_e_coeffs::PolarizationEnergyCoefficients;
    basis_set_settings::BasisSetSettings;
end

function Base.copy(mol::Molecule)
    return Molecule(
        mol.name,
        copy(mol.atoms)
    );
end

function Base.copy(at::Atom)
    return Atom(
        at.atomic_number,
        at.coordinates,
        at.core_electrons,
        at.valence_electrons,
        at.polarization_coefficient,
        at.electron_cloud_shells
    );
end

function Base.display(at::Atom)
    print(@sprintf "%26s  " "Element:");
    print(@sprintf "%12s" get_element_symbol(at.atomic_number));
    print("\n");

    print(@sprintf "%26s  " "Atomic Number:");
    print(@sprintf "%12d" at.atomic_number);
    print("\n");

    print(@sprintf "%26s  " "Core Electrons:");
    print(@sprintf "%12d" at.core_electrons);
    print("\n");

    print(@sprintf "%26s  " "Valence Electrons:");
    print(@sprintf "%12d" at.valence_electrons);
    print("\n");

    print(@sprintf "%26s  " "Electron Shells:");
    print(@sprintf "%12d" length(at.electron_cloud_shells));
    print("\n");

    print(@sprintf "%26s  " "Polarization Coefficient:");
    print(@sprintf "%12.6f " at.polarization_coefficient);
    print("\n");

    print(@sprintf "%26s  " "Coordinates:");
    print(@sprintf "%12.6f " at.coordinates[1]);
    print(@sprintf "%12.6f " at.coordinates[2]);
    print(@sprintf "%12.6f " at.coordinates[3]);
    print("\n");

    print("\n");
end

function Base.display(mol::Molecule)
    function center(s::AbstractString, width::Int)
        pad = width - length(s)
        pad <= 0 && return s
        lpad(rpad(s, length(s) + pad ÷ 2), width)
    end

    function print_delimeter(width::Int)
        for i in 1:width
            print("-");
        end
        print(" ");
    end

    print("\n");

    print_delimeter(14);
    print_delimeter(14);
    print_delimeter(14);
    print_delimeter(14);
    print_delimeter(20);
    print_delimeter(62);
    print("\n");

    print(center("",14)*" ");
    print(center("Atomic",14)*" ");
    print(center("Core",14)*" ");
    print(center("Valence",14)*" ");
    print(center("Polarization",20)*" ");
    print(center("Coordinates [Bohr]",62)*" ");
    print("\n");
    
    print(center("Element",14)*" ");
    print(center("Number",14)*" ");
    print(center("Electrons",14)*" ");
    print(center("Electrons",14)*" ");
    print(center("Coefficient",20)*" ");
    print(center("x",20)*" ");
    print(center("y",20)*" ");
    print(center("z",20)*" ");
    print("\n");

    print_delimeter(14);
    print_delimeter(14);
    print_delimeter(14);
    print_delimeter(14);
    print_delimeter(20);
    print_delimeter(20);
    print_delimeter(20);
    print_delimeter(20);
    print("\n");

    for atom in mol.atoms
        print(@sprintf "%14s " get_element_symbol(atom.atomic_number));
        print(@sprintf "%14d " atom.atomic_number);
        print(@sprintf "%14d " atom.core_electrons);
        print(@sprintf "%14d " atom.valence_electrons);
        print(@sprintf "%20.6f " atom.polarization_coefficient);
        print(@sprintf "%20.6f " atom.coordinates[1]);
        print(@sprintf "%20.6f " atom.coordinates[2]);
        print(@sprintf "%20.6f " atom.coordinates[3]);
        print("\n");
    end

    print_delimeter(14);
    print_delimeter(14);
    print_delimeter(14);
    print_delimeter(14);
    print_delimeter(20);
    print_delimeter(20);
    print_delimeter(20);
    print_delimeter(20);
    print("\n");
end

function Molecule()
    name = "Kryptonite";
    atoms = Vector{Atom}();
    return Molecule(name,atoms);
end

function MolecularSystem()
    name = "Kryptonite++";
    molecules = Vector{Molecule}();
    charge = 0;
    energy = 0;
    chemical_potential = 0;

    return MolecularSystem(name, molecules, charge, 
        energy, chemical_potential);
end

function number_of_atoms(molecule::Molecule)
    # Returns the number of atoms in the molecule.
    return length(molecule.atoms);
end

function atom_atomic_number(molecule::Molecule, i::Int)
    # Returns the atomic number of the iᵗʰ atom.
    return molecule.atoms[i].atomic_number;
end

function atom_core_electrons(molecule::Molecule, i::Int)
    # Returns the number of core electrons.
    return molecule.atoms[i].core_electrons;
end

function atom_valence_electrons(molecule::Molecule, i::Int)
    # Returns the number of valence electrons.
    return molecule.atoms[i].valence_electrons;
end

function atom_position(molecule::Molecule, i::Int)
    # Returns the position vector of the iᵗʰ atom.
    return molecule.atoms[i].position;
end

function atom_polarization_coeff(molecule::Molecule, i::Int)
    # Returns the polarization coefficient of the iᵗʰ atom.
    return molecule.atoms[i].polarization_coefficient;
end

function set_atom_pol_coeff!(molecule::Molecule, ζ::Number, i::Int)
    # Sets the atom iᵗʰ atom's polarization coefficient and its corresponding 
    # electron cloud.
    molecule.atoms[i].polarization_coefficient = ζ;
end

function set_atom_partial_charge!(molecule::Molecule, ρ::Number, i::Int)
    # Sets the atom iᵗʰ atom's polarization coefficient and its corresponding 
    # electron cloud so that they match the desired partial charge.
    z_eff = atom_valence_electrons(molecule,i);
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

function set_atom_coordinates!(molecule::Molecule, r::Vector, i::Int)
    # Sets the coordinates of the iᵗʰ atom and its corresponding electron cloud.
    molecule.atoms[i].coordinates[:] = r[:];
end

function erf_arg_cutoff()
    return 20.0;
end

function atoms_dist_cutoff()
    return 1.0E-6;
end

function get_element_symbol(atomic_number::Int)
    symbols = [
        "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  
        "Ne", "Na", "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar"
    ];

    if 1 <= atomic_number <= 18
        return symbols[atomic_number]
    else
        return "Invalid atomic number. Please enter a number from 1 to 18."
    end
end

function get_atomic_number(symbol::String)
    symbols = Dict(
         "H" =>  1, "He" => 2, "Li" => 3, "Be" => 4,   "B" => 5,   "C" => 6,
         "N" =>  7,  "O" => 8,  "F" => 9, "Ne" => 10, "Na" => 11, "Mg" => 12,
        "Al" => 13, "Si" => 14, "P" => 15, "S" => 16, "Cl" => 17, "Ar" => 18
    );

    symbol = strip(symbol);
    if haskey(symbols, symbol)
        return symbols[symbol]
    else
        return "Invalid element symbol. Please enter a symbol from H to Ar."
    end
end
