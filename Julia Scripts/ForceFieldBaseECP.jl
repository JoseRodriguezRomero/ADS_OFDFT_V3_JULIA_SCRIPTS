using Printf;

include("ForceFieldBase.jl")

function get_element_symbol(atomic_number::Int)
    symbols = [
        "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  
        "Ne", "Na", "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar"
    ]

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
    )

    symbol = strip(symbol);
    if haskey(symbols, symbol)
        return symbols[symbol]
    else
        return "Invalid element symbol. Please enter a symbol from H to Ar."
    end
end

function save_1b_coeffs(coeffs::EmpiricalCoefficients, 
    mappings::EmpiricalCoefficientMappings, type::String)
    base_dir = "ECP XC Coeffs/"*type*"/One Body/";
    xc_coeffs_a_id = base_dir*"xc_coeffs_A.txt";
    xc_coeffs_c_id = base_dir*"xc_coeffs_C.txt";
    xc_coeffs_e_id = base_dir*"xc_coeffs_E.txt";
    xc_coeffs_f_id = base_dir*"xc_coeffs_F.txt";

    xc_a_1b = coeffs.xc_a_1b;
    xc_c_1b = coeffs.xc_c_1b;
    xc_e_1b = coeffs.xc_e_1b;
    xc_f_1b = coeffs.xc_f_1b;

    mappings_1b = mappings.coeff_1b_map;
    max_atomic_number = coeffs.max_atomic_number;

    function save_to_file(file_name::String, coeffs::Matrix{Number})
        file_id = open(file_name,"w");
        for i in 1:max_atomic_number
            write(file_id, @sprintf "%18.8E " coeffs[mappings_1b[i],1]);
            write(file_id, @sprintf "%18.8E " coeffs[mappings_1b[i],2]);
            write(file_id, @sprintf "%5s " get_element_symbol(i));
            write(file_id, @sprintf "\n");
        end

        close(file_id);
    end

    save_to_file(xc_coeffs_a_id,xc_a_1b);
    save_to_file(xc_coeffs_c_id,xc_c_1b);
    save_to_file(xc_coeffs_e_id,xc_e_1b);
    save_to_file(xc_coeffs_f_id,xc_f_1b);
end

function save_1b_pol_e_coeffs(simulation::SimulationSystem)
    save_1b_coeffs(simulation.pol_e_xc_coeffs,
        simulation.coeff_mappings,"Polarization");
end

function save_1b_tot_e_coeffs(simulation::SimulationSystem)
    save_1b_coeffs(simulation.tot_e_xc_coeffs,
        simulation.coeff_mappings,"Energy");
end

function save_2b_coeffs(coeffs::EmpiricalCoefficients,
    mappings::EmpiricalCoefficientMappings, type::String)
    base_dir = "ECP XC Coeffs/"*type*"/Two Body/"
    xc_coeffs_a_id = base_dir*"xc_coeffs_A.txt";
    xc_coeffs_b_id = base_dir*"xc_coeffs_B.txt";
    xc_coeffs_c_id = base_dir*"xc_coeffs_C.txt";
    xc_coeffs_d_id = base_dir*"xc_coeffs_D.txt";

    xc_a_2b = coeffs.xc_a_2b;
    xc_b_2b = coeffs.xc_b_2b;
    xc_c_2b = coeffs.xc_c_2b;
    xc_d_2b = coeffs.xc_d_2b;

    mappings_2b = mappings.coeff_2b_map;
    max_atomic_number = coeffs.max_atomic_number;

    function save_to_file(file_name::String, coeffs::Matrix{Number})
        file_id = open(file_name,"w");
        for i in 1:max_atomic_number
            for j in i:max_atomic_number
                ij = mappings_2b[(i,j)];
                write(file_id, @sprintf "%18.8E " coeffs[ij,1]);
                write(file_id, @sprintf "%18.8E " coeffs[ij,2]);
                write(file_id, @sprintf "%18.8E " coeffs[ij,3]);
                write(file_id, @sprintf "%5s " get_element_symbol(i));
                write(file_id, @sprintf "%5s " get_element_symbol(j));
                write(file_id, @sprintf "\n");
            end
        end

        close(file_id);
    end

    save_to_file(xc_coeffs_a_id,xc_a_2b);
    save_to_file(xc_coeffs_b_id,xc_b_2b);
    save_to_file(xc_coeffs_c_id,xc_c_2b);
    save_to_file(xc_coeffs_d_id,xc_d_2b);
end

function save_2b_pol_e_coeffs(simulation::SimulationSystem)
    save_2b_coeffs(simulation.pol_e_xc_coeffs,
        simulation.coeff_mappings,"Polarization");
end

function save_2b_tot_e_coeffs(simulation::SimulationSystem)
    save_2b_coeffs(simulation.tot_e_xc_coeffs,
        simulation.coeff_mappings,"Energy");
end

function save_fitted_coeffs(simulation::SimulationSystem)
    save_1b_pol_e_coeffs(simulation);
    save_2b_pol_e_coeffs(simulation);

    save_1b_tot_e_coeffs(simulation);
    save_2b_tot_e_coeffs(simulation);
end

function load_1b_coeffs(type::String)
    base_dir = "ECP XC Coeffs/"*type*"/One Body/";
    xc_coeffs_a_id = base_dir*"xc_coeffs_A.txt";
    xc_coeffs_c_id = base_dir*"xc_coeffs_C.txt";
    xc_coeffs_e_id = base_dir*"xc_coeffs_E.txt";
    xc_coeffs_f_id = base_dir*"xc_coeffs_F.txt";

    function read_file_data(file_name::String)
        lines = readlines(file_name);
        coeffs = zeros(Float64,length(lines),2);
        for i in eachindex(lines)
            line_splitted = split(lines[i]);
            coeffs[i,1] = parse(Float64,line_splitted[1]);
            coeffs[i,2] = parse(Float64,line_splitted[2]);
        end

        return coeffs;
    end

    xc_a_1b = read_file_data(xc_coeffs_a_id);
    xc_c_1b = read_file_data(xc_coeffs_c_id);
    xc_e_1b = read_file_data(xc_coeffs_e_id);
    xc_f_1b = read_file_data(xc_coeffs_f_id);

    return xc_a_1b, xc_c_1b, xc_e_1b, xc_f_1b;
end

function load_1b_pol_e_coeffs!(simulation::SimulationSystem)
    xc_a_1b, xc_c_1b, xc_e_1b, xc_f_1b = load_1b_coeffs("Polarization");
    
    simulation.pol_e_xc_coeffs.xc_a_1b = xc_a_1b;
    simulation.pol_e_xc_coeffs.xc_c_1b = xc_c_1b;
    simulation.pol_e_xc_coeffs.xc_e_1b = xc_e_1b;
    simulation.pol_e_xc_coeffs.xc_f_1b = xc_f_1b;
end

function load_1b_tot_e_coeffs!(simulation::SimulationSystem)
    xc_a_1b, xc_c_1b, xc_e_1b, xc_f_1b = load_1b_coeffs("Energy");

    simulation.tot_e_xc_coeffs.xc_a_1b = xc_a_1b;
    simulation.tot_e_xc_coeffs.xc_c_1b = xc_c_1b;
    simulation.tot_e_xc_coeffs.xc_e_1b = xc_e_1b;
    simulation.tot_e_xc_coeffs.xc_f_1b = xc_f_1b;
end

function load_2b_coeffs(type::String)
    base_dir = "ECP XC Coeffs/"*type*"/Two Body/";
    xc_coeffs_a_id = base_dir*"xc_coeffs_A.txt";
    xc_coeffs_b_id = base_dir*"xc_coeffs_B.txt";
    xc_coeffs_c_id = base_dir*"xc_coeffs_C.txt";
    xc_coeffs_d_id = base_dir*"xc_coeffs_D.txt";

    function read_file_data(file_name::String)
        lines = readlines(file_name);
        coeffs = zeros(Float64,length(lines),3);
        for i in eachindex(lines)
            line_splitted = split(lines[i]);
            coeffs[i,1] = parse(Float64,line_splitted[1]);
            coeffs[i,2] = parse(Float64,line_splitted[2]);
            coeffs[i,3] = parse(Float64,line_splitted[3]);
        end

        return coeffs;
    end

    xc_a_2b = read_file_data(xc_coeffs_a_id);
    xc_b_2b = read_file_data(xc_coeffs_b_id);
    xc_c_2b = read_file_data(xc_coeffs_c_id);
    xc_d_2b = read_file_data(xc_coeffs_d_id);

    return xc_a_2b, xc_b_2b, xc_c_2b, xc_d_2b;
end

function load_2b_pol_e_coeffs!(simulation::SimulationSystem)
    xc_a_2b, xc_b_2b, xc_c_2b, xc_d_2b = load_2b_coeffs("Polarization");

    simulation.pol_e_xc_coeffs.xc_a_2b = xc_a_2b;
    simulation.pol_e_xc_coeffs.xc_b_2b = xc_b_2b;
    simulation.pol_e_xc_coeffs.xc_c_2b = xc_c_2b;
    simulation.pol_e_xc_coeffs.xc_d_2b = xc_d_2b;
end

function load_2b_tot_e_coeffs!(simulation::SimulationSystem)
    xc_a_2b, xc_b_2b, xc_c_2b, xc_d_2b = load_2b_coeffs("Energy");

    simulation.tot_e_xc_coeffs.xc_a_2b = xc_a_2b;
    simulation.tot_e_xc_coeffs.xc_b_2b = xc_b_2b;
    simulation.tot_e_xc_coeffs.xc_c_2b = xc_c_2b;
    simulation.tot_e_xc_coeffs.xc_d_2b = xc_d_2b;
end

function load_fitted_coeffs!(simulation::SimulationSystem)
    load_1b_pol_e_coeffs!(simulation);
    load_2b_pol_e_coeffs!(simulation);

    load_1b_tot_e_coeffs!(simulation);
    load_2b_tot_e_coeffs!(simulation);
end

function thomas_fermi_ke(Z::Int, element_cloud_coeffs::Matrix)
    # Calculates the Thomas-Fermi kinetic energy of the atom whose atomic 
    # number is specified in the argument of this function.
    tf_coeff = (3.0/10.0) * ((3*(π^2))^(2.0/3.0));

    num_clouds = Int((size(element_cloud_coeffs)[2])/2);
    gauss_c = element_cloud_coeffs[Z,1:2:end];
    gauss_λ = element_cloud_coeffs[Z,2:2:end];

    function foo(r,p)
        # Auxiliary function to integrate the density in spherical coordinates.
        ret_val = 0.0;

        for i in 1:num_clouds
            aux_val = (gauss_λ[i]/π)^(3.0/2.0);
            aux_val *= gauss_c[i]*exp(-gauss_λ[i]*(r^2));
            ret_val += aux_val;
        end

        return 4*π*(r^2)*(ret_val^(5.0/3.0));
    end

    prob = IntegralProblem(foo,(0.0,300.0));
    sol = solve(prob, HCubatureJL(), reltol = 1.0E-8, abstol = 1.0E-8);
    return tf_coeff*(sol[1]);
end

function von_weizsacker_ke(Z::Int, element_cloud_coeffs::Matrix)
    # Calculates the von Weizsäcker kinetic energy of the atom whose atomic 
    # number is specified in the argument of this function.
    vw_coeff = (1.0/8.0);

    num_clouds = Int((size(element_cloud_coeffs)[2])/2);
    gauss_c = element_cloud_coeffs[Z,1:2:end];
    gauss_λ = element_cloud_coeffs[Z,2:2:end];

    function foo(r,p)
        # Auxiliary function to integrate the density in spherical coordinates.
        ρ = 0.0;
        for i in 1:num_clouds
            aux_val = (gauss_λ[i]/π)^(3.0/2.0);
            aux_val *= gauss_c[i]*exp(-gauss_λ[i]*(r^2));
            ρ += aux_val;
        end

        norm_∇ρ = 0.0;
        for i in 1:num_clouds
            aux_λ = gauss_λ[i];
            aux_val = gauss_c[i]*((gauss_λ[i])^(5.0/2.0));
            aux_val *= exp(-aux_λ*(r^2));
            norm_∇ρ += (2.0*r)/(π^(3.0/2.0))*aux_val
        end

        ret_val = 4*π*(r^2)*((norm_∇ρ^2)/ρ);

        if isnan(ret_val)
            return 0.0;
        end

        return ret_val;
    end

    prob = IntegralProblem(foo,(0.0,300.0));
    sol = solve(prob, HCubatureJL(), reltol = 1.0E-8, abstol = 1.0E-8);
    return vw_coeff*(sol[1]);
end

function load_ecp_basis_set()
    max_atomic_number = 18;
    clouds_per_atom = 3;

    file_name = "../Reference Data/atom_density/ecp_coeffs.txt";
    file_id = open(file_name, "r");
    lines = readlines(file_id);
    element_cloud_coeffs = zeros(Float64, 18, 6);

    for i in 1:18
        line_splitted = split(lines[i+1]);
        for j in 1:6
            element_cloud_coeffs[i, j] = parse(Float64, line_splitted[j+1]);
        end
    end

    close(file_id);

    thomas_fermi_kes = zeros(Float64, max_atomic_number);
    von_weizsacker_kes = zeros(Float64, max_atomic_number);

    for Z in 1:max_atomic_number
        thomas_fermi_kes[Z] = thomas_fermi_ke(Z,element_cloud_coeffs);
        von_weizsacker_kes[Z] = von_weizsacker_ke(Z,element_cloud_coeffs);
    end

    return BasisSetSettings(clouds_per_atom,max_atomic_number,
        element_cloud_coeffs,thomas_fermi_kes,von_weizsacker_kes);
end

function initialize_ecp_simulation_environment()
    basis_set_settings = load_ecp_basis_set();

    max_atomic_number = basis_set_settings.max_atomic_number;
    coeff_1b_combinations = max_atomic_number;
    coeff_2b_combinations = 
        round(Int, coeff_1b_combinations*(coeff_1b_combinations+1)/2);

    coeff_1b_map = Dict{Int,Int}();
    inv_coeff_1b_map = Dict{Int,Int}();

    for i in 1:max_atomic_number
        coeff_1b_map[i] = i;
        inv_coeff_1b_map[i] = i;
    end

    coeff_2b_map = Dict{Tuple{Int,Int},Int}();
    inv_coeff_2b_map = Dict{Int,Tuple{Int,Int}}();

    num_xc_coeffs = 0;
    for i in 1:max_atomic_number
        for j in i:max_atomic_number
            num_xc_coeffs += 1;
            coeff_2b_map[(i,j)] = num_xc_coeffs;
            coeff_2b_map[(j,i)] = num_xc_coeffs;
            inv_coeff_2b_map[num_xc_coeffs] = (i,j);
        end
    end    

    coeff_mappings = EmpiricalCoefficientMappings(max_atomic_number,
        coeff_1b_combinations,coeff_2b_combinations,coeff_1b_map,
        inv_coeff_1b_map,coeff_2b_map,inv_coeff_2b_map)

    xc_a_1b = zeros(Float64,coeff_1b_combinations,2);
    xc_c_1b = zeros(Float64,coeff_1b_combinations,2);
    xc_e_1b = zeros(Float64,coeff_1b_combinations,2);
    xc_f_1b = zeros(Float64,coeff_1b_combinations,2);

    xc_a_2b = zeros(Float64,coeff_2b_combinations,3);
    xc_b_2b = zeros(Float64,coeff_2b_combinations,3);
    xc_c_2b = zeros(Float64,coeff_2b_combinations,3);
    xc_d_2b = zeros(Float64,coeff_2b_combinations,3);

    tot_e_xc_coeffs = EmpiricalCoefficients(max_atomic_number,xc_a_1b,xc_c_1b,
        xc_e_1b,xc_f_1b,xc_a_2b,xc_b_2b,xc_c_2b,xc_d_2b);
    pol_e_xc_coeffs = EmpiricalCoefficients(max_atomic_number,xc_a_1b,xc_c_1b,
        xc_e_1b,xc_f_1b,xc_a_2b,xc_b_2b,xc_c_2b,xc_d_2b);

    cast_types = SimulationCastTypes(0,1,2);

    system = MolecularSystem();
    simulation = SimulationSystem(system,pol_e_xc_coeffs,tot_e_xc_coeffs,
        coeff_mappings,basis_set_settings,cast_types);
    load_fitted_coeffs!(simulation);

    return simulation;
end

function reset_xc_coeffs()
    # Use this with caution, it will set all the coefficients in the saved 
    # files to zero.
    simulation = initialize_ecp_simulation_environment();

    function set_coeffs_to_zero!(coeffs::EmpiricalCoefficients)
        coeffs.xc_a_1b .= 0;
        coeffs.xc_c_1b .= 0;
        coeffs.xc_e_1b .= 0;
        coeffs.xc_f_1b .= 0;

        coeffs.xc_a_2b .= 0;
        coeffs.xc_b_2b .= 0;
        coeffs.xc_c_2b .= 0;
        coeffs.xc_d_2b .= 0;
    end

    set_coeffs_to_zero!(simulation.pol_e_xc_coeffs);
    set_coeffs_to_zero!(simulation.tot_e_xc_coeffs);

    save_fitted_coeffs(simulation);
end

function make_atom_system(Z::Int,charge::Int)
    # Makes a molecule object with a single atom whose atomic number is Z.
    simulation = initialize_ecp_simulation_environment();
    simulation.system.charge = charge;

    element_cloud_coeffs = simulation.basis_set_settings.element_cloud_coeffs;
    atom = Molecule();

    atom.atoms_data = zeros(Float64,1,6);
    atom.cloud_data = zeros(Float64,3,6);

    z_eff = round(sum(element_cloud_coeffs[Z,1:2:end]));
    atom.atoms_data[1,5] = Z;
    atom.atoms_data[1,4] = z_eff;
    atom.atoms_data[1,6] = (z_eff - charge) / z_eff;

    atom.cloud_data[1:3,4] .= element_cloud_coeffs[Z,1:2:end];
    atom.cloud_data[1:3,5] .= element_cloud_coeffs[Z, 2:2:end];
    atom.cloud_data[1:3,6] .= (z_eff - charge) / z_eff;

    simulation.system.molecules = [atom];
    return simulation;
end

function make_diatomic_system(Z1::Int,Z2::Int,d::Number,charge::Int)
    # Makes a diatomic molecule with a diatomic separation d in Bohr with the 
    # specified charge.
    simulation = initialize_ecp_simulation_environment();
    simulation.system.charge = charge;

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

    simulation.system.molecules = [molecule];
    return simulation;
end

function make_diatomic_system(Z::Int,d::Number,charge::Int)
    # Makes a homonuclear diatomic molecule with a diatomic separation d in 
    # Bohr with the specified charge.
    return make_diatomic_system(Z,Z,d,charge);
end
