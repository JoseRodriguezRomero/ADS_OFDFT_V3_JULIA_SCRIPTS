include("ForceFieldBase.jl")

function save_1b_coeffs(coeffs::Dict{Int,Number},
    max_atomic_number::Int, type::String, name::String)
    base_dir = "Fitted Coefficients/"*type*"/One Body/";
    file_name = base_dir*name*".txt";

    file_id = open(file_name,"w");
    for i in 1:max_atomic_number
        if haskey(coeffs,i)
            write(file_id, @sprintf "%18.8E " coeffs[i]);
        else
            write(file_id, @sprintf "%18.8E " 0.0);
        end
        
        write(file_id, @sprintf "%5s " get_element_symbol(i));
        write(file_id, @sprintf "\n");
    end

    close(file_id);
    return;
end

function save_1b_coeffs(coeffs::PolarizationEnergyCoefficients, 
    max_atomic_number::Int)
    xc_a_1b = coeffs.xc_coeffs.xc_a_1b;
    xc_b_1b = coeffs.xc_coeffs.xc_b_1b;
    xc_c_1b = coeffs.xc_coeffs.xc_c_1b;
    xc_d_1b = coeffs.xc_coeffs.xc_d_1b;
    ke_e_1b = coeffs.ke_coeffs.ke_e_1b;
    ke_f_1b = coeffs.ke_coeffs.ke_f_1b;

    type = "Polarization";
    save_1b_coeffs(xc_a_1b,max_atomic_number,type,"xc_coeffs_A");
    save_1b_coeffs(xc_b_1b,max_atomic_number,type,"xc_coeffs_B");
    save_1b_coeffs(xc_c_1b,max_atomic_number,type,"xc_coeffs_C");
    save_1b_coeffs(xc_d_1b,max_atomic_number,type,"xc_coeffs_D");
    save_1b_coeffs(ke_e_1b,max_atomic_number,type,"ke_coeffs_E");
    save_1b_coeffs(ke_f_1b,max_atomic_number,type,"ke_coeffs_F");
    
    return;
end

function save_1b_coeffs(coeffs::TotalEnergyCoefficients, 
    max_atomic_number::Int)
    xc_a_1b = coeffs.xc_coeffs.xc_a_1b;
    xc_b_1b = coeffs.xc_coeffs.xc_b_1b;
    xc_c_1b = coeffs.xc_coeffs.xc_c_1b;
    xc_d_1b = coeffs.xc_coeffs.xc_d_1b;
    ke_e_1b = coeffs.ke_coeffs.ke_e_1b;
    ke_f_1b = coeffs.ke_coeffs.ke_f_1b;

    type = "Energy";
    save_1b_coeffs(xc_a_1b,max_atomic_number,type,"xc_coeffs_A");
    save_1b_coeffs(xc_b_1b,max_atomic_number,type,"xc_coeffs_B");
    save_1b_coeffs(xc_c_1b,max_atomic_number,type,"xc_coeffs_C");
    save_1b_coeffs(xc_d_1b,max_atomic_number,type,"xc_coeffs_D");
    save_1b_coeffs(ke_e_1b,max_atomic_number,type,"ke_coeffs_E");
    save_1b_coeffs(ke_f_1b,max_atomic_number,type,"ke_coeffs_F");
    
    return;
end

function save_1b_pol_e_coeffs(simulation::SimulationSystem)
    coeffs = simulation.pol_e_coeffs;
    max_atomic_number = simulation.basis_set_settings.max_atomic_number;

    save_1b_coeffs(coeffs,max_atomic_number);
    return;
end

function save_1b_tot_e_coeffs(simulation::SimulationSystem)
    coeffs = simulation.tot_e_coeffs;
    max_atomic_number = simulation.basis_set_settings.max_atomic_number;

    save_1b_coeffs(coeffs,max_atomic_number);
    return;
end

function save_2b_coeffs(coeffs::Dict{Tuple{Int,Int},Number},
    max_atomic_number::Int, type::String, name::String)

    base_dir = "Fitted Coefficients/"*type*"/Two Body/";
    file_name = base_dir*name*".txt";

    file_id = open(file_name,"w");
    for i in 1:max_atomic_number
        for j in i:max_atomic_number
            if haskey(coeffs,(i,j))
                write(file_id, @sprintf "%18.8E " coeffs[(i,j)]);
            else
                write(file_id, @sprintf "%18.8E " 0.0);
            end
            
            write(file_id, @sprintf "%5s " get_element_symbol(i));
            write(file_id, @sprintf "%5s " get_element_symbol(j));
            write(file_id, @sprintf "\n");
        end
    end

    close(file_id);
    return;
end

function save_2b_coeffs(coeffs::PolarizationEnergyCoefficients, 
    max_atomic_number::Int)
    xc_a_2b = coeffs.xc_coeffs.xc_a_2b;
    xc_b_2b = coeffs.xc_coeffs.xc_b_2b;
    xc_c_2b = coeffs.xc_coeffs.xc_c_2b;
    xc_d_2b = coeffs.xc_coeffs.xc_d_2b;

    type = "Polarization";
    save_2b_coeffs(xc_a_2b,max_atomic_number,type,"xc_coeffs_A");
    save_2b_coeffs(xc_b_2b,max_atomic_number,type,"xc_coeffs_B");
    save_2b_coeffs(xc_c_2b,max_atomic_number,type,"xc_coeffs_C");
    save_2b_coeffs(xc_d_2b,max_atomic_number,type,"xc_coeffs_D");
    
    return;
end

function save_2b_coeffs(coeffs::TotalEnergyCoefficients, 
    max_atomic_number::Int)
    xc_a_2b = coeffs.xc_coeffs.xc_a_2b;
    xc_b_2b = coeffs.xc_coeffs.xc_b_2b;
    xc_c_2b = coeffs.xc_coeffs.xc_c_2b;
    xc_d_2b = coeffs.xc_coeffs.xc_d_2b;

    morse_depth = coeffs.non_polarizable_coeffs.depth;
    morse_stiffness_parameter = 
        coeffs.non_polarizable_coeffs.stiffness_parameter;
    morse_equilibrium_distance = 
        coeffs.non_polarizable_coeffs.equilibrium_distance;

    type = "Energy";
    save_2b_coeffs(xc_a_2b,max_atomic_number,type,"xc_coeffs_A");
    save_2b_coeffs(xc_b_2b,max_atomic_number,type,"xc_coeffs_B");
    save_2b_coeffs(xc_c_2b,max_atomic_number,type,"xc_coeffs_C");
    save_2b_coeffs(xc_d_2b,max_atomic_number,type,"xc_coeffs_D");

    save_2b_coeffs(morse_depth,max_atomic_number,
        type,"non_polarizable_morse_depth");
    save_2b_coeffs(morse_stiffness_parameter,max_atomic_number,
        type,"non_polarizable_morse_stiffness_parameter");
    save_2b_coeffs(morse_equilibrium_distance,max_atomic_number,
        type,"non_polarizable_morse_equilibrium_distance");
    
    return;
end

function save_2b_pol_e_coeffs(simulation::SimulationSystem)
    coeffs = simulation.pol_e_coeffs;
    max_atomic_number = simulation.basis_set_settings.max_atomic_number;

    save_2b_coeffs(coeffs,max_atomic_number);
    return;
end

function save_2b_tot_e_coeffs(simulation::SimulationSystem)
    coeffs = simulation.tot_e_coeffs;
    max_atomic_number = simulation.basis_set_settings.max_atomic_number;

    save_2b_coeffs(coeffs,max_atomic_number);
    return;
end

function load_1b_coeffs(type::String, name::String)
    base_dir = "Fitted Coefficients/"*type*"/One Body/";
    file_name = base_dir*name*".txt";

    coeffs = Dict{Int,Number}();
    if !isfile(file_name)
        return coeffs;
    end

    lines = readlines(file_name);

    for i in eachindex(lines)
        line_splitted = string.(split(lines[i]));
        elem_symbol = line_splitted[2];
        atomic_number = get_atomic_number(elem_symbol);

        coeffs[atomic_number] = parse(Float64,line_splitted[1]);
    end

    return coeffs;
end

function load_1b_pol_e_coeffs!(simulation::SimulationSystem)
    type = "Polarization";
    xc_a_1b = load_1b_coeffs(type,"xc_coeffs_A");
    xc_b_1b = load_1b_coeffs(type,"xc_coeffs_B");
    xc_c_1b = load_1b_coeffs(type,"xc_coeffs_C");
    xc_d_1b = load_1b_coeffs(type,"xc_coeffs_D");
    ke_e_1b = load_1b_coeffs(type,"ke_coeffs_E");
    ke_f_1b = load_1b_coeffs(type,"ke_coeffs_F");
    
    simulation.pol_e_coeffs.xc_coeffs.xc_a_1b = xc_a_1b;
    simulation.pol_e_coeffs.xc_coeffs.xc_b_1b = xc_b_1b;
    simulation.pol_e_coeffs.xc_coeffs.xc_c_1b = xc_c_1b;
    simulation.pol_e_coeffs.xc_coeffs.xc_d_1b = xc_d_1b;
    simulation.pol_e_coeffs.ke_coeffs.ke_e_1b = ke_e_1b;
    simulation.pol_e_coeffs.ke_coeffs.ke_f_1b = ke_f_1b;

    return;
end

function load_1b_tot_e_coeffs!(simulation::SimulationSystem)
    type = "Energy";
    xc_a_1b = load_1b_coeffs(type,"xc_coeffs_A");
    xc_b_1b = load_1b_coeffs(type,"xc_coeffs_B");
    xc_c_1b = load_1b_coeffs(type,"xc_coeffs_C");
    xc_d_1b = load_1b_coeffs(type,"xc_coeffs_D");
    ke_e_1b = load_1b_coeffs(type,"ke_coeffs_E");
    ke_f_1b = load_1b_coeffs(type,"ke_coeffs_F");
    
    simulation.tot_e_coeffs.xc_coeffs.xc_a_1b = xc_a_1b;
    simulation.tot_e_coeffs.xc_coeffs.xc_b_1b = xc_b_1b;
    simulation.tot_e_coeffs.xc_coeffs.xc_c_1b = xc_c_1b;
    simulation.tot_e_coeffs.xc_coeffs.xc_d_1b = xc_d_1b;
    simulation.tot_e_coeffs.ke_coeffs.ke_e_1b = ke_e_1b;
    simulation.tot_e_coeffs.ke_coeffs.ke_f_1b = ke_f_1b;

    return;
end

function load_2b_coeffs(type::String, name::String)
    base_dir = "Fitted Coefficients/"*type*"/Two Body/";
    file_name = base_dir*name*".txt";

    coeffs = Dict{Tuple{Int,Int},Number}();
    if !isfile(file_name)
        return coeffs;
    end

    lines = readlines(file_name);
    
    for i in eachindex(lines)
        line_splitted = string.(split(lines[i]));

        elem_symbol_1 = line_splitted[end-1];
        elem_symbol_2 = line_splitted[end-0];
        atomic_number_1 = get_atomic_number(elem_symbol_1);
        atomic_number_2 = get_atomic_number(elem_symbol_2);
        dict_key_1 = (atomic_number_1,atomic_number_2);
        dict_key_2 = (atomic_number_2,atomic_number_1);

        coeffs[dict_key_1] = parse(Float64,line_splitted[1]);
        coeffs[dict_key_2] = parse(Float64,line_splitted[1]);
    end

    return coeffs;
end

function load_2b_pol_e_coeffs!(simulation::SimulationSystem)
    type = "Polarization";
    xc_a_2b = load_2b_coeffs(type,"xc_coeffs_A");
    xc_b_2b = load_2b_coeffs(type,"xc_coeffs_B");
    xc_c_2b = load_2b_coeffs(type,"xc_coeffs_C");
    xc_d_2b = load_2b_coeffs(type,"xc_coeffs_D");

    simulation.pol_e_coeffs.xc_coeffs.xc_a_2b = xc_a_2b;
    simulation.pol_e_coeffs.xc_coeffs.xc_b_2b = xc_b_2b;
    simulation.pol_e_coeffs.xc_coeffs.xc_c_2b = xc_c_2b;
    simulation.pol_e_coeffs.xc_coeffs.xc_d_2b = xc_d_2b;

    return;
end

function load_2b_tot_e_coeffs!(simulation::SimulationSystem)
    type = "Energy";
    xc_a_2b = load_2b_coeffs(type,"xc_coeffs_A");
    xc_b_2b = load_2b_coeffs(type,"xc_coeffs_B");
    xc_c_2b = load_2b_coeffs(type,"xc_coeffs_C");
    xc_d_2b = load_2b_coeffs(type,"xc_coeffs_D");

    non_polarizable_morse_depth = load_2b_coeffs(type,
        "non_polarizable_morse_depth");
    non_polarizable_morse_stiffness_parameter = load_2b_coeffs(type,
        "non_polarizable_morse_stiffness_parameter");
    non_polarizable_morse_equilibrium_distance = load_2b_coeffs(type,
        "non_polarizable_morse_equilibrium_distance");

    simulation.tot_e_coeffs.xc_coeffs.xc_a_2b = xc_a_2b;
    simulation.tot_e_coeffs.xc_coeffs.xc_b_2b = xc_b_2b;
    simulation.tot_e_coeffs.xc_coeffs.xc_c_2b = xc_c_2b;
    simulation.tot_e_coeffs.xc_coeffs.xc_d_2b = xc_d_2b;

    simulation.tot_e_coeffs.non_polarizable_coeffs.depth = 
        non_polarizable_morse_depth;
    simulation.tot_e_coeffs.non_polarizable_coeffs.stiffness_parameter = 
        non_polarizable_morse_stiffness_parameter;
    simulation.tot_e_coeffs.non_polarizable_coeffs.equilibrium_distance = 
        non_polarizable_morse_equilibrium_distance;

    return;
end

function load_fitted_coeffs!(simulation::SimulationSystem)
    load_1b_pol_e_coeffs!(simulation);
    load_2b_pol_e_coeffs!(simulation);

    load_1b_tot_e_coeffs!(simulation);
    load_2b_tot_e_coeffs!(simulation);

    return;
end

function save_fitted_coeffs(simulation::SimulationSystem)
    save_1b_pol_e_coeffs(simulation);
    save_2b_pol_e_coeffs(simulation);

    save_1b_tot_e_coeffs(simulation);
    save_2b_tot_e_coeffs(simulation);

    return;
end

function reset_fitted_coeffs()
    max_atomic_number = 18;

    xc_a_1b = Dict{Int,Number}();
    xc_b_1b = Dict{Int,Number}();
    xc_c_1b = Dict{Int,Number}();
    xc_d_1b = Dict{Int,Number}();
    ke_e_1b = Dict{Int,Number}();
    ke_f_1b = Dict{Int,Number}();

    xc_a_2b = Dict{Tuple{Int,Int},Number}();
    xc_b_2b = Dict{Tuple{Int,Int},Number}();
    xc_c_2b = Dict{Tuple{Int,Int},Number}();
    xc_d_2b = Dict{Tuple{Int,Int},Number}();

    morse_depth = Dict{Int,Number}();
    morse_stiffness_parameter = Dict{Int,Number}();
    morse_equilibrium_distance = Dict{Int,Number}();

    xc_coeffs = EmpiricalXCCoefficients(
        xc_a_1b,xc_b_1b,xc_c_1b,xc_d_1b,xc_a_2b,xc_b_2b,xc_c_2b,xc_d_2b);
    ke_coeffs = EmpiricalKECoefficients(ke_e_1b,ke_f_1b);
    non_polarizable_coeffs = EmpiricalMorseCoefficients(morse_depth,
        morse_stiffness_parameter,morse_equilibrium_distance);

    tot_e_coeffs = TotalEnergyCoefficients(
        xc_coeffs,ke_coeffs,non_polarizable_coeffs);
    pol_e_coeffs = PolarizationEnergyCoefficients(xc_coeffs,ke_coeffs);
    
    save_1b_coeffs(tot_e_coeffs,max_atomic_number);
    save_2b_coeffs(tot_e_coeffs,max_atomic_number);

    save_1b_coeffs(pol_e_coeffs,max_atomic_number);
    save_2b_coeffs(pol_e_coeffs,max_atomic_number);

    return;
end
