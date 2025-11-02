include("ForceFieldBase.jl")

function save_1b_coeffs(
    xc_a_1b::Dict{Int,Number}, xc_b_1b::Dict{Int,Number},
    xc_c_1b::Dict{Int,Number}, xc_d_1b::Dict{Int,Number}, 
    ke_e_1b::Dict{Int,Number}, ke_f_1b::Dict{Int,Number}, 
    max_atomic_number::Int, type::String)
    
    base_dir = "XC Coeffs/"*type*"/One Body/";
    xc_coeffs_a_id = base_dir*"xc_coeffs_A.txt";
    xc_coeffs_b_id = base_dir*"xc_coeffs_B.txt";
    xc_coeffs_c_id = base_dir*"xc_coeffs_C.txt";
    xc_coeffs_d_id = base_dir*"xc_coeffs_D.txt";
    ke_coeffs_e_id = base_dir*"ke_coeffs_E.txt";
    ke_coeffs_f_id = base_dir*"ke_coeffs_F.txt";

    function save_to_file(file_name::String, coeffs::Dict{Int,Number})
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
    end

    save_to_file(xc_coeffs_a_id,xc_a_1b);
    save_to_file(xc_coeffs_b_id,xc_b_1b);
    save_to_file(xc_coeffs_c_id,xc_c_1b);
    save_to_file(xc_coeffs_d_id,xc_d_1b);
    save_to_file(ke_coeffs_e_id,ke_e_1b);
    save_to_file(ke_coeffs_f_id,ke_f_1b);

    return;
end

function save_1b_coeffs(coeffs::TotalEnergyCoefficients)
    xc_a_1b = coeffs.tot_e_xc_coeffs.xc_a_1b;
    xc_b_1b = coeffs.tot_e_xc_coeffs.xc_b_1b;
    xc_c_1b = coeffs.tot_e_xc_coeffs.xc_c_1b;
    xc_d_1b = coeffs.tot_e_xc_coeffs.xc_d_1b;
    ke_e_1b = coeffs.tot_e_ke_coeffs.ke_e_1b;
    ke_f_1b = coeffs.tot_e_ke_coeffs.ke_f_1b;

    max_atomic_number = coeffs.max_atomic_number;

    type = "Energy";
    save_1b_coeffs(
        xc_a_1b,xc_b_1b,xc_c_1b,xc_d_1b,
        ke_e_1b,ke_f_1b,max_atomic_number,type);
    
    return;
end

function save_1b_coeffs(coeffs::PolarizationEnergyCoefficients)
    xc_a_1b = coeffs.pol_e_xc_coeffs.xc_a_1b;
    xc_b_1b = coeffs.pol_e_xc_coeffs.xc_b_1b;
    xc_c_1b = coeffs.pol_e_xc_coeffs.xc_c_1b;
    xc_d_1b = coeffs.pol_e_xc_coeffs.xc_d_1b;
    ke_e_1b = coeffs.pol_e_ke_coeffs.ke_e_1b;
    ke_f_1b = coeffs.pol_e_ke_coeffs.ke_f_1b;

    max_atomic_number = coeffs.max_atomic_number;

    type = "Polarization";
    save_1b_coeffs(
        xc_a_1b,xc_b_1b,xc_c_1b,xc_d_1b,
        ke_e_1b,ke_f_1b,max_atomic_number,type);
    
    return;
end

function save_1b_pol_e_coeffs(simulation::SimulationSystem)
    save_1b_coeffs(simulation.pol_e_coeffs);
    return;
end

function save_1b_tot_e_coeffs(simulation::SimulationSystem)
    save_1b_coeffs(simulation.tot_e_coeffs);
    return;
end

function save_2b_coeffs(
    xc_a_2b::Dict{Tuple{Int,Int},XCCoeff2B}, xc_b_2b::Dict{Tuple{Int,Int},XCCoeff2B},
    xc_c_2b::Dict{Tuple{Int,Int},XCCoeff2B}, xc_d_2b::Dict{Tuple{Int,Int},XCCoeff2B},
    max_atomic_number::Int, type::String)

    base_dir = "XC Coeffs/"*type*"/Two Body/";
    xc_coeffs_a_id = base_dir*"xc_coeffs_A.txt";
    xc_coeffs_b_id = base_dir*"xc_coeffs_B.txt";
    xc_coeffs_c_id = base_dir*"xc_coeffs_C.txt";
    xc_coeffs_d_id = base_dir*"xc_coeffs_D.txt";

    function save_to_file(file_name::String, coeffs::Dict{Tuple{Int,Int},XCCoeff2B})
        file_id = open(file_name,"w");
        for i in 1:max_atomic_number
            for j in i:max_atomic_number
                if haskey(coeffs,(i,j))
                    write(file_id, @sprintf "%18.8E " coeffs[(i,j)].m);
                    write(file_id, @sprintf "%18.8E " coeffs[(i,j)].b);
                else
                    write(file_id, @sprintf "%18.8E " 0.0);
                    write(file_id, @sprintf "%18.8E " 0.0);
                end
                
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

    return;
end

function save_2b_coeffs(
    morse_2b::Dict{Tuple{Int,Int},MorsePotentialCoefficients},
    max_atomic_number::Int, type::String)

    base_dir = "XC Coeffs/"*type*"/Two Body/";
    morse_coeffs_id = base_dir*"morse_coeffs.txt";

    function save_to_file(file_name::String, coeffs::Dict{Tuple{Int,Int},MorsePotentialCoefficients})
        file_id = open(file_name,"w");
        for i in 1:max_atomic_number
            for j in i:max_atomic_number
                if haskey(coeffs,(i,j))
                    write(file_id, @sprintf "%18.8E " coeffs[(i,j)].depth);
                    write(file_id, @sprintf "%18.8E " coeffs[(i,j)].exponential_decay);
                    write(file_id, @sprintf "%18.8E " coeffs[(i,j)].equilibrium_distance);
                else
                    write(file_id, @sprintf "%18.8E " 0.0);
                    write(file_id, @sprintf "%18.8E " 0.0);
                    write(file_id, @sprintf "%18.8E " 0.0);
                end
                
                write(file_id, @sprintf "%5s " get_element_symbol(i));
                write(file_id, @sprintf "%5s " get_element_symbol(j));
                write(file_id, @sprintf "\n");
            end
        end

        close(file_id);
    end

    save_to_file(morse_coeffs_id,morse_2b);

    return;
end

function save_2b_coeffs(coeffs::TotalEnergyCoefficients)
    xc_a_2b = coeffs.tot_e_xc_coeffs.xc_a_2b;
    xc_b_2b = coeffs.tot_e_xc_coeffs.xc_b_2b;
    xc_c_2b = coeffs.tot_e_xc_coeffs.xc_c_2b;
    xc_d_2b = coeffs.tot_e_xc_coeffs.xc_d_2b;

    morse_2b = coeffs.tot_e_static_coeffs.morse_2b

    max_atomic_number = coeffs.max_atomic_number;

    type = "Energy";
    save_2b_coeffs(xc_a_2b,xc_b_2b,xc_c_2b,xc_d_2b,max_atomic_number,type);
    save_2b_coeffs(morse_2b,max_atomic_number,type);
    
    return;
end

function save_2b_coeffs(coeffs::PolarizationEnergyCoefficients)
    xc_a_2b = coeffs.pol_e_xc_coeffs.xc_a_2b;
    xc_b_2b = coeffs.pol_e_xc_coeffs.xc_b_2b;
    xc_c_2b = coeffs.pol_e_xc_coeffs.xc_c_2b;
    xc_d_2b = coeffs.pol_e_xc_coeffs.xc_d_2b;

    max_atomic_number = coeffs.max_atomic_number;

    type = "Polarization";
    save_2b_coeffs(xc_a_2b,xc_b_2b,xc_c_2b,xc_d_2b,max_atomic_number,type);
    
    return;
end

function save_2b_pol_e_coeffs(simulation::SimulationSystem)
    save_2b_coeffs(simulation.pol_e_coeffs);
    return;
end

function save_2b_tot_e_coeffs(simulation::SimulationSystem)
    save_2b_coeffs(simulation.tot_e_coeffs);
    return;
end

function load_1b_xc_coeffs(type::String)
    base_dir = "XC Coeffs/"*type*"/One Body/";
    xc_coeffs_a_id = base_dir*"xc_coeffs_A.txt";
    xc_coeffs_b_id = base_dir*"xc_coeffs_B.txt";
    xc_coeffs_c_id = base_dir*"xc_coeffs_C.txt";
    xc_coeffs_d_id = base_dir*"xc_coeffs_D.txt";

    function read_file(file_name::String)
        lines = readlines(file_name);
        coeffs = Dict{Int,Number}();

        for i in eachindex(lines)
            line_splitted = string.(split(lines[i]));
            elem_symbol = line_splitted[2];
            atomic_number = get_atomic_number(elem_symbol);

            coeffs[atomic_number] = parse(Float64,line_splitted[1]);
        end

        return coeffs;
    end

    xc_a_1b = read_file(xc_coeffs_a_id);
    xc_b_1b = read_file(xc_coeffs_b_id);
    xc_c_1b = read_file(xc_coeffs_c_id);
    xc_d_1b = read_file(xc_coeffs_d_id);

    return xc_a_1b, xc_b_1b, xc_c_1b, xc_d_1b;
end

function load_1b_ke_coeffs(type::String)
    base_dir = "XC Coeffs/"*type*"/One Body/";
    ke_coeffs_e_id = base_dir*"ke_coeffs_E.txt";
    ke_coeffs_f_id = base_dir*"ke_coeffs_F.txt";

    function read_file(file_name::String)
        lines = readlines(file_name);
        coeffs = Dict{Int,Number}();

        for i in eachindex(lines)
            line_splitted = string.(split(lines[i]));
            elem_symbol = line_splitted[2];
            atomic_number = get_atomic_number(elem_symbol);

            coeffs[atomic_number] = parse(Float64,line_splitted[1]);
        end

        return coeffs;
    end

    ke_e_1b = read_file(ke_coeffs_e_id);
    ke_f_1b = read_file(ke_coeffs_f_id);

    return ke_e_1b, ke_f_1b;
end

function load_1b_pol_e_coeffs!(simulation::SimulationSystem)
    xc_a_1b, xc_b_1b, xc_c_1b, xc_d_1b = load_1b_xc_coeffs("Polarization");
    ke_e_1b, ke_f_1b = load_1b_ke_coeffs("Polarization");
    
    simulation.pol_e_coeffs.pol_e_xc_coeffs.xc_a_1b = xc_a_1b;
    simulation.pol_e_coeffs.pol_e_xc_coeffs.xc_b_1b = xc_b_1b;
    simulation.pol_e_coeffs.pol_e_xc_coeffs.xc_c_1b = xc_c_1b;
    simulation.pol_e_coeffs.pol_e_xc_coeffs.xc_d_1b = xc_d_1b;
    simulation.pol_e_coeffs.pol_e_ke_coeffs.ke_e_1b = ke_e_1b;
    simulation.pol_e_coeffs.pol_e_ke_coeffs.ke_f_1b = ke_f_1b;

    return;
end

function load_1b_tot_e_coeffs!(simulation::SimulationSystem)
    xc_a_1b, xc_b_1b, xc_c_1b, xc_d_1b = load_1b_xc_coeffs("Energy");
    ke_e_1b, ke_f_1b = load_1b_ke_coeffs("Energy");
    
    simulation.tot_e_coeffs.tot_e_xc_coeffs.xc_a_1b = xc_a_1b;
    simulation.tot_e_coeffs.tot_e_xc_coeffs.xc_b_1b = xc_b_1b;
    simulation.tot_e_coeffs.tot_e_xc_coeffs.xc_c_1b = xc_c_1b;
    simulation.tot_e_coeffs.tot_e_xc_coeffs.xc_d_1b = xc_d_1b;
    simulation.tot_e_coeffs.tot_e_ke_coeffs.ke_e_1b = ke_e_1b;
    simulation.tot_e_coeffs.tot_e_ke_coeffs.ke_f_1b = ke_f_1b;

    return;
end

function load_2b_xc_coeffs(type::String)
    base_dir = "XC Coeffs/"*type*"/Two Body/";
    xc_coeffs_a_id = base_dir*"xc_coeffs_A.txt";
    xc_coeffs_b_id = base_dir*"xc_coeffs_B.txt";
    xc_coeffs_c_id = base_dir*"xc_coeffs_C.txt";
    xc_coeffs_d_id = base_dir*"xc_coeffs_D.txt";

    function read_xc_file(file_name::String)
        lines = readlines(file_name);
        coeffs = Dict{Tuple{Int,Int},XCCoeff2B}();

        for i in eachindex(lines)
            line_splitted = string.(split(lines[i]));

            elem_symbol_1 = line_splitted[end-1];
            elem_symbol_2 = line_splitted[end-0];
            atomic_number_1 = get_atomic_number(elem_symbol_1);
            atomic_number_2 = get_atomic_number(elem_symbol_2);
            dict_key_1 = (atomic_number_1,atomic_number_2);
            dict_key_2 = (atomic_number_2,atomic_number_1);

            coeff = XCCoeff2B(
                parse(Float64,line_splitted[1]),
                parse(Float64,line_splitted[2])
            )

            coeffs[dict_key_1] = coeff;
            coeffs[dict_key_2] = coeff;
        end

        return coeffs;
    end

    xc_a_2b = read_xc_file(xc_coeffs_a_id);
    xc_b_2b = read_xc_file(xc_coeffs_b_id);
    xc_c_2b = read_xc_file(xc_coeffs_c_id);
    xc_d_2b = read_xc_file(xc_coeffs_d_id);

    return xc_a_2b, xc_b_2b, xc_c_2b, xc_d_2b;
end

function load_2b_morse_coeffs()
    base_dir = "XC Coeffs/Energy/Two Body/";
    morse_coeffs_id = base_dir*"morse_coeffs.txt";

    function read_xc_file(file_name::String)
        lines = readlines(file_name);
        coeffs = Dict{Tuple{Int,Int},MorsePotentialCoefficients}();

        for i in eachindex(lines)
            line_splitted = string.(split(lines[i]));

            elem_symbol_1 = line_splitted[end-1];
            elem_symbol_2 = line_splitted[end-0];
            atomic_number_1 = get_atomic_number(elem_symbol_1);
            atomic_number_2 = get_atomic_number(elem_symbol_2);
            dict_key_1 = (atomic_number_1,atomic_number_2);
            dict_key_2 = (atomic_number_2,atomic_number_1);

            coeff = MorsePotentialCoefficients(0.0,0.0,0.0);
            coeff.depth = parse(Float64,line_splitted[1]);
            coeff.exponential_decay = parse(Float64,line_splitted[2]);
            coeff.equilibrium_distance = parse(Float64,line_splitted[3]);

            coeffs[dict_key_1] = coeff;
            coeffs[dict_key_2] = coeff;
        end

        return coeffs;
    end

    morse_2b = read_xc_file(morse_coeffs_id);

    return morse_2b;
end

function load_2b_pol_e_coeffs!(simulation::SimulationSystem)
    xc_a_2b, xc_b_2b, xc_c_2b, xc_d_2b = load_2b_xc_coeffs("Polarization");

    simulation.pol_e_coeffs.pol_e_xc_coeffs.xc_a_2b = xc_a_2b;
    simulation.pol_e_coeffs.pol_e_xc_coeffs.xc_b_2b = xc_b_2b;
    simulation.pol_e_coeffs.pol_e_xc_coeffs.xc_c_2b = xc_c_2b;
    simulation.pol_e_coeffs.pol_e_xc_coeffs.xc_d_2b = xc_d_2b;

    return;
end

function load_2b_tot_e_coeffs!(simulation::SimulationSystem)
    xc_a_2b, xc_b_2b, xc_c_2b, xc_d_2b = load_2b_xc_coeffs("Energy");
    morse_2b = load_2b_morse_coeffs();

    simulation.tot_e_coeffs.tot_e_xc_coeffs.xc_a_2b = xc_a_2b;
    simulation.tot_e_coeffs.tot_e_xc_coeffs.xc_b_2b = xc_b_2b;
    simulation.tot_e_coeffs.tot_e_xc_coeffs.xc_c_2b = xc_c_2b;
    simulation.tot_e_coeffs.tot_e_xc_coeffs.xc_d_2b = xc_d_2b;
    simulation.tot_e_coeffs.tot_e_static_coeffs.morse_2b = morse_2b;

    return;
end

function load_fitted_coeffs!(simulation::SimulationSystem)
    load_1b_pol_e_coeffs!(simulation);
    load_2b_pol_e_coeffs!(simulation);

    load_1b_tot_e_coeffs!(simulation);
    load_2b_tot_e_coeffs!(simulation);

    simulation.pol_e_coeffs.max_atomic_number = 
        length(simulation.pol_e_coeffs.pol_e_xc_coeffs.xc_a_1b);

    simulation.tot_e_coeffs.max_atomic_number = 
        length(simulation.tot_e_coeffs.tot_e_xc_coeffs.xc_a_1b);

    return;
end

function save_fitted_coeffs(simulation::SimulationSystem)
    save_1b_coeffs(simulation.tot_e_coeffs);
    save_2b_coeffs(simulation.tot_e_coeffs);
    
    save_1b_coeffs(simulation.pol_e_coeffs);
    save_2b_coeffs(simulation.pol_e_coeffs);

    return;
end

function get_xc_2b_coeff(xc_coeffs::XCCoeff2B, d::Number)
    return xc_coeffs.b + d * xc_coeffs.m;
end

function reset_fitted_coeffs()
    max_atomic_number = 18;

    xc_a_1b = Dict{Int,Number}();
    xc_b_1b = Dict{Int,Number}();
    xc_c_1b = Dict{Int,Number}();
    xc_d_1b = Dict{Int,Number}();
    ke_e_1b = Dict{Int,Number}();
    ke_f_1b = Dict{Int,Number}();

    xc_a_2b = Dict{Tuple{Int,Int},XCCoeff2B}();
    xc_b_2b = Dict{Tuple{Int,Int},XCCoeff2B}();
    xc_c_2b = Dict{Tuple{Int,Int},XCCoeff2B}();
    xc_d_2b = Dict{Tuple{Int,Int},XCCoeff2B}();
    morse_2b = Dict{Tuple{Int,Int},MorsePotentialCoefficients}();

    tot_e_xc_coeffs = EmpiricalXCCoefficients(
        xc_a_1b,xc_b_1b,xc_c_1b,xc_d_1b,xc_a_2b,xc_b_2b,xc_c_2b,xc_d_2b);
    tot_e_ke_coeffs = EmpiricalKECoefficients(ke_e_1b,ke_f_1b);
    tot_e_static_coeffs = EmpiricalStaticCoefficients(morse_2b);

    pol_e_xc_coeffs = EmpiricalXCCoefficients(
        xc_a_1b,xc_b_1b,xc_c_1b,xc_d_1b,xc_a_2b,xc_b_2b,xc_c_2b,xc_d_2b);
    pol_e_ke_coeffs = EmpiricalKECoefficients(ke_e_1b,ke_f_1b);

    default_tot_e_coeffs = TotalEnergyCoefficients(max_atomic_number,
        tot_e_xc_coeffs,tot_e_ke_coeffs,tot_e_static_coeffs);
    default_pol_e_coeffs = PolarizationEnergyCoefficients(max_atomic_number,
        pol_e_xc_coeffs,pol_e_ke_coeffs);
    
    save_1b_coeffs(default_tot_e_coeffs);
    save_2b_coeffs(default_tot_e_coeffs);

    save_1b_coeffs(default_pol_e_coeffs);
    save_2b_coeffs(default_pol_e_coeffs);

    return;
end
