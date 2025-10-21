include("ForceFieldBase.jl")

function save_1b_coeffs(coeffs::EmpiricalXCCoefficients, type::String)
    base_dir = "XC Coeffs/"*type*"/One Body/";
    xc_coeffs_a_id = base_dir*"xc_coeffs_A.txt";
    xc_coeffs_b_id = base_dir*"xc_coeffs_B.txt";
    xc_coeffs_c_id = base_dir*"xc_coeffs_C.txt";
    xc_coeffs_d_id = base_dir*"xc_coeffs_D.txt";
    xc_coeffs_e_id = base_dir*"xc_coeffs_E.txt";
    xc_coeffs_f_id = base_dir*"xc_coeffs_F.txt";

    xc_a_1b = coeffs.xc_a_1b;
    xc_b_1b = coeffs.xc_b_1b;
    xc_c_1b = coeffs.xc_c_1b;
    xc_d_1b = coeffs.xc_d_1b;
    xc_e_1b = coeffs.xc_e_1b;
    xc_f_1b = coeffs.xc_f_1b;

    max_atomic_number = coeffs.max_atomic_number;

    function save_to_file(file_name::String, coeffs::Dict{Int, XCCoefficient})
        file_id = open(file_name,"w");
        for i in 1:max_atomic_number
            if haskey(coeffs,i)
                write(file_id, @sprintf "%18.8E " coeffs[i].m);
                write(file_id, @sprintf "%18.8E " coeffs[i].b);
            else
                write(file_id, @sprintf "%18.8E " 0.0);
                write(file_id, @sprintf "%18.8E " 0.0);
            end
            
            write(file_id, @sprintf "%5s " get_element_symbol(i));
            write(file_id, @sprintf "\n");
        end

        close(file_id);
    end

    function save_to_file(file_name::String, coeffs::Dict{Int, Number})
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
    save_to_file(xc_coeffs_e_id,xc_e_1b);
    save_to_file(xc_coeffs_f_id,xc_f_1b);

    return;
end

function save_1b_pol_e_coeffs(simulation::SimulationSystem)
    save_1b_coeffs(simulation.pol_e_xc_coeffs,"Polarization");

    return;
end

function save_1b_tot_e_coeffs(simulation::SimulationSystem)
    save_1b_coeffs(simulation.tot_e_xc_coeffs,"Energy");

    return;
end

function save_2b_coeffs(coeffs::EmpiricalXCCoefficients, type::String)
    base_dir = "XC Coeffs/"*type*"/Two Body/";
    xc_coeffs_a_id = base_dir*"xc_coeffs_A.txt";
    xc_coeffs_b_id = base_dir*"xc_coeffs_B.txt";
    xc_coeffs_c_id = base_dir*"xc_coeffs_C.txt";
    xc_coeffs_d_id = base_dir*"xc_coeffs_D.txt";

    xc_a_2b = coeffs.xc_a_2b;
    xc_b_2b = coeffs.xc_b_2b;
    xc_c_2b = coeffs.xc_c_2b;
    xc_d_2b = coeffs.xc_d_2b;

    max_atomic_number = coeffs.max_atomic_number;

    function save_to_file(file_name::String, coeffs::Dict{Tuple{Int, Int}, XCCoefficient})
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

function save_2b_pol_e_coeffs(simulation::SimulationSystem)
    save_2b_coeffs(simulation.pol_e_xc_coeffs,"Polarization");

    return;
end

function save_2b_tot_e_coeffs(simulation::SimulationSystem)
    save_2b_coeffs(simulation.tot_e_xc_coeffs,"Energy");

    return;
end

function load_1b_coeffs(type::String)
    base_dir = "XC Coeffs/"*type*"/One Body/";
    xc_coeffs_a_id = base_dir*"xc_coeffs_A.txt";
    xc_coeffs_b_id = base_dir*"xc_coeffs_B.txt";
    xc_coeffs_c_id = base_dir*"xc_coeffs_C.txt";
    xc_coeffs_d_id = base_dir*"xc_coeffs_D.txt";
    xc_coeffs_e_id = base_dir*"xc_coeffs_E.txt";
    xc_coeffs_f_id = base_dir*"xc_coeffs_F.txt";

    function read_xc_file(file_name::String)
        lines = readlines(file_name);
        coeffs = Dict{Int,XCCoefficient}();

        for i in eachindex(lines)
            line_splitted = string.(split(lines[i]));
            elem_symbol = line_splitted[end];
            atomic_number = get_atomic_number(elem_symbol);

            m = parse(Float64,line_splitted[1]);
            b = parse(Float64,line_splitted[2]);
            coeffs[atomic_number] = XCCoefficient(m,b);
        end

        return coeffs;
    end

    function read_ke_file(file_name::String)
        lines = readlines(file_name);
        coeffs = Dict{Int,Number}();

        for i in eachindex(lines)
            line_splitted = string.(split(lines[i]));
            elem_symbol = line_splitted[end];
            atomic_number = get_atomic_number(elem_symbol);

            coeffs[atomic_number] = parse(Float64,line_splitted[1]);
        end

        return coeffs;
    end

    xc_a_1b = read_xc_file(xc_coeffs_a_id);
    xc_b_1b = read_xc_file(xc_coeffs_b_id);
    xc_c_1b = read_xc_file(xc_coeffs_c_id);
    xc_d_1b = read_xc_file(xc_coeffs_d_id);
    xc_e_1b = read_ke_file(xc_coeffs_e_id);
    xc_f_1b = read_ke_file(xc_coeffs_f_id);

    return xc_a_1b, xc_b_1b, xc_c_1b, xc_d_1b, xc_e_1b, xc_f_1b;
end

function load_1b_pol_e_coeffs!(simulation::SimulationSystem)
    xc_a_1b, xc_b_1b, xc_c_1b, xc_d_1b, xc_e_1b, xc_f_1b =
        load_1b_coeffs("Polarization");
    
    simulation.pol_e_xc_coeffs.xc_a_1b = xc_a_1b;
    simulation.pol_e_xc_coeffs.xc_b_1b = xc_b_1b;
    simulation.pol_e_xc_coeffs.xc_c_1b = xc_c_1b;
    simulation.pol_e_xc_coeffs.xc_d_1b = xc_d_1b;
    simulation.pol_e_xc_coeffs.xc_e_1b = xc_e_1b;
    simulation.pol_e_xc_coeffs.xc_f_1b = xc_f_1b;

    return;
end

function load_1b_tot_e_coeffs!(simulation::SimulationSystem)
    xc_a_1b, xc_b_1b, xc_c_1b, xc_d_1b, xc_e_1b, xc_f_1b =  
        load_1b_coeffs("Energy");

    simulation.tot_e_xc_coeffs.xc_a_1b = xc_a_1b;
    simulation.tot_e_xc_coeffs.xc_b_1b = xc_b_1b;
    simulation.tot_e_xc_coeffs.xc_c_1b = xc_c_1b;
    simulation.tot_e_xc_coeffs.xc_d_1b = xc_d_1b;
    simulation.tot_e_xc_coeffs.xc_e_1b = xc_e_1b;
    simulation.tot_e_xc_coeffs.xc_f_1b = xc_f_1b;

    return;
end

function load_2b_coeffs(type::String)
    base_dir = "XC Coeffs/"*type*"/Two Body/";
    xc_coeffs_a_id = base_dir*"xc_coeffs_A.txt";
    xc_coeffs_b_id = base_dir*"xc_coeffs_B.txt";
    xc_coeffs_c_id = base_dir*"xc_coeffs_C.txt";
    xc_coeffs_d_id = base_dir*"xc_coeffs_D.txt";

    function read_xc_file(file_name::String)
        lines = readlines(file_name);
        coeffs = Dict{Tuple{Int,Int},XCCoefficient}();

        for i in eachindex(lines)
            line_splitted = string.(split(lines[i]));

            elem_symbol_1 = line_splitted[end-1];
            elem_symbol_2 = line_splitted[end-0];
            atomic_number_1 = get_atomic_number(elem_symbol_1);
            atomic_number_2 = get_atomic_number(elem_symbol_2);
            dict_key_1 = (atomic_number_1,atomic_number_2);
            dict_key_2 = (atomic_number_2,atomic_number_1);

            m = parse(Float64,line_splitted[1]);
            b = parse(Float64,line_splitted[2]);

            coeff = XCCoefficient(m,b);
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

function load_2b_pol_e_coeffs!(simulation::SimulationSystem)
    xc_a_2b, xc_b_2b, xc_c_2b, xc_d_2b = load_2b_coeffs("Polarization");

    simulation.pol_e_xc_coeffs.xc_a_2b = xc_a_2b;
    simulation.pol_e_xc_coeffs.xc_b_2b = xc_b_2b;
    simulation.pol_e_xc_coeffs.xc_c_2b = xc_c_2b;
    simulation.pol_e_xc_coeffs.xc_d_2b = xc_d_2b;

    return;
end

function load_2b_tot_e_coeffs!(simulation::SimulationSystem)
    xc_a_2b, xc_b_2b, xc_c_2b, xc_d_2b = load_2b_coeffs("Energy");

    simulation.tot_e_xc_coeffs.xc_a_2b = xc_a_2b;
    simulation.tot_e_xc_coeffs.xc_b_2b = xc_b_2b;
    simulation.tot_e_xc_coeffs.xc_c_2b = xc_c_2b;
    simulation.tot_e_xc_coeffs.xc_d_2b = xc_d_2b;

    return;
end

function load_fitted_coeffs!(simulation::SimulationSystem)
    load_1b_pol_e_coeffs!(simulation);
    load_2b_pol_e_coeffs!(simulation);

    load_1b_tot_e_coeffs!(simulation);
    load_2b_tot_e_coeffs!(simulation);

    simulation.pol_e_xc_coeffs.max_atomic_number = 
        length(simulation.pol_e_xc_coeffs.xc_a_1b);
    simulation.tot_e_xc_coeffs.max_atomic_number = 
        length(simulation.tot_e_xc_coeffs.xc_a_1b);

    return;
end

function save_fitted_coeffs(simulation::SimulationSystem)
    save_1b_pol_e_coeffs(simulation);
    save_2b_pol_e_coeffs(simulation);

    save_1b_tot_e_coeffs(simulation);
    save_2b_tot_e_coeffs(simulation);

    return;
end

function get_xc_coeff(d::Number, xc_coeff::XCCoefficient)
    # Calculates the XC coefficient for the interatomic separation d in Bohr.
    return  xc_coeff.m * exp(d * xc_coeff.b);
end

function reset_xc_coeffs()
    max_atomic_number = 18;

    xc_a_1b = Dict{Int, XCCoefficient}();
    xc_b_1b = Dict{Int, XCCoefficient}();
    xc_c_1b = Dict{Int, XCCoefficient}();
    xc_d_1b = Dict{Int, XCCoefficient}();
    xc_e_1b = Dict{Int, Number}();
    xc_f_1b = Dict{Int, Number}();

    xc_a_2b = Dict{Tuple{Int,Int}, XCCoefficient}();
    xc_b_2b = Dict{Tuple{Int,Int}, XCCoefficient}();
    xc_c_2b = Dict{Tuple{Int,Int}, XCCoefficient}();
    xc_d_2b = Dict{Tuple{Int,Int}, XCCoefficient}();

    default_xc_coeffs = EmpiricalXCCoefficients(max_atomic_number,
        xc_a_1b,xc_b_1b,xc_c_1b,xc_d_1b,xc_e_1b,xc_f_1b,
        xc_a_2b, xc_b_2b, xc_c_2b, xc_d_2b);

    save_1b_coeffs(default_xc_coeffs,"Polarization");
    save_1b_coeffs(default_xc_coeffs,"Energy");

    save_2b_coeffs(default_xc_coeffs,"Polarization");
    save_2b_coeffs(default_xc_coeffs,"Energy");

    return;
end
