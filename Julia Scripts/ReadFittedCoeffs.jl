using Printf;

include("ForceFieldBase.jl")

function GenericSaveCoeffs(all_coeffs::Vector{Matrix{Float64}}, 
    one_body_file_names::Vector{String}, two_body_file_names::Vector{String})
    global num_1b_coeffs, num_2b_coeffs, num_coeffs;
    global num_elems, num_2b_combs;

    coeffs_1b = all_coeffs[1];
    for i in eachindex(one_body_file_names)
        file_name = one_body_file_names[i];
        fileID = open(file_name,"w");

        print_coeffs = coeffs_1b[i:num_1b_coeffs:end,:];
        num_rows = size(print_coeffs)[1];
        num_cols = size(print_coeffs)[2];
        
        for ii in 1:num_rows
            for jj in 1:num_cols
                coeff = print_coeffs[ii,jj];
                write(fileID, @sprintf "%18.8E " coeff);
            end

            # Print the key in the vector
            write(fileID, @sprintf "%5s " getElementSymbol(ii));

            write(fileID, "\n");
        end

        close(fileID);
    end

    coeffs_2b = all_coeffs[2];
    for i in eachindex(two_body_file_names)
        file_name = two_body_file_names[i];
        fileID = open(file_name,"w");

        print_coeffs = coeffs_2b[i:num_2b_coeffs:end,:];
        num_rows = size(print_coeffs)[1];
        num_cols = size(print_coeffs)[2];

        for ii in 1:num_rows
            for jj in 1:num_cols
                coeff = print_coeffs[ii,jj];
                if occursin("force_coeffs", file_name)
                    write(fileID, @sprintf "%18.8E " abs(coeff));
                else
                    write(fileID, @sprintf "%18.8E " coeff);
                end
            end

            # Print the key in the hash-map
            which_key = inv_coeff_map_2b[ii];
            write(fileID, @sprintf "%5s " getElementSymbol(which_key[1]));
            write(fileID, @sprintf "%5s " getElementSymbol(which_key[2]));

            write(fileID, "\n");
        end

        close(fileID);
    end
end

function SaveCoeffs()
    global num_1b_coeffs, num_2b_coeffs, num_coeffs;
    global num_elems, num_2b_combs;
    global all_coeffs;

    one_body_file_names = [
        "XC Coeffs/Energy/One Body/xc_coeffs_A.txt",
        "XC Coeffs/Energy/One Body/xc_coeffs_C.txt",
        "XC Coeffs/Energy/One Body/xc_coeffs_E.txt",
        "XC Coeffs/Energy/One Body/xc_coeffs_F.txt"
    ];

    two_body_file_names = [
        "XC Coeffs/Energy/Two Body/xc_coeffs_A.txt",
        "XC Coeffs/Energy/Two Body/xc_coeffs_B.txt",
        "XC Coeffs/Energy/Two Body/xc_coeffs_C.txt",
        "XC Coeffs/Energy/Two Body/xc_coeffs_D.txt"
    ];

    GenericSaveCoeffs(all_coeffs, one_body_file_names, two_body_file_names);
end

function SavePolCoeffs()
    global num_1b_coeffs, num_2b_coeffs, num_coeffs;
    global num_elems, num_2b_combs;
    global pol_all_coeffs;

    one_body_file_names = [
        "XC Coeffs/Polarization/One Body/xc_coeffs_A.txt",
        "XC Coeffs/Polarization/One Body/xc_coeffs_C.txt",
        "XC Coeffs/Polarization/One Body/xc_coeffs_E.txt",
        "XC Coeffs/Polarization/One Body/xc_coeffs_F.txt"
    ];

    two_body_file_names = [
        "XC Coeffs/Polarization/Two Body/xc_coeffs_A.txt",
        "XC Coeffs/Polarization/Two Body/xc_coeffs_B.txt",
        "XC Coeffs/Polarization/Two Body/xc_coeffs_C.txt",
        "XC Coeffs/Polarization/Two Body/xc_coeffs_D.txt"
    ];

    GenericSaveCoeffs(pol_all_coeffs, one_body_file_names, two_body_file_names);
end

function SaveBothCoeffs()
    SaveCoeffs();
    SavePolCoeffs();
end

function GenericLoadCoeffs(one_body_file_names::Vector{String},
    two_body_file_names::Vector{String})
    global num_1b_coeffs, num_2b_coeffs;
    global num_elems, num_2b_combs;

    coeffs_1b = zeros(Float64, num_1b_coeffs*num_elems, 2);
    for i in eachindex(one_body_file_names)
        file_name = one_body_file_names[i];

        if !isfile(file_name)
            continue;
        end

        fileID = open(file_name,"r");
        lines = readlines(fileID);
        coeffs = zeros(Float64, length(lines),2);
        for i in eachindex(lines)
            line_splitted = split(lines[i]);
            for j in 1:2
                coeffs[i,j] = parse(Float64, line_splitted[j]);
            end
        end
        close(fileID);

        coeffs_1b[i:num_1b_coeffs:end,:] = coeffs;
    end

    coeffs_2b = zeros(Float64, num_2b_coeffs*num_2b_combs, 3);
    for i in eachindex(two_body_file_names)
        file_name = two_body_file_names[i];

        if !isfile(file_name)
            continue;
        end

        fileID = open(file_name,"r");
        lines = readlines(fileID);
        coeffs = zeros(Float64, length(lines), 3);
        for i in eachindex(lines)
            line_splitted = split(lines[i]);
            for j in 1:3
                coeffs[i,j] = parse(Float64, line_splitted[j]);
            end
        end
        close(fileID);

        coeffs_2b[i:num_2b_coeffs:end,:] = coeffs;
    end

    return [coeffs_1b, coeffs_2b];
end

function LoadCoeffs()
    global num_1b_coeffs, num_2b_coeffs;
    global num_elems, num_2b_combs;

    one_body_file_names = [
        "XC Coeffs/Energy/One Body/xc_coeffs_A.txt",
        "XC Coeffs/Energy/One Body/xc_coeffs_C.txt",
        "XC Coeffs/Energy/One Body/xc_coeffs_E.txt",
        "XC Coeffs/Energy/One Body/xc_coeffs_F.txt"
    ];

    two_body_file_names = [
        "XC Coeffs/Energy/Two Body/xc_coeffs_A.txt",
        "XC Coeffs/Energy/Two Body/xc_coeffs_B.txt",
        "XC Coeffs/Energy/Two Body/xc_coeffs_C.txt",
        "XC Coeffs/Energy/Two Body/xc_coeffs_D.txt"
    ];

    return GenericLoadCoeffs(one_body_file_names, two_body_file_names);
end

function LoadPolCoeffs()
    global num_1b_coeffs, num_2b_coeffs;
    global num_elems, num_2b_combs;

    one_body_file_names = [
        "XC Coeffs/Polarization/One Body/xc_coeffs_A.txt",
        "XC Coeffs/Polarization/One Body/xc_coeffs_C.txt",
        "XC Coeffs/Polarization/One Body/xc_coeffs_E.txt",
        "XC Coeffs/Polarization/One Body/xc_coeffs_F.txt"
    ];

    two_body_file_names = [
        "XC Coeffs/Polarization/Two Body/xc_coeffs_A.txt",
        "XC Coeffs/Polarization/Two Body/xc_coeffs_B.txt",
        "XC Coeffs/Polarization/Two Body/xc_coeffs_C.txt",
        "XC Coeffs/Polarization/Two Body/xc_coeffs_D.txt"
    ];

    return GenericLoadCoeffs(one_body_file_names, two_body_file_names);
end

function ResetCoeffs()
    global num_elems, num_1b_combs, num_2b_combs;
    global all_coeffs = Vector{Matrix{Float64}}();
    resize!(all_coeffs,2);

    num_1b_combs = num_elems;
    num_2b_combs = round(Int,num_elems*(num_elems+1)/2.0);

    all_coeffs[1] = zeros(Float64,num_1b_coeffs*num_1b_combs,2);
    all_coeffs[2] = zeros(Float64,num_2b_coeffs*num_2b_combs,3);

    SaveCoeffs();
end

function ResetPolCoeffs()
    global num_elems, num_1b_combs, num_2b_combs;
    global pol_all_coeffs = Vector{Matrix{Float64}}();
    resize!(pol_all_coeffs,2);

    num_1b_combs = num_elems;
    num_2b_combs = round(Int,num_elems*(num_elems+1)/2.0);

    pol_all_coeffs[1] = zeros(Float64,num_1b_coeffs*num_1b_combs,2);
    pol_all_coeffs[2] = zeros(Float64,num_2b_coeffs*num_2b_combs,3);

    SavePolCoeffs();
end

function ResetBothCoeffs()
    ResetCoeffs();
    ResetPolCoeffs();
end

function LoadBothCoeffs()
    all_coeffs = LoadCoeffs();
    pol_all_coeffs = LoadPolCoeffs();
    return all_coeffs, pol_all_coeffs;
end
