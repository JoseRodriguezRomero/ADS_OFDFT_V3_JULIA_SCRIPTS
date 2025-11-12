using Plots, GridInterpolations, Glob;
using Optim, Interpolations;
using Base.Threads;
using Printf;

function read_gaussian_cube_file(file_name::String)
    # Reads the Gaussian formatted cubefile and returns a three dimensional 
    # meshes for the Cartesian coordinates and theirs respective values.
    lines = readlines(file_name);

    line_splitted = split(lines[3]);
    x0 = parse(Float64,line_splitted[2]);
    y0 = parse(Float64,line_splitted[3]);
    z0 = parse(Float64,line_splitted[4]);

    function parse_vector_line(line::String)
        local line_splitted = split(line);
        N = parse(Int,line_splitted[1]);

        dr = zeros(Float64,3);
        dr[1] = parse(Float64,line_splitted[2]);
        dr[2] = parse(Float64,line_splitted[3]);
        dr[3] = parse(Float64,line_splitted[4]);

        return N, dr;
    end

    Ni, dri = parse_vector_line(lines[4]);
    Nj, drj = parse_vector_line(lines[5]);
    Nk, drk = parse_vector_line(lines[6]);

    U = zeros(Float64,Ni,Nj,Nk);
    X = zeros(Float64,Ni,Nj,Nk) .+ x0;
    Y = zeros(Float64,Ni,Nj,Nk) .+ y0;
    Z = zeros(Float64,Ni,Nj,Nk) .+ z0;

    i = 0;
    for line in lines[8:end]
        line_splitted = split(line);
        for j in eachindex(line_splitted)
            U[i+j] = parse(Float64,line_splitted[j]);
        end

        i += length(line_splitted);
    end

    for i in 1:Ni
        for j in 1:Nj
            for k in 1:Nk
                X[i,j,k] += (i-1)*dri[1];
                Y[i,j,k] += (i-1)*dri[2];
                Z[i,j,k] += (i-1)*dri[3];

                X[i,j,k] += (j-1)*drj[1];
                Y[i,j,k] += (j-1)*drj[2];
                Z[i,j,k] += (j-1)*drj[3];

                X[i,j,k] += (k-1)*drk[1];
                Y[i,j,k] += (k-1)*drk[2];
                Z[i,j,k] += (k-1)*drk[3];
            end
        end
    end

    return U, X, Y, Z;
end

function element_symbol(Z::Int)
    symbols = [
        "H",  "He",
        "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne",
        "Na", "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar"
    ];

    return symbols[Z];
end

function get_element_file_names(Z::Int)
    base_dir = "../Psi4 Files/" * element_symbol(Z) * "/";
    files_x = glob("*.cube",base_dir*"x/");
    files_y = glob("*.cube",base_dir*"y/");
    files_z = glob("*.cube",base_dir*"z/");
    files = vcat(files_x,files_y,files_z);

    row_1_a_names = ["_a_1_"];
    row_1_b_names = ["_b_1_"];
    row_1_names = vcat(row_1_a_names,row_1_b_names);

    row_2_a_names = ["_a_2_","_a_3_","_a_4_","_a_5_"];
    row_2_b_names = ["_b_2_","_b_3_","_b_4_","_b_5_"];
    row_2_names = vcat(row_2_a_names,row_2_b_names);

    row_3_a_names = ["_a_6_","_a_7_","_a_8_","_a_9_"];
    row_3_b_names = ["_b_6_","_b_7_","_b_8_","_b_9_"];
    row_3_names = vcat(row_3_a_names,row_3_b_names);

    all_shell_file_names = [row_1_names, row_2_names, row_3_names];

    file_names = Vector{Vector{String}}();
    resize!(file_names,3)
    for i in 1:3
        file_names[i] = Vector{String}();
    end

    for file_name in files
        for i in eachindex(all_shell_file_names)
            shell_file_names = all_shell_file_names[i];

            for shell_file_name in shell_file_names
                if contains(file_name, shell_file_name)
                    file_names[i] = vcat(file_names[i], file_name);
                end
            end
        end
    end

    return file_names;
end

function parse_shell_densities(Z::Int)
    file_names = get_element_file_names(Z);

    _, X, Y, Z = read_gaussian_cube_file(file_names[1][1]);
    num_rows = length(X[:]);
    I0 = ceil(Int, num_rows/2.0);

    min_x = minimum(X[:]);
    max_x = maximum(X[:]);

    min_y = minimum(Y[:]);
    max_y = maximum(Y[:]);

    min_z = minimum(Z[:]);
    max_z = maximum(Z[:]);

    r0 = minimum([min_x, min_y, min_z]);
    r1 = maximum([max_x, max_y, max_z]);

    dr = (r1 - r0) / (num_rows - 1);

    radial_densities = zeros(Float64,num_rows, 4);
    radial_densities[:,1] = collect(r0:dr:r1);

    for i in eachindex(file_names)
        for file_name in file_names[i]
            U, _, _, _ = read_gaussian_cube_file(file_name);

            radial_densities[:,i+1] += (1.0/3.0) .* U[:] .^ 2;
        end
    end

    return radial_densities[I0:end,:];
end

function num_gaussians()
    return factorial(num_fit_gaussians());
end

function num_fit_gaussians()
    return 3;
end

function fit_gaussians(data_x::Vector, data_y::Vector, num_electrons::Int)
    # Fits three s-type Gaussian basis-functions to the data.
    
    peak_value = maximum((data_x .^ 2) .* data_y);

    for i in eachindex(data_x)
        comp_value = (data_x[i] ^ 2) * data_y[i];
        if comp_value > 0.25 * peak_value;
            data_y = data_y[i:end];
            data_x = data_x[i:end];
            break;
        end
    end

    function cost_func(params::Vector)
        λ = abs.(params[1:2:end]);
        c = params[2:2:end];

        n_threads = Threads.nthreads();
        ret_val = zeros(typeof(params[1]),n_threads);
        @threads for thread_id in 1:n_threads
            for i in thread_id:n_threads:length(data_x)
                r = data_x[i];
                ρ = data_y[i];

                model_ρ = 0;
                for k in eachindex(c)
                    model_ρ += c[k] * exp(- λ[k] * r^2);
                end

                model_ρ = model_ρ ^ 2;
                model_diff = (r^2)*(ρ - model_ρ);
                ret_val[thread_id] += model_diff^2;
            end
        end

        ret_val = sum(ret_val);
        ret_val /= length(data_x);

        aux_int = 0.0;
        for i in eachindex(c)
            for j in eachindex(c)
                int_ij = c[i] * c[j];
                int_ij *= (π / (λ[i] + λ[j])) ^ (3/2);
                aux_int += int_ij;
            end
        end

        ret_val += (aux_int - num_electrons)^2;
        return ret_val;
    end

    params_0 = 5.0 .* (rand(Float64,2*num_fit_gaussians()) .- 0.5);

    sol = Optim.optimize(cost_func,params_0, NelderMead(),
        Optim.Options(show_trace=true));
    params_0 = Optim.minimizer(sol);

    sol = Optim.optimize(cost_func,params_0, LBFGS(), autodiff = :forward,
        Optim.Options(show_trace=true));
    params_0 = Optim.minimizer(sol);

    old_params_0 = copy(params_0);
    for _ in 1:20
        params_0 = 5.0 .* (rand(Float64,2*num_fit_gaussians()) .- 0.5);

        sol = Optim.optimize(cost_func,params_0, NelderMead(),
            Optim.Options(show_trace=true));
        params_0 = Optim.minimizer(sol);

        sol = Optim.optimize(cost_func,params_0, LBFGS(), autodiff = :forward,
            Optim.Options(show_trace=true));
        params_0 = Optim.minimizer(sol);

        if cost_func(params_0) < cost_func(old_params_0)
            old_params_0 = copy(params_0);
        end
    end

    A = abs.(old_params_0[1:2:end]);
    B = old_params_0[2:2:end];
    ret_params = zeros(Float64,2*num_gaussians());

    ret_params[1]  = 2 * A[1];
    ret_params[3]  = 2 * A[2];
    ret_params[5]  = 2 * A[3];
    ret_params[7]  = A[1] + A[2];
    ret_params[9]  = A[1] + A[3];
    ret_params[11] = A[2] + A[3];

    ret_params[2]  = (B[1] ^ 2) * ((π / (2 * A[1])) ^ (3/2));
    ret_params[4]  = (B[2] ^ 2) * ((π / (2 * A[2])) ^ (3/2));
    ret_params[6]  = (B[3] ^ 2) * ((π / (2 * A[3])) ^ (3/2));
    ret_params[8]  = 2 * (B[1] * B[2]) * ((π / (A[1] + A[2])) ^ (3/2));
    ret_params[10] = 2 * (B[1] * B[3]) * ((π / (A[1] + A[3])) ^ (3/2));
    ret_params[12] = 2 * (B[2] * B[3]) * ((π / (A[2] + A[3])) ^ (3/2));

    return ret_params;
end

function eval_fitted_gaussians(data_x::Vector, gaussian_params::Vector)
    λ = gaussian_params[1:2:end];
    c = gaussian_params[2:2:end];

    model_ρ = zeros(Float64,length(data_x));

    for i in eachindex(data_x)
        ρ = 0;
        r = data_x[i];
        for k in eachindex(c)
            ρ += c[k] * ((λ[k]/π)^(3.0/2.0)) * exp(-λ[k]*(r^2));
        end

        model_ρ[i] = ρ;
    end

    return model_ρ;
end

function label_power(power::Int)
    if power == 1
        return "¹";
    elseif power == 2
        return "²";
    elseif power == 3
        return "³";
    elseif power == 4
        return "⁴";
    elseif power == 5
        return "⁵";
    elseif power == 6
        return "⁶";
    end
end

function save_fitted_parameters(file_name::String, params::Vector)
    λ = params[1:2:end];
    c = params[2:2:end];

    file_id = open(file_name,"w");
    for i in eachindex(λ)
        write(file_id, @sprintf "%18.8E " λ[i]);
        write(file_id, @sprintf "%18.8E " c[i]);
        write(file_id, "\n")
    end

    close(file_id);
end

function plot_density(Z::Int)
    densities = parse_shell_densities(Z);

    r = densities[:,1];
    shell_1 = densities[:,2];
    shell_2 = densities[:,3];
    shell_3 = densities[:,4];

    p = scatter(r,(r.^2).*shell_1, marker=true, color=palette(:auto)[1], 
        markersize=2);
    model_r0 = -3.0;
    model_r1 = 3.0;
    model_dr = (model_r1-model_r0)/1000.0;
    model_r = collect(model_r0:model_dr:model_r1);
    model_r = 10.0 .^ model_r;

    shell_1_fit_sum = 0.0;
    shell_2_fit_sum = 0.0;
    shell_3_fit_sum = 0.0;

    if sum(shell_1) > 1.0E-6
        num_electrons = minimum([Z,2]);
        fitted_gaussians = fit_gaussians(r,shell_1,num_electrons);
        shell_1_fit_sum = sum(fitted_gaussians[2:2:end]);

        file_name = "../Psi4 Files/"*element_symbol(Z)*"/shell_1.txt";
        save_fitted_parameters(file_name,fitted_gaussians);

        model_shell_1 = eval_fitted_gaussians(model_r,fitted_gaussians);
        model_shell_1 .*= model_r .^ 2;

        plot_label = "1s"*label_power(num_electrons);

        p = scatter(r,(r.^2).*shell_1, marker=true, color=palette(:auto)[1], 
            markersize=2, label=plot_label*" (B3LYP/aug-cc-pVTZ)");
        plot!(model_r,model_shell_1, color=palette(:auto)[1], 
            label=plot_label*" (fitted gaussians)");
    end

    if sum(shell_2) > 1.0E-6
        num_electrons = minimum([8,Z-2]);
        fitted_gaussians = fit_gaussians(r,shell_2,num_electrons);
        shell_2_fit_sum = sum(fitted_gaussians[2:2:end]);

        file_name = "../Psi4 Files/"*element_symbol(Z)*"/shell_2.txt";
        save_fitted_parameters(file_name,fitted_gaussians);

        if num_electrons < 3
            plot_label = "2s"*label_power(num_electrons);
        else
            plot_label = "2s"*label_power(2)*" + ";
            plot_label *= "2p"*label_power(num_electrons-2);
        end

        model_shell_2 = eval_fitted_gaussians(model_r,fitted_gaussians);
        model_shell_2 .*= model_r .^ 2;

        scatter!(r,(r.^2).*shell_2, marker=true, color=palette(:auto)[2], 
            markersize=2, label=plot_label*" (B3LYP/aug-cc-pVTZ)");
        plot!(model_r,model_shell_2, color=palette(:auto)[2], 
            label=plot_label*" (fitted gaussians)");
    end

    if sum(shell_3) > 1.0E-6
        num_electrons = minimum([8,Z-10]);
        fitted_gaussians = fit_gaussians(r,shell_3,num_electrons);
        shell_3_fit_sum = sum(fitted_gaussians[2:2:end]);

        file_name = "../Psi4 Files/"*element_symbol(Z)*"/shell_3.txt";
        save_fitted_parameters(file_name,fitted_gaussians);

        if num_electrons < 3
            plot_label = "3s"*label_power(num_electrons);
        else
            plot_label = "3s"*label_power(2)*" + ";
            plot_label *= "3p"*label_power(num_electrons-2);
        end

        model_shell_3 = eval_fitted_gaussians(model_r,fitted_gaussians);
        model_shell_3 .*= model_r .^ 2;

        scatter!(r,(r.^2).*shell_3, marker=true, color=palette(:auto)[3], 
            markersize=2, label=plot_label*" (B3LYP/aug-cc-pVTZ)");
        plot!(model_r,model_shell_3, color=palette(:auto)[3], 
            label=plot_label*" (fitted gaussians)");
    end

    plot!(xaxis = :log10);
    plot!(xlims = [1.0E-3,1.0E3]);
    plot!(xlims = [1.0E-3,1.0E3]);

    plot!(framestyle=:box);

    plot!(ylabel="r² ρ(r)");
    plot!(xlabel="r [Bohr]");

    if Z == 1
        r = -3:0.01:3;
        r = 10.0.^r;

        function hydrogen_density(r)
            return exp(-2*r)/π;
        end

        plot!(r,(r.^2) .* hydrogen_density.(r), label="exact");
    end

    display([shell_1_fit_sum,shell_2_fit_sum,shell_3_fit_sum]);

    file_name = "../Psi4 Files/"*element_symbol(Z)*"/fitted_densities.pdf";
    savefig(file_name);

    return p;
end
