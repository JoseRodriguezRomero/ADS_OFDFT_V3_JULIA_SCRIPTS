using Optim, Plots, Evolutionary, Distributed;
using Printf, LinearAlgebra;
using LaTeXStrings, Latexify, Measures;
using Base.Threads;

include("FitHomonuclearCoeffs_General.jl")

which_atomic_numbers = [[1,6],[1,8],[6,8],[7,8]];
function FitCoeffsFoo(Z1::Integer, Z2::Integer)
    global pol_all_coeffs;
    neutral_data, cation_data, anion_data = readAllSanitizedData(Z1,Z2);
    all_data = vcat(neutral_data, cation_data, anion_data);

    num_vars = 3*num_2b_coeffs;
    aux_X = zeros(Float64,num_vars);

    all_molecs = Vector{Molecule}();
    resize!(all_molecs,length(all_data));
    for i in eachindex(all_data)
        all_molecs[i] = MakeMoleculeFromParsedFile(all_data[i]);
    end

    function CostFunc(aux_X::Vector)
        global pol_all_coeffs;
        aux_type = typeof(aux_X[1]);
        aux_pol_all_coeffs = Vector{Matrix}();
        resize!(aux_pol_all_coeffs,2);

        aux_pol_all_coeffs[1] = aux_type.(pol_all_coeffs[1]);
        aux_pol_all_coeffs[2] = aux_type.(pol_all_coeffs[2]);
        SetFittedCoeffs!(aux_pol_all_coeffs,Z1,Z2,aux_X);

        Z1_eff = all_molecs[1].atoms_data[1,4];
        Z2_eff = all_molecs[1].atoms_data[2,4];

        n_threads = Threads.nthreads();
        ret_val = zeros(aux_type,n_threads);
        @threads for thread_id in 1:n_threads
            for i in thread_id:n_threads:(length(all_molecs))
                _, _, aux_m, aux_y = 
                    MatPolarizeMolecules(all_molecs[i],aux_pol_all_coeffs);

                aux_x = zeros(aux_type,3);
                aux_x[1] = all_molecs[i].cloud_data[1,6];
                aux_x[2] = all_molecs[i].cloud_data[4,6];
                aux_x[3] = all_molecs[i].chem_μ;

                aux_diff = (aux_m * aux_x) - aux_y;
                ret_val[thread_id] += dot(aux_diff,aux_diff);
            end
        end

        return sum(ret_val);
    end

    sol = Optim.optimize(CostFunc, aux_X[:], LBFGS(), autodiff=:forward,
        Optim.Options(show_trace=true,iterations=8000));
    aux_X = Optim.minimizer(sol);

    # Save the currently fitted coefficients
    SetFittedCoeffs!(pol_all_coeffs,Z2,Z1,aux_X);

    # Save the newly fitted coefficients
    SavePolCoeffs();
end

for atomic_numbers in which_atomic_numbers
    Z1 = atomic_numbers[1];
    Z2 = atomic_numbers[2];

    # Load the previously fitted coefficients
    global pol_all_coeffs = LoadPolCoeffs();

    FitCoeffsFoo(Z1,Z2);
end

function TestResultChemμ(Z1::Integer, Z2::Integer)
    pol_all_coeffs = LoadPolCoeffs();
    neutral_data, cation_data, anion_data = readAllSanitizedData(Z1,Z2);
    all_data = vcat(neutral_data, cation_data, anion_data);

    model_chem_μ = zeros(Float64,length(all_data));
    dft_chem_μ = zeros(Float64,length(all_data));

    function GetData(data::Vector{ParsedFile})
        model_chem_μ = zeros(Float64,length(data));
        dft_chem_μ = zeros(Float64,length(data));

        n_threads = Threads.nthreads();
        @threads for thread_id in 1:n_threads
            for i in thread_id:n_threads:length(data)
                mol = MakeMoleculeFromParsedFile(data[i]);
                dft_chem_μ[i] = mol.chem_μ;
            
                PolarizeMolecules!(mol,pol_all_coeffs);
                model_chem_μ[i] = mol.chem_μ;
            end
        end

        return model_chem_μ, dft_chem_μ;
    end

    neutral_model_chem_μ, neutral_dft_chem_μ = GetData(neutral_data);
    cation_model_chem_μ, cation_dft_chem_μ = GetData(cation_data);
    anion_model_chem_μ, anion_dft_chem_μ = GetData(anion_data);

    model_chem_μ = vcat(neutral_model_chem_μ,cation_model_chem_μ,anion_model_chem_μ);
    dft_chem_μ = vcat(neutral_dft_chem_μ,cation_dft_chem_μ,anion_dft_chem_μ);

    mean_dft_chem_μ = sum(dft_chem_μ) / length(dft_chem_μ);
    SS_res = sum((model_chem_μ - dft_chem_μ).^2.0);
    SS_tot = sum((dft_chem_μ .- mean_dft_chem_μ).^2.0);
    R² = 1.0 - SS_res / SS_tot;

    min_x = minimum([minimum(model_chem_μ), minimum(dft_chem_μ)]);
    max_x = maximum([maximum(model_chem_μ), maximum(dft_chem_μ)]);

    elem_symbol_1 = getElementSymbol(Z1);
    elem_symbol_2 = getElementSymbol(Z2);

    neutral_label = elem_symbol_1*" + "*elem_symbol_2;
    cation_label = "("*elem_symbol_1*" + "*elem_symbol_2*")⁺";
    anion_label = "("*elem_symbol_1*" + "*elem_symbol_2*")⁻";

    Pts = [-5,5];
    p = plot(Pts,Pts,label=false);
    scatter!(neutral_dft_chem_μ, neutral_model_chem_μ, label=neutral_label);
    scatter!(cation_dft_chem_μ, cation_model_chem_μ, label=cation_label);
    scatter!(anion_dft_chem_μ, anion_model_chem_μ, label=anion_label);
    plot!(framestyle = :box);
    plot!(legend=:topleft);

    return p, R²;
end

function CompareChemμ()
    y_label_all = L"$\tilde{\mu} \quad \mathrm{(This \ Work)}$";
    x_label_all = L"$\tilde{\mu} \quad (\mathrm{KS{-}DFT})$";

    R²_rel_pos_x = (5.0/8.0);
    R²_rel_pos_y = (1.0/8.0);

    # HC Plots
    pHC, R² = TestResultChemμ(1,6);
    aux_lims = [-1.0,0.6];
    aux_ticks = -1.0:0.4:0.6;
    plot!(xlims=aux_lims);
    plot!(ylims=aux_lims);
    plot!(xticks=aux_ticks);
    plot!(yticks=aux_ticks);
    plot!(ylabel=y_label_all);
    plot!(xlabel=x_label_all);
    plot!(left_margin=6mm, bottom_margin=8mm);

    l_x_pos = aux_lims[1] + R²_rel_pos_x*(aux_lims[2] - aux_lims[1]);
    l_y_pos = aux_lims[1] + R²_rel_pos_y*(aux_lims[2] - aux_lims[1]);
    R² = round(R²,digits=5);
    annotate!(l_x_pos, l_y_pos, text("R² = "*string(R²), :center, 10));

    # HO Plots
    pHO, R² = TestResultChemμ(1,8);
    aux_lims = [-1.0,0.6];
    aux_ticks = -1.0:0.4:0.6;
    plot!(xlims=aux_lims);
    plot!(ylims=aux_lims);
    plot!(xticks=aux_ticks);
    plot!(yticks=aux_ticks);
    plot!(xlabel=x_label_all);

    l_x_pos = aux_lims[1] + R²_rel_pos_x*(aux_lims[2] - aux_lims[1]);
    l_y_pos = aux_lims[1] + R²_rel_pos_y*(aux_lims[2] - aux_lims[1]);
    R² = round(R²,digits=5);
    annotate!(l_x_pos, l_y_pos, text("R² = "*string(R²), :center, 10));

    # CO Plots
    pCO, R² = TestResultChemμ(6,8);
    aux_lims = [-1.0,0.6];
    aux_ticks = -1.0:0.4:0.6;
    plot!(xlims=aux_lims);
    plot!(ylims=aux_lims);
    plot!(xticks=aux_ticks);
    plot!(yticks=aux_ticks);
    plot!(xlabel=x_label_all);

    l_x_pos = aux_lims[1] + R²_rel_pos_x*(aux_lims[2] - aux_lims[1]);
    l_y_pos = aux_lims[1] + R²_rel_pos_y*(aux_lims[2] - aux_lims[1]);
    R² = round(R²,digits=5);
    annotate!(l_x_pos, l_y_pos, text("R² = "*string(R²), :center, 10));

    # NO Plots
    pNO, R² = TestResultChemμ(7,8);
    aux_lims = [-1.0,0.6];
    aux_ticks = -1.0:0.4:0.6;
    plot!(xlims=aux_lims);
    plot!(ylims=aux_lims);
    plot!(xticks=aux_ticks);
    plot!(yticks=aux_ticks);
    plot!(xlabel=x_label_all);

    l_x_pos = aux_lims[1] + R²_rel_pos_x*(aux_lims[2] - aux_lims[1]);
    l_y_pos = aux_lims[1] + R²_rel_pos_y*(aux_lims[2] - aux_lims[1]);
    R² = round(R²,digits=5);
    annotate!(l_x_pos, l_y_pos, text("R² = "*string(R²), :center, 10));

    p = plot(pHC,pHO, pCO, pNO, layout=(1,4), size = (1100,235));
    savefig("HeteronuclearFitComps.pdf");

    return p;
end

function Foo(Z1::Integer, Z2::Integer)
    pol_all_coeffs = LoadPolCoeffs();
    neutral_data, cation_data, anion_data = readAllSanitizedData(Z1,Z2);

    function GetPlotData(data::Vector{ParsedFile})
        atoms_dist = zeros(Float64,length(data));
        charges_dft = zeros(Float64,length(data));
        charges_model = zeros(Float64,length(data));

        n_threads = Threads.nthreads();
        @threads for thread_id in 1:n_threads
            mol = MakeMoleculeFromParsedFile(data[1]);
            Z1_eff = mol.atoms_data[1,4];

            for i in thread_id:n_threads:length(data)
                mol = MakeMoleculeFromParsedFile(data[i]);
                # charges_dft[i] = Z1_eff * (1 - mol.cloud_data[1,6]);
                charges_dft[i] = mol.chem_μ;

                mol.cloud_data[:,6] .= 1.0;
                mol.chem_μ = 0.0;

                PolarizeMolecules!(mol,pol_all_coeffs);
                # charges_model[i] = Z1_eff * (1 - mol.cloud_data[1,6]);
                charges_model[i] = mol.chem_μ;

                atoms_dist[i] = mol.atoms_data[2,3];
            end
        end

        return atoms_dist, charges_dft, charges_model;
    end
    
    neutral_atoms_dist, neutral_charges_dft, neutral_charges_model = 
        GetPlotData(neutral_data);
    cation_atoms_dist, cation_charges_dft, cation_charges_model = 
        GetPlotData(cation_data);
    anion_atoms_dist, anion_charges_dft, anion_charges_model = 
        GetPlotData(anion_data);

    p1 = plot(neutral_atoms_dist, neutral_charges_model, 
        label="This Work");
    plot!(neutral_atoms_dist, neutral_charges_dft,
        label="B3LYP/aug-cc-pVTZ");
    plot!(title="neutral");

    p2 = plot(cation_atoms_dist, cation_charges_model, 
        label="This Work");
    plot!(cation_atoms_dist, cation_charges_dft,
        label="B3LYP/aug-cc-pVTZ");
    plot!(title="anion");

    p3 = plot(anion_atoms_dist, anion_charges_model, 
        label="This Work");
    plot!(anion_atoms_dist, anion_charges_dft,
        label="B3LYP/aug-cc-pVTZ");
    plot!(title="anion");

    p = plot(p1,p2,p3,layout=(3,1), size = (450,550));
    return p;
end

CompareChemμ()
