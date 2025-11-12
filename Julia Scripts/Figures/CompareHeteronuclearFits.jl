using LaTeXStrings, Latexify, Measures;
using Base.Threads;
using Plots;

include("../FitCoeffs_General.jl")

function test_result_ΔE(Z1::Int, Z2::Int)
    neutral_data, cation_data, anion_data = 
        read_all_sanitized_data(Z1,Z2);

    function get_data(data::Vector{ParsedFile})
        at1 = make_monoatomic_system(Z1,0);
        at2 = make_monoatomic_system(Z2,0);

        dft_atom_energies, _, _ = get_reference_atom_total_energy();

        at1_dft_energy = dft_atom_energies[Z1];
        at2_dft_energy = dft_atom_energies[Z2];
        dft_ref_energy = at1_dft_energy + at2_dft_energy;

        model_ref_energy = total_energy(at1) + total_energy(at2);

        dft_ΔE = zeros(Float64,length(data));
        model_ΔE = zeros(Float64,length(data));

        n_threads = Threads.nthreads();
        @threads for thread_id in 1:n_threads
            simulation = make_system_from_parsed_file(data[1]);

            for i in thread_id:n_threads:length(data)
                set_diatomic_system_to_parsed_file!(simulation,data[i]);

                dft_ΔE[i] = data[i].total_energy;
                dft_ΔE[i] -= dft_ref_energy;

                polarize_molecules!(simulation);
                model_ΔE[i] = total_energy(simulation);
                model_ΔE[i] -= model_ref_energy;
            end
        end

        return dft_ΔE, model_ΔE;
    end

    neutral_dft_ΔE, neutral_model_ΔE = get_data(neutral_data);
    cation_dft_ΔE, cation_model_ΔE = get_data(cation_data);
    anion_dft_ΔE, anion_model_ΔE = get_data(anion_data);

    model_ref_E = minimum(neutral_model_ΔE);
    neutral_model_ΔE .-= model_ref_E;
    cation_model_ΔE .-= model_ref_E;
    anion_model_ΔE .-= model_ref_E;

    dft_ref_E = minimum(neutral_dft_ΔE);
    neutral_dft_ΔE .-= dft_ref_E;
    cation_dft_ΔE .-= dft_ref_E;
    anion_dft_ΔE .-= dft_ref_E;

    model_ΔE = vcat(neutral_model_ΔE,cation_model_ΔE,anion_model_ΔE);
    dft_ΔE = vcat(neutral_dft_ΔE,cation_dft_ΔE,anion_dft_ΔE);

    mean_dft_ΔE = sum(dft_ΔE) / length(dft_ΔE);
    SS_res = sum((model_ΔE - dft_ΔE).^2.0);
    SS_tot = sum((dft_ΔE .- mean_dft_ΔE).^2.0);
    R² = 1.0 - SS_res / SS_tot;

    elem_symbol1 = get_element_symbol(Z1);
    elem_symbol2 = get_element_symbol(Z2);

    neutral_label = elem_symbol1*" + "*elem_symbol2;
    cation_label = "("*elem_symbol1*" + "*elem_symbol2*")⁺";
    anion_label = "("*elem_symbol1*" + "*elem_symbol2*")⁻";

    Pts = [-100,100];
    p = plot(Pts,Pts,label=false, color=palette(:default)[4]);
    scatter!(neutral_dft_ΔE, neutral_model_ΔE, 
        label=neutral_label, color=palette(:default)[1]);
    scatter!(cation_dft_ΔE, cation_model_ΔE, 
        label=cation_label, color=palette(:default)[2]);
    scatter!(anion_dft_ΔE, anion_model_ΔE, 
        label=anion_label, color=palette(:default)[3]);
    plot!(framestyle = :box);
    plot!(legend=:topleft);

    return p, R²;
end

function test_result_chemical_potential(Z1::Int, Z2::Int)
    neutral_data, cation_data, anion_data = 
        read_all_sanitized_data(Z1,Z2,true);

    function get_data(data::Vector{ParsedFile})
        model_chem_μ = zeros(Float64,length(data));
        dft_chem_μ = zeros(Float64,length(data));

        n_threads = Threads.nthreads();
        simulation = Vector{SimulationSystem}();
        resize!(simulation,n_threads);
        for thread_id in 1:n_threads
            simulation[thread_id] = make_system_from_parsed_file(data[1]);
        end

        @threads for thread_id in 1:n_threads
            for i in thread_id:n_threads:length(data)
                # dft_μ = (data[i].HOMO_energy + data[i].LUMO_energy)/2.0;
                dft_μ = data[i].chemical_potential;

                set_diatomic_system_to_parsed_file!(
                    simulation[thread_id],data[i]);
                polarize_molecules!(simulation[thread_id]);
                model_μ = simulation[thread_id].system.chemical_potential;

                model_chem_μ[i] = model_μ;
                dft_chem_μ[i] = dft_μ;
            end
        end

        return model_chem_μ, dft_chem_μ;
    end

    neutral_model_chem_μ, neutral_dft_chem_μ = get_data(neutral_data);
    cation_model_chem_μ, cation_dft_chem_μ = get_data(cation_data);
    anion_model_chem_μ, anion_dft_chem_μ = get_data(anion_data);

    model_chem_μ = vcat(
            neutral_model_chem_μ,cation_model_chem_μ,anion_model_chem_μ);
    dft_chem_μ = vcat(
            neutral_dft_chem_μ,cation_dft_chem_μ,anion_dft_chem_μ);

    mean_dft_chem_μ = sum(dft_chem_μ) / length(dft_chem_μ);
    SS_res = sum((model_chem_μ - dft_chem_μ).^2.0);
    SS_tot = sum((dft_chem_μ .- mean_dft_chem_μ).^2.0);
    R² = 1.0 - SS_res / SS_tot;

    elem_symbol1 = get_element_symbol(Z1);
    elem_symbol2 = get_element_symbol(Z2);

    neutral_label = elem_symbol1*" + "*elem_symbol2;
    cation_label = "("*elem_symbol1*" + "*elem_symbol2*")⁺";
    anion_label = "("*elem_symbol1*" + "*elem_symbol2*")⁻";

    Pts = [-2,2];
    p = plot(Pts,Pts,label=false, color=palette(:default)[4]);
    scatter!(neutral_dft_chem_μ, neutral_model_chem_μ, 
        label=neutral_label, color=palette(:default)[1]);
    scatter!(cation_dft_chem_μ, cation_model_chem_μ, 
        label=cation_label, color=palette(:default)[2]);
    scatter!(anion_dft_chem_μ, anion_model_chem_μ, 
        label=anion_label, color=palette(:default)[3]);
    plot!(framestyle = :box);
    plot!(legend=:topleft);

    return p, R²;
end

function compare_data()
    y_label_all = L"$\Delta E \quad \mathrm{(This \ Work)}$";
    x_label_all = L"$\Delta E \quad (\mathrm{KS{-}DFT})$";

    R²_rel_pos_x = (14.0/15.0);
    R²_rel_pos_y = (1.0/8.0);

    # Energy Differences
    # HC Plots
    eHC, R² = test_result_ΔE(1,6);
    aux_lims = [-0.3,1.7];
    aux_ticks = -0.3:0.5:1.7;
    plot!(xlims=aux_lims);
    plot!(ylims=aux_lims);
    plot!(xticks=aux_ticks);
    plot!(yticks=aux_ticks);
    plot!(ylabel=y_label_all);
    plot!(xlabel=x_label_all);
    plot!(left_margin=6mm, bottom_margin=8mm);

    l_x_pos = aux_lims[1] + R²_rel_pos_x*(aux_lims[2] - aux_lims[1]);
    l_y_pos = aux_lims[1] + R²_rel_pos_y*(aux_lims[2] - aux_lims[1]);
    R² = @sprintf "%.5lf" round(R²,digits=5);
    annotate!(l_x_pos, l_y_pos, text("R² = "*R², :right, 10));

    # HO Plots
    eHO, R² = test_result_ΔE(1,8);
    aux_lims = [-0.3,1.7];
    aux_ticks = -0.3:0.5:1.7;
    plot!(xlims=aux_lims);
    plot!(ylims=aux_lims);
    plot!(xticks=aux_ticks);
    plot!(xlabel=x_label_all);
    plot!(yticks=aux_ticks);

    l_x_pos = aux_lims[1] + R²_rel_pos_x*(aux_lims[2] - aux_lims[1]);
    l_y_pos = aux_lims[1] + R²_rel_pos_y*(aux_lims[2] - aux_lims[1]);
    R² = @sprintf "%.5lf" round(R²,digits=5);
    annotate!(l_x_pos, l_y_pos, text("R² = "*R², :right, 10));

    # CO Plots
    eCO, R² = test_result_ΔE(6,8);  
    aux_lims = [-0.4,2.4];
    aux_ticks = -0.4:0.7:2.4;
    plot!(xlims=aux_lims);
    plot!(ylims=aux_lims);
    plot!(xticks=aux_ticks);
    plot!(xlabel=x_label_all);
    plot!(yticks=aux_ticks);

    l_x_pos = aux_lims[1] + R²_rel_pos_x*(aux_lims[2] - aux_lims[1]);
    l_y_pos = aux_lims[1] + R²_rel_pos_y*(aux_lims[2] - aux_lims[1]);
    R² = @sprintf "%.5lf" round(R²,digits=5);
    annotate!(l_x_pos, l_y_pos, text("R² = "*R², :right, 10));

    # NO Plots
    eNO, R² = test_result_ΔE(7,8);
    aux_lims = [-0.2,1.0];
    aux_ticks = -0.2:0.3:1.0;
    plot!(xlims=aux_lims);
    plot!(ylims=aux_lims);
    plot!(xticks=aux_ticks);
    plot!(xlabel=x_label_all);
    plot!(yticks=aux_ticks);

    l_x_pos = aux_lims[1] + R²_rel_pos_x*(aux_lims[2] - aux_lims[1]);
    l_y_pos = aux_lims[1] + R²_rel_pos_y*(aux_lims[2] - aux_lims[1]);
    R² = @sprintf "%.5lf" round(R²,digits=5);
    annotate!(l_x_pos, l_y_pos, text("R² = "*R², :right, 10));

    # Chemical Potentials
    y_label_all = L"$\tilde{\mu} \quad \mathrm{(This \ Work)}$";
    x_label_all = L"$\tilde{\mu} \quad (\mathrm{KS{-}DFT})$";

    # HC Plots
    pHC, R² = test_result_chemical_potential(1,6);
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
    R² = @sprintf "%.5lf" round(R²,digits=5);
    annotate!(l_x_pos, l_y_pos, text("R² = "*R², :right, 10));

    # HO Plots
    pHO, R² = test_result_chemical_potential(1,8);
    aux_lims = [-1.2,0.4];
    aux_ticks = -1.2:0.4:0.4;
    plot!(xlims=aux_lims);
    plot!(ylims=aux_lims);
    plot!(xticks=aux_ticks);
    plot!(yticks=aux_ticks);
    plot!(xlabel=x_label_all);

    l_x_pos = aux_lims[1] + R²_rel_pos_x*(aux_lims[2] - aux_lims[1]);
    l_y_pos = aux_lims[1] + R²_rel_pos_y*(aux_lims[2] - aux_lims[1]);
    R² = @sprintf "%.5lf" round(R²,digits=5);
    annotate!(l_x_pos, l_y_pos, text("R² = "*R², :right, 10));

    # CO Plots
    pCO, R² = test_result_chemical_potential(6,8);
    aux_lims = [-1.0,0.6];
    aux_ticks = -1.0:0.4:0.6;
    plot!(xlims=aux_lims);
    plot!(ylims=aux_lims);
    plot!(xticks=aux_ticks);
    plot!(yticks=aux_ticks);
    plot!(xlabel=x_label_all);

    l_x_pos = aux_lims[1] + R²_rel_pos_x*(aux_lims[2] - aux_lims[1]);
    l_y_pos = aux_lims[1] + R²_rel_pos_y*(aux_lims[2] - aux_lims[1]);
    R² = @sprintf "%.5lf" round(R²,digits=5);
    annotate!(l_x_pos, l_y_pos, text("R² = "*R², :right, 10));

    # NO Plots
    pNO, R² = test_result_chemical_potential(7,8);
    aux_lims = [-1.0,0.6];
    aux_ticks = -1.0:0.4:0.6;
    plot!(xlims=aux_lims);
    plot!(ylims=aux_lims);
    plot!(xticks=aux_ticks);
    plot!(yticks=aux_ticks);
    plot!(xlabel=x_label_all);

    l_x_pos = aux_lims[1] + R²_rel_pos_x*(aux_lims[2] - aux_lims[1]);
    l_y_pos = aux_lims[1] + R²_rel_pos_y*(aux_lims[2] - aux_lims[1]);
    R² = @sprintf "%.5lf" round(R²,digits=5);
    annotate!(l_x_pos, l_y_pos, text("R² = "*R², :right, 10));

    # Join all plots
    p = plot(eHC,eHO,eCO,eNO,pHC,pHO,pCO,pNO,layout=(2,4), size=(1100,470));
    savefig("Figures/HeteronuclearFitComps.pdf");

    return p;
end

function test_result_ΔE2(Z1::Int, Z2::Int)
    neutral_data, cation_data, anion_data = 
        read_all_sanitized_data(Z1,Z2);

    function get_data(data::Vector{ParsedFile})
        dft_r = zeros(Float64,length(data));
        dft_ΔE = zeros(Float64,length(data));

        bohr_to_angstrom = 0.529177;
        angstrom_to_bohr = 1.88973;

        for i in eachindex(data)
            dft_ΔE[i] = data[i].total_energy;
            dft_r[i] = bohr_to_angstrom*data[i].atomic_separation;
        end

        r0 = 0.0;
        r1 = angstrom_to_bohr*6.5;

        model_r = collect(r0:0.01:r1);
        model_ΔE = zeros(Float64,length(model_r));

        n_threads = Threads.nthreads();
        @threads for thread_id in 1:n_threads
            simulation = make_system_from_parsed_file(data[end]);
            r1 = zeros(Float64,3);
            r2 = zeros(Float64,3);
            set_atom_coordinates!(simulation.system.molecules[1],r1,1);
            set_atom_coordinates!(simulation.system.molecules[1],r2,2);

            for i in thread_id:n_threads:length(model_r)
                d = angstrom_to_bohr * model_r[i];
                r2[3] = d;

                set_atom_coordinates!(simulation.system.molecules[1],r1,1);
                set_atom_coordinates!(simulation.system.molecules[1],r2,2);
                polarize_molecules!(simulation);

                model_ΔE[i] = total_energy(simulation);
            end
        end

        hartree_to_ev = 27.2114;
        model_ΔE .*= hartree_to_ev;
        dft_ΔE .*= hartree_to_ev;

        return dft_ΔE, model_ΔE, dft_r, model_r;
    end

    neutral_dft_ΔE, neutral_model_ΔE, neutral_dft_r, neutral_model_r = 
        get_data(neutral_data);
    cation_dft_ΔE, cation_model_ΔE, cation_dft_r, cation_model_r = 
        get_data(cation_data);
    anion_dft_ΔE, anion_model_ΔE, anion_dft_r, anion_model_r = 
        get_data(anion_data);

    i_min_dft_e = argmin(neutral_dft_ΔE)
    min_dft_e = neutral_dft_ΔE[i_min_dft_e];
    min_dft_r = neutral_dft_r[i_min_dft_e];
    neutral_dft_ΔE .-= min_dft_e;
    cation_dft_ΔE .-= min_dft_e;
    anion_dft_ΔE .-= min_dft_e;

    hartree_to_ev = 27.2114;
    neutral_mol = make_system_from_parsed_file(neutral_data[i_min_dft_e]);
    min_model_e = total_energy(neutral_mol);
    min_model_e *= hartree_to_ev;
    
    neutral_model_ΔE .-= min_model_e;
    cation_model_ΔE .-= min_model_e;
    anion_model_ΔE .-= min_model_e;

    elem_symbol1 = get_element_symbol(Z1);
    elem_symbol2 = get_element_symbol(Z2);
    neutral_label = elem_symbol1*" + "*elem_symbol2;
    cation_label = "("*elem_symbol1*" + "*elem_symbol2*")⁺";
    anion_label = "("*elem_symbol1*" + "*elem_symbol2*")⁻";

    neutral_color = palette(:default)[1];
    cation_color = palette(:default)[2];
    anion_color = palette(:default)[3];

    p = plot(neutral_model_r, neutral_model_ΔE, color = neutral_color,
        label=neutral_label, linewidth = 2);
    scatter!(neutral_dft_r, neutral_dft_ΔE, color = neutral_color,
        label=false, markersize=3);

    plot!(cation_model_r, cation_model_ΔE, color = cation_color,
        label = cation_label, linewidth = 2);
    scatter!(cation_dft_r, cation_dft_ΔE, color = cation_color, 
        label=false, markersize=3);

    plot!(anion_model_r, anion_model_ΔE, color = anion_color,
        label = anion_label, linewidth = 2);
    scatter!(anion_dft_r, anion_dft_ΔE, color = anion_color, 
        label=false, markersize=3);

    plot!(framestyle = :box);
    return p;
end

function comp_heteronuclear_scan(pair1::Tuple{Int,Int}, pair2::Tuple{Int,Int})
    x_max = 6;

    p1 = test_result_ΔE2(pair1[1],pair1[2]);
    plot!(ylabel=L"$\Delta E \quad \mathrm{[eV]}$");
    plot!(ylims=[-5,55],yticks=-5:15:55);
    plot!(xlims=[0.0,x_max],xticks=(0:2:6,[]));
    plot!(legend = :topright)

    p2 = test_result_ΔE2(pair2[1],pair2[2]);
    plot!(xlabel=L"$d \quad \mathrm{[\AA ngstrom]}$");
    plot!(ylabel=L"$\Delta E \quad \mathrm{[eV]}$");
    plot!(ylims=[-5,55],yticks=-5:15:55);
    plot!(xlims=[0.0,x_max],xticks=0:2:6);
    plot!(legend = :topright)

    p = plot(p1,p2, layout=(2,1), size = (500,380));
    savefig("Figures/HeteronuclearScanComp.pdf");

    return p;
end

function plot_chemical_potential_scan(Z1::Int, Z2::Int, Q::Int)
    sim = make_diatomic_system(Z1,Z2,0.0,Q);

    r = collect(0.0:0.01:5);
    u = zeros(Float64,length(r));

    for i in eachindex(r)
        sim.system.molecules[1].atoms[1].coordinates[3] = r[i];
        aux_m, aux_y = polarization_matrix_problem(sim);
    
        u[i] = (aux_m \ aux_y)[end];
    end

    bohr_to_angstrom = 0.529177;
    r .*= bohr_to_angstrom;

    p = plot(r,u,label=false);
    plot!(xlabel="Atomic Separation [Å]");
    plot!(ylabel="Chemical Potential [Hartree]");
    plot!(framestyle = :box);
    plot!(xlims=[r[1],r[end]]);

    return p;
end

function plot_partial_charges(Z1::Int, Z2::Int, Q::Int)
    sim = make_diatomic_system(Z1,Z2,0.0,Q);

    ref_data = Vector{ParsedFile}();
    neutral_data, cation_data, anion_data = 
        read_all_sanitized_data(Z1,Z2);
    if Q == 0
        ref_data = neutral_data;
    elseif Q == 1
        ref_data = cation_data;
    else
        ref_data = anion_data;
    end

    r0 = ref_data[1].atomic_separation * 0.95;
    r1 = ref_data[end].atomic_separation * 1.05;

    r = collect(r0:0.01:r1);
    q1 = zeros(Float64,length(r));
    q2 = zeros(Float64,length(r));

    r_dft = zeros(Float64,length(ref_data));
    q1_dft = zeros(Float64,length(ref_data));
    q2_dft = zeros(Float64,length(ref_data));

    for i in eachindex(r)
        sim.system.molecules[1].atoms[1].coordinates[3] = r[i];
        aux_m, aux_y = polarization_matrix_problem(sim);
    
        aux_x = aux_m \ aux_y;

        q1[i] = aux_x[1];
        q2[i] = aux_x[2];
    end

    for i in eachindex(ref_data)
        r_dft[i] = ref_data[i].atomic_separation;
        q1_dft[i] = ref_data[i].partial_charge_1;
        q2_dft[i] = ref_data[i].partial_charge_2;
    end

    z1_eff = sim.system.molecules[1].atoms[1].valence_electrons;
    z2_eff = sim.system.molecules[1].atoms[2].valence_electrons;

    q1 = (1.0 .- q1) .* z1_eff;
    q2 = (1.0 .- q2) .* z2_eff;

    bohr_to_angstrom = 0.529177;
    r .*= bohr_to_angstrom;
    r_dft .*= bohr_to_angstrom;

    p = plot(r,q1,label=get_element_symbol(Z1),
        color=palette(:auto)[1]);
    plot!(r,q2,label=get_element_symbol(Z2),
        color=palette(:auto)[2]);
    plot!(xlabel="Atomic Separation [Å]");
    plot!(ylabel="Partial Charge");
    plot!(framestyle = :box);
    plot!(xlims=[r[1],r[end]]);

    scatter!(r_dft,q1_dft,label=false,color=palette(:auto)[1]);
    scatter!(r_dft,q2_dft,label=false,color=palette(:auto)[2]);

    return p;
end

function plot_compare_partial_charges(Z1::Int, Z2::Int, Q::Int)
    sim = make_diatomic_system(Z1,Z2,0.0,Q);

    ref_data = Vector{ParsedFile}();
    neutral_data, cation_data, anion_data = 
        read_all_sanitized_data(Z1,Z2);
    if Q == 0
        ref_data = neutral_data;
    elseif Q == 1
        ref_data = cation_data;
    else
        ref_data = anion_data;
    end

    q1 = zeros(Float64,length(ref_data));
    q2 = zeros(Float64,length(ref_data));

    q1_dft = zeros(Float64,length(ref_data));
    q2_dft = zeros(Float64,length(ref_data));

    for i in eachindex(ref_data)
        set_diatomic_system_to_parsed_file!(sim,ref_data[i]);
        polarize_molecules!(sim);


        q1_dft[i] = ref_data[i].partial_charge_1;
        q2_dft[i] = ref_data[i].partial_charge_2;

        q1[i] = sim.system.molecules[1].atoms[1].polarization_coefficient;
        q2[i] = sim.system.molecules[1].atoms[2].polarization_coefficient;
    end

    z1_eff = sim.system.molecules[1].atoms[1].valence_electrons;
    z2_eff = sim.system.molecules[1].atoms[2].valence_electrons;

    # q1 = (1.0 .- q1) .* z1_eff;
    # q2 = (1.0 .- q2) .* z2_eff;

    q1 *= z1_eff;
    q2 *= z2_eff;

    q1_dft = z1_eff .- q1_dft;
    q2_dft = z2_eff .- q2_dft;

    q12_dft = vcat(q1_dft,q2_dft);
    q12_model = vcat(q1,q2);

    min_x = minimum(vcat(q12_dft,q12_model));
    max_x = maximum(vcat(q12_dft,q12_model));

    p = plot([min_x,max_x],[min_x,max_x],label=false)
    scatter!(q12_dft,q12_model);

    return p;
end

comp_heteronuclear_scan((7,8),(6,8));
compare_data()
