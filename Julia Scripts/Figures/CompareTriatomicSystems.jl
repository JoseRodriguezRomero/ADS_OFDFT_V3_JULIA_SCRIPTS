using Measures, Plots;

include("../FitCoeffs_General.jl");
pyplot();

function plot_data_triatomic_system(Z1::Int, Z2::Int, 
    θ1::Real, θ2::Real, r1::Real, r2::Real)
    simulation = initialize_simulation_environment();
    reference_atoms = simulation.basis_set_settings.reference_atoms;

    hartree_to_ev = 27.2114;
    angstrom_to_bohr = 1.88973;

    all_angles = collect(θ1:((θ2-θ1)/200.0):θ2);
    all_separations = collect(r1:((r2-r1)/200.0):r2);
    
    energy = zeros(Float64,length(all_separations),length(all_angles));
    
    min_e = 1.0E8;
    curr_min = (0,0);

    for i in eachindex(all_separations)
        for j in eachindex(all_angles)
            separation = all_separations[i] * angstrom_to_bohr;
            angle = all_angles[j];

            molecule = make_triatomic_molecule(
                reference_atoms,Z2,Z1,Z1,separation,separation,angle);
            simulation.system.molecules = [molecule];

            polarize_molecules!(simulation);
            energy[i,j] = total_energy(simulation);

            if energy[i,j] < min_e
                min_e = energy[i,j];
                curr_min = (i,j);
            end
        end
    end

    println(@sprintf "Bond length: %12.3lf" all_separations[curr_min[1]]);
    println(@sprintf "Bond angle : %12.3lf" all_angles[curr_min[2]]);

    energy .-= minimum(energy[:]);
    energy .*= hartree_to_ev;

    return all_angles, all_separations, energy;
end

function plot_h2o_comparison()
    H2O_angle = 104.5;
    H2O_bond_lenght = 0.9854;

    r1 = 0.75;
    r2 = 1.75;

    θ1 = 25.0;
    θ2 = 180.0;

    Z1 = 1;
    Z2 = 8;

    all_angles, all_separations, energy = 
        plot_data_triatomic_system(Z1,Z2,θ1,θ2,r1,r2);

    x_label = get_element_symbol(Z1) * "-" * get_element_symbol(Z2) *
        "-" * get_element_symbol(Z1) * " angle [degrees]";
    y_label = get_element_symbol(Z2) * "-" * get_element_symbol(Z1) *
        " length\n[Ångstrom]"
    c_label = "ΔE [eV]";

    energy[energy .> 12] .= 12;

    p = contourf(all_angles, all_separations, energy,
        xlabel=x_label, ylabel=y_label,levels=14,
        clims=(0,12),colorbar_title=c_label, color=:turbo);
    
    contour!(all_angles, all_separations, energy,
        xlabel=x_label,levels=14,clims=(0,12),
        linewidth = 1.0,linestyle=:solid,linecolor=[:black]);

    scatter!([H2O_angle],[H2O_bond_lenght],label=false,
        color=:white,markersize=6);
    plot!(framestyle = :box);
    plot!(xticks=45:45:180);
    plot!(size = (500,180));
    plot!(xlims=[θ1,θ2],ylims=[r1,r2]);
    plot!(left_margin=2mm, bottom_margin=1mm, right_margin=-6mm);
    # savefig("Figures/H2O_energy_scan.pdf");

    contourf!(colorbar_ticks=(0:4:12,string.(0:4:12)));

    return p;
end

function plot_co2_comparison()
    CO2_angle = 180.0;
    CO2_bond_lenght = 1.19;

    r1 = 1.0;
    r2 = 1.8;

    θ1 = 45.0;
    θ2 = 225.0;

    Z1 = 8;
    Z2 = 6;

    all_angles, all_separations, energy = 
        plot_data_triatomic_system(Z1,Z2,θ1,θ2,r1,r2);

    x_label = get_element_symbol(Z1) * "-" * get_element_symbol(Z2) *
        "-" * get_element_symbol(Z1) * " angle [degrees]";
    y_label = get_element_symbol(Z2) * "-" * get_element_symbol(Z1) *
        " length\n[Ångstrom]"
    c_label = "ΔE [eV]";

    energy[energy .> 24] .= 24;

    p = contourf(all_angles, all_separations, energy,
        xlabel=x_label, ylabel=y_label,levels=14,
        clims=(0,24),colorbar_title=c_label,color=:turbo);

    contour!(all_angles, all_separations, energy,
        xlabel=x_label, ylabel=y_label,levels=14,
        clims=(0,24),linewidth = 1.0,linestyle=:solid,
        linecolor=[:black]);

    scatter!([CO2_angle],[CO2_bond_lenght],
        label=false,color=:white,markersize=6);
    plot!(framestyle = :box);
    plot!(xticks=45:45:225);
    plot!(size = (500,180));
    plot!(xlims=[θ1,θ2],ylims=[r1,r2]);
    plot!(left_margin=2mm, bottom_margin=1mm, right_margin=-6mm);
    contourf!(clims=(0,20));
    plot!(colorbar_ticks = 0:6:24);
    # savefig("Figures/CO2_energy_scan.pdf");

    contourf!(colorbar_ticks=(0:6:24,string.(0:6:24)));

    return p;
end

function plot_triatomic_comparisons()
    p1 = plot_h2o_comparison();
    p2 = plot_co2_comparison();

    p = plot(p1,p2,layout=(2,1));
    plot!(size = (500,310));

    savefig("Figures/TriatomicMoleculesScans.pdf");

    return p;
end
