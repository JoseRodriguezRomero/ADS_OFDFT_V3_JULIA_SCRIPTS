using Glob;

include("ForceFieldBase.jl");

function thomas_fermi_ke(electron_cloud::ElectronCloud)
    # Calculates the Thomas-Fermi kinetic energy of the atom whose 
    # atomic number is specified in the argument of this function.
    tf_coeff = (3.0/10.0) * ((3*(π^2))^(2.0/3.0));

    num_clouds = length(electron_cloud.basis_function_amplitude);
    gauss_c = electron_cloud.basis_function_amplitude;
    gauss_λ = electron_cloud.basis_function_decay;

    function thomas_fermi_ke_integrand(r,p)
        # Auxiliary function to integrate the density in spherical coordinates.
        ret_val = 0.0;

        for i in 1:num_clouds
            aux_val = (gauss_λ[i]/π)^(3.0/2.0);
            aux_val *= gauss_c[i]*exp(-gauss_λ[i]*(r^2));
            ret_val += aux_val;
        end

        return 4*π*(r^2)*(abs(ret_val)^(5.0/3.0));
    end

    prob = IntegralProblem(thomas_fermi_ke_integrand,(0.0,300.0));
    sol = solve(prob, HCubatureJL(), reltol = 1.0E-8, abstol = 1.0E-8);
    return tf_coeff*(sol[1]);
end

function von_weizsacker_ke(electron_cloud::ElectronCloud)
    # Calculates the von Weizsäcker kinetic energy of the atom whose atomic 
    # number is specified in the argument of this function.
    vw_coeff = (1.0/8.0);

    num_clouds = length(electron_cloud.basis_function_amplitude);
    gauss_c = electron_cloud.basis_function_amplitude;
    gauss_λ = electron_cloud.basis_function_decay;

    function von_weizsacker_ke_integrand(r,p)
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

    prob = IntegralProblem(von_weizsacker_ke_integrand,(0.0,300.0));
    sol = solve(prob, HCubatureJL(), reltol = 1.0E-8, abstol = 1.0E-8);
    return vw_coeff*(sol[1]);
end

function load_isolated_atoms_chemical_potential()
    max_elements = 18;

    neutral_energies = zeros(Float64,max_elements);
    cation_energies = zeros(Float64,max_elements);
    anion_energies = zeros(Float64,max_elements);

    function read_energies(file_name::String)
        lines = readlines(file_name);
        energies = zeros(Float64,length(lines));

        for i in eachindex(lines)
            energies[i] = parse(Float64,lines[i]);
        end

        return energies;
    end

    neutral_energies = read_energies("Atom Energies/Neutral/elem_total_energy.txt");
    cation_energies = read_energies("Atom Energies/Cation/elem_total_energy.txt");
    anion_energies = read_energies("Atom Energies/Anion/elem_total_energy.txt");

    atoms_μ0 = Dict{Int,Float64}();
    for i in 1:max_elements
        e_neutral = neutral_energies[i];
        e_cation = cation_energies[i];
        e_anion = anion_energies[i];

        M = zeros(Float64,3,3);
        Y = zeros(Float64,3);

        M[1,:] = [ (0.0)^2,  (0.0)^1, 1.0];
        M[2,:] = [ (1.0)^2,  (1.0)^1, 1.0];
        M[3,:] = [(-1.0)^2, (-1.0)^1, 1.0];

        Y[:] = [e_neutral, e_cation, e_anion];
        X = - (M \ Y)[:];

        q_neutral = 0.0;
        atom_μ0 = 2 * X[1] * q_neutral + X[2];
        atoms_μ0[i] = atom_μ0;
    end

    return atoms_μ0;
end

function load_basis_set(compute_ke::Bool = false)
    max_atomic_number = 18;

    function make_reference_atom()
        atomic_number = 1;
        coordinates = zeros(Float64,3);
        core_electrons = 0;
        valence_electrons = 0;
        polarization_coefficient = 1.0;
        electron_cloud_shells = Vector{ElectronCloud}();

        return Atom(atomic_number,coordinates,core_electrons,valence_electrons,
            polarization_coefficient,electron_cloud_shells);
    end

    function make_reference_electron_cloud()
        basis_function_amplitude = zeros(Float64,0);
        basis_function_decay = zeros(Float64,0);
        thomas_fermi_ke = 0.0;
        von_weizsacker_ke = 0.0;

        return ElectronCloud(basis_function_amplitude, basis_function_decay,
            thomas_fermi_ke, von_weizsacker_ke);
    end

    function electrons_in_cloud(electron_cloud::ElectronCloud)
        num_electrons = sum(electron_cloud.basis_function_amplitude);
        return round(Int,num_electrons);
    end

    function electrons_in_atom(atom::Atom)
        num_electrons = 0;
        for electron_cloud in atom.electron_cloud_shells
            num_electrons += electrons_in_cloud(electron_cloud);
        end

        return num_electrons;
    end

    function valence_electrons_in_atom(atom::Atom)
        return electrons_in_cloud(atom.electron_cloud_shells[end]);
    end

    function core_electrons_in_atom(atom::Atom)
        return electrons_in_atom(atom) - valence_electrons_in_atom(atom);
    end

    reference_atoms = Dict{Int,Atom}();
    for i in 1:max_atomic_number
        elem_symbol = get_element_symbol(i);
        reference_atom = make_reference_atom();

        base_dir = "Basis Set/"*elem_symbol*"/";

        max_shells = 3;
        for num_shell in 1:max_shells
            file_name = base_dir * "shell_" * string(num_shell) * ".txt";
            
            if !isfile(file_name)
                break;
            end

            lines = readlines(file_name);
            refere_electron_cloud = make_reference_electron_cloud();

            for line in lines
                line_splitted = split(line);
                guassian_λ = parse(Float64,line_splitted[1]);
                guassian_c = parse(Float64,line_splitted[2]);

                append!(refere_electron_cloud.basis_function_decay,
                    guassian_λ);
                append!(refere_electron_cloud.basis_function_amplitude, 
                    guassian_c);
            end

            tf_file_name = 
                base_dir * "shell_" * string(num_shell) * "_tf_ke.txt";
            vw_file_name = 
                base_dir * "shell_" * string(num_shell) * "_vw_ke.txt";

            if compute_ke == true
                refere_electron_cloud.thomas_fermi_ke = 
                    thomas_fermi_ke(refere_electron_cloud);
                refere_electron_cloud.von_weizsacker_ke = 
                    von_weizsacker_ke(refere_electron_cloud);

                tf_file_id = open(tf_file_name,"w");
                vw_file_id = open(vw_file_name,"w");

                tf_ke = refere_electron_cloud.thomas_fermi_ke;
                vw_ke = refere_electron_cloud.von_weizsacker_ke;

                write(tf_file_id, @sprintf "%.12E" tf_ke);
                write(vw_file_id, @sprintf "%.12E" vw_ke);

                close(tf_file_id);
                close(vw_file_id);
            else
                tf_lines = readlines(tf_file_name);
                vw_lines = readlines(vw_file_name);

                tf_ke = parse(Float64,tf_lines[1]);
                vw_ke = parse(Float64,vw_lines[1]);

                refere_electron_cloud.thomas_fermi_ke = tf_ke;
                refere_electron_cloud.von_weizsacker_ke = vw_ke;
            end

            append!(reference_atom.electron_cloud_shells,
                [refere_electron_cloud]);
        end

        reference_atom.atomic_number = i;
        reference_atom.core_electrons = core_electrons_in_atom(reference_atom);
        reference_atom.valence_electrons = valence_electrons_in_atom(reference_atom);

        reference_atoms[i] = reference_atom;
    end

    atoms_μ0 = load_isolated_atoms_chemical_potential();

    return BasisSetSettings(max_atomic_number,reference_atoms,atoms_μ0);
end
