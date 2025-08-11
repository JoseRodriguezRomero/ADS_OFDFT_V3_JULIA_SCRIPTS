using Optim;
using Integrals;
using SpecialFunctions;
using LinearAlgebra, GenericLinearAlgebra;

mutable struct Molecule
    charge::Number;
    energy::Number;
    chem_μ::Number;
    name::String;
    atoms_data::Matrix;
    cloud_data::Matrix;
end

function Base.copy(mol::Molecule)
    return Molecule(
        copy(mol.charge),
        copy(mol.energy),
        copy(mol.chem_μ),
        mol.name,
        copy(mol.atoms_data),
        copy(mol.cloud_data)
    )
end

function Molecule(type = Float64)
    charge = type(0.0);
    energy = type(0.0);
    chem_μ = type(0.0);
    name = "Kryptonite";
    atoms_data = zeros(type,0,6);
    cloud_data = zeros(type,0,6);
    return Molecule(charge,energy,chem_μ,name,atoms_data,cloud_data);
end

function ReadMolecule(molec_data::String, type = Float64)
    fileID = open(molec_data,"r");
    molecule = Molecule(type);

    readline(fileID);

    while true
        line = readline(fileID);
        line_splitted = split(line);

        if length(line_splitted) < 4
            break;
        end

        x = type(parse(Float64,line_splitted[1]));
        y = type(parse(Float64,line_splitted[2]));
        z = type(parse(Float64,line_splitted[3]));
        q = type(parse(Float64,line_splitted[4]));
        Z = type(parse(Float64,line_splitted[5]));
        
        new_atom = zeros(type,1,5);
        new_atom[1] = x;
        new_atom[2] = y;
        new_atom[3] = z;
        new_atom[4] = q;    # Effective Atomic Number 
        new_atom[5] = Z;    # Atomic Number

        molecule.atoms_data = [molecule.atoms_data; new_atom];
    end

    readline(fileID);
    while true
        line = readline(fileID);
        line_splitted = split(line);

        if length(line_splitted) < 5
            break;
        end

        x = type(parse(Float64,line_splitted[1]));
        y = type(parse(Float64,line_splitted[2]));
        z = type(parse(Float64,line_splitted[3]));
        c = type(parse(Float64,line_splitted[4]));
        λ = type(parse(Float64,line_splitted[5]));

        new_cloud = zeros(type,1,6)
        new_cloud[1] = x;
        new_cloud[2] = y;
        new_cloud[3] = z;
        new_cloud[4] = c;         # Magnitude
        new_cloud[5] = λ;         # Decay coefficient
        new_cloud[6] = type(1.0); # Polarization coefficient

        molecule.cloud_data = [molecule.cloud_data; new_cloud];
    end

    close(fileID);
    return molecule;
end

function CenterAtAtomIndex!(molecule::Molecule,index::Int)
    r0 = molecule.atoms_data[index,1:3];
    molecule.atoms_data[:,1:3] .-= r0';
    molecule.cloud_data[:,1:3] .-= r0';
end

function MoveAndRotateMolec!(molecule::Molecule,angs_and_disp::AbstractVector)
    θx = angs_and_disp[1];
    θy = angs_and_disp[2];
    θz = angs_and_disp[3];
    dr = angs_and_disp[4:6];

    aux_type = typeof(angs_and_disp[1,1]);

    rot_x = zeros(aux_type,3,3);
    rot_x[1,1] = 1.0;
    rot_x[2,2] = cos(θx);
    rot_x[2,3] = -sin(θx);
    rot_x[3,2] = sin(θx);
    rot_x[3,3] = cos(θx);

    rot_y = zeros(aux_type,3,3);
    rot_y[1,1] = cos(θy);
    rot_y[1,3] = sin(θy);
    rot_y[2,2] = 1.0;
    rot_y[3,1] = -sin(θy);
    rot_y[3,3] = cos(θy);

    rot_z = zeros(aux_type,3,3);
    rot_z[1,1] = cos(θz);
    rot_z[1,2] = -sin(θz);
    rot_z[2,1] = sin(θz);
    rot_z[2,2] = cos(θz);
    rot_z[3,3] = 1.0;

    rot = (rot_x*rot_y*rot_z)';

    molecule.atoms_data[:,1:3] = molecule.atoms_data[:,1:3]*rot;
    molecule.atoms_data[:,1:3] .+= dr';

    molecule.cloud_data[:,1:3] = molecule.cloud_data[:,1:3]*rot;
    molecule.cloud_data[:,1:3] .+= dr';
end

num_elems = 18;

num_1b_combs = num_elems;
num_2b_combs = round(Int,num_elems*(num_elems+1)/2.0);

num_1b_coeffs = 4;
num_2b_coeffs = 4;

coeff_map_2b = Dict{Tuple{Int,Int},Int}();
inv_coeff_map_2b = Dict{Int,Tuple{Int,Int}}();

num_xc_coeffs = 0;
for i in 1:num_elems
    global num_xc_coeffs;
    for j in i:num_elems
        num_xc_coeffs += 1;
        coeff_map_2b[(i,j)] = num_xc_coeffs;
        coeff_map_2b[(j,i)] = num_xc_coeffs;
        inv_coeff_map_2b[num_xc_coeffs] = (i,j);
    end
end

num_coeffs = num_elems*num_1b_coeffs + num_2b_combs*num_2b_coeffs;

factorials = gamma.(BigFloat.(collect(1:25) .+ 1, precision = 256));
factorials = Float64.(1.0 ./ factorials);

function ReadCloudCoeffs()
    # Reads the coefficients used to fit the electron density of all elements 
    # in the first three rows of the periodic table using a linear combination 
    # of s-type Gaussian basis-functions.
    file_name = "../Reference Data/atom_density/ecp_coeffs.txt";
    fileID = open(file_name, "r");
    lines = readlines(fileID);
    cloud_data = zeros(Float64, 18, 6);

    for i in 1:18
        line_splitted = split(lines[i+1]);
        for j in 1:6
            cloud_data[i, j] = parse(Float64, line_splitted[j+1]);
        end
    end

    close(fileID);
    return cloud_data;
end

function MakeAtom(Z::Integer,charge::Integer)
    # Makes a molecule object with a single atom whose atomic number is Z.
    all_cloud_coeffs = ReadCloudCoeffs();
    atom = Molecule();
    atom.charge = charge;

    atom.atoms_data = zeros(Float64,1,5);
    atom.cloud_data = zeros(Float64,3,6);

    Zeff = round(sum(all_cloud_coeffs[Z,1:2:end]));
    atom.atoms_data[1,5] = Z;
    atom.atoms_data[1,4] = Zeff;

    atom.cloud_data[1:3,4] .= all_cloud_coeffs[Z,1:2:end];
    atom.cloud_data[1:3,5] .= all_cloud_coeffs[Z, 2:2:end];
    atom.cloud_data[1:3,6] .= (Zeff - charge) / Zeff;

    return atom;
end

function ThomasFermiUnpolIntegral(Z::Integer)
    # Calculates the Thomas-Fermi kinetic energy of the atom whose atomic 
    # number is specified in the argument of this function.
    #
    # NOTICE: This function is slow and it is meant to be cached in a table 
    # for all elements in the periodic table to avoid bottlenecks.
    cloud_coeffs = ReadCloudCoeffs();
    TF_coeff = (3.0/10.0) * ((3*(π^2))^(2.0/3.0));

    num_clouds = Int((size(cloud_coeffs)[2])/2);
    gauss_c = cloud_coeffs[Z,1:2:end];
    gauss_λ = cloud_coeffs[Z,2:2:end];

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
    return TF_coeff*(sol[1]);
end

function vonWeizsackerUnpolIntegral(Z::Integer)
    # Calculates the von Weizsäcker kinetic energy of the atom whose atomic 
    # number is specified in the argument of this function.
    #
    # NOTICE: This function is slow and it is meant to be cached in a table 
    # for all elements in the periodic table to avoid bottlenecks.
    cloud_coeffs = ReadCloudCoeffs();
    vW_coeff = (1.0/8.0);

    num_clouds = Int((size(cloud_coeffs)[2])/2);
    gauss_c = cloud_coeffs[Z,1:2:end];
    gauss_λ = cloud_coeffs[Z,2:2:end];

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
    return vW_coeff*(sol[1]);
end

unpol_TF_KE = ThomasFermiUnpolIntegral.(1:num_elems);
unpol_vW_KE = vonWeizsackerUnpolIntegral.(1:num_elems);

function NaiveCoulombInt(λ::Real, d::Real)
    if (abs(d) < 1.0E-6)
        return 2*sqrt(λ/π);
    end

    return erf(d*sqrt(λ))/d;
end 

function XC_Sph(λ::Real, d::Real)
    if (abs(d) < 1.0E-6)
        return -((4*λ^(3/2)*exp(-λ))/sqrt(π));
    end

    return (1/d)*sqrt(λ/π)*((1-exp(4*λ*d))/(exp(λ*(d+1)^2)));
end

function XC_Cyl(λ::Real, d::Real)
    return (1/d)*sqrt(λ/π)*((((exp(2*λ*d)-1)^2)/(exp(λ*(d+1)^2))) + 
        ((2*(exp(-λ)-1))/(exp(λ*d^2))));
end

function OverlapFuncInt(λ::Real, d::Real)
    return exp(-λ*(d^2))*((λ/π)^(3/2));
end

function GetOneBodyCoeffs(one_body_coeffs::AbstractMatrix)
    # Splits the matrix of coefficients one_body_coeffs (with dimensions N×2) 
    # to smaller chuncks containing all the same kind of coefficient 
    # contributions to the total and polarization energies.
    global num_elems, num_2b_combs, num_1b_coeffs, num_2b_coeffs;

    tot_num_1b_coeffs = num_1b_coeffs*num_elems;

    XCa_1b = one_body_coeffs[1:num_1b_coeffs:tot_num_1b_coeffs,:];  # EE sph
    XCc_1b = one_body_coeffs[2:num_1b_coeffs:tot_num_1b_coeffs,:];  # EN sph
    XCe_1b = one_body_coeffs[3:num_1b_coeffs:tot_num_1b_coeffs,:];  # TF KE
    XCf_1b = one_body_coeffs[4:num_1b_coeffs:tot_num_1b_coeffs,:];  # vW KE

    return XCa_1b, XCc_1b, XCe_1b, XCf_1b;
end

function GetTwoBodyCoeffs(two_body_coeffs::AbstractMatrix)
    # Splits the matrix of coefficients two_body_coeffs (with dimensions N×3) 
    # to smaller chuncks containing all the same kind of coefficient 
    # contributions to the total and polarization energies.
    global num_elems, num_1b_coeffs, num_2b_coeffs;

    XCa_2b = two_body_coeffs[1:num_2b_coeffs:end,:];    # EE sph
    XCb_2b = two_body_coeffs[2:num_2b_coeffs:end,:];    # EE cyl
    XCc_2b = two_body_coeffs[3:num_2b_coeffs:end,:];    # EN sph
    XCd_2b = two_body_coeffs[4:num_2b_coeffs:end,:];    # EN cyl

    return XCa_2b, XCb_2b, XCc_2b, XCd_2b;
end

function CombXCFunc1B(coeffs::AbstractMatrix, Z::Integer, ζ::Number)
    # Calculates the XC coefficient based on the polarization coefficients.
    # (One body cases)
    b = coeffs[Z,1];
    m = coeffs[Z,2];

    return m*ζ + b;
end

function CombXCFuncPol1B(coeffs::AbstractMatrix, Z::Integer, ζ::Number)
    # Calculates the XC coefficient based on the polarization coefficients.
    # (One body cases)
    #
    # NOTICE: Use this when taking the functional derivative of the atomic 
    # contribution of the atom with atomic number Z and polarization 
    # coefficient ζ.
    b = coeffs[Z,1];
    m = coeffs[Z,2];

    return 2*m*ζ + b;
end

function CombXCFunc2B(coeffs::AbstractMatrix, Z1::Integer, Z2::Integer,
    ζ1::Number, ζ2::Number)
    # Calculates the XC coefficient based on the polarization coefficients.
    # (One body cases)
    b = coeffs[coeff_map_2b[(Z1,Z2)],1];
    m1 = coeffs[coeff_map_2b[(Z1,Z2)],2];
    m2 = coeffs[coeff_map_2b[(Z1,Z2)],3];

    if (Z1 < Z2)
        return m1*ζ1 + m2*ζ2 + b;
    else
        return m1*ζ2 + m2*ζ1 + b;
    end
end

function CombXCFuncPol2B(coeffs::AbstractMatrix, Z1::Integer, Z2::Integer,
    ζ1::Number, ζ2::Number)
    # Calculates the XC coefficient based on the polarization coefficients.
    # (One body cases)
    b = coeffs[coeff_map_2b[(Z1,Z2)],1];
    m1 = coeffs[coeff_map_2b[(Z1,Z2)],2];
    m2 = coeffs[coeff_map_2b[(Z1,Z2)],3];

    if (Z1 < Z2)
        return 2*m1*ζ1 + m2*ζ2 + b;
    else
        return m1*ζ2 + 2*m2*ζ1 + b;
    end
end

function MatPolarizeMolecules(molecule::Molecule, 
    all_coeffs::AbstractVector, type_cast::Int = 0)
    # Calculates the polarization coefficients of all atoms in the system.
    return MatPolarizeMolecules([molecule],all_coeffs,type_cast);
end

function MatPolarizeMolecules(molecules::AbstractVector{Molecule}, 
    all_coeffs::AbstractVector, type_cast::Int = 0)
    # Calculates the polarization coefficients of all atoms in the system.
    global num_elems, factorials;
    global unpol_TF_KE, unpol_vW_KE;

    if (type_cast == 0)
        # Cast to the type of the XC coefficients.
        aux_type = typeof(all_coeffs[1][1]);
    else
        # Cast to the type in the data matrices of the atoms and electron 
        # clouds of all the molecules.
        aux_type = typeof(molecules[1].cloud_data[1,1]);
    end

    XCa_1b, XCc_1b, XCe_1b, XCf_1b = GetOneBodyCoeffs(all_coeffs[1]);
    XCa_2b, XCb_2b, XCc_2b, XCd_2b = GetTwoBodyCoeffs(all_coeffs[2]);

    clouds_per_atom = 0;

    tot_num_atoms = 0;
    tot_num_clouds = 0;
    atom_ind_base = [];
    for molecule in molecules
        if isempty(atom_ind_base)
            push!(atom_ind_base,1);
            tot_num_atoms += length(molecule.atoms_data[:,1]);
            tot_num_clouds += length(molecule.cloud_data[:,1]);
        else
            aux_num_atoms = size(molecule.atoms_data)[1];
            aux_num_clouds = size(molecule.cloud_data)[1];
            
            push!(atom_ind_base, tot_num_atoms + 1);
            tot_num_atoms += aux_num_atoms;
            tot_num_clouds += aux_num_clouds;
        end
    end

    num_vars = tot_num_atoms + 1;
    aux_M = zeros(aux_type, num_vars, num_vars);
    aux_Y = zeros(aux_type, num_vars, 1);

    for ii in eachindex(molecules)
        molecule1 = molecules[ii];
        num_atoms1 = size(molecule1.atoms_data)[1];
        num_clouds1 = size(molecule1.cloud_data)[1];
        clouds_per_atom = ceil(Int,num_clouds1 / num_atoms1);

        # CHEMICAL POTENTIAL
        # Lagrange multiplier to keep the number of electrons constant in the 
        # new fictitious electron density.
        ind1 = atom_ind_base[ii];
        ind2 = ind1 + num_atoms1 - 1;

        aux_M[ind1:ind2,end] = -molecule1.atoms_data[:,4];
        aux_M[end,ind1:ind2] = molecule1.atoms_data[:,4];
        aux_Y[end] += round(Int, sum(molecule1.atoms_data[:,4]));
        aux_Y[end] -= molecule1.charge;

        # Kinetic Energy Contributions
        for i in 1:num_atoms1
            Z1 = round(Int, molecule1.atoms_data[i,5]);
            ζ1 = molecule1.cloud_data[clouds_per_atom*(i-1)+1,6];

            ii0 = atom_ind_base[ii] + i - 1;

            TF_fit_coeff = 0.0;
            TF_fit_coeff += (5.0/3.0)*XCe_1b[Z1,1]*(ζ1^(2.0/3.0));
            TF_fit_coeff += (8.0/3.0)*XCe_1b[Z1,2]*(ζ1^(5.0/3.0));

            vW_fit_coeff = 0.0;
            vW_fit_coeff += 1.0*XCf_1b[Z1,1]*(ζ1^(0.0));
            vW_fit_coeff += 2.0*XCf_1b[Z1,2]*(ζ1^(1.0));

            aux_M[ii0,ii0] += TF_fit_coeff*unpol_TF_KE[Z1];
            aux_M[ii0,ii0] += vW_fit_coeff*unpol_vW_KE[Z1];
        end

        # Fill the rest of the matrix
        for jj in eachindex(molecules)
            molecule2 = molecules[jj];
            num_atoms2 = size(molecule2.atoms_data)[1];
            num_clouds2 = size(molecule2.cloud_data)[1];

            # cloud-nuclei
            for i in 1:num_clouds1
                for j in 1:num_atoms2
                    Z1 = round(Int, 
                        molecule1.atoms_data[ceil(Int,i/clouds_per_atom),5]);
                    Z2 = round(Int, molecule2.atoms_data[j,5]);

                    c1 = molecule1.cloud_data[i,4];
                    λ1 = molecule1.cloud_data[i,5];
                    q2 = molecule2.atoms_data[j,4];

                    ζ1 = molecule1.cloud_data[i,6];
                    ζ2 = molecule2.cloud_data[j*clouds_per_atom,6];
                    
                    ii0 = atom_ind_base[ii] + ceil(Int,i/clouds_per_atom) - 1;

                    aux_dist = molecule1.cloud_data[i,1:3];
                    aux_dist -= molecule2.atoms_data[j,1:3];
                    aux_dist = norm(aux_dist);

                    λ = aux_type(λ1);
                    aux_dist = aux_type(aux_dist);

                    # Naive contribution
                    aux_Y[ii0] += c1*q2*NaiveCoulombInt(λ,aux_dist);
                    
                    # Get the XC coefficients
                    xc_coeff_1 = aux_type(0.0);
                    xc_coeff_2 = aux_type(0.0);

                    if (aux_dist > 1.0E-6)
                        xc_coeff_1 = CombXCFuncPol2B(XCc_2b,Z1,Z2,ζ1,ζ2);
                        xc_coeff_2 = CombXCFuncPol2B(XCd_2b,Z1,Z2,ζ1,ζ2);
                        
                        # XC contributions
                        aux_Y[ii0] -= xc_coeff_1*c1*q2*XC_Sph(λ,aux_dist);
                        aux_Y[ii0] -= xc_coeff_2*c1*q2*XC_Cyl(λ,aux_dist);
                    else
                        xc_coeff_1 = CombXCFuncPol1B(XCc_1b,Z1,ζ1);

                        # XC contributions
                        aux_Y[ii0] -= xc_coeff_1*c1*q2*XC_Sph(λ,aux_dist);
                    end
                end
            end

            # cloud-cloud
            for i in 1:num_clouds1
                for j in 1:num_clouds2
                    Z1 = round(Int, 
                        molecule1.atoms_data[ceil(Int,i/clouds_per_atom),5]);
                    Z2 = round(Int, 
                        molecule2.atoms_data[ceil(Int,j/clouds_per_atom),5]);

                    ζ1 = molecule1.cloud_data[i,6];
                    ζ2 = molecule2.cloud_data[j,6];

                    c1 = molecule1.cloud_data[i,4];
                    λ1 = molecule1.cloud_data[i,5];
                    c2 = molecule2.cloud_data[j,4];
                    λ2 = molecule2.cloud_data[j,5];

                    ii0 = atom_ind_base[ii] + ceil(Int,i/clouds_per_atom) - 1;
                    jj0 = atom_ind_base[jj] + ceil(Int,j/clouds_per_atom) - 1;

                    aux_dist = molecule1.cloud_data[i,1:3];
                    aux_dist -= molecule2.cloud_data[j,1:3];
                    aux_dist = norm(aux_dist);

                    λ = aux_type((λ1 * λ2) / (λ1 + λ2));
                    aux_dist = aux_type(aux_dist);

                    # Naive contribution
                    aux_M[ii0,jj0] += c1*c2*NaiveCoulombInt(λ,aux_dist);

                    # Get the XC coefficients
                    xc_coeff_1 = aux_type(0.0);
                    xc_coeff_2 = aux_type(0.0);

                    if (aux_dist > 1.0E-6)
                        xc_coeff_1 = CombXCFuncPol2B(XCa_2b,Z1,Z2,ζ1,ζ2);
                        xc_coeff_2 = CombXCFuncPol2B(XCb_2b,Z1,Z2,ζ1,ζ2);

                        # XC contributions
                        aux_M[ii0,jj0] += xc_coeff_1*c1*c2*XC_Sph(λ,aux_dist);
                        aux_M[ii0,jj0] += xc_coeff_2*c1*c2*XC_Cyl(λ,aux_dist);
                    else
                        xc_coeff_1 = CombXCFuncPol1B(XCa_1b,Z1,ζ1);
                        
                        # XC contributions
                        aux_M[ii0,jj0] += xc_coeff_1*c1*c2*XC_Sph(λ,aux_dist);
                    end
                end
            end
        end
    end

    return atom_ind_base, clouds_per_atom, aux_M, aux_Y;
end

function SCFMinFunc(molecule::Molecule, all_coeffs::AbstractVector)
    return SCFMinFunc([molecule],all_coeffs);
end

function SCFMinFunc(molecules::AbstractVector{Molecule}, 
    all_coeffs::AbstractVector)
    tot_num_atoms = 0;
    tot_num_clouds = 0;
    atom_ind_base = [];
    for molecule in molecules
        if isempty(atom_ind_base)
            push!(atom_ind_base,1);
            tot_num_atoms += length(molecule.atoms_data[:,1]);
            tot_num_clouds += length(molecule.cloud_data[:,1]);
        else
            aux_num_atoms = size(molecule.atoms_data)[1];
            aux_num_clouds = size(molecule.cloud_data)[1];
            
            push!(atom_ind_base, tot_num_atoms + 1);
            tot_num_atoms += aux_num_atoms;
            tot_num_clouds += aux_num_clouds;
        end
    end

    num_vars = tot_num_atoms + 1;

    mol_1_num_atoms = size(molecules[1].atoms_data)[1];
    mol_1_num_clouds = size(molecules[1].cloud_data)[1];
    clouds_per_atom = Int(mol_1_num_clouds / mol_1_num_atoms);

    # Make the initial guess of the polarization coefficients and the 
    # chemical potential.
    # x0 = ones(Float64,num_vars);
    # x0[end] = 0.0;

    x0 = zeros(Float64,0);
    for molecule in molecules
        abs_ζ = abs.(molecule.cloud_data[1:clouds_per_atom:end,6]);
        x0 = vcat(x0, abs_ζ);
    end
    x0 = sqrt.(x0);
    x0 = vcat(x0,molecules[1].chem_μ);

    needs_casting = true;
    cpy_molecules = deepcopy(molecules);
    for i in eachindex(cpy_molecules)
        cpy_molecules[i].cloud_data[:,6] .= 1.0;
    end

    function MinFunc(aux_X::AbstractVector)
        # Error function that if the polarization coefficients and chemical 
        # potential are exact will always be equal to zero. This function 
        # always returns non-negative values.

        aux_type = typeof(aux_X[1]);
        if (needs_casting == true)
            needs_casting = false;

            for ii in eachindex(cpy_molecules)
                cpy_molecules[ii].cloud_data = 
                    aux_type.(cpy_molecules[ii].cloud_data);

                cpy_molecules[ii].atoms_data = 
                    aux_type.(cpy_molecules[ii].atoms_data);

                cpy_molecules[ii].chem_μ = 
                    aux_type(cpy_molecules[ii].chem_μ);
            end
        end

        # This ensures that the polarization coefficients are limited 
        # to non-negative numbers.
        aux_X[1:(end-1)] .^= 2.0;

        # Copy the trial polarization coefficients to the deepcopied 
        # vector of molecules.
        for ii in eachindex(cpy_molecules)
            ii0 = atom_ind_base[ii];
            num_atoms = size(cpy_molecules[ii].atoms_data)[1];
            cpy_molecules[ii].chem_μ = aux_X[end];

            ζ = aux_X[ii0:(ii0+num_atoms-1)];
            for i in 1:clouds_per_atom
                cpy_molecules[ii].cloud_data[i:clouds_per_atom:end,6] = ζ;
            end
        end

        # Get the system of equations associated with making the  
        # polarization Lagrangian stationary.
        _, _, aux_M, aux_Y = MatPolarizeMolecules(cpy_molecules, all_coeffs, 1);
        ret_vec = aux_M * aux_X - aux_Y;

        return dot(ret_vec,ret_vec);
    end

    cost_trace = zeros(Float64,0);

    # Callback function
    function trace_callback(state)
        push!(cost_trace, state.value);
        return false;
    end

    sol = Optim.optimize(MinFunc, x0, LBFGS(), autodiff=:forward,
        Optim.Options(show_trace=false, callback=trace_callback,
        iterations=8000));
    sol = Optim.minimizer(sol);
    sol[1:(end-1)] .^= 2.0;

    return cost_trace, sol;
end

function PolarizeMolecules!(molecule::Molecule, all_coeffs::AbstractVector)
    # Calculates the polarization coefficients of all atoms in the system.
    return PolarizeMolecules!([molecule], all_coeffs);
end

function PolarizeMolecules!(molecules::AbstractVector{Molecule}, 
    all_coeffs::AbstractVector)

    _, minimizer = SCFMinFunc(molecules,all_coeffs);

    tot_num_atoms = 0;
    tot_num_clouds = 0;
    atom_ind_base = [];
    for molecule in molecules
        if isempty(atom_ind_base)
            push!(atom_ind_base,1);
            tot_num_atoms += length(molecule.atoms_data[:,1]);
            tot_num_clouds += length(molecule.cloud_data[:,1]);
        else
            aux_num_atoms = size(molecule.atoms_data)[1];
            aux_num_clouds = size(molecule.cloud_data)[1];
            
            push!(atom_ind_base, tot_num_atoms + 1);
            tot_num_atoms += aux_num_atoms;
            tot_num_clouds += aux_num_clouds;
        end
    end

    for i in eachindex(molecules)
        i0 = atom_ind_base[i];
        num_atoms = size(molecules[i].atoms_data)[1];
        num_clouds = size(molecules[i].cloud_data)[1];
        clouds_per_atom = round(Int,num_clouds/num_atoms);

        molecules[i].chem_μ = minimizer[end];

        for ii in 1:num_atoms
            ζ = minimizer[i0+ii-1];

            for jj in 1:clouds_per_atom
                jj0 = (ii-1)*clouds_per_atom + jj;
                molecules[i].cloud_data[jj0,6] = ζ;
            end
        end
    end
end

function NaiveEnergyFromDensity(molecule::Molecule, all_coeffs::AbstractVector)
    # Returns the surface energy from the fitted electron clouds and nuclei
    # positions using the polarized naive model.
    return NaiveEnergyFromDensity([molecule], all_coeffs);
end

function NaiveEnergyFromDensity(molecules::AbstractVector{Molecule},
    all_coeffs::AbstractVector)
    # Returns the surface energy from the fitted electron clouds and nuclei
    # positions using the polarized naive model.
    global num_elems, factorials;

    num_molecules = length(molecules);
    
    aux_type = typeof(all_coeffs[1][1]);
    XCa_1b, XCc_1b, XCe_1b, XCf_1b = GetOneBodyCoeffs(all_coeffs[1]);
    XCa_2b, XCb_2b, XCc_2b, XCd_2b = GetTwoBodyCoeffs(all_coeffs[2]);

    energy = aux_type(0.0);
    kinetic_energy = aux_type(0.0);

    for ii in 1:num_molecules
        molecule1 = molecules[ii];
        num_atoms1 = size(molecule1.atoms_data)[1];
        num_clouds1 = size(molecule1.cloud_data)[1];
        clouds_per_atom = ceil(Int,num_clouds1/num_atoms1);

        # Kinetic Energy
        for i in 1:num_atoms1
            Z1 = round(Int,molecule1.atoms_data[i,5]);
            ζ1 = molecule1.cloud_data[i*clouds_per_atom,6];

            TF_fit_coeff = 0.0;
            TF_fit_coeff += XCe_1b[Z1,1]*(ζ1^(5.0/3.0));
            TF_fit_coeff += XCe_1b[Z1,2]*(ζ1^(8.0/3.0));

            vW_fit_coeff = 0.0;
            vW_fit_coeff += XCf_1b[Z1,1]*(ζ1^(1.0));
            vW_fit_coeff += XCf_1b[Z1,2]*(ζ1^(2.0));

            kinetic_energy += TF_fit_coeff*unpol_TF_KE[Z1];
            kinetic_energy += vW_fit_coeff*unpol_vW_KE[Z1];
        end

        # Intramolecular interactions
        # nuclei-nuclei
        for i in 1:num_atoms1
            for j in (i+1):num_atoms1
                q1 = molecule1.atoms_data[i,4];
                q2 = molecule1.atoms_data[j,4];

                aux_dist = molecule1.atoms_data[i,1:3];
                aux_dist -= molecule1.atoms_data[j,1:3];
                aux_dist = norm(aux_dist);

                energy += q1*q2*(1.0/aux_dist);
            end
        end

        # nuclei-cloud
        for i in 1:num_atoms1
            for j in 1:num_clouds1
                q1 = molecule1.atoms_data[i,4];
                c2 = molecule1.cloud_data[j,4];
                λ2 = molecule1.cloud_data[j,5];

                aux_dist = molecule1.atoms_data[i,1:3];
                aux_dist -= molecule1.cloud_data[j,1:3];
                aux_dist = norm(aux_dist);
                
                ζ2 = molecule1.cloud_data[j,6];

                # polarization
                c2 *= ζ2;

                energy -= q1*c2*NaiveCoulombInt(λ2,aux_dist);
            end
        end

        # cloud-cloud
        for i in 1:num_clouds1
            for j in i:num_clouds1
                c1 = molecule1.cloud_data[i,4];
                λ1 = molecule1.cloud_data[i,5];
                c2 = molecule1.cloud_data[j,4];
                λ2 = molecule1.cloud_data[j,5];

                ζ1 = molecule1.cloud_data[i,6];
                ζ2 = molecule1.cloud_data[j,6];

                ii0 = ceil(Int,i/clouds_per_atom);
                jj0 = ceil(Int,j/clouds_per_atom);

                aux_dist = molecule1.cloud_data[i,1:3];
                aux_dist -= molecule1.cloud_data[j,1:3];
                aux_dist = norm(aux_dist);
                
                λ = (λ1 * λ2) / (λ1 + λ2);

                if (ii0 == jj0) && (i == j)
                    c1 *= 0.5;
                end

                # polarization
                c1 *= ζ1;
                c2 *= ζ2;

                # Naive Model
                energy += c1*c2*NaiveCoulombInt(λ,aux_dist);
            end
        end

        # Intermolecular interactions
        for jj in (ii+1):num_molecules
            molecule2 = molecules[jj];
            num_atoms2 = size(molecule2.atoms_data)[1];
            num_clouds2 = size(molecule2.cloud_data)[1];

            # nuclei-nuclei
            for i in 1:num_atoms1
                for j in 1:num_atoms2
                    q1 = molecule1.atoms_data[i,4];
                    q2 = molecule2.atoms_data[j,4];
        
                    aux_dist = molecule1.atoms_data[i,1:3];
                    aux_dist -= molecule2.atoms_data[j,1:3];
                    aux_dist = norm(aux_dist);
        
                    energy += q1*q2*(1.0/aux_dist);
                end
            end

            # nuclei-cloud
            for i in 1:num_atoms1
                for j in 1:num_clouds2
                    q1 = molecule1.atoms_data[i,4];
                    c2 = molecule2.cloud_data[j,4];
                    λ2 = molecule2.cloud_data[j,5];

                    ζ2 = molecule2.cloud_data[j,6];

                    # polarization
                    c2 *= ζ2;

                    aux_dist = molecule1.atoms_data[i,1:3];
                    aux_dist -= molecule2.cloud_data[j,1:3];
                    aux_dist = norm(aux_dist);

                    energy -= q1*c2*NaiveCoulombInt(λ2,aux_dist);
                end
            end

            for i in 1:num_clouds1
                for j in 1:num_atoms2
                    Z1 = round(Int, 
                        molecule1.atoms_data[ceil(Int,i/clouds_per_atom),5]);
                    Z2 = round(Int, molecule2.atoms_data[j,5]);

                    c1 = molecule1.cloud_data[i,4];
                    λ1 = molecule1.cloud_data[i,5];
                    q2 = molecule2.atoms_data[j,4];

                    ζ1 = molecule1.cloud_data[i,6];

                    # polarization
                    c1 *= ζ1;

                    aux_dist = molecule1.cloud_data[i,1:3];
                    aux_dist -= molecule2.atoms_data[j,1:3];
                    aux_dist = norm(aux_dist);

                    energy -= c1*q2*NaiveCoulombInt(λ1,aux_dist);
                end
            end

            # cloud-cloud
            for i in 1:num_clouds1
                for j in 1:num_clouds2
                    Z1 = round(Int, 
                        molecule1.atoms_data[ceil(Int,i/clouds_per_atom),5]);
                    Z2 = round(Int, 
                        molecule2.atoms_data[ceil(Int,j/clouds_per_atom),5]);

                    c1 = molecule1.cloud_data[i,4];
                    λ1 = molecule1.cloud_data[i,5];
                    c2 = molecule2.cloud_data[j,4];
                    λ2 = molecule2.cloud_data[j,5];

                    ζ1 = molecule1.cloud_data[i,6];
                    ζ2 = molecule2.cloud_data[j,6];

                    λ = (λ1 * λ2) / (λ1 + λ2);

                    aux_dist = molecule1.cloud_data[i,1:3];
                    aux_dist -= molecule2.cloud_data[j,1:3];
                    aux_dist = norm(aux_dist);

                    # polarization
                    c1 *= ζ1;
                    c2 *= ζ2;

                    # Naive Model
                    energy += c1*c2*NaiveCoulombInt(λ,aux_dist);
                end
            end
        end
    end

    return energy + kinetic_energy;
end

function XCEnergyFromDensity(molecule::Molecule, all_coeffs::AbstractVector)
    # Returns the surface XC energy from the fitted polarized electron clouds 
    # and nuclei positions.
    return XCEnergyFromDensity([molecule], all_coeffs);
end

function XCEnergyFromDensity(molecules::AbstractVector{Molecule}, 
    all_coeffs::AbstractVector)
    # Returns the surface XC energy from the fitted polarized electron clouds 
    # and nuclei positions.
    global num_elems, factorials;
    
    aux_type = typeof(all_coeffs[1][1]);
    XCa_1b, XCc_1b, XCe_1b, XCf_1b = GetOneBodyCoeffs(all_coeffs[1]);
    XCa_2b, XCb_2b, XCc_2b, XCd_2b = GetTwoBodyCoeffs(all_coeffs[2]);

    num_molecules = length(molecules);
    energy = aux_type(0.0);

    for ii in 1:num_molecules
        molecule1 = molecules[ii];
        num_atoms1 = size(molecule1.atoms_data)[1];
        num_clouds1 = size(molecule1.cloud_data)[1];
        clouds_per_atom = ceil(Int,num_clouds1/num_atoms1);

        # Intramolecular interactions
        # nuclei-cloud 
        for i in 1:num_atoms1
            for j in 1:num_clouds1
                Z1 = round(Int,molecule1.atoms_data[i,5]);
                Z2 = round(Int,
                    molecule1.atoms_data[ceil(Int,j/clouds_per_atom),5]);

                q1 = molecule1.atoms_data[i,4];
                c2 = molecule1.cloud_data[j,4];
                λ2 = molecule1.cloud_data[j,5];

                ζ1 = molecule1.cloud_data[i*clouds_per_atom,6];
                ζ2 = molecule1.cloud_data[j,6];

                aux_dist = molecule1.atoms_data[i,1:3];
                aux_dist -= molecule1.cloud_data[j,1:3];
                aux_dist = norm(aux_dist);

                # polarization
                c2 *= ζ2;

                # XC contributions
                if (aux_dist > 1.0E-6)
                    xc_coeff_1 = CombXCFunc2B(XCc_2b,Z1,Z2,ζ1,ζ2);
                    xc_coeff_2 = CombXCFunc2B(XCd_2b,Z1,Z2,ζ1,ζ2);
                    
                    energy += xc_coeff_1*q1*c2*XC_Sph(λ2,aux_dist);
                    energy += xc_coeff_2*q1*c2*XC_Cyl(λ2,aux_dist);
                else
                    xc_coeff_1 = CombXCFunc1B(XCc_1b,Z1,ζ1);

                    energy += xc_coeff_1*q1*c2*XC_Sph(λ2,aux_dist);
                end
            end
        end

        # cloud-cloud
        for i in 1:num_clouds1
            for j in i:num_clouds1
                Z1 = round(Int,
                    molecule1.atoms_data[ceil(Int,i/clouds_per_atom),5]);
                Z2 = round(Int,
                    molecule1.atoms_data[ceil(Int,j/clouds_per_atom),5]);

                c1 = molecule1.cloud_data[i,4];
                λ1 = molecule1.cloud_data[i,5];
                c2 = molecule1.cloud_data[j,4];
                λ2 = molecule1.cloud_data[j,5];

                ζ1 = molecule1.cloud_data[i,6];
                ζ2 = molecule1.cloud_data[j,6];

                ii0 = ceil(Int,i/clouds_per_atom);
                jj0 = ceil(Int,j/clouds_per_atom);

                aux_dist = molecule1.cloud_data[i,1:3];
                aux_dist -= molecule1.cloud_data[j,1:3];
                aux_dist = norm(aux_dist);

                λ = (λ1 * λ2) / (λ1 + λ2);

                if (ii0 == jj0) && (i == j)
                    c1 *= 0.5;
                end

                # polarization
                c1 *= ζ1;
                c2 *= ζ2;

                # XC contributions
                if (aux_dist > 1.0E-6)
                    xc_coeff_1 = CombXCFunc2B(XCa_2b,Z1,Z2,ζ1,ζ2);
                    xc_coeff_2 = CombXCFunc2B(XCb_2b,Z1,Z2,ζ1,ζ2);
                    
                    energy += xc_coeff_1*c1*c2*XC_Sph(λ,aux_dist);
                    energy += xc_coeff_2*c1*c2*XC_Cyl(λ,aux_dist);
                else
                    xc_coeff_1 = CombXCFunc1B(XCa_1b,Z1,ζ1);

                    energy += xc_coeff_1*c1*c2*XC_Sph(λ,aux_dist);
                end
            end
        end

        # Intermolecular interactions
        for jj in (ii+1):num_molecules
            molecule2 = molecules[jj];
            num_atoms2 = size(molecule2.atoms_data)[1];
            num_clouds2 = size(molecule2.cloud_data)[1];

            # nuclei-cloud 
            for i in 1:num_atoms1
                for j in 1:num_clouds2
                    Z1 = round(Int,molecule1.atoms_data[i,5]);
                    Z2 = round(Int,
                        molecule2.atoms_data[ceil(Int,j/clouds_per_atom),5]);

                    q1 = molecule1.atoms_data[i,4];
                    c2 = molecule2.cloud_data[j,4];
                    λ2 = molecule2.cloud_data[j,5];

                    ζ1 = molecule1.cloud_data[i*clouds_per_atom,6];
                    ζ2 = molecule2.cloud_data[j,6];

                    # polarization
                    c2 *= ζ2;

                    aux_dist = molecule1.atoms_data[i,1:3];
                    aux_dist -= molecule2.cloud_data[j,1:3];
                    aux_dist = norm(aux_dist);    

                    # XC contributions
                    xc_coeff_1 = CombXCFunc2B(XCc_2b,Z1,Z2,ζ1,ζ2);
                    xc_coeff_2 = CombXCFunc2B(XCd_2b,Z1,Z2,ζ1,ζ2);
                    
                    energy += xc_coeff_1*q1*c2*XC_Sph(λ2,aux_dist);
                    energy += xc_coeff_2*q1*c2*XC_Cyl(λ2,aux_dist);
                end
            end

            for i in 1:num_clouds1
                for j in 1:num_atoms2
                    Z1 = round(Int,
                        molecule1.atoms_data[ceil(Int,i/clouds_per_atom),5]);
                    Z2 = round(Int,molecule2.atoms_data[j,5]);

                    c1 = molecule1.cloud_data[i,4];
                    λ1 = molecule1.cloud_data[i,5];
                    q2 = molecule2.atoms_data[j,4];

                    ζ1 = molecule1.cloud_data[i,6];
                    ζ2 = molecule2.cloud_data[j*clouds_per_atom,6];

                    # polarization
                    c1 *= ζ1;

                    aux_dist = molecule1.cloud_data[i,1:3];
                    aux_dist -= molecule2.atoms_data[j,1:3];
                    aux_dist = norm(aux_dist);

                    # XC contributions
                    xc_coeff_1 = CombXCFunc2B(XCc_2b,Z1,Z2,ζ1,ζ2);
                    xc_coeff_2 = CombXCFunc2B(XCd_2b,Z1,Z2,ζ1,ζ2);
                    
                    energy += xc_coeff_1*c1*q2*XC_Sph(λ1,aux_dist);
                    energy += xc_coeff_2*c1*q2*XC_Cyl(λ1,aux_dist);
                end
            end

            # cloud-cloud
            for i in 1:num_clouds1
                for j in 1:num_clouds2
                    Z1 = round(Int,
                        molecule2.atoms_data[ceil(Int,j/clouds_per_atom),5]);
                    Z2 = round(Int,
                        molecule1.atoms_data[ceil(Int,i/clouds_per_atom),5]);
                    
                    c1 = molecule1.cloud_data[i,4];
                    λ1 = molecule1.cloud_data[i,5];
                    c2 = molecule2.cloud_data[j,4];
                    λ2 = molecule2.cloud_data[j,5];

                    ζ1 = molecule1.cloud_data[i,6];
                    ζ2 = molecule2.cloud_data[j,6];

                    λ = (λ1 * λ2) / (λ1 + λ2);

                    aux_dist = molecule1.cloud_data[i,1:3];
                    aux_dist -= molecule2.cloud_data[j,1:3];
                    aux_dist = norm(aux_dist);

                    # polarization
                    c1 *= ζ1;
                    c2 *= ζ2;

                    # XC contributions
                    xc_coeff_1 = CombXCFunc2B(XCa_2b,Z1,Z2,ζ1,ζ2);
                    xc_coeff_2 = CombXCFunc2B(XCb_2b,Z1,Z2,ζ1,ζ2);
                    
                    energy += xc_coeff_1*c1*c2*XC_Sph(λ,aux_dist);
                    energy += xc_coeff_2*c1*c2*XC_Cyl(λ,aux_dist);
                end
            end
        end
    end

    return energy;
end

function GetPartialCharges(mol::Molecule)
    # Prints the partial charges of the molecule accounting for its polarization 
    # coefficients.
    num_clouds = size(mol.cloud_data)[1];
    num_atoms = size(mol.atoms_data)[1];
    clouds_per_atom = round(Int, num_clouds / num_atoms);
    
    aux_type = eltype(mol.atoms_data[1,1]);
    charges = zeros(aux_type, num_atoms);
        
    for i in 1:num_atoms
        Z = mol.atoms_data[i,4];
        charge = Z;

        for j in 1:clouds_per_atom
            cloud_index = (i-1)*clouds_per_atom + j;
            c = mol.cloud_data[cloud_index,4];
            ζ = mol.cloud_data[cloud_index,6];
            charge -= c * ζ;
        end

        charges[i] = charge;
    end

    return charges;
end
