using Plots, LaTeXStrings, Measures;
using Polynomials, SpecialPolynomials;
using SpecialFunctions;

function xc_sph_trunc(K::Integer, λ::Real, d::Real)
    # Old definition
    aux_ret = 0.0;

    for k in 2:(2*K+1)
        aux_val = -2*(sqrt(λ)*(1im))^k;
        aux_val /= d*sqrt(π)*factorial(big(k));
        aux_val *= basis(Hermite,k-1)(sqrt(λ)*d);
        aux_val *= exp(-λ*(d^2))
        aux_ret += real(aux_val);
    end

    return aux_ret;
end

function xc_cyl_trunc(K::Integer, λ::Real, d::Real)
    # Old definition
    aux_ret = 0.0;

    for k in 3:(2*K+2)
        aux_val = 2*(sqrt(λ)*(1im))^(k-1);
        aux_val /= (d^2)*sqrt(π)*(factorial(big(k)));
        aux_val *= basis(Hermite,k-2)(sqrt(λ)*d);
        aux_val *= exp(-λ*(d^2))

        aux_ret += real(aux_val);
    end

    for k in 3:(2*K+2)
        aux_val = 2*(sqrt(λ)*(1im))^k;
        aux_val /= d*sqrt(π)*(factorial(big(k)));
        aux_val *= basis(Hermite,k-1)(sqrt(λ)*d);
        aux_val *= exp(-λ*(d^2));

        aux_ret += imag(aux_val);
    end

    return aux_ret;
end

function xc_sph(λ::Real, d::Real)
    # Old definition
    return real((erf(d*sqrt(λ)+(1im)*sqrt(λ))+erf(d*sqrt(λ)-
        (1im)*sqrt(λ)))/(2*d)-erf(d*sqrt(λ))/d);
end

function xc_cyl(λ::Real, d::Real)
    # Old definition
    return real(exp(-λ*d^2)/(d*sqrt(π))*((exp(λ)*sin(2*d*λ))/(d*sqrt(λ))-
        2*sqrt(λ))-(erf(sqrt(λ)*(d+1im))+erf(sqrt(λ)*(d-1im))-
        2*erf(sqrt(λ)*d))/(2*d^2));
end

r = collect(0:0.001:3.0);
e_sphe_t10_l8_λ8 = xc_sph_trunc.(10,8,r);
e_sphe_t20_l8_λ8 = xc_sph_trunc.(20,8,r);
e_sphe_exact_l8_λ8 = xc_sph.(8,r);

e_cyl_t10_l8_λ8 = xc_cyl_trunc.(10,8,r);
e_cyl_t20_l8_λ8 = xc_cyl_trunc.(20,8,r);
e_cyl_exact_l8_λ8 = xc_cyl.(8,r);

p1 = plot(r, e_sphe_t10_l8_λ8, label="Truncated Sum (K = 10)", linewidth = 2);
plot!(r, e_sphe_t20_l8_λ8, label="Truncated Sum (K = 20)", linewidth = 2);
plot!(r, e_sphe_exact_l8_λ8, label="Exact", linewidth = 2,
    linestyle = :dot);

plot!(xlims=[0,3], framestyle = :box);
plot!(xticks=(0:1:3,[]));
plot!(ylims=[-2,4]);
plot!(yticks=-2:2:4);
plot!(ylabel=L"$\mathrm{XC}_\mathrm{Sph}^\mathrm{EN} (d)$");

l_x_pos = 2.5;
l_y_pos = -2 + (4 - (-2)) * (1.0/6.0);
annotate!(l_x_pos, l_y_pos, text(L"$\lambda = 8$", :center, 10));

p2 = plot(r, e_cyl_t10_l8_λ8, label="Truncated Sum (K = 10)", linewidth = 2);
plot!(r, e_cyl_t20_l8_λ8, label="Truncated Sum (K = 20)", linewidth = 2);
plot!(r, e_cyl_exact_l8_λ8, label="Exact", linewidth = 2,
    linestyle = :dot);

plot!(xlims=[0.0,3], framestyle = :box);
plot!(xlabel=L"d");
plot!(ylabel=L"$\mathrm{XC}_\mathrm{Cyl}^\mathrm{EN} (d)$");
plot!(xticks=0:1:3);
plot!(ylims=[-2,4]);
plot!(yticks=-2:2:4);

l_x_pos = 2.5;
l_y_pos = -2 + (4 - (-2)) * (1.0/6.0);
annotate!(l_x_pos, l_y_pos, text(L"$\lambda = 8$", :center, 10));

plot(p1,p2,layout=(2,1), size = (500, 400))
savefig("Figures/XC_exp-poly_comp.pdf");

function XC_Sph_Trunc(K::Integer, λ::Real, d::Real)
    # New definition
    aux_ret = 0.0;

    for k in 1:K
        aux_val = 2*(λ^k)*(basis(Hermite,2*k-1)(sqrt(λ)*d));
        aux_val /= d*sqrt(π)*factorial(big(2*k-1));
        aux_val *= exp(-λ*(d^2));
        aux_ret -= aux_val;
    end

    return aux_ret;
end

function XC_Sph(λ::Real, d::Real)
    return (1/d)*sqrt(λ/π)*((1-exp(4*λ*d))/(exp(λ*(d+1)^2)));
end

function XC_Cyl_Trunc(K::Integer, λ::Real, d::Real)
    # New definition
    aux_ret = 0.0;

    for k in 1:K
        aux_val = 2*(λ^(k+0.5))*(basis(Hermite,2*k)(sqrt(λ)*d));
        aux_val /= d*sqrt(π)*factorial(big(2*k));
        aux_val *= exp(-λ*(d^2));
        aux_ret += aux_val;
    end

    return aux_ret;
end

function XC_Cyl(λ::Real, d::Real)
    return (1/d)*sqrt(λ/π)*((((exp(2*λ*d)-1)^2)/(exp(λ*(d+1)^2))) + 
        ((2*(exp(-λ)-1))/(exp(λ*d^2))));
end

λ = 8;

r = collect(0:0.01:3);
xc_sph_K10 = -XC_Sph_Trunc.(10, λ, r);
xc_sph_K20 = -XC_Sph_Trunc.(20, λ, r);
xc_sph_exact = -XC_Sph.(λ, r);

p1 = plot(r, xc_sph_K10, label="Truncated Sum (K = 10)", linewidth = 2);
plot!(r, xc_sph_K20, label="Truncated Sum (K = 20)", linewidth = 2);
plot!(r, xc_sph_exact, label="Exact", linewidth = 2, linestyle=:dot);
plot!(ylims=[-1,2],xlims=[0,3]);
plot!(xticks=(0:1:3,[]), framestyle = :box);
plot!(ylabel=L"$- \mathrm{XC}_\mathrm{Sph}^\mathrm{EN} (d)$");

l_x_pos = 2.5;
l_y_pos = -1 + (2 - (-1)) * (1.0/6.0);
annotate!(l_x_pos, l_y_pos, text(L"$\lambda = 8$", :center, 10));

r = collect(0:0.01:3);
xc_cyl_K10 = -XC_Cyl_Trunc.(10, λ, r);
xc_cyl_K20 = -XC_Cyl_Trunc.(20, λ, r);
xc_cyl_exact = -XC_Cyl.(λ, r);

p2 = plot(r, xc_cyl_K10, label="Truncated Sum (K = 10)", linewidth = 2);
plot!(r, xc_cyl_K20, label="Truncated Sum (K = 20)", linewidth = 2);
plot!(r, xc_cyl_exact, label="Exact", linewidth = 2, linestyle=:dot);
plot!(ylims=[-3,6],xlims=[0,3]);
plot!(xticks=0:1:3, framestyle = :box);
plot!(xlabel=L"d");
plot!(ylabel=L"$- \mathrm{XC}_\mathrm{Cyl}^\mathrm{EN} (d)$");

l_x_pos = 2.5;
l_y_pos = -3 + (6 - (-3)) * (1.0/6.0);
annotate!(l_x_pos, l_y_pos, text(L"$\lambda = 8$", :center, 10));

plot(p1,p2,layout=(2,1), size = (500, 400))
savefig("Figures/XC_exact_comp.pdf")
