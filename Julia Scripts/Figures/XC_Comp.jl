using Plots, LaTeXStrings, Measures;
using Polynomials, SpecialPolynomials;
using SpecialFunctions, ForwardDiff;

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

λ = 7;

r = collect(0:0.001:3.0);
e_sphe_t10 = xc_sph_trunc.(10,λ,r);
e_sphe_t20 = xc_sph_trunc.(20,λ,r);
e_sphe_exact = xc_sph.(λ,r);

e_cyl_t10 = xc_cyl_trunc.(10,λ,r);
e_cyl_t20 = xc_cyl_trunc.(20,λ,r);
e_cyl_exact = xc_cyl.(λ,r);

p1 = plot(r, e_sphe_t10, label="Truncated Sum (K = 10)", linewidth = 2);
plot!(r, e_sphe_t20, label="Truncated Sum (K = 20)", linewidth = 2);
plot!(r, e_sphe_exact, label="Exact (K = ∞)", linewidth = 2, linestyle = :dot);

plot!(xlims=[0,3], framestyle = :box);
plot!(xlabel=L"d");
plot!(ylims=[-2,4]);
plot!(yticks=-2:2:4);
plot!(ylabel=L"$\mathrm{XC}_\mathrm{Sph}^\mathrm{EN} (d)$");

l_x_pos = 2.5;
l_y_pos = -2 + (4 - (-2)) * (1.0/6.0);
annotate!(l_x_pos, l_y_pos, text(L"$\lambda = 7$", :center, 10));

p2 = plot(r, e_cyl_t10, label="Truncated Sum (K = 10)", linewidth = 2);
plot!(r, e_cyl_t20, label="Truncated Sum (K = 20)", linewidth = 2);
plot!(r, e_cyl_exact, label="Exact (K = ∞)", linewidth = 2, linestyle = :dot);

plot!(xlims=[0.0,3], framestyle = :box);
plot!(xlabel=L"d");
plot!(ylabel=L"$\mathrm{XC}_\mathrm{Cyl}^\mathrm{EN} (d)$");
plot!(xticks=0:1:3);
plot!(ylims=[-2,4]);
plot!(yticks=-2:2:4);

l_x_pos = 2.5;
l_y_pos = -2 + (4 - (-2)) * (1.0/6.0);
annotate!(l_x_pos, l_y_pos, text(L"$\lambda = 7$", :center, 10));

plot(p1,p2,layout=(1,2), size = (1100, 260), 
    left_margin = [7mm 7mm], bottom_margin = [7mm 7mm])
savefig("Figures/XC_exact_comp_old.pdf");

function XC_Sph_Trunc(K::Integer, λ::Real, d::Real)
    # New definition
    aux_ret = 0.0;

    for k in 1:K
        aux_val = 2*(λ^k)*(basis(Hermite,2*k-1)(sqrt(λ)*d));
        aux_val /= sqrt(π)*factorial(big(2*k-1));
        aux_val *= exp(-λ*(d^2));
        aux_ret -= aux_val;
    end

    for k in 1:K
        aux_val = 2*(λ^(k+0.5))*(basis(Hermite,2*k)(sqrt(λ)*d));
        aux_val /= sqrt(π)*factorial(big(2*k));
        aux_val *= exp(-λ*(d^2));
        aux_ret += aux_val;
    end

    return aux_ret;
end

function XC_Sph(λ::Real, d::Real)
    return 2*sqrt(λ/π)*(exp(-2*λ*d-λ)-1)/exp(λ*d^2);
end

function XC_Cyl_Trunc(K::Integer, λ::Real, d::Real)
    # New definition
    function foo(d::Number)
        return XC_Sph_Trunc(K,λ,d);
    end

    return ForwardDiff.derivative(foo,d);
end

function XC_Cyl(λ::Real, d::Real)
    return (4*exp(-(d^2*λ)-2*d*λ-λ)*(exp(2*d*λ+λ)*d-d-1)*λ^(3/2))/sqrt(π);
end

λ = 7;

r = collect(0.0:0.01:3);
xc_sph_K10 = -XC_Sph_Trunc.(10, λ, r);
xc_sph_K20 = -XC_Sph_Trunc.(20, λ, r);
xc_sph_exact = -XC_Sph.(λ, r);

p1 = plot(r, xc_sph_K10, label="Truncated Sum (K = 10)", linewidth = 2);
plot!(r, xc_sph_K20, label="Truncated Sum (K = 20)", linewidth = 2);
plot!(r, xc_sph_exact, label="Exact (K = ∞)", linewidth = 2, linestyle=:dot);

plot!(ylims=[-3,6],yticks=-3:3:6,xlims=[0,3]);
plot!(framestyle = :box);
plot!(xlabel=L"d");
plot!(ylabel=L"$- \mathrm{XC}_\mathrm{Sph}^\mathrm{EN} (d)$");

l_x_pos = 2.5;
l_y_pos = -3 + (6 - (-3)) * (1.0/6.0);
annotate!(l_x_pos, l_y_pos, text(L"$\lambda = 7$", :center, 10));

r = collect(0.0:0.01:3);
xc_cyl_K10 = XC_Cyl_Trunc.(10, λ, r);
xc_cyl_K20 = XC_Cyl_Trunc.(20, λ, r);
xc_cyl_exact = XC_Cyl.(λ, r);

p2 = plot(r, xc_cyl_K10, label="Truncated Sum (K = 10)", linewidth = 2);
plot!(r, xc_cyl_K20, label="Truncated Sum (K = 20)", linewidth = 2);
plot!(r, xc_cyl_exact, label="Exact (K = ∞)", linewidth = 2, linestyle=:dot);

plot!(ylims=[-5,10],xlims=[0,3]);
plot!(xticks=0:1:3, yticks=-5:5:10, framestyle = :box);
plot!(xlabel=L"d");
plot!(ylabel=L"$\mathrm{XC}_\mathrm{Cyl}^\mathrm{EN} (d)$");

l_x_pos = 2.5;
l_y_pos = -5 + (10 - (-5)) * (1.0/6.0);
annotate!(l_x_pos, l_y_pos, text(L"$\lambda = 7$", :center, 10));

plot(p1,p2,layout=(1,2), size = (1100, 260), 
    left_margin = [7mm 7mm], bottom_margin = [7mm 7mm])
savefig("Figures/XC_exact_comp_new.pdf")
