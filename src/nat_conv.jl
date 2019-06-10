# Define types first
struct substance
    name::String
    density::Function
    enthalpy::Function
    heat_capacity::Function
    viscosity::Function
    thermal_conductivity::Function
    alt_heat_capacity::Function
    visc_kin::Function
    exp_coeff::Function
    therm_diff::Function
    Pr::Function
end

function substance(name, x::Vector{Function})
    dens, enth, Cp, visc, thermCond = x
    function alt_heat_capacity(p, T)
        return derivative(T_ -> enth(p, T_), T)
    end
    function visc_kin(p, T)
        return visc(p, T)/dens(p, T)
    end
    function exp_coeff(p, T)
        dRdT = derivative(T_ -> dens(p, T_), T)
        return (-1/dens(p, T))*dRdT
    end
    function thermDiff(p, T)
        if Cp == nothing
            Cp_ = alt_heat_capacity(p, T)
        else
            Cp_ = Cp(p, T)
        end
        return thermCond(p, T) / (dens(p,T)*Cp_)
    end
    function Pr(p, T)
        if Cp == nothing
            Cp_ = alt_heat_capacity(p, T)
        else
            Cp_ = Cp(p, T)
        end
        return Cp_ * visc(p, T) / thermCond(p, T)
    end
    return substance(name, dens, enth, Cp, visc, thermCond, alt_heat_capacity, visc_kin, exp_coeff, thermDiff, Pr)
end

@enum hx_type conc_parallel=0 conc_counter=1 st_oneshell=2 st_twoshell=3 Cr_0=4

function ϵ_Cp(NTU, Cr)
    top = 1 - exp(-NTU *(1+Cr))
    return top/(1+Cr)
end

function ϵ_Cc(NTU, Cr)

    if Cr==1
        return NTU/(1+NTU)
    end

    top = 1 - exp(-NTU *(1+Cr))
    bottom = 1 - Cr*exp(-NTU *(1+Cr))
    return top/bottom

end

function ϵ_St₁(NTU, Cr)
    term = exp(-NTU * √(1+Cr^2))
    return 2*(1 + Cr + √(1+Cr^2) * (1+term)/(1-term))
end

function ϵ_Stₙ(NTU, Cr, n)
    term = ((1 - ϵ_St₁(NTU, Cr)*Cr) / (1 - ϵ_St₁(NTU, Cr)))^n
    return (term - 1)/(term - Cr)
end

function ϵ_Cr0(NTU, Cr=0)
    return 1 - exp(-NTU)
end

ϵ_functions = Dict(conc_parallel => ϵ_Cp, conc_counter => ϵ_Cc, st_oneshell => ϵ_St₁, st_twoshell => ϵ_Stₙ, Cr_0 => ϵ_Cr0)

function ϵ(NTU, Cr, hx::hx_type)
    return ϵ_functions[hx](NTU, Cr)
end
