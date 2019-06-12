
# using Calculus
# define types first
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
    stens::Function
    Pr::Function
end

function substance(name, x::Vector{Function})
    dens, enth, Cp, visc, thermCond, stens = x
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
    return substance(name, dens, enth, Cp, visc, thermCond, alt_heat_capacity, visc_kin, exp_coeff, thermDiff, stens, Pr)
end

function na_enthalpy(p, T)
    """This function takes the temperature of liquid sodium in K and returns the enthalpy in kJ/kg.
    The reference state is taken as solid sodium at 298 K. The correlation takes this form:
    H(l, T) - H(s, 298.15) = 365.77 + 1.6582 T - 4.2395 × 10^-4 * T^2 + 1.4847 × 10^-7 * T^3 + 2992.6 * T^-1
    This correlation is recommended by a report prepared by Argonne National Lab that can be found here:
    http://www.ne.anl.gov/eda/ANL-RE-95-2.pdf
    The correlation is valid between 371 and 2000 K."""
    a = -365.77
    b = 1.6582
    c = 4.2395E-4
    d = 1.4847E-7
    e = 2992.6
    return (a + b * T - c * T ^ 2 + d * T ^ 3 + e / T) * 1000
end

function na_Cp(p, T)
    """This function takes the temperature of liquid sodium and returns the heat
    capacity at constant pressure in kJ/kg*K. The Argonne report gives a very
    complicated correlation for Cp but notes that it's almost equal to the
    derivative of enthalpy with temperature up to about 1800 K. So I used the derivative here."""
    a = -365.77
    b = 1.6582
    c = 4.2395E-4
    d = 1.4847E-7
    e = 2992.6
    return (b - 2*c * T + 3 * d * T ^ 2 - e / (T^2)) * 1000
    #return derivative(enthalpy, T, dx = 1E-6)
end

function na_density(p, T)
    """This function takes the temperature of liquid sodium in K and returns the density in kg/m^3.
    The correlation takes this form:
    rho = rho_c + f * (1 - Tr) + g * (1 - Tr)^h
    where rho_c is the critical density of sodium, Tr is the reduced temperature T/T_c,
    and f g and h are constants given in the function.
    This correlation is recommended by a report prepared by Argonne National Lab that can be found here:
    http://www.ne.anl.gov/eda/ANL-RE-95-2.pdf
    The correlation is valid between 371 and ~2000 K."""
    rho_c = 219
    f = 275.32
    g = 511.58
    h = 0.5
    Tc = 2503.7
    Tr = T / Tc
    return rho_c + f * (1 - Tr) + g * (1 - Tr) ^ h
end

function na_thermal_conductivity(p, T)
    """This function takes the temperature of liquid sodium in K and returns the thermal conductivity in W/m*K.
    The reference state is taken as solid sodium at 298 K. The correlation takes this form:
    k = 124.67 - 0.11381 * T + 5.5226 × 10^-5 * T^2 - 1.1842 × 10^-8 * T^3
    This correlation is recommended by a report prepared by Argonne National Lab that can be found here:
    http://www.ne.anl.gov/eda/ANL-RE-95-2.pdf
    The correlation is valid between 371 and 2000 K."""
    a = 124.67
    b = 0.11381
    c = 5.5226E-5
    d = 1.1842E-8
    return a - b * T + c * T ^ 2 - d * T ^ 3
end

function na_dynamic_viscosity(p, T)
    """This function takes the temperature of liquid sodium in K and returns the dynamic viscosity in Pa * s.
    The reference state is taken as solid sodium at 298 K. The correlation takes this form:
    ln eta = - 6.4406 - 0.3958 ln (T) + 556.835 / T
    This correlation is recommended by a report prepared by Argonne National Lab that can be found here:
    http://www.ne.anl.gov/eda/ANL-RE-95-2.pdf
    The correlation is valid between 371 and 2500 K."""
    a = -6.4406
    b = -0.3958
    c = 556.835
    return exp(a + b * log(T) + c / T)
end

