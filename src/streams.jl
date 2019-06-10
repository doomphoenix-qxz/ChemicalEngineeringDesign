struct stream{T_<:Real}
    subst_::substance
    ṁ::T_
    V̇::T_
    p::T_
    T::T_
    ρ::T_
    Ĥ::T_
    H::T_
    Cp::T_
    μ::T_
    k::T_
    ν::T_
    α::T_
    Pr::T_
    β::T_
end

function stream(subst_::substance, ṁ, p, T)
    ρ = subst_.density(p, T)
    V̇ = ṁ / ρ
    Ĥ = subst_.enthalpy(p, T)
    H = Ĥ * ṁ
    Cp = subst_.heat_capacity(p, T)
    μ = subst_.viscosity(p, T)
    k = subst_.thermal_conductivity(p, T)
    ν = subst_.visc_kin(p, T)
    α = subst_.therm_diff(p, T)
    Pr = subst_.Pr(p, T)
    β = subst_.exp_coeff(p, T)
    return stream(subst_, ṁ, V̇, p, T, ρ, Ĥ, H, Cp, μ, k, ν, α, Pr, β)
end

# function stream(subst_::String, ṁ, p, T)
#     return stream(stream_substances[subst_], ṁ, p, T)
# end


struct flowstream{T<:Real}
    tpipe::pipes.pipe
    tstream::streams.stream
    Re::T
    Pr::T
    f::T
    ΔP::T
    Δh::T
end

"""
Flowstream: Basically takes in a stream and calculates relevant parameters for
its flow through a pipe, including the Reynolds and Prandtl numbers, the
friction factor, and pressure/head losses. Assumes negligible change in fluid
properties (notably density and viscosity) down the length of the pipe.

Uses well-known equations that can be found in (or derived from) any good book
on fluid mechanics.
"""
function flowstream(tpipe::pipes.pipe, tstream::streams.stream)
    d = tpipe.diameter
    Re = tstream.ṁ * d / (tstream.μ * tpipe.flowarea)
    f = pipes.getf(Re, tpipe)
    vel = tstream.V̇ / tpipe.flowarea
    Δh = (f*tpipe.length / d) * vel^2 / (2*9.81) + tpipe.elevation_change
    ΔP = Δh * (9.81 * tstream.ρ)
    return flowstream(tpipe, tstream, Re, tstream.Pr, f, ΔP, Δh)
end
