struct Stream{T_<:Real}
    subst_::Substance
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

function Stream(subst_::Substance, ṁ, p, T)
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
    return Stream(subst_, ṁ, V̇, p, T, ρ, Ĥ, H, Cp, μ, k, ν, α, Pr, β)
end

# function Stream(subst_::String, ṁ, p, T)
#     return Stream(Stream_substances[subst_], ṁ, p, T)
# end


struct FlowStream{T<:Real}
    tpipe::Pipe
    tstream::Stream
    Re::T
    Pr::T
    f::T
    ΔP::T
    Δh::T
end

"""
FlowStream: Basically takes in a Stream and calculates relevant parameters for
its flow through a pipe, including the Reynolds and Prandtl numbers, the
friction factor, and pressure/head losses. Assumes negligible change in fluid
properties (notably density and viscosity) down the length of the pipe.

Uses well-known equations that can be found in (or derived from) any good book
on fluid mechanics.
"""
function FlowStream(tpipe::Pipe, tStream::Stream)
    d = tpipe.diameter
    Re = tStream.ṁ * d / (tStream.μ * tpipe.flowarea)
    f = getf(Re, tpipe)
    vel = tStream.V̇ / tpipe.flowarea
    Δh = (f*tpipe.length / d) * vel^2 / (2*9.81) + tpipe.elevation_change
    ΔP = Δh * (9.81 * tStream.ρ)
    return FlowStream(tpipe, tStream, Re, tStream.Pr, f, ΔP, Δh)
end
