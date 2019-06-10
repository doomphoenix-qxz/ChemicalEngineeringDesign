# meos.jl
################################################################################
# Note: The following is based on the nomenclature in the paper
# "Revised Release on the IAPWS Formulation 1995 for the Thermodynamic
# Properties of Ordinary Water Substance for General and Scientific Use"
# by IAPWS, the International Association for the Properties of Water and Steam.
#
# It is more or less as follows:
#
# ϕ : Dimensionless Helmholtz energy as a fucntion of density and temperature
# δ : Reduced density (= ρ/ρc)
# τ : Inverse reduced temperature (= Tc/T)
# n, d, t, α, β, γ, ϵ, A, B, C, D : these are parameters used in fitting the
# equation of state to data. As far as we're concerned, they are constants with
# no real physical significance.
#
# Note that although the code was written by consulting the IAPWS paper, it is
# applicable to any multiparameter equation of state with the same form (notably
# Span and Wagner's equation of state for carbon dioxide, which is also provided
# by this package; see R. Span and W. Wagner, J. Phys. Chem. Ref. Data Vol 25,
# No. 6, 1996 Titled "A New Equation of State for Carbon Dioxide Covering the
# Fluid Region from the Triple-Point Temperature to 1100 K at Pressures up to
# 800 MPa")
################################################################################

include("substances.jl")

__precompile__()
module meos
#import ForwardDiff
import Calculus
import pr
#using Optim
using Roots.fzero
struct multiparameos{T<:AbstractFloat}
    Tc::T
    ρc::T
    pc::T
    ω::T
    R::T

    #parameters for ideal gas portion
    n₀::Vector{T}
    γ₀::Vector{T}

    # parameters for residual term 1
    n₁::Vector{T}
    d₁::Vector{T}
    t₁::Vector{T}

    # parameters for residual term 2
    n₂::Vector{T}
    d₂::Vector{T}
    t₂::Vector{T}
    c₂::Vector{T}

    # parameters for residual term 3
    n₃::Vector{T}
    d₃::Vector{T}
    t₃::Vector{T}
    α₃::Vector{T}
    β₃::Vector{T}
    γ₃::Vector{T}
    ϵ₃::Vector{T}

    # parameters for residual term 4
    n₄::Vector{T}
    a₄::Vector{T}
    b₄::Vector{T}
    β₄::Vector{T}
    A₄::Vector{T}
    B₄::Vector{T}
    C₄::Vector{T}
    D₄::Vector{T}

    # Functions needed to get density guess values
    psat_guess::Function
    PR_::pr.pengRobinson

end

function Σ(expr)
    return sum(eval(@. expr))
end

function ϕ₀(eos::multiparameos, δ, τ)
    first = log(δ) + eos.n₀[1] + eos.n₀[2]*τ + eos.n₀[3]*log(τ)
    end_ = sum(eos.n₀[4:8].*log.(1-exp.(-eos.γ₀.*τ)))
    return first + end_
end

function res_1(eos::multiparameos, δ, τ)
    return sum(eos.n₁ .* (δ .^ eos.d₁) .* (τ .^ eos.t₁))
end

function res_2(eos::multiparameos, δ, τ)
    return sum(eos.n₂ .* (δ .^ eos.d₂) .* (τ .^ eos.t₂) .* exp.(-δ .^ eos.c₂))
end

function res_3(eos::multiparameos, δ, τ)
    return sum(eos.n₃ .* (δ .^ eos.d₃) .* (τ .^ eos.t₃) .* exp.(-eos.α₃ .* (δ - eos.ϵ₃).^2 - eos.β₃ .* (τ - eos.γ₃).^2))
end

function θ(eos::multiparameos, δ, τ)
    return sum((1-τ) + eos.A₄.*((δ - 1).^2).^(1 ./(2 .* eos.β₄)))
end

function Δ(eos::multiparameos, δ, τ)
    return sum(θ(eos, δ, τ).^2 + eos.B₄ .* ((δ - 1).^2) .* eos.a₄)
end

function Ψ(eos::multiparameos, δ, τ)
    return sum(exp(-eos.C₄ .* (δ - 1)^2 - eos.D₄ .* (τ - 1)^2))
end

function res_4(eos::multiparameos, δ, τ)
    return sum(eos.n₄ .* Δ(eos,δ,τ).^eos.b₄ .* δ * Ψ(eos,δ,τ))
end

function ϕr(eos::multiparameos, δ, τ)
    return res_1(eos, δ, τ)+res_2(eos, δ, τ)+res_3(eos, δ, τ)+res_4(eos, δ, τ)
end

function ϕ₀_δ(eos::multiparameos, δ, τ)
    f(δ_) = ϕ₀(eos, δ_, τ)
    return Calculus.derivative(f, δ)
end
function ϕ₀_δδ(eos::multiparameos, δ, τ)
    f(δ_) = ϕ₀_δ(eos, δ_, τ)
    return Calculus.derivative(f, δ)
end
function ϕ₀_τ(eos::multiparameos, δ, τ)
    f(τ_) = ϕ₀(eos, δ, τ_)
    return Calculus.derivative(f, τ)
end
function ϕ₀_ττ(eos::multiparameos, δ, τ)
    f(τ_) = ϕ₀_τ(eos, δ, τ_)
    return Calculus.derivative(f, τ)
end
function ϕ₀_δτ(eos::multiparameos, δ, τ)
    return 0.0
end

function ϕr_δ(eos::multiparameos, δ, τ)
    f(δ_) = ϕr(eos, δ_, τ)
    return Calculus.derivative(f, δ)
end
function ϕr_δδ(eos::multiparameos, δ, τ)
    f(δ_) = ϕr_δ(eos, δ_, τ)
    return Calculus.derivative(f, δ)
end
function ϕr_τ(eos::multiparameos, δ, τ)
    f(τ_) = ϕr(eos, δ, τ_)
    return Calculus.derivative(f, τ)
end
function ϕr_ττ(eos::multiparameos, δ, τ)
    f(τ_) = ϕr_τ(eos, δ, τ_)
    return Calculus.derivative(f, τ)
end
function ϕr_δτ(eos::multiparameos, δ, τ)
    f(τ_) = ϕr_δ(eos, δ, τ_)
    return Calculus.derivative(f, τ)
end

function p(eos::multiparameos, ρ, T)
    δ = ρ/eos.ρc
    τ = eos.Tc/T
    p_nondim = 1 + δ*ϕr_δ(eos, δ, τ)
    return p_nondim * ρ*eos.R*T
end

"""
This function is arguably the most critical piece of the whole equation of
state. It allows all other functions to be expressed in terms of pressure and
temperature rather than density and temperature.
"""
function ρ_from_p(eos::multiparameos, p_, T)
  pr_ = p_/eos.pc
  Tr = T/eos.Tc
  if abs(pr_-1) < 0.05 || abs(Tr-1) < 0.05
    # Interpolated guess value near critical point
    ρg = eos.ρc + (2*(p_-eos.pc)/eos.pc - (T-eos.Tc)/eos.Tc)
  elseif T > eos.Tc
    # Guess value for supercritical fluid
    ρg = eos.ρc
  elseif p_ > psat(eos, T)
    # Guess value for liquid region
    ρg = 1/pr.vl(eos.PR_,p_, T)
  else
    # Guess value for vapor region
    ρg = 1/pr.vv(eos.PR_,p_, T)
  end
  return ρ_(eos, p_, T, ρg, 2)
end

function ρ_(eos::multiparameos, p_, T, ρg, order_=0)
    function solveit(ρ_guess)
        return abs(p_ - p(eos, ρ_guess, T))
    end
    ans = 0
    try
      a = fzero(solveit, ρg, order=order_)
      ans = a
    catch ex
      lastval = float(ex.reason[18:end])
      #@printf "\nLast value: %f" lastval
      pcalc = p(eos, lastval,T)
      #@printf "\np calculated at %f with last value of ρ" pcalc
      #@printf "\nError: %f" abs(p_ - pcalc)
      ans = lastval
    end
    #print(ans)
    return ans
end

function psat(eos::multiparameos, T)
  if T > eos.Tc
    throw(DomainError("Can't have a vapor pressure above the critical point!"))
  end

  guess = eos.psat_guess(T)
  function solveit(pg)
    ρl_g = 1/pr.vl(eos.PR_,pg, T)
    ρv_g = 1/pr.vv(eos.PR_,pg, T)
    Gl = G_(eos, ρl_g, T)
    Gv = G_(eos, ρv_g, T)
    return Gl-Gv
  end

  return fzero(solveit, guess)
end

function u(eos::multiparameos, p_, T)
    ρ = ρ_from_p(eos, p_, T)
    δ = ρ/eos.ρc
    τ = eos.Tc/T
    u_nondim = τ * (ϕ₀_τ(eos, δ, τ) + ϕr_τ(eos, δ, τ))
    return u_nondim * T*eos.R
end

function s(eos::multiparameos, p_, T)
    ρ = ρ_from_p(eos, p_, T)
    δ = ρ/eos.ρc
    τ = eos.Tc/T
    s_nondim = τ * (ϕ₀_τ(eos, δ, τ) + ϕr_τ(eos, δ, τ)) - (ϕ₀(eos, δ, τ) + ϕr(eos, δ, τ))
    return s_nondim * eos.R
end

function H(eos::multiparameos, p_, T)
    ρ = ρ_from_p(eos, p_, T)
    δ = ρ/eos.ρc
    τ = eos.Tc/T
    H_nondim = 1 + τ * (ϕ₀_τ(eos, δ, τ) + ϕr_τ(eos, δ, τ)) + δ*ϕr_δ(eos, δ, τ)
    return H_nondim * T*eos.R
end

function G(eos::multiparameos, p_, T)
    ρ = ρ_from_p(eos, p_, T)
    δ = ρ/eos.ρc
    τ = eos.Tc/T
    G_nondim = 1 + (ϕ₀(eos, δ, τ) + ϕr(eos, δ, τ)) + δ*ϕr_δ(eos, δ, τ)
    return G_nondim * T*eos.R
end

function G_(eos::multiparameos, ρ, T)
    δ = ρ/eos.ρc
    τ = eos.Tc/T
    G_nondim = 1 + (ϕ₀(eos, δ, τ) + ϕr(eos, δ, τ)) + δ*ϕr_δ(eos, δ, τ)
    return G_nondim * T*eos.R
end

function Cp(eos::multiparameos, p_, T)
    ρ = ρ_from_p(eos, p_, T)
    δ = ρ/eos.ρc
    τ = eos.Tc/T
    part1 = -τ^2 * (ϕ₀_ττ(eos, δ, τ) + ϕr_ττ(eos, δ, τ))
    p2_top = (1 + δ*ϕr_δ(eos, δ, τ) - δ*τ*ϕr_δτ(eos, δ, τ))^2
    p2_bot = 1 + 2δ*ϕr_δ(eos, δ, τ) + δ^2*ϕr_δδ(eos, δ, τ)
    Cp_nondim = part1 + p2_top/p2_bot
    return Cp_nondim * eos.R
end

end # ends module meos

module IAPWS
using meos
using pr.pengRobinson
function gen_iapws()
    h2o_n₀ = [-8.3204464837497, 6.6832105275932, 3.00632, 0.012436, 0.97315, 1.27950,
              0.96956, 0.24873]
    h2o_γ₀ = [1.28728967, 3.53734222, 7.74073708, 9.24437796, 27.5075105]

    h2o_n₁ = [0.12533547935523e-1,
              0.78957634722828e1,
              -0.87803203303561e1,
              0.31802509345418,
              -0.26145533859358,
              -0.78199751687981e-2,
              0.88089493102134e-2]
    h2o_d₁ = [1, 1, 1, 2, 2, 3, 4]
    h2o_t₁ = [-0.5, 0.875, 1, 0.5, 0.75, 0.375, 1]

    h2o_n₂ = [-0.66856572307965,
              0.20433810950965,
              -0.66212605039687e-4,
              -0.19232721156002,
              -0.25709043003438,
              0.16074868486251,
              -0.40092828925807e-1,
              0.39343422603254e-6,
              -0.75941377088144e-5,
              0.56250979351888e-3,
              -0.15608652257135e-4,
              0.11537996422951e-8,
              0.36582165144204e-6,
              -0.13251180074668e-11,
              -0.62639586912454e-9,
              -0.10793600908932,
              0.17611491008752e-1,
              0.22132295167546,
              -0.40247669763528,
              0.58083399985759,
              0.49969146990806e-2,
              -0.31358700712549e-1,
              -0.74315929710341,
              0.47807329915480,
              0.20527940895948e-1,
              -0.13636435110343,
              0.14180634400617e-1,
              0.83326504880713e-2,
              -0.29052336009585e-1,
              0.38615085574206e-1,
              -0.20393486513704e-1,
              -0.16554050063734e-2,
              0.19955571979541e-2,
              0.15870308324157e-3,
              -0.16388568342530e-4,
              0.43613615723811e-1,
              0.34994005463765e-1,
              -0.76788197844621e-1,
              0.22446277332006e-1,
              -0.62689710414685e-4,
              -0.55711118565645e-9,
              -0.19905718354408,
              0.31777497330738,
              -0.11841182425981]
    h2o_c₂ = [1,1,1,1,1, 1,1,1,1,1, 1,1,1,1,1, 2,2,2,2,2,2,2,2, 2,2,2,2,2,2,2,2,2,
              2,2,2, 3,3,3,3, 4,6,6,6,6]
    h2o_d₂ = [1,1,1,2,2,3,4,4,5,7,9,10,11,13,15,1,2,2,2,3,4,4,4,5,6,6,7,9,9,9,9,9,
              10,10,12,3,4,4,5,14,3,6,6,6]
    h2o_t₂ = [4,6,12,1,5,4,2,13,9,3,4,11,4,13,1,7,1,9,10,10,3,7,10,10,6,10,10,1,2,3,
              4,8,6,9,8,16,22,23,23,10,50,44,46,50]

    h2o_d₃ = [3, 3, 3]
    h2o_t₃ = [0,1,4]
    h2o_n₃ = [-0.31306260323435e2, 0.31546140237781e2, -0.25213154341695e4]
    h2o_α₃ = [20,20,20]
    h2o_β₃ = [150,150,150]
    h2o_γ₃ = [1.21,1.21,1.25]
    h2o_ϵ₃ = [1,1,1]

    h2o_a₄ = [3.5,3.5]
    h2o_b₄ = [0.85,0.95]
    h2o_β₄ = [0.3,0.3]
    h2o_n₄ = [-0.14874640856724, 0.31806110878444]
    h2o_A₄ = [0.32,0.32]
    h2o_B₄ = [0.2,0.2]
    h2o_C₄ = [28,32]
    h2o_D₄ = [700,800]

    Tc = 647.096
    pc = 22.064e6
    ρc = 322
    ω = 0.344
    R = 0.46151805*1000
    function psat_dippr(T)
      return exp(73.649 - 7258.2/T - 7.3037*log(T) + 4.1653e-6T^2)
    end
    pr_thing = pengRobinson(pc, Tc, ω)

    return meos.multiparameos(Tc,ρc,pc,R,ω,h2o_n₀,h2o_γ₀,h2o_n₁,h2o_d₁,h2o_t₁,
      h2o_n₂,h2o_d₂,h2o_t₂,h2o_c₂,h2o_n₃,h2o_d₃,h2o_t₃,h2o_α₃,h2o_β₃,h2o_γ₃,
      h2o_ϵ₃,h2o_n₄,h2o_a₄,h2o_b₄,h2o_β₄,h2o_A₄,h2o_B₄,h2o_C₄,h2o_D₄, psat_dippr, pr_thing)
end

end # ends module IAPWS

module co2
import meos
using pr.pengRobinson
using subst.substance
function gen_co2()
    co2_n₀= [ 8.37304456, -3.70454304,  2.5       ,  1.99427042,  0.62105248,
        0.41195293,  1.04028922,  0.08327678]

    co2_γ₀ = [3.15163,   6.1119 ,   6.77708,
            11.32384,  27.08792]

    co2_n₁ = [ 0.38856823,  2.93854759, -5.58671885, -0.767532  ,  0.31729004,
            0.54803316,  0.12279411]
    co2_d₁ = [ 1.,  1.,  1.,  1.,  2.,  2.,  3.]
    co2_t₁ = [ 0.  ,  0.75,  1.  ,  2.  ,  0.75,  2.  ,  0.75]
    co2_n₂ = [  2.16589615e+00,   1.58417351e+00,  -2.31327054e-01,
             5.81169164e-02,  -5.53691372e-01,   4.89466159e-01,
            -2.42757398e-02,   6.24947905e-02,  -1.21758602e-01,
            -3.70556853e-01,  -1.67758797e-02,  -1.19607366e-01,
            -4.56193625e-02,   3.56127893e-02,  -7.44277271e-03,
            -1.73957049e-03,  -2.18101213e-02,   2.43321666e-02,
            -3.47401334e-02,   1.43387158e-01,  -1.34919691e-01,
            -2.31512251e-02,   1.23631255e-02,   2.10583220e-03,
            -3.39585190e-04,   5.59936518e-03,  -3.03351181e-04]

    co2_d₂ = [  1.,   2.,   4.,   5.,   5.,   5.,   6.,   6.,   6.,   1.,   1.,
             4.,   4.,   4.,   7.,   8.,   2.,   3.,   3.,   5.,   5.,   6.,
             7.,   8.,  10.,   4.,   8.]

    co2_c₂ = [ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  2.,  2.,  2.,  2.,
            2.,  2.,  2.,  3.,  3.,  3.,  4.,  4.,  4.,  4.,  4.,  4.,  5.,  6.]

    co2_t₂ = [  1.5,   1.5,   2.5,   0. ,   1.5,   2. ,   0. ,   1. ,   2. ,
             3. ,   6. ,   3. ,   6. ,   8. ,   6. ,   0. ,   7. ,  12. ,
            16. ,  22. ,  24. ,  16. ,  24. ,   8. ,   2. ,  28. ,  14. ]

    co2_n₃ = [  -213.65488688,  26641.56914927, -24027.21220456,  -283.41603424,
              212.472844  ]
    co2_n₄ = [-0.66642277,  0.72608632,  0.05506867]
    co2_d₃ = [2.,2.,2.,3.,3.]
    co2_t₃ = [1.,0.,1.,3.,3.]
    co2_α₃ = [25.,25.,25.,15.,20.]
    co2_β₃ = [325.,300.,300.,275.,275.]
    co2_γ₃ = [1.16,1.19,1.19,1.25,1.22]
    co2_ϵ₃ = [1.,1.,1.,1.,1.]

    co2_a₄ = [3.5,3.5,3.0]
    co2_b₄ = [0.875,0.925,0.875]
    co2_β₄ = [0.3,0.3,0.3]
    co2_A₄ = [0.7,0.7,0.7]
    co2_B₄ = [0.3,0.3,1.0]
    co2_C₄ = [10.0,10.0,12.5]
    co2_D₄ = [275.,275.,275.]

    Tc = 304.1282
    ρc = 467.6
    pc = 7.38e6
    R = 188.9241
    ω = 0.228
    function psat_dippr(T)
      return exp(47.0169 - 2839/T - 3.86388*log(T) + 2.81115e-16T^6)
    end
    pr_thing = pengRobinson(pc, Tc, ω)
    return meos.multiparameos(Tc,ρc,pc,ω,R,co2_n₀,co2_γ₀,co2_n₁,co2_d₁,co2_t₁,
      co2_n₂,co2_d₂,co2_t₂,co2_c₂,co2_n₃,co2_d₃,co2_t₃,co2_α₃,co2_β₃,co2_γ₃,
      co2_ϵ₃,co2_n₄,co2_a₄,co2_b₄,co2_β₄,co2_A₄,co2_B₄,co2_C₄,co2_D₄,psat_dippr,pr_thing)
end

function co2_viscosity(eos::meos.multiparameos, p, T)

  i1 = [0,1,2,3,4]
  a = [0.235156,
       -0.491266,
       5.211155e-2,
       5.347906e-2,
       -1.537102e-2]
  i2 = [1,2,6,8,8]
  j = [0,0,3,0,1]
  d = [0.407119e-2,
       0.7198037e-4,
       0.211697e-16,
       0.2971072e-22,
       -0.1627888e-22]
  ρ = meos.ρ_from_p(eos, p, T)
  Tr = T / 251.196
  rhor = ρ/eos.ρc
  top = 1.00697 * T ^ 0.5
  C_star = 0
  for i in i1
    C_star += a[i+1] * log(Tr) ^ i
  end
  bottom = exp(C_star)
  eta0 = top / bottom
  etar = 0
  for (num, i) in enumerate(i2)
    etar += d[num] * (rhor ^ i) / (Tr ^ j[num])
  end
  eta = eta0 + etar
  return eta * 1e-6
end

function co2_thermalConductivity(eos::meos.multiparameos, p, T)
  ik1 = [1,2,3,4,5,6,7,8,9,10]
  g = [0.0,0.0,1.5,0.0,1.0,1.5,1.5,1.5,3.5,5.5]
  h = [1.0,5.0,1.0,1.0,2.0,0.0,5.0,9.0,0.0,0.0]
  n = [7.69857587,
       0.159885811,
       1.56918621,
       -6.73400790,
       16.3890156,
       3.69415242,
       22.3205514,
       66.1420950,
       -0.171779133,
       0.00433043347,
       0.775547504]
  a = [3.0,
       6.70697,
       0.94604,
       0.30,
       0.30,
       0.39751,
       0.33791,
       0.77963,
       0.79857,
       0.9,
       0.02,
       0.2]
  ρ = meos.ρ_from_p(eos, p, T)
  rhor = ρ / eos.ρc
  Tr = T / eos.Tc
  thermal_cond = 0
  for i in ik1
    if i != 0
      index = i
      term = n[index] * Tr ^ g[index] * rhor ^ h[index]
      if i >= 4
        term *= exp(-5 * rhor ^ 2)
      end
      thermal_cond += term
    end
  end
  alpha = 1 - a[10] * acosh(1 + a[11] * ((1 - Tr) ^ 2) ^ a[12])
  crit_correction = rhor * exp((-rhor ^ a[1]) / a[1] - (a[2] * (Tr - 1)) ^ 2 - (a[3] * (rhor - 1)) ^ 2) /
                      ((((1 - 1 / Tr) + a[4] * ((rhor - 1) ^ 2) ^ (1 / (2 * a[5]))) ^ 2) ^ a[6]
                       + ((a[7] * (rhor - alpha)) ^ 2) ^ a[8]) ^ a[9]
  nc = n[9]
  thermal_cond += nc * crit_correction
  kc = 4.81314e-3
  return thermal_cond * kc  # Should return units of W / m *K
end

function init_co2()
  general_eos = gen_co2()
  function co2_density(p, T)
    return meos.ρ_from_p(general_eos, p, T)
  end
  function co2_enthalpy(p, T)
    return meos.H(general_eos, p, T)
  end
  function co2_Cp(p, T)
    return meos.Cp(general_eos, p, T)
  end
  function co2_μ(p, T)
    return co2_viscosity(general_eos, p, T)
  end
  function co2_k(p, T)
    return co2_thermalConductivity(general_eos, p, T)
  end
  return substance("carbon dioxide", [co2_density, co2_enthalpy, co2_Cp,
                                      co2_μ, co2_k])
end

end #ends module co2
