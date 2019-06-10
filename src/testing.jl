#test multiparameter_eos.jl
include("parameterized_eos.jl")


function test_meos(eos::multiparameos, δ, τ)
    print("Testing module meos:\n")
    print("δ = "*string(δ))
    print("\nτ = "*string(τ))
    print("\nϕ₀ = "*string(ϕ₀(eos, δ, τ)))
    print("\nϕ₀_δ = "*string(ϕ₀_δ(eos, δ, τ)))
    print("\nϕ₀_δδ = "*string(ϕ₀_δδ(eos, δ, τ)))
    print("\nϕ₀_τ = "*string(ϕ₀_τ(eos, δ, τ)))
    print("\nϕ₀_ττ = "*string(ϕ₀_ττ(eos, δ, τ)))
    print("\nϕ₀_δτ = "*string(ϕ₀_δτ(eos, δ, τ)))
    print("\nϕr = "*string(ϕr(eos, δ, τ)))
    print("\nϕr_δ = "*string(ϕr_δ(eos, δ, τ)))
    print("\nϕr_δδ = "*string(ϕr_δδ(eos, δ, τ)))
    print("\nϕr_τ = "*string(ϕr_τ(eos, δ, τ)))
    print("\nϕr_ττ = "*string(ϕr_ττ(eos, δ, τ)))
    print("\nϕr_δτ = "*string(ϕr_δτ(eos, δ, τ)))
    print("\n\nNew test conditions:")
    ρ = 467.5483
    T = 305
    @printf "\nT = %f K, ρ = %f kg/m³" T ρ

    p_ = 7.5e6
    print("\np = "*string(p_))
    p_ = p(eos, ρ, T)
    print("\np = "*string(p_))
    print("\nSolving in reverse, with p = "string(p_)*" and T="*string(T))
    #print(ρ_from_p(eos, p_, 300))
    @printf "\n%f kg/m³" ρ_from_p(eos, p_, T)
end

#eos = gen_iapws()
#Tc = 647.096
#ρc = 322
#δ₁ = 838.025/ρc
#τ₁ = Tc/500
#meos.test_peos(eos, δ₁, τ₁)
#meos.ϕ₀(eos, δ₁, τ₁)

function test_co2_transport()
  @printf "Testing thermal conductivity equation\n"
  pspace = [8, 9, 10, 15, 20]
  t1 = 600
  pspace2 = [0.1, 0.1, 0.1, 7, 15, 50, 75]
  t2space = [220, 300, 800, 304, 220, 300, 800]

  for p in pspace
    @printf "Temperature is 600 K. Pressure is %f MPa.\n" p
    rho = meos.ρ_from_p(eos, p * 10 ^ 6, t1)
    @printf "Density is %f kg/m^3 " rho
    @printf "Thermal cond. is %f W / m * K" co2_thermalConductivity(eos, p * 10 ^ 6, t1)
  end

  @printf "Testing viscosity equation"
  for (i, p) in enumerate(pspace2)
    t2 = t2space[i]
    @printf "Temperature is %f K. Pressure is %f MPa." t2 p
    rho = meos.ρ_from_p(eos, p*1e6, t2)
    @printf "Density is %f kg/m^3" rho
    @printf "Viscosity is %f micropoise" co2_viscosity(eos, p*1e6, t2)
  end
end
#test_co2_transport()
#Tc = 304.1282
#ρc = 467.6
#δ₁ = 838.025/ρc
#τ₁ = Tc/500
#meos.test_peos(eos, δ₁, τ₁)

#Test pipes.jl
import pipes

function main_()
    pipe1 = pipes.pipe(1, 20e-3, 4e-3, 0.002e-3, 0.2)
    Re1 = 1e7
    ans = pipes.getf(Re1, pipe1)
    print("Friction factor is "*string(ans))
end

main_()
