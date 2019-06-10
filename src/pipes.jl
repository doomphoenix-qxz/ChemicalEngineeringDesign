#jl
struct pipe{T<:Real}
    length::T
    diameter::T
    radius::T
    flowarea::T
    thickness::T
    outer_diameter::T
    outer_radius::T
    totalarea::T
    roughness::T
    elevation_change::T
    therm_cond::Function
end

const valid_pipe_diams = [0.125, 0.25, 0.375, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5,
                    3.0, 3.5, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 16.0]

function pipe(l, d, t, ϵ, z, kfunc)
    # Assumes a circular pipe of diameter d
    return pipe(l,d,d/2,pi*(d/2)^2,t,d+2t,(d+2t)/2,pi*((d+2t)/2)^2,ϵ,z,kfunc)
end

mutable struct pipeProps{T<:Real}
    dᵢ::T
    dₒ::T
    thickness::T
    rᵢ::T
    rₒ::T
    units::String
end

function pipeProps(d, t, units="in")
  di = d - 2t
  ro = d/2
  ri = di/2
  return pipeProps(di, d, t, ri, ro, units)
end

const conversions = Dict(["cm","m"] => 1e-2,
                   ["m","cm"] => 100,
                   ["in","cm"] => 2.54,
                   ["cm","in"] => 1/2.54,
                   ["in","m"] => 2.54e-2,
                   ["m","in"] => 1/2.54e-2,
                   ["m","m"] => 1.0,
                   ["cm","cm"] => 1.0,
                   ["in","in"] => 1.0)

function convert(p::pipeProps, newUnits::String)
    key_ = [p.units, newUnits]
    conversion_factor = conversions[key_]
    p.dᵢ *= conversion_factor
    p.dₒ *= conversion_factor
    p.rᵢ *= conversion_factor
    p.rₒ *= conversion_factor
    p.units = newUnits
end

const schedule40_props = [pipeProps(0.405, 0.068), pipeProps(0.540, 0.088),
                    pipeProps(0.675, 0.091), pipeProps(0.840, 0.109),
                    pipeProps(1.050, 0.113), pipeProps(1.315, 0.133),
                    pipeProps(1.660, 0.140), pipeProps(1.900, 0.145),
                    pipeProps(2.375, 0.154), pipeProps(2.875, 0.203),
                    pipeProps(3.500, 0.216), pipeProps(4.000, 0.226),
                    pipeProps(4.500, 0.237), pipeProps(5.563, 0.258),
                    pipeProps(6.625, 0.280), pipeProps(8.625, 0.322),
                    pipeProps(10.750, 0.365), pipeProps(12.750, 0.406),
                    pipeProps(16.000, 0.500)]
# Note: These outer diameters and wall thicknesses were obtained at this web
# site: http://www.engineeringtoolbox.com/steel-dimensions-d_43.html

for props in schedule40_props
    convert(props, "m")
end

const props40_by_nominal_diameter = Dict(zip(valid_pipe_diams, schedule40_props))

const schedule_props = Dict("40" => props40_by_nominal_diameter)


function pipe(schedule_info::Tuple, mat::subst.material, l, Δz)
  # schedule_info is a Tuple: the schedule number comes first and the nominal
  # diameter comes second
  props = schedule_props[schedule_info[1]][schedule_info[2]]
  return pipe(l,props.dᵢ,props.dᵢ/2,pi*(props.rᵢ)^2,props.thickness,props.dₒ,
            props.rₒ,pi*(props.rₒ)^2,mat.roughness,Δz,mat.thermal_conductivity)
end

function pipe(schedule_info::Tuple, mat::String, l, Δz)
  return pipe(schedule_info, pipe_materials[mat], l, Δz)
end

function colebrook(Re, p::pipe, f_g)
    if f_g < 0
        return f_g * 1e5
    else
        return 10^(-1 / (2 * √(f_g))) - (p.roughness / (3.7 * p.diameter) + 2.51 / (Re * √(f_g)))
    end
end

function getf(Re, p::pipe)
    function dummysolve!(f_guess, fvec)
        fvec[1] = colebrook(Re, p, f_guess[1])
    end
    return NLsolve.nlsolve(dummysolve!, [0.001]).zero[1]

end

"""
Struct material
Used for solid materials in pipes. So far includes just name, thermal
conductivity, and surface roughness.
"""
struct material
    name::String
    thermal_conductivity::Function
    roughness::Real
end

function stainless_steel_thermal_conductivity(T)
    # Source: http://www.mace.manchester.ac.uk/project/research/structures/strucfire/materialInFire/Steel/StainlessSteel/thermalProperties.htm
    return 14.6 + 1.27e-2 *(T-273.15)
end

function copper_thermal_conductivity(T)
    # Source: http://www-ferp.ucsd.edu/LIB/PROPS/PANOS/cu.html
    return 14.6 + 1.27e-2 *(T-273.15)
end

const ss_roughness_by_finish = Dict("2D" => 1E-6,
                              "2B" => 0.5E-6,
                              "2R" => 0.2E-6,
                              "BA" => 0.2E-6,
                              "2BB" => 0.1E-6)
# These surface roughness values are in m and were obtained from the following
# corporate site: http://www.outokumpu.com/en/products-properties/
# more-stainless/stainless-steel-surface-finishes/cold-rolled-finishes/Pages/default.aspx
# Note that the given values are the maximum roughness values in the range on
# the website, except for finish 2BB, which didn't have a roughness value
# so I took the lowest value for finish 2B since the compnay says this:
# "Due to the fact that surface roughness of the finish 2BB is lower than that
# of 2B, some modifications on lubrication during forming might be needed."

const ss_mat = material("Stainless Steel", stainless_steel_thermal_conductivity,
                  ss_roughness_by_finish["2B"])

const cu_roughness = 0.03e-3
# Copper roughness.
# Source: http://www.pressure-drop.com/Online-Calculator/rauh.html

const cu_mat = material("Copper", copper_thermal_conductivity, cu_roughness)
