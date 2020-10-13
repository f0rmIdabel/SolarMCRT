using HDF5
using Unitful
import PhysicalConstants.CODATA2018: σ_e

struct Atmosphere
   # Dimensions
   x::Array{<:Unitful.Length, 1}
   y::Array{<:Unitful.Length, 1}
   z::Array{<:Unitful.Length, 1}

   # Local properties
   velocity_x::Array{<:Unitful.Velocity, 3}
   velocity_y::Array{<:Unitful.Velocity, 3}
   velocity_z::Array{<:Unitful.Velocity, 3}
   temperature::Array{<:Unitful.Temperature, 3}
   χ_continuum::Array{<:Unitful.Quantity{<:Real, Unitful.𝐋^(-1)}, 3}
   ε_continuum::Array{Real, 3}
end

"""
   function get_atmosphere_data(atmos_data,
                                rh_output)

Reads atmosphere parameters and reworks them to fit simulation.
"""
function get_atmosphere_data(atmos_data,
                             rh_output)

    path_atmos = "../../basement/MScProject/Atmospheres/"*atmos_data
    path_RH = "../../basement/MScProject/Atmospheres/"*rh_output

    # ===========================================================
    # READ ATMOSPHERE FILE
    # ===========================================================
    x = h5read(path_atmos, "x")[450:end]u"m"
    y = -h5read(path_atmos, "y")[450:end]u"m" # get pos values
    z = h5read(path_atmos, "z")[:,1]u"m"

    velocity_x = h5read(path_atmos, "velocity_x")[:,:,:,1]u"m/s"
    velocity_y = h5read(path_atmos, "velocity_y")[:,:,:,1]u"m/s"
    velocity_z = h5read(path_atmos, "velocity_z")[:,:,:,1]u"m/s"

    temperature = h5read(path_atmos, "temperature")[:,450:end,450:end,1]u"K"
    electron_density = h5read(path_atmos, "electron_density")[:,450:end,450:end,1]u"m^-3"

    # ===========================================================
    # READ RH OUTPUT
    # ===========================================================
    χ_absorption = h5read(path_RH, "chi_continuum")[1,:,450:end,450:end]u"m^-1"  #WL
    # ===========================================================
    # RE-WORK PARAMETERS
    # ===========================================================
    # Add endpoints for box calculations
    x = push!(x, 2*x[end] - x[end-1])
    y = push!(y, 2*y[end] - y[end-1])
    z = push!(z, 2*z[end] - z[end-1])

    # Calculate epsilon and chi                                       #WL
    χ_thomson = σ_e*electron_density
    χ_continuum = χ_absorption .+ χ_thomson
    ε_continuum = χ_absorption ./ χ_continuum

    # Transpose all 3D-space arrays,(k,i,j) -> (i,j,k)
    velocity_x = permutedims(velocity_x, [2,3,1])
    velocity_y = permutedims(velocity_y, [2,3,1])
    velocity_z = permutedims(velocity_z, [2,3,1])
    temperature = permutedims(temperature, [2,3,1])
    χ_continuum = permutedims(χ_continuum, [2,3,1])
    ε_continuum = permutedims(ε_continuum, [2,3,1])

    return x, y, z, velocity_x, velocity_y, velocity_z,
           temperature, χ_continuum, ε_continuum
end

#
