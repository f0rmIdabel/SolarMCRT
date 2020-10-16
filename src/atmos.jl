using HDF5
using Unitful
import PhysicalConstants.CODATA2018: σ_e

struct Atmosphere
   # Dimensions
   # (nx+1), (ny+1), (nz+1)
   x::Array{<:Unitful.Length, 1}
   y::Array{<:Unitful.Length, 1}
   z::Array{<:Unitful.Length, 1}

   # Local properties
   # (nx, ny, nz)
   velocity_x::Array{<:Unitful.Velocity, 3}
   velocity_y::Array{<:Unitful.Velocity, 3}
   velocity_z::Array{<:Unitful.Velocity, 3}
   temperature::Array{<:Unitful.Temperature, 3}

   # (λ, nx, ny, nz)
   χ_continuum::Array{<:Unitful.Quantity{<:Real, Unitful.𝐋^(-1)}, 3} ###################
   ε_continuum::Array{Real, 3} #########################################################
end

"""
   function get_atmosphere_data(atmos_data,
                                rh_output)

Reads atmosphere parameters and reworks them to fit simulation.
"""
function get_atmosphere_data(atmos_data,
                             rh_output)

    path_atmos = "../../../basement/MScProject/Atmospheres/"*atmos_data
    path_rh = "../../../basement/MScProject/Atmospheres/"*rh_output
    # ===========================================================
    # READ ATMOSPHERE FILE
    # ===========================================================
    atmos = h5open(path_atmos, "r")

    x = read(atmos, "x")u"m"
    y = read(atmos, "y")u"m"
    z = read(atmos, "z")[:,1]u"m"

    velocity_x = read(atmos, "velocity_x")[:,:,:,1]u"m/s"
    velocity_y = read(atmos, "velocity_y")[:,:,:,1]u"m/s"
    velocity_z = read(atmos, "velocity_z")[:,:,:,1]u"m/s"

    temperature = read(atmos, "temperature")[:,:,:,1]u"K"
    electron_density = read(atmos, "electron_density")[:,:,:,1]u"m^-3"

    close(atmos)
    # ===========================================================
    # READ RH OUTPUT
    # ===========================================================
    rh = h5open(path_rh, "r")

    χ_absorption = read(rh, "chi_continuum")[1,:,:,:]u"m^-1"  ############################

    close(rh)
    # ===========================================================
    # RE-WORK PARAMETERS
    # ===========================================================

    # Calculate epsilon and chi
    χ_thomson = σ_e*electron_density                    ##########################################
    χ_continuum = χ_absorption .+ χ_thomson
    ε_continuum = χ_absorption ./ χ_continuum

    # Transpose all 3D-space arrays,(k,i,j) -> (i,j,k)
    velocity_x = permutedims(velocity_x, [2,3,1])
    velocity_y = permutedims(velocity_y, [2,3,1])
    velocity_z = permutedims(velocity_z, [2,3,1])
    temperature = permutedims(temperature, [2,3,1])
    χ_continuum = permutedims(χ_continuum, [2,3,1]) #########################################
    ε_continuum = permutedims(ε_continuum, [2,3,1]) #########################################

    # Make sure x and y are increasing and z decreasing
    if x[1] > x[end]
        x = reverse(x)
        velocity_x = velocity_x[end:-1:1,:,:]
        velocity_y = velocity_y[end:-1:1,:,:]
        velocity_z = velocity_z[end:-1:1,:,:]
        temperature = temperature[end:-1:1,:,:]
        χ_continuum = χ_continuum[end:-1:1,:,:]
        ε_continuum = ε_continuum[end:-1:1,:,:]
    end

    if y[1] > y[end]
        y = reverse(y)
        velocity_x = velocity_x[:,end:-1:1,:]
        velocity_y = velocity_y[:,end:-1:1,:]
        velocity_z = velocity_z[:,end:-1:1,:]
        temperature = temperature[:,end:-1:1,:]
        χ_continuum = χ_continuum[:,end:-1:1,:]
        ε_continuum = ε_continuum[:,end:-1:1,:]
    end

    if z[1] < z[end]
        z = reverse(z)
        velocity_x = velocity_x[:,:,end:-1:1]
        velocity_y = velocity_y[:,:,end:-1:1]
        velocity_z = velocity_z[:,:,end:-1:1]
        temperature = temperature[:,:,end:-1:1]
        χ_continuum = χ_continuum[:,:,end:-1:1]
        ε_continuum = ε_continuum[:,:,end:-1:1]
    end

    nx, ny, nz = size(χ_continuum)
    # Add endpoints for box calculations
    x = push!(x, 2*x[end] - x[end-1])
    y = push!(y, 2*y[end] - y[end-1])
    z = push!(z, 2*z[end] - z[end-1])

    return x, y, z, velocity_x, velocity_y, velocity_z,
           temperature, χ_continuum, ε_continuum
end
