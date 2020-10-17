using HDF5
using Unitful
using UnitfulRecipes
using Transparency
import PhysicalConstants.CODATA2018: c_0
@derived_dimension NumberDensity Unitful.𝐋^-3
@derived_dimension PerLength Unitful.𝐋^-1

struct Atmosphere
   # Dimensions - position of box edges
   # (nx+1), (ny+1), (nz+1)
   x::Array{<:Unitful.Length, 1} # Increasing
   y::Array{<:Unitful.Length, 1} # Increasing
   z::Array{<:Unitful.Length, 1} # Decreasing

   # Local box properties
   # (nx, ny, nz)
   velocity_x::Array{<:Unitful.Velocity, 3}
   velocity_y::Array{<:Unitful.Velocity, 3}
   velocity_z::Array{<:Unitful.Velocity, 3}
   temperature::Array{<:Unitful.Temperature, 3}

   # (λ, nx, ny, nz)
   χ::Array{<:PerLength, 3}# 4 ##########################################
   ε::Array{Real, 3} # 4 ########################################################
end

"""
From Tiago
"""
function α_abs(λ::Unitful.Length,
               temperature::Unitful.Temperature,
               electron_density::NumberDensity,
               h_ground_density::NumberDensity,
               proton_density::NumberDensity)

    α = Transparency.hminus_ff_stilley(λ, temperature, h_ground_density, electron_density)
    α += Transparency.hminus_bf_geltman(λ, temperature, h_ground_density, electron_density)
    α += hydrogenic_ff(c_0 / λ, temperature, electron_density, proton_density, 1)
    α += h2plus_ff(λ, temperature, h_ground_density, proton_density)
    α += h2plus_bf(λ, temperature, h_ground_density, proton_density)
    return α
end

"""
From Tiago
"""
function α_scatt(λ::Unitful.Length,
                 electron_density::NumberDensity,
                 h_ground_density::NumberDensity)

    α = thomson(electron_density)
    α += rayleigh_h(λ, h_ground_density)
    return α
end

"""
   function get_atmosphere_data(atmos_data,
                                rh_output)

Reads atmosphere parameters and reworks them to fit simulation.
"""
function get_atmosphere_data(atmos_data, λ)

    path_atmos = "../../../basement/MScProject/Atmospheres/"*atmos_data
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
    hydrogen_populations = read(atmos, "hydrogen_populations")[:,:,:,:,1]u"m^-3"
    close(atmos)
    # ===========================================================
    # RE-WORK PARAMETERS
    # ===========================================================

    # Calculate epsilon and chi
    hydrogen_density = sum(hydrogen_populations, dims = 4)
    ionised_hydrogen_density = hydrogen_populations[:,:,:,end]
    neutral_hydrogen_density = hydrogen_density .- ionised_hydrogen_density

    χ_abs = α_abs.(λ, temperature, electron_density, neutral_hydrogen_density, ionised_hydrogen_density)[:,:,:,1]
    χ_scatt = α_scatt.(λ, electron_density, neutral_hydrogen_density)[:,:,:,1]
    χ = χ_abs .+ χ_scatt
    ε = χ_abs ./ χ

    # Transpose all 3D-space arrays,(k,i,j) -> (i,j,k)
    velocity_x = permutedims(velocity_x, [2,3,1])
    velocity_y = permutedims(velocity_y, [2,3,1])
    velocity_z = permutedims(velocity_z, [2,3,1])
    temperature = permutedims(temperature, [2,3,1])
    χ = permutedims(χ, [2,3,1]) #########################################
    ε = permutedims(ε, [2,3,1]) #########################################

    # Make sure x and y are increasing and z decreasing

    if x[1] > x[end]
        x = reverse(x)
        velocity_x = velocity_x[end:-1:1,:,:]
        velocity_y = velocity_y[end:-1:1,:,:]
        velocity_z = velocity_z[end:-1:1,:,:]
        temperature = temperature[end:-1:1,:,:]
        χ = χ[end:-1:1,:,:]
        ε = ε[end:-1:1,:,:]
    end

    if y[1] > y[end]
        y = reverse(y)
        velocity_x = velocity_x[:,end:-1:1,:]
        velocity_y = velocity_y[:,end:-1:1,:]
        velocity_z = velocity_z[:,end:-1:1,:]
        temperature = temperature[:,end:-1:1,:]
        χ = χ[:,end:-1:1,:]
        ε = ε[:,end:-1:1,:]
    end

    if z[1] < z[end]
        z = reverse(z)
        velocity_x = velocity_x[:,:,end:-1:1]
        velocity_y = velocity_y[:,:,end:-1:1]
        velocity_z = velocity_z[:,:,end:-1:1]
        temperature = temperature[:,:,end:-1:1]
        χ = χ[:,:,end:-1:1]
        ε = ε[:,:,end:-1:1]
    end

    # Add endpoints for box calculations
    x = push!(x, 2*x[end] - x[end-1])
    y = push!(y, 2*y[end] - y[end-1])
    z = push!(z, 2*z[end] - z[end-1])

    return x, y, z, velocity_x, velocity_y, velocity_z,
           temperature, χ, ε
end
