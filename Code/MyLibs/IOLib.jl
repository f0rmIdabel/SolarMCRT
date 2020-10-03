using Printf
using HDF5
using Unitful
import PhysicalConstants.CODATA2018: σ_e
# @derived_dimension PerLength Unitful.𝐋^-1 # can't make this work


function get_atmosphere_data(path_atmos, path_RH)

    # Read atmosphere file
    #------------------------------------------------------------
    x = h5read(path_atmos, "x")u"m"
    y = h5read(path_atmos, "y")u"m" # !Beware of negative values!
    z = h5read(path_atmos, "z")[:,1]u"m"

    velocity_x = h5read(path_atmos, "velocity_x")[:,:,:,1]u"m/s"
    velocity_y = h5read(path_atmos, "velocity_y")[:,:,:,1]u"m/s"
    velocity_z = h5read(path_atmos, "velocity_z")[:,:,:,1]u"m/s"

    temperature = h5read(path_atmos, "temperature")[:,:,:,1]u"K"
    electron_density = h5read(path_atmos, "electron_density")[:,:,:,1]u"m^-3"

    # Read RH output
    #---------------------------------------------------------------
    χ_absorption = h5read(path_RH, "chi_continuum")[1,:,:,:]u"m^-1"

    # Re-work parameters
    #---------------------------------------------------------------

    # Add endpoints
    x = push!(x, 2*x[end] - x[end-1])
    y = push!(y, 2*y[end] - y[end-1])
    z = push!(z, 2*z[end] - z[end-1])

    # Calculate epsilon and chi
    χ_thomson = σ_e*electron_density
    χ_continuum = χ_absorption .+ χ_thomson
    ε_continuum = χ_absorption ./ χ_continuum

    return x, y, z, velocity_x, velocity_y, velocity_z,
           temperature, χ_continuum, ε_continuum
end


function write_results_to_file(current_time, τ_max, packets, destroyed, escaped, scatterings,
                               elapsed_time, meanJ, minJ, maxJ, filename="Results.txt")

    f = open("/mn/stornext/u3/idarhan/SolarMCRT/Results/"*filename, "a")
    results = string(current_time)*(@sprintf("%13.1f%13.1e%20.5e%18g%22g%22.1f%22.1f%17g%17g\n",
                                             τ_max, packets, destroyed, escaped, scatterings,
                                             elapsed_time, meanJ, minJ, maxJ))
    write(f, results)
end
