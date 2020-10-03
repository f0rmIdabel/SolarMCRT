include("MCRT.jl")
include("atmos.jl")
include("MyLibs/IOLib.jl")
include("MyLibs/plotLib.jl")
include("MyLibs/physLib.jl")

using Dates
using Printf

function main(max_scatterings = 1e10)

    # Initialise atmosphere
    parameters = get_atmosphere_data("/mn/stornext/u3/idarhan/basement/MScProject/Atmospheres/bifrost_cb24bih_s385_fullv.ncdf",
                                     "/mn/stornext/u3/idarhan/basement/MScProject/Atmospheres/output_ray.hdf5")
    atmosphere = Atmosphere(parameters...)

    # Initialise wavelengths
    wavelength = 500u"nm"

    # Get user input: τ_max, packets
    println("Choose maximum optical depth:")
    τ_max = readline()
    τ_max = parse(Float64, τ_max)

    println("Choose number of photon packets:")
    packets = readline()
    packets = parse(Float64, packets)

    # Current time
    current_time = now()

    # Run and time simulation
    simulation = @timed simulate(atmosphere, wavelength,
                                 max_scatterings, τ_max, packets)

    surface, destroyed, escaped, scatterings, J = simulation.value
    elapsed_time = simulation.time

    # Evaluate field above boundary
    meanJ, minJ, maxJ = field_above_boundary(atmosphere.z,
                                             atmopshere.chi_continuum, J, τ_max)
    # Signal to noise ratio
    SNR = sqrt(maximum(surface))

    # Print results
    println("Packets: ", packets)
    println("Destroyed: ", destroyed)
    println("Escaped: ", escaped)
    println("Scatterings: ", scatterings)
    println(@sprintf("Mean/Med/Min/Max field: %.1f / %g / %g", meanJ,  minJ, maxJ))
    println("S/N: ", SNR)

    # Append results to file
    write_results_to_file(current_time, τ_max, packets, destroyed,
                          escaped, scatterings, elapsed_time, meanJ, minJ, maxJ)

    # Plot surface intensity
    surface_intensity(surface, τ_max, packets)
    escape_direction(surface, τ_max, packets)
end

main()
