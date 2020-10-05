include("MCRT.jl")
include("atmos.jl")
include("MyLibs/plotLib.jl")
include("MyLibs/physLib.jl")
using Dates

function main(max_scatterings = 1e10)

    println("Reading atmosphere...")

    # Load atmosphere data
    parameters = get_atmosphere_data("bifrost_cb24bih_s385_fullv.ncdf",
                                     "output_ray.hdf5")
    atmosphere = Atmosphere(parameters...)

    # Choose wavelengths
    wavelength = 500u"nm"

    # Get user input: τ_max, packets
    println("Choose maximum optical depth:")
    τ_max = readline()
    τ_max = parse(Float64, τ_max)

    println("Choose number of photon packets:")
    packets = readline()
    packets = parse(Float64, packets)

    # Number of threads used, >export JULIA_NUM_THREADS=4
    threads = Threads.nthreads()

    # Current time
    current_time = now()

    # Run and time simulation
    simulation = @timed simulate(atmosphere, wavelength,
                                 max_scatterings, τ_max, packets)

    surface, destroyed, escaped, scatterings, J = simulation.value
    elapsed_time = simulation.time

    # Evaluate field above boundary
    meanJ, minJ, maxJ = field_above_boundary(atmosphere.z,
                                             atmosphere.χ_continuum,
                                             J, τ_max)
    # Signal to noise ratio
    SNR = sqrt(maximum(surface))

    # Print results for quick check
    println("Sneak peak:\n----------- ")
    println("Threads: ", threads)
    println("Packets: ", packets) #check actual value
    println("Destroyed: ", destroyed)
    println("Escaped: ", escaped)
    println("Scatterings: ", scatterings)
    println(@sprintf("Mean/Med/Min/Max field: %.1f / %d / %d", meanJ,  minJ, maxJ))
    println("S/N: ", SNR)

    # Append results to file
    """f = open("../Results/Results.txt", "a")
    results = (@sprintf("%23s%13.1f%13.1e%20d%18d%22d%22.1f%22.1f%17d%17d\n",
                         string(current_time), τ_max, packets, destroyed, escaped, scatterings,
                         elapsed_time, meanJ, minJ, maxJ))
    write(f, results)"""

    # Plot surface intensity
    surface_intensity(surface, τ_max, packets)
    escape_direction(surface, τ_max, packets)
end

main()
