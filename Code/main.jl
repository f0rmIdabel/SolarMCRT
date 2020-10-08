include("MCRT.jl")
include("atmos.jl")
include("MyLibs/plotLib.jl")
include("MyLibs/IOLib.jl")
using Dates

function main(max_scatterings = 1e10)
    println("\nMCRT calculation in the solar atmosphere")
    println("----------------------------------------")
    println("\n--Reading atmosphere...")

    # Load atmosphere data
    parameters = get_atmosphere_data("bifrost_cb24bih_s385_fullv.ncdf",
                                     "output_ray.hdf5")
    atmosphere = Atmosphere(parameters...)

    # Choose wavelengths
    wavelength = 500u"nm"

    # Get user input: τ_max, # packets
    print("--Choose maximum optical depth: ")
    τ_max = readline()
    τ_max = parse(Float64, τ_max)

    print("--Choose number of photon packets: ")
    target_packets = readline()
    target_packets = parse(Float64, target_packets)

    # Number of threads used, >export JULIA_NUM_THREADS=4
    threads = Threads.nthreads()

    # Current time
    current_time = string(now())

    # Run and time simulation
    simulation = @timed simulate(atmosphere, wavelength,
                                 max_scatterings, τ_max, target_packets)

    packet_data, J_data, surface_intensity = simulation.value
    elapsed_time = simulation.time

    # Signal to noise ratio
    SNR = sqrt(maximum(surface_intensity))

    # Print results for quick check
    quick_print(threads, packet_data, J_data[2:4], SNR)

    # Append data to file
    write_results_to_file(current_time, threads, elapsed_time,
                          τ_max, packet_data, J_data[2:4], SNR)

    print("--Plotting stuff...")
    # Plot surface intensity and escape direction distribution
    plot_surface_intensity(surface_intensity, τ_max, packet_data[1])
    plot_escape_direction(surface_intensity, τ_max, packet_data[1])
    traverse_field_gif(J_data[1], 504)
end

main()
