include("mcRT.jl")
include("atmos.jl")
include("MyLibs/plotLib.jl")
include("MyLibs/ioLib.jl")
using Dates

function main(max_scatterings = 1e10)
    println("\n","-"^42,"\n MCRT calculation in the solar atmosphere \n","-"^42)
    println("\n--Reading atmosphere model...")

    # ==================================================================
    # LOAD ATMOSPHERE DATA
    # ==================================================================
    parameters = get_atmosphere_data("bifrost_cb24bih_s385_fullv.ncdf",
                                     "output_ray.hdf5")
    atmosphere = Atmosphere(parameters...)

    # ==================================================================
    # CHOOSE PARAMETERS
    # ==================================================================
    # Choose wavelengths
    wavelength = 499.86u"nm"

    # Get user input: τ_max, # packets
    print("--Choose maximum optical depth: ")
    τ_max = readline()
    τ_max = parse(Float64, τ_max)

    print("--Choose number of photon packets: ")
    target_packets = readline()
    target_packets = parse(Float64, target_packets)

    # ==================================================================
    # SIMULATION
    # ==================================================================
    threads = Threads.nthreads()
    current_time = string(now())

    # Run and time simulation
    simulation = @timed simulate(atmosphere, wavelength,
                                 max_scatterings, τ_max, target_packets)

    packet_data, J_data, surface_intensity = simulation.value
    elapsed_time = simulation.time

    # Signal to noise ratio
    SNR = sqrt(maximum(surface_intensity))

    # ==================================================================
    # WRITE RESULTS
    # ==================================================================
    # Print results for quick check
    quick_print(packet_data, J_data[3:5], SNR)

    # Append data to file
    write_results_to_file(current_time, threads, elapsed_time,
                          τ_max, packet_data, J_data[3:5], SNR)

    # ==================================================================
    # PLOTTING
    # ==================================================================
    println("\n--Plotting stuff...")
    plot_surface_intensity(surface_intensity, τ_max, packet_data[1], :[:,:])
    plot_escape_direction(surface_intensity, τ_max, packet_data[1])
    traverse_field_gif(log10.(J_data[1]), atmosphere.x, atmosphere.z)
    println("--Finished successfully\n")
end


main()
