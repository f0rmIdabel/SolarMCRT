using Printf

function quick_print(threads::Int64, packet_data::Array{Int64,1},
                     J_data::Array{Any,1}, SNR::Float64)

    packets, destroyed, escaped, scatterings = packet_data
    mean_J, min_J, max_J = J_data

    println("\nSneak peak:\n----------- ")
    println("Threads: ", threads)
    println("Packets: ", packets) #check actual value
    println("Destroyed: ", destroyed/packets)
    println("Escaped: ", escaped/packets)
    println("Scatterings: ", scatterings/packets)
    println(@sprintf("Mean/Med/Min/Max field: %.1f / %d / %d",
            mean_J,  min_J, max_J))
    println("S/N: ", SNR)
end

function write_results_to_file(time::String, threads::Int64, elapsed_time::Float64,
                               τ_max::Real, packet_data::Array{Int64, 1},
                               J_data::Array{Any, 1}, SNR::Real)

    packets, destroyed, escaped, scatterings = packet_data
    mean_J, min_J, max_J = J_data

    f = open("/mn/stornext/u3/idarhan/SolarMCRT/Results/Results.txt", "a")
    results = (@sprintf("%23s%13d%19.1f%13.1f%20d%20d%20d%20d%22.1f%22d%22d%13.1f\n",
                         time, threads, elapsed_time, τ_max,
                         packets, destroyed, escaped, scatterings,
                         mean_J, min_J, max_J, SNR))
    write(f, results)
    close(f)
end
