using Plots
using Random
using Dates
using Printf
using TickTock

import MCRT
import Bifrost


function main(atmosphere::Module)
    
    #atmosphere = Bifrost
    max_scatterings = 10_000_000_000 
    
    # Get user input 
    println("Choose maximum optical depth:")
    tau_max = readline()
    tau_max = parse(Float64, tau_max)

    println("Choose number of photon packets:")
    packets = readline()
    packets = parse(Float64, packets)
    
    # Set up rng with seed for reproducability
    current_time = now()
    seed = Int(floor(datetime2unix(current_time)))
    rng = MersenneTwister(seed)
    
    # Time function -- but use something else like @time 
    tick()
    # Run simulation
    surface, destroyed, escaped, scatterings, J = MCRT.simulate(atmosphere, max_scatterings, 
                                                             tau_max, packets, rng)

    elapsed_time = tok()
    # Print results
    println("Packets: ", packets) 
    println("Destroyed: ", destroyed)
    println("Escaped: ", escaped)
    println("Scatterings: ", scatterings)
    
    # Append results to file
    f = open("/mn/stornext/u3/idarhan/SolarMCRT/Results/Results.txt", "a")
    results = string(current_time)*(@sprintf("%13.1f%13.1e%20g%18g%22g%22.1f\n", tau_max, packets, destroyed, escaped, scatterings, elapsed_time))
    write(f, results)
    
    # To avoid ssh display problems 
    ENV["GKSwstype"]="nul"

    heatmap(1:size(surface,1), 1:size(surface,2), surface, c=:grays, title = "")
    fig = @sprintf("/mn/stornext/u3/idarhan/SolarMCRT/Results/Plots/Surface/bf_%.1f_%.0e", tau_max, packets)
    png(fig)
    
end

main(Bifrost)