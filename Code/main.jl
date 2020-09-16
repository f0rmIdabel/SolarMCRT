using Plots
using Random
using Dates

import MCRT
import Bifrost
import MyPlots
import IOLib

function main(atmosphere::Module, max_scatterings = 10_000_000_000 )
    
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
    
    # Run and time simulation
    simulation = @timed MCRT.simulate(atmosphere, max_scatterings, tau_max, packets, rng)
    surface, destroyed, escaped, scatterings, J = simulation.value
    elapsed_time = simulation.time

    # Print results
    println("Packets: ", packets) 
    println("Destroyed: ", destroyed)
    println("Escaped: ", escaped)
    println("Scatterings: ", scatterings)
    
    # Append results to file
    IOLib.write_results_to_file(current_time, tau_max, packets, destroyed, escaped, scatterings, elapsed_time)
    
    # Plot surface intensity
    MyPlots.surface_intensity(surface, tau_max, packets)
    
    
end

main(Bifrost)