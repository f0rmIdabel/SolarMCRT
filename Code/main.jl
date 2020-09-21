using Random
using Dates
using Printf
using Statistics

import MCRT
import Bifrost
import MyPlots
import PhysLib
import IOLib

function main(atmosphere::Module, max_scatterings = 10_000_000_000 )

    # Get user input
    println("Choose maximum optical depth:")
    τ_max = readline()
    τ_max = parse(Float64, τ_max)

    println("Choose number of photon packets:")
    packets = readline()
    packets = parse(Float64, packets)

    # Set up rng with seed for reproducability
    current_time = now()
    seed = Int(floor(datetime2unix(current_time)))
    rng = MersenneTwister(seed)

    # Run and time simulation
    simulation = @timed MCRT.simulate(atmosphere, max_scatterings, τ_max, packets, rng)
    surface, destroyed, escaped, scatterings, J = simulation.value
    elapsed_time = simulation.time

    # Evaluate field
    meanJ = mean(J)
    medJ  = median(J)
    minJ  = minimum(J)
    maxJ  = maximum(J)

    # Print results
    println("Packets: ", packets)
    println("Destroyed: ", destroyed)
    println("Escaped: ", escaped)
    println("Scatterings: ", scatterings)
    println(@sprintf("Mean/Med/Min/Max field: %.1f / %.1f / %g / %g", meanJ,  medJ, minJ, maxJ))

    # Append results to file
    IOLib.write_results_to_file(current_time, τ_max, packets, destroyed,
                                escaped, scatterings, elapsed_time, meanJ,  medJ, minJ, maxJ)

    # Plot surface intensity
    MyPlots.surface_intensity(surface, τ_max, packets)
    MyPlots.escape_direction(surface, τ_max, packets)
end

main(Bifrost)
