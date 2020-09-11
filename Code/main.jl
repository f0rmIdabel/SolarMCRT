import MCRT

using Plots

function main_Bifrost()

    # 
    max_scatterings = 10_000_000_000 
    scale_emission = 1e4 * 2
    tau_max = 20

    # 2-4 minutes runtime
    #1e7, 1e-3
    #1e6, 1e-2
    #1e5, 1e-1
    #1e4, 1e0
    
    atmosphere = Bifrost

    hits, destroyed, escaped, total_packets, avg_scatterings = MCRT.simulate(atmosphere, max_scatterings, scale_emission, tau_max)

    println("Packets: ", total_packets)
    println("Destroyed: ", destroyed/total_packets)
    println("Avg scatterings: ", avg_scatterings)
    println("Escaped: ", escaped/total_packets)
    println("Total hits:", sum(hits))

    heatmap(1:size(hits,1), 1:size(hits,2), hits, c=:grays)
    savefig("/mn/stornext/u3/idarhan/SolarMCRT/Plots/Bifrost_hits.png")
end

import Bifrost
main_Bifrost()

function main_FALC()
    max_scatterings = 10_000_000_000
    scale_emission = 1e6
    tau_max = 20.0

    atmosphere = FALC

    hits, destroyed, escaped, total_packets, avg_scatterings = MCRT.simulate(atmosphere, max_scatterings, scale_emission, tau_max)

    println("Packets: ", total_packets)
    println("Destroyed: ", destroyed/total_packets)
    println("Avg scatterings: ", avg_scatterings)
    println("Escaped: ", escaped/total_packets)
    println("Total hits:", sum(hits))

    heatmap(1:size(hits,1), 1:size(hits,2), hits, c=:grays)
    savefig("/mn/stornext/u3/idarhan/SolarMCRT/Plots/FALC_hits.png")
end

#import FALC
#main_FALC()