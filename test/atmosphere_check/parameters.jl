include("../../src/radiation.jl")
import Plots
import Statistics

function parameters()

      # ==================================================================
      # LOAD WAVELENGTHS
      # ==================================================================
      print("\n--Loading wavelengths......................")
      λ = get_λ()[1]
      println(@sprintf(" %d wavelength(s) loaded.", length(λ)))

      # ==================================================================
      # LOAD ATMOSPHERE DATA AND CALCULATE BOUNDARY
      # ==================================================================
      print("--Loading atmosphere data..................")
      atmosphere_parameters = collect_atmosphere_data(λ)
      atmosphere = Atmosphere(atmosphere_parameters...)
      println(@sprintf(" Atmosphere loaded with tau_max = %.1f.", get_τ_max()))

      # ==================================================================
      # LOAD RADIATION DATA AND CALCULATE # PACKETS
      # ==================================================================
      print("--Loading radiation data...................")
      radiation_parameters = collect_radiation_data(atmosphere, λ)
      radiation = Radiation(λ, radiation_parameters...)
      println(@sprintf(" Radiation loaded with %.2e packets.", sum(radiation.S)))

      # ==================================================================
      # LOAD PARAMETERS
      # ==================================================================

      x = ustrip(atmosphere.x)
      y = ustrip(atmosphere.y)
      z = ustrip(atmosphere.z)[1:end-1]
      χ = ustrip(atmosphere.χ)
      ε = ustrip(atmosphere.ε)
      boundary = ustrip(atmosphere.boundary)
      T = ustrip(atmosphere.temperature)
      packets = ustrip(radiation.S)

      # ==================================================================
      # AVG PARAMETERS
      # ==================================================================
      print("--Calculate averages.......................")
      mean_χ = average_column(χ)
      mean_ε = average_column(ε)
      mean_T = average_column(T)
      mean_packets = average_column(packets)
      println(" Averages calculated")

      # ==================================================================
      # PLOT PARAMETERS
      # ==================================================================

      ENV["GKSwstype"]="nul"

      p1 = Plots.plot(mean_χ, z, xlabel = "̄χ", ylabel = "z", xscale=:log10)
      p2 = Plots.plot(mean_ε, z, xlabel = "̄ε", ylabel = "z")
      p3 = Plots.plot(mean_T, z, xlabel = "temperature", ylabel = "z" )
      p4 = Plots.plot(mean_packets, z, xlabel = "packets", ylabel = "z")
      Plots.plot(p1, p2, p3, p4, layout = (2, 2), legend = false)
      Plots.png("parameters")

end

function average_column(array)
      Statistics.mean(array, dims=[2,3])[:,1,1]
end

parameters()
