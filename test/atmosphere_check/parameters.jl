include("../../src/radiation.jl")
using Plots
using UnitfulRecipes

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
      # PARAMETERS
      # ==================================================================

      x = ustrip(atmosphere.x)
      y = ustrip(atmosphere.y)
      z = ustrip(atmosphere.z)
      χ = ustrip(atmosphere.χ)
      ε = ustrip(atmosphere.ε)
      boundary = ustrip(atmosphere.boundary)
      temperature = ustrip(atmosphere.temperature)
      S = ustrip(radiation.S)

      ENV["GKSwstype"]="nul"

      plot(z[1:end-1], χ[:,100,100])
      png("chi")

      plot(z[1:end-1], ε[:,100,100])
      png("epsilon")

      #plot(z[1:end-1], temperature[:,100,100])
      plot(z[1:end-1], S[:,100,100])
      png("temperature")



end

parameters()
