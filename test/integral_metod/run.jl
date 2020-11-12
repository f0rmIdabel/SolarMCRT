include("feautrier.jl")

function run()
    println("\n", "="^83, "\n", " "^20,
              "FEAUTRIER RT IN A SOLAR ATMOSPHERE MODEL",
              "\n", "="^83)

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
      atmosphere_parameters = collect_atmosphere_data(λ, false, false)
      atmosphere = Atmosphere(atmosphere_parameters...)
      println(" Full atmosphere loaded.")

      # ==================================================================
      # FEAUTRIER CALCULATION
      # ==================================================================
      print("--Starting Feautrier.......................")
      feautrier(atmosphere, λ, 3, 4)
      println(" Feautrier finished.")
end

run()
