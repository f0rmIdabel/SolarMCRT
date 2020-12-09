include("../src/mcrt.jl")

function run()
    println("\n", "="^83, "\n", " "^30,
            "SOLAR ATMOSPHERE MCRT",
            "\n", "="^83)

    # ==================================================================
    # LOAD WAVELENGTHS
    # ==================================================================
    print("\n--Loading wavelengths......................")
    λs = get_λ()
    println(@sprintf(" %d wavelength(s) loaded.", length(λs)))

    # ==================================================================
    # LOOP OVER WAVELENGTHS
    # ==================================================================

    for λ in λs
        println("\n--λ = " * string(λ) * ".............................")
        # ==================================================================
        # LOAD ATMOSPHERE DATA AND CALCULATE BOUNDARY
        # ==================================================================
        print("--Loading atmosphere data..................")
        atmosphere_parameters = collect_atmosphere_data(λ)
        atmosphere = Atmosphere(atmosphere_parameters...)
        println(@sprintf(" Atmosphere loaded with tau_max = %.1f.", get_cut_off()))

        # ==================================================================
        # LOAD RADIATION DATA AND CALCULATE # PACKETS
        # ==================================================================
        print("--Loading radiation data...................")
        radiation_parameters = collect_radiation_data(atmosphere, λ)
        radiation = Radiation(λ, radiation_parameters...)
        println(@sprintf(" Radiation loaded with %.2e packets.", sum(radiation.S)))

        # ==================================================================
        # SIMULATION
        # ==================================================================
        mcrt(atmosphere, radiation)
    end

end

run()
