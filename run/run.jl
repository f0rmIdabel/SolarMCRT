include("../src/mcrt.jl")

function run()
    println("\n", "="^83, "\n", " "^30,
            "SOLAR ATMOSPHERE MCRT",
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
    atmosphere_parameters = collect_atmosphere_data(λ, false)
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
    # SIMULATION
    # ==================================================================
    mcrt(atmosphere, radiation)
end

run()
