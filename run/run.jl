include("../src/mcrt.jl")

function run()
    println("\n", "="^83, "\n", " "^30,
            "SOLAR ATMOSPHERE MCRT",
            "\n", "="^83)

    # ==================================================================
    # LOAD ATMOSPHERE DATA
    # ==================================================================
    print("--Loading atmosphere data..................")
    atmosphere_parameters = collect_atmosphere_data()
    atmosphere = Atmosphere(atmosphere_parameters...)
    println("Atmosphere loaded with dimensions ", size(atmosphere.temperature), ".")

    mode = get_mode()

    if mode == "wave"

        # ==================================================================
        # LOAD WAVELENGTH
        # ==================================================================

        print("--Loading wavelength....................")
        λ = get_test_λ()
        println(λ)
        println("Wavelength λ = ", λ, " nm loaded.")

        # ==================================================================
        # LOAD RADIATION DATA
        # ==================================================================
        print("--Loading radiation data...................")
        radiation_parameters = collect_radiation_data(atmosphere, λ)
        radiation = Radiation(radiation_parameters...)
        write_to_file(radition)
        println(@sprintf(" Radiation loaded with %.2e packets.", sum(radiation.S)))

        # ==================================================================
        # SIMULATION
        # ==================================================================
        mcrt(atmosphere, radiation)

    else if mode == "atom"

        # ==================================================================
        # LOAD ATOM
        # ==================================================================
        atom_parameters = collect_atom_data()
        atom = AtomicLine(collect_atom_data()...)

        # ==================================================================
        # LOAD RADIATION DATA
        # ==================================================================
        print("--Loading radiation data...................")
        radiation_parameters = collect_radiation_data(atmosphere, atom)
        radiation = Radiation(radiation_parameters...)
        write_to_file(radition)
        println(@sprintf(" Radiation loaded with %.2e packets.", sum(radiation.S)))

        # ==================================================================
        # SIMULATION
        # ==================================================================
        mcrt(atmosphere, radiation)

    end

    # ==================================================================
    # LOOP OVER WAVELENGTHS
    # ==================================================================

    for λ in λs
        println("\n--λ = " * string(λ) * ".............................")


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
