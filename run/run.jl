include("../src/mcrt.jl")

function run()
    println("\n", "="^83, "\n", " "^30,
            "SOLAR ATMOSPHERE MCRT",
            "\n", "="^83, "\n")

    # ==================================================================
    # LOAD ATMOSPHERE DATA
    # ==================================================================
    print("--Loading atmosphere data..................")
    atmosphere_parameters = collect_atmosphere_data()
    atmosphere = Atmosphere(atmosphere_parameters...)
    println("Atmosphere loaded with dimensions ", size(atmosphere.temperature), ".")

    mode = get_mode()

    if mode == "test"

        # ==================================================================
        # LOAD WAVELENGTH
        # ==================================================================

        print("--Loading wavelength.......................")
        λ = get_test_λ()
        println("Wavelength λ = ", λ, " loaded.")


        # ==================================================================
        # LOAD RADIATION DATA
        # ==================================================================
        print("--Loading radiation data...................")
        radiation_parameters = collect_radiation_data(atmosphere, λ)
        radiation = Radiation(radiation_parameters...)
        println(@sprintf("Radiation loaded with %.2e packets.", sum(radiation.packets)))

        # ==================================================================
        # SIMULATION
        # ==================================================================
        mcrt(atmosphere, radiation)

    elseif mode == "atom"

        # ==================================================================
        # LOAD ATOM
        # ==================================================================
        print("--Loading atom.............................")
        atom_parameters = collect_atom_data()
        atom = AtomicLine(collect_atom_data()...)
        nλ_bb, nλ_bf = get_nλ()
        nλ = nλ_bf*2 + nλ_bb
        println("Atom loaded with ", nλ, " wavelengths.")

        # ==================================================================
        # LOAD INITIAL ATOM POPULATIONS
        # ==================================================================
        new_atom_populations = atmosphere.hydrogen_populations[:,:,:,1:2]
        converged = false

        while !converged
            atom_populations = new_atom_populations
            # ==================================================================
            # LOAD RADIATION DATA
            # ==================================================================
            print("--Loading radiation data...................")
            radiation_parameters = collect_radiation_data(atmosphere, atom, atom_populations)
            radiation = Radiation(radiation_parameters...)
            write_to_file(radiation)
            println(@sprintf("Radiation loaded with %.2e packets.", sum(radiation.packets[1,:,:,:])))

            # ==================================================================
            # SIMULATION
            # ==================================================================
            mcrt(atmosphere, radiation)

            new_atom_population = get_revised_populations(atom, J)
            converged = check_conversion(atom_populations, new_atom_populations)
        end

    end
end

run()
