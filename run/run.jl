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
        write_to_file(radiation)
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

        max_iterations = 10

        new_populations =  collect_initial_populations(atmosphere.hydrogen_populations)
        converged_populations = false
        error = Array{Float64,2}(undef, max_iterations, nλ)

        # ==================================================================
        # CALCULATE RADIATION PROPERTIES AND RUN MCRT UNTIL POP CONVERGE
        # ==================================================================

        for n=1:1#max_iterations
            println("\nITERATION ", n)
            println("="^83)
            populations = new_populations
            # ==================================================================
            # LOAD RADIATION DATA WITH CURRENT POPULATIONS
            # ==================================================================
            print("--Loading radiation data...................")
            radiation_parameters = collect_radiation_data(atmosphere, atom, populations)
            radiation = Radiation(radiation_parameters...)
            write_to_file(radiation)
            println(@sprintf("Radiation loaded with %.2e packets.", sum(radiation.packets[1,:,:,:])))

            # ==================================================================
            # SIMULATION
            # ==================================================================
            mcrt(atmosphere, radiation)

            """# ==================================================================
            # CHECK IF POPULATIONS CONVERGE
            # ==================================================================
            new_population = get_revised_populations(atom)
            converged = check_converged(populations, new_populations, error, n)

            if converged
                println("--Convergence at iteration n = ", n, ". Population-iteration finished.")
                break
            end"""

        end
    end
end

run()
