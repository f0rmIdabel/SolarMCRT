include("../src/mcrt.jl")

function run()
    println("\n", "="^91, "\n", " "^34,
            "SOLAR ATMOSPHERE MCRT",
            "\n", "="^91, "\n")

    # ==================================================================
    # LOAD ATMOSPHERE DATA
    # ==================================================================
    print("--Loading atmosphere data..................")
    atmosphere_parameters = collect_atmosphere_data()
    atmosphere = Atmosphere(atmosphere_parameters...)
    println("Atmosphere loaded with dimensions ", size(atmosphere.temperature), ".")

    if test_mode()

        # ==================================================================
        # LOAD WAVELENGTH
        # ==================================================================
        print("--Loading wavelength.......................")
        λ = get_background_λ()
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

    else

        # ==================================================================
        # LOAD ATOM
        # ==================================================================
        print("--Loading atom.............................")
        atom = Atom(collect_atom_data()...)
        println("Atom loaded with ", atom.nλ_bb + 2*atom.nλ_bf, " wavelengths.")

        # ==================================================================
        # LOAD INITIAL ATOM POPULATIONS
        # ==================================================================
        new_populations =  collect_initial_populations(atmosphere, atom)
        converged_populations = false

        #error = Array{Float64,2}(undef, max_iterations, nλ)
        # ==================================================================
        # CALCULATE RADIATION PROPERTIES AND RUN MCRT UNTIL POP CONVERGE
        # ==================================================================

        for n=1:1#max_iterations
            println("\n  ITERATION ", n)
            println("="^91)
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
