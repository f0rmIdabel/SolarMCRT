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
        radiation = RadiationBackground(radiation_parameters...)
        write_to_file(radiation)
        println(@sprintf("Radiation loaded with %.2e packets.", sum(radiation.packets)))

        # ==================================================================
        # SIMULATION
        # ==================================================================
        mcrt(atmosphere, radiation)

    else
        # =======================================================================
        # LOAD ATOM
        # =======================================================================
        print("--Loading atom.............................")
        atom_parameters = collect_atom_data(atmosphere)
        atom = Atom(atom_parameters...)
        write_to_file(atom.λ)
        println("Atom loaded with ", atom.nλ_bb + 2*atom.nλ_bf, " wavelengths.")

        # =======================================================================
        # LOAD INITIAL POPULATIONS
        # =======================================================================
        print("--Loading initial populations..............")
        populations = collect_initial_populations()
        println("Initial populations loaded.")

        # =======================================================================
        # CALCULATE INITIAL TRANSITION RATES
        # =======================================================================
        print("--Loading transition rates.................")
        Bλ = blackbody_lambda(atom.λ, atmosphere.temperature)
        rate_parameters = calculate_transition_rates(atom, atmosphere, Bλ)
        rates = TransitionRates(rate_parameters...)
        println("Rates loaded.")

        # =======================================================================
        # RUN MCRT UNTIL POPULATIONS CONVERGE
        # =======================================================================
        converged_populations = false
        max_iterations = get_max_iterations()

        for n=1:max_iterations
            println("\n  ITERATION ", n, "\n", "="^91)
            # ==================================================================
            # LOAD RADIATION DATA WITH CURRENT POPULATIONS
            # ==================================================================
            print("--Loading radiation data...................")
            radiation_parameters = collect_radiation_data(atmosphere, atom, rates)
            radiation = Radiation(radiation_parameters...)
            write_to_file(radiation)
            println(@sprintf("Radiation loaded with %.2e packets.", sum(radiation.packets[1,:,:,:])))

            # ==================================================================
            # SIMULATION
            # ==================================================================
            mcrt(atmosphere, radiation)

            # =======================================================================
            # CALCULATE NEW TRANSITION RATES
            # =======================================================================
            print("--Update transition rates..............")
            Jλ = get_Jλ()
            rate_parameters = calculate_transition_rates(atom, atmosphere, Jλ)
            rates = TransitionRates(rate_parameters...)
            println("Rates updated.")

            # ==================================================================
            # CALCULATE NEW POPULATIONS
            # ==================================================================
            print("--Update populations................")
            new_populations = get_revised_populations(atom, rates)
            write_to_file(new_populations)
            println("Populations updated.")

            # ==================================================================
            # CHECK POPULATION CONVERGENCE
            # ==================================================================
            converged = check_population_convergence(populations, new_populations)
            populations = copy(new_population)

            if converged
                println("--Convergence at iteration n = ", n, ". Population-iteration finished.")
                break
            else
                println("--No convergence. Error = ", get_error(), ".")
            end
        end
    end
end

run()
