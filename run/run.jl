include("../src/mcrt.jl")
include("../src/populations.jl")

function run()
    println("\n", "="^91, "\n", " "^34,
            "SOLAR ATMOSPHERE MCRT",
            "\n", "="^91, "\n")

    # =============================================================================
    # LOAD ATMOSPHERE DATA
    # =============================================================================
    print("--Loading atmosphere data..................")
    atmosphere_parameters = collect_atmosphere_data()
    atmosphere = Atmosphere(atmosphere_parameters...)
    println("Atmosphere loaded with dimensions ", size(atmosphere.temperature), ".")

    if test_mode()
        # =============================================================================
        # INITIALISE OUTPUT FILE
        # =============================================================================
        output_path = get_output_path()
        max_scatterings = get_max_scatterings()
        target_packets = get_target_packets()
        cut_off = get_cut_off()

        # =============================================================================
        # LOAD WAVELENGTH
        # =============================================================================
        print("--Loading wavelength.......................")
        λ = get_background_λ()
        println("Wavelength λ = ", λ, " loaded.")

        # =============================================================================
        # LOAD RADIATION DATA
        # =============================================================================
        print("--Loading radiation data...................")
        radiation_parameters = collect_radiation_data(atmosphere, λ, cut_off, target_packets)
        radiation = RadiationBackground(radiation_parameters...)
        println(@sprintf("Radiation loaded with %.2e packets.", sum(radiation.packets)))

        # =============================================================================
        # SIMULATION
        # =============================================================================
        mcrt(atmosphere, radiation, max_scatterings, output_path)
        write_to_file(radiation, output_path)

        # =============================================================================
        # END OF TEST MODE
        # =============================================================================
    else
        # =============================================================================
        # INITIALISE OUTPUT FILE
        # =============================================================================
        output_path = get_output_path()
        max_iterations = get_max_iterations()
        max_scatterings = get_max_scatterings()
        target_packets = get_target_packets()
        cut_off = get_cut_off()

        population_distribution = get_population_distribution()
        # =============================================================================
        # LOAD ATOM
        # =============================================================================
        print("--Loading atom.............................")
        atom_parameters = collect_atom_data(atmosphere)
        atom = Atom(atom_parameters...)
        println("Atom loaded with ", atom.nλ_bb + 2*atom.nλ_bf, " wavelengths.")

        # =============================================================================
        # LOAD INITIAL POPULATIONS
        # =============================================================================
        print("--Loading initial populations..............")
        populations = collect_initial_populations(atom, atmosphere.temperature,
                                                        atmosphere.electron_density)
        println("Initial ", population_distribution, "-populations loaded.")

        # =============================================================================
        # CALCULATE INITIAL TRANSITION RATES
        # =============================================================================
        print("--Loading initial transition rates.........")
        Bλ = blackbody_lambda(atom.λ, atmosphere.temperature)
        rate_parameters = calculate_transition_rates(atom, atmosphere.temperature,
                                                           atmosphere.electron_density, Bλ)
        rates = TransitionRates(rate_parameters...)
        println("Initial transition rates loaded.")

        # =============================================================================
        # RUN MCRT UNTIL POPULATIONS CONVERGE
        # =============================================================================
        converged_populations = false

        for n=1:max_iterations
            println("\n  ITERATION ", n, "\n", "="^91)
            # =============================================================================
            # LOAD RADIATION DATA WITH CURRENT POPULATIONS
            # =============================================================================
            print("--Loading radiation data...................")
            radiation_parameters = collect_radiation_data(atmosphere, atom, rates, populations, cut_off, target_packets)
            radiation = Radiation(radiation_parameters...)
            println(@sprintf("Radiation loaded with %.2e packets per λ.", sum(radiation.packets[1,:,:,:])))

            # =============================================================================
            # SIMULATION
            # =============================================================================
            mcrt(atmosphere, radiation, atom, max_scatterings, output_path)

            # =============================================================================
            # CALCULATE NEW TRANSITION RATES
            # =============================================================================
            print("\n--Update transition rates..................")
            Jλ = get_Jλ(output_path, radiation.intensity_per_packet)
            rate_parameters = calculate_transition_rates(atom, atmosphere.temperature,
                                                               atmosphere.electron_density, Jλ)
            rates = TransitionRates(rate_parameters...)
            println("Transition rates updated.")

            # =============================================================================
            # CALCULATE NEW POPULATIONS
            # =============================================================================
            print("--Update populations.......................")
            new_populations = get_revised_populations(rates, atom.density)
            println("Populations updated.")

            # =============================================================================
            # CHECK POPULATION CONVERGENCE
            # =============================================================================
            converged, error = check_population_convergence(populations, new_populations, output_path)
            populations = copy(new_populations)

            if converged
                println(@sprintf("--Convergence at iteration n = %d. Error = %.1e.\n", n, error))
                write_to_file(new_populations, output_path)
                break
            else
                println(@sprintf("\n--No convergence. Error = %.1e.\n", error))
            end

            # =============================================================================
            # END OF ITERATION
            # =============================================================================
        end

        write_to_file(atom, output_path)
        write_to_file(radiation, output_path)

        # =============================================================================
        # END OF ATOM MODE
        # =============================================================================
    end
end

run()

function how_much_data()
