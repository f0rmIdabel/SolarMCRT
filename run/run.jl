include("../src/mcrt.jl")
include("../src/populations.jl")
include("../test/check_parameters.jl")

function run()
    println("\n", "="^91, "\n", " "^34,
            "SOLAR ATMOSPHERE MCRT",
            "\n", "="^91, "\n")

    # =============================================================================
    # LOAD ATMOSPHERE DATA
    # =============================================================================
    print("--Loading atmosphere data..................")
    atmosphere_parameters = collect_atmosphere()
    atmosphere = Atmosphere(atmosphere_parameters...)
    atmosphere_size = size(atmosphere.temperature)
    println("Atmosphere loaded with dimensions ", atmosphere_size, ".")

    if background_mode()
        # =============================================================================
        # READ CONFIG FILE
        # =============================================================================
        output_path = get_output_path()
        max_scatterings = get_max_scatterings()
        target_packets = get_target_packets()
        cut_off = get_cut_off()

        # =============================================================================
        # LOAD WAVELENGTH
        # =============================================================================
        print("--Loading wavelength.......................")
        λ = [get_background_λ()]
        nλ = length(λ)
        println("Wavelength λ = ", λ[1], " loaded.")

        # =============================================================================
        # LOAD RADIATION DATA
        # =============================================================================
        print("--Loading radiation data...................")
        radiation_parameters = collect_background_radiation(atmosphere, λ, cut_off, target_packets)
        radiation = RadiationContinuum(radiation_parameters...)
        println(@sprintf("Radiation loaded with %.2e packets.", sum(radiation.packets)))

        # =============================================================================
        # CREATE OUTPUT FILE
        # =============================================================================
        print("--Initialise output file...................")
        create_output_file(output_path, nλ, atmosphere_size)
        write_to_file(λ, output_path)
        write_to_file(radiation, output_path)
        println(@sprintf("%.1f GBs of data initialised.", how_much_data(nλ, atmosphere_size)))

        # =============================================================================
        # SIMULATION
        # =============================================================================
        mcrt_continuum(atmosphere, radiation, λ, max_scatterings, 1, 0, 1, output_path)

        # =============================================================================
        # END OF TEST MODE
        # =============================================================================
    else
        # =============================================================================
        # READ CONFIG FILE
        # =============================================================================
        output_path = get_output_path()
        max_iterations = get_max_iterations()
        max_scatterings = get_max_scatterings()
        target_packets_bf = get_target_packets_bf()
        target_packets_bb = get_target_packets_bb()
        cut_off = get_cut_off()
        population_distribution = get_population_distribution()
        write_rates = get_write_rates()

        # =============================================================================
        # LOAD ATOM
        # =============================================================================
        print("--Loading atom.............................")
        atom_parameters = collect_atom_data(atmosphere)
        atom = Atom(atom_parameters...)
        println("Atom loaded with ", atom.nλ , " wavelengths.")

        n_levels = atom.n_levels
        n_lines = atom.n_lines
        n_transitions = n_levels + n_lines

        # =============================================================================
        # LOAD INITIAL POPULATIONS
        # =============================================================================
        print("--Loading initial populations..............")
        populations = collect_initial_populations(atmosphere, atom)
        println("Initial ", population_distribution, "-populations loaded.")

        # =============================================================================
        # CALCULATE INITIAL TRANSITION RATES
        # =============================================================================
        print("--Loading initial transition rates.........")

        Bλ = []

        for t=1:n_transitions
            append!(Bλ, [blackbody_lambda(atom.λ[t], atmosphere.temperature)])
        end

        rate_parameters = calculate_transition_rates(atmosphere, atom, Bλ)

        rates = TransitionRates(rate_parameters...)
        println("Initial transition rates loaded.")

        # =============================================================================
        # CREATE OUTPUT FILE
        # =============================================================================
        print("--Initialise output file...................")
        create_output_file(output_path, max_iterations, atom.nλ, n_levels, atmosphere_size, write_rates)
        write_to_file(atom.λ, output_path)
        write_to_file(populations, 0, output_path)
        if write_rates; write_to_file(rates, 0, output_path); end
        println(@sprintf("%.1f GBs of data initialised.", how_much_data(atom.nλ, atmosphere_size, max_iterations, write_rates)))

        # =============================================================================
        # RUN MCRT UNTIL POPULATIONS CONVERGE
        # =============================================================================
        converged_populations = false

        for n=1:max_iterations
            println("\n  ITERATION ", n, "\n", "="^91, "\n", "="^91)

            nλ0 = 0
            for level=1:n_levels

                println("\n  TRANSITION ", level, " -> c", "\n", "="^91)

                λ = atom.λ[level]

                # =============================================================================
                # LOAD RADIATION DATA WITH CURRENT POPULATIONS
                # =============================================================================
                print("--Loading radiation data...................")
                radiation_parameters = collect_bf_radiation(atmosphere, atom, level, rates, populations, cut_off, target_packets_bf[level])
                radiation = RadiationContinuum(radiation_parameters...)
                write_to_file(radiation, n, nλ0, output_path)
                println(@sprintf("Radiation loaded with %.2e packets per λ.", sum(radiation.packets[1,:,:,:])))

                # =============================================================================
                # SIMULATION
                # =============================================================================
                mcrt_continuum(atmosphere, radiation, λ, max_scatterings, n, nλ0, atom.nλ, output_path)
                nλ0 += length(λ)
            end

            line_number = 0
            for l=1:n_levels-1
                for u=l+1:n_levels

                    println("\n  TRANSITION ", l, " -> ", u, "\n", "="^91)

                    line_number += 1
                    λ = atom.λ[n_levels + line_number]

                    line_parameters = collect_line_data(atmosphere, atom, u, l)
                    line = Line(line_parameters...)

                    # =============================================================================
                    # LOAD RADIATION DATA WITH CURRENT POPULATIONS
                    # =============================================================================
                    print("--Loading radiation data...................")
                    radiation_parameters = collect_bb_radiation(atmosphere, λ, line, rates, populations, cut_off, target_packets_bb[line_number])
                    radiation = RadiationLine(radiation_parameters...)
                    write_to_file(radiation, n, nλ0, output_path)
                    println(@sprintf("Radiation loaded with %.2e packets per λ.", sum(radiation.packets[1,:,:,:])))

                    # =============================================================================
                    # SIMULATION
                    # =============================================================================
                    mcrt_line(atmosphere, radiation, λ, line, max_scatterings, n, nλ0, atom.nλ, output_path)
                    nλ0 += length(λ)
                end
            end

            # =============================================================================
            # CALCULATE NEW TRANSITION RATES
            # =============================================================================
            print("\n--Update transition rates..................")
            Jλ = get_Jλ(output_path, n, atom.λ)
            rate_parameters = calculate_transition_rates(atmosphere, atom, Jλ)
            rates = TransitionRates(rate_parameters...)
            if write_rates; write_to_file(rates, n, output_path); end
            println("Transition rates updated.")

            # =============================================================================
            # CALCULATE NEW POPULATIONS
            # =============================================================================
            print("--Update populations.......................")
            new_populations = get_revised_populations(rates, atom.density)
            write_to_file(new_populations, n, output_path)
            println("Populations updated.")

            # =============================================================================
            # CHECK FOR UNVALID VARIABLES
            # =============================================================================
            #check_field(Jλ)
            #check_rates(rates, atmosphere_size)
            #check_populations(new_populations, atmosphere_size)

            # =============================================================================
            # CHECK POPULATION CONVERGENCE
            # =============================================================================
            converged, error = check_population_convergence(populations, new_populations)
            populations = copy(new_populations)

            if converged
                println(@sprintf("--Convergence at iteration n = %d. Error = %.1e.\n", n, error))
                cut_output_file(output_path, n, write_rates)
                break
            else
                println(@sprintf("\n--No convergence. Error = %.1e.\n", error))
            end

            # =============================================================================
            # END OF ITERATION
            # =============================================================================
        end

        # =============================================================================
        # END OF ATOM MODE
        # =============================================================================
    end

end

run()
