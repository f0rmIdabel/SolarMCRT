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
        # =======================================================================
        # LOAD ATOM
        # =======================================================================
        print("--Loading atom.............................")
        atom_parameters = collect_atom_data(atmosphere)
        atom = Atom(atom_parameters...)
        println("Atom loaded with ", atom.nλ_bb + 2*atom.nλ_bf, " wavelengths.")

        # =======================================================================
        # CALCULATE RADIATION PROPERTIES AND RUN MCRT UNTIL POPULATIONS CONVERGE
        # =======================================================================
        atom_density = sum(atom.populations, dims=4)[:,:,:,1]
        LTE_populations = LTE_populations(atom, atom_density, atmosphere.temperature, atmosphere.electron_density)
        tr = calculate_transition_rates(atom, LTE_populations, atomsphere.temperature, atmosphere.electron_density)
        transitionRates = TransitionRates(tr...)
        converged_populations = false

        for n=1:1#max_iterations
            println("\n  ITERATION ", n, "\n", "="^91)
            # ==================================================================
            # LOAD RADIATION DATA WITH CURRENT POPULATIONS
            # ==================================================================
            print("--Loading radiation data...................")
            radiation_parameters = collect_radiation_data(atmosphere, atom, transitionRates, populations)
            radiation = Radiation(radiation_parameters...)
            write_to_file(radiation)
            println(@sprintf("Radiation loaded with %.2e packets.",
                    sum(radiation.packets[1,:,:,:])))

            # ==================================================================
            # SIMULATION
            # ==================================================================
            mcrt(atmosphere, radiation)

            # ==================================================================
            # CHECK IF POPULATIONS CONVERGE
            # ==================================================================
            tr = calculate_transition_rates(atom, LTE_populations, atomsphere.temperature, atmosphere.electron_density)
            rates = TransitionRates(tr...)
            new_population = get_revised_populations(atom, LTE_populations, atmosphere.temperature)
            converged = check_converged(atom.populations, new_populations, error, n)

            if converged
                println("--Convergence at iteration n = ", n, ". Population-iteration finished.")
                break
            else
                # Update populations
                atom.populations = copy(new_population)
                atom.αlc = α_line_const(atom.line, new_populations[:,:,:,1], new_populations[:,:,:,2])
            end

        end
    end
end

run()
