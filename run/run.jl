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
        # LOAD WAVELENGTH
        # =============================================================================
        print("--Loading wavelength.......................")
        λ = get_background_λ()
        println("Wavelength λ = ", λ, " loaded.")

        # =============================================================================
        # LOAD RADIATION DATA
        # =============================================================================
        print("--Loading radiation data...................")
        radiation_parameters = collect_radiation_data(atmosphere, λ)
        radiation = RadiationBackground(radiation_parameters...)
        write_to_file(radiation) # creates new file
        println(@sprintf("Radiation loaded with %.2e packets.", sum(radiation.packets)))

        # =============================================================================
        # SIMULATION
        # =============================================================================
        mcrt(atmosphere, radiation)

    else
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
        populations = collect_initial_populations()
        println("Initial populations loaded.")

        # =============================================================================
        # CALCULATE INITIAL TRANSITION RATES
        # =============================================================================
        print("--Loading initial transition rates.........")
        Bλ = blackbody_lambda(atom.λ, atmosphere.temperature)
        rate_parameters = calculate_transition_rates(atom, atmosphere, populations, Bλ)
        rates = TransitionRates(rate_parameters...)
        println("Initial transition rates loaded.")

        # =============================================================================
        # RUN MCRT UNTIL POPULATIONS CONVERGE
        # =============================================================================
        converged_populations = false
        max_iterations = get_max_iterations()

        for n=1:max_iterations
            println("\n  ITERATION ", n, "\n", "="^91)
            # =============================================================================
            # LOAD RADIATION DATA WITH CURRENT POPULATIONS
            # =============================================================================
            print("--Loading radiation data...................")
            radiation_parameters = collect_radiation_data(atmosphere, atom, rates, populations)
            radiation = Radiation(radiation_parameters...)
            write_to_file(radiation) # creates new file
            write_to_file(atom.λ)
            println(@sprintf("Radiation loaded with %.2e packets per λ.", sum(radiation.packets[1,:,:,:])))

            # =============================================================================
            # SIMULATION
            # =============================================================================
            mcrt(atmosphere, radiation, atom)

            # =============================================================================
            # CALCULATE NEW TRANSITION RATES
            # =============================================================================
            print("\n--Update transition rates..................")
            Jλ = get_Jλ()
            rate_parameters = calculate_transition_rates(atom, atmosphere, populations, Jλ)
            rates = TransitionRates(rate_parameters...)
            println("Transition rates updated.")

            # =============================================================================
            # CALCULATE NEW POPULATIONS
            # =============================================================================
            print("--Update populations.......................")
            new_populations = get_revised_populations(atom, rates, populations)
            write_to_file(new_populations)
            println("Populations updated.")

            # =============================================================================
            # CHECK POPULATION CONVERGENCE
            # =============================================================================
            converged = check_population_convergence(populations, new_populations, n)
            populations = copy(new_populations)

            if converged
                println("--Convergence at iteration n = ", n, ". Error = ", get_error(n),"\n")
                break
            else
                println("--No convergence. Error = ", get_error(n), ".\n")
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

#run()


function fix()
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
    populations = collect_initial_populations()
    println("Initial populations loaded.")

    original = 1
    h5open("../../../basement/MScProject/Atmospheres/raw/bifrost_qs006023_s525_quarter.hdf5", "r") do file
        original = read(file, "hydrogen_populations")[:,:,:,:,1]u"m^-3"
    end

    pop = LTE_populations(atom, original, atmosphere.temperature, atmosphere.electron_density)

    h5open("../../../basement/MScProject/Atmospheres/bifrost_qs006023_s525_quarter_reworked.hdf5", "r+") do file
        write(file, "hydrogen_populations", ustrip.(pop))
    end

end

fix()
