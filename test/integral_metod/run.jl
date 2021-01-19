include("lambda_iteration.jl")


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
        println(@sprintf("Radiation loaded with %.2e packets.", sum(radiation.S)))

        # ==================================================================
        # FEAUTRIER CALCULATION
        # ==================================================================
        print("--Starting λ-iteration.....................")
        lambda_iteration(atmosphere, radiation)
        println(" λ-iteration finished.")

    # Needs massive fix
    elseif mode == "atom"

        # ==================================================================
        # LOAD ATOM
        # ==================================================================
        atom_parameters = collect_atom_data()
        atom = AtomicLine(collect_atom_data()...)

        # ==================================================================
        # LOAD INITIAL ATOM POPULATIONS
        # ==================================================================
        #populations = collect_initial_populations()

        hydrogen_populations = atmosphere.hydrogen_populations
        nz, ny, nz, nl = size(hydrogen_populations)
        atom_populations = Array{Float64, 4}(undef, nz, nx, ny, 3)
        atom_populations[:,:,:,1:2] = atmosphere.hydrogen_populations[:,:,:,1:2]
        atom_populations[:,:,:,end] = atmosphere.hydrogen_populations[:,:,:,end]

        # ==================================================================
        # LOAD RADIATION DATA
        # ==================================================================
        print("--Loading radiation data...................")
        radiation_parameters = collect_radiation_data(atmosphere, atom, populations)
        radiation = Radiation(radiation_parameters...)
        write_to_file(radition)
        println(@sprintf(" Radiation loaded with %.2e packets.", sum(radiation.S)))

        # ==================================================================
        # FEAUTRIER CALCULATION
        # ==================================================================
        print("--Starting λ-iteration.....................")
        lambda_iteration(atmosphere, radiation)
        println(" λ-iteration finished.")

    end
end


run()
