include("../../src/radiation.jl")
import Plots
import Statistics
using Test

function test()
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
                println("--Convergence at iteration n = ", n, ".\n")
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

run()

function test_atmosphere(atmosphere::Atmosphere)

    z = atmosphere.z
    x = atmosphere.x
    y = atmosphere.y
    T = atmosphere.temperature
    v = atmosphere.velocity
    vz = atmosphere.velocity_z
    electron_density = atmosphere.electron_density
    hydrogen_populations = atmosphere.hydrogen_populations

    # ===========================================================
    # CHECK DIMENSIONS
    # ===========================================================
    size(T) = nz, nx, ny

    @assert size(T) == size(v) == size(electron_density) == size(vz) == size(v)
    @assert length(v[1,1,1]) == 3
    @assert size(hydrogen_populations) == (nz, nx, ny, 3)
    @assert length(z) == nz+1
    @assert length(x) == nx+1
    @assert length(y) == ny+1

    # ===========================================================
    # CHECK UNITS
    # ===========================================================
    @test all(unit.(z) .== unit(1u"m"))
    @test all(unit.(x) .== unit(1u"m"))
    @test all(unit.(y) .== unit(1u"m"))
    @test all(unit.(T) .== unit(1u"K"))
    @test all(unit.(v) .== unit(1u"m/s"))
    @test all(unit.(vz) .== unit(1u"m/s"))
    @test all(unit.(electron_density) .== unit(1u"m^-3"))
    @test all(unit.(hydrogen_populations) .== unit(1u"m^-3"))

    # ===========================================================
    # NO NEGAITVE VALUES
    # ===========================================================
    @test all(T .>= 0u"K")
    @test all(electron_density .>= 0u"m^-3")
    @test all(hydrogen_populations .>= 0u"m^-3")

    # ===========================================================
    # DECREASING Z, INCREASING X AND Y
    # ===========================================================
    dz = z[2:end] .- z[1:end-1]
    dx = x[2:end] .- x[1:end-1]
    dy = y[2:end] .- y[1:end-1]

    @test all(dz .<= 0u"m")
    @test all(dx .>= 0u"m")
    @test all(dy .>= 0u"m")
end


function test_radiation(atom::Atom)

    h5open("../out/output.h5", "w") do file
        α_continuum = read(file, "extinction_continuum")u"m^-1"
        ε_continuum = read(file, "destruction_continuum")
        ε_line = read(file, "destruction_line")
        α_line_constant = read(file, "extinction_line_constant")
        packets = read(file, "packets")
        boundary = read(file, "boundary")
        λ = read(file, "wavelength")
        intensity_per_packet = read(file, "intensity_per_packet")u"kW / m^2 / sr / nm"
    end

    λ0 = atom.line.λ0
    damping_constant = atom.damping_constant
    doppler_width = atom.doppler_width
    nλ_bb = atom.nλ_bb

    nλ,nz,nx,ny = size(α_continuum)
    α_line = Array{Unitful.PerLength, 4}(undef,nλ,nz,nx,ny)
    for l=1:nλ_bb
        α_line[l,:,:,:] = line_extinction.(λ, λ0, doppler_width, damping_constant, α_line_constant)
    end
end


function check_parameters(atmosphere::Atmosphere, radiation::Radiation)

      # ==================================================================
      # LOAD ATMOSPHERE PARAMETERS
      # ==================================================================
      z = atmosphere.z
      x = atmosphere.x
      y = atmosphere.y
      T = atmosphere.temperature

      λ = radiation.λ
      χ = radiation.χ
      ε = radiation.ε
      boundary = radiation.boundary
      packets = radiation.packets
      intensity_per_packet = radiation.intensity_per_packet

      nλ, nz, nx, ny = size(χ)




      # ==================================================================
      # AVG PARAMETERS
      # ==================================================================
      print("--Calculate averages.......................")

      mean_χ = Array{Float64, 2}(undef, nλ, nz)
      mean_ε = Array{Float64, 2}(undef, nλ, nz)
      mean_packets = Array{Float64, 2}(undef, nλ, nz)
      mean_boundary = Array{Float64, 1}(undef, nλ)
      mean_T = average_column(T)

      for l=1:nλ
          mean_χ[l,:] = average_column(χ[l,:,:,:])
          mean_ε[l,:] = average_column(ε[l,:,:,:])
          mean_packets[l,:] = average_column(packets[l,:,:,:])
          mean_boundary[l] = Statistics.mean(boundary)
      end

      println(" Averages calculated")

      # ==================================================================
      # PLOT PARAMETERS
      # ==================================================================
      plot(mean_χ[1,:], z[1:end-1], xlabel = "̄χ", ylabel = "z", xscale=:log10)

      """ENV["GKSwstype"]="nul"

      p1 = Plots.plot(mean_χ, z, xlabel = "̄χ", ylabel = "z", xscale=:log10)
      p2 = Plots.plot(mean_ε, z, xlabel = "̄ε", ylabel = "z")
      p3 = Plots.plot(mean_T, z, xlabel = "temperature", ylabel = "z" )
      p4 = Plots.plot(mean_packets, z, xlabel = "packets", ylabel = "z")
      Plots.plot(p1, p2, p3, p4, layout = (2, 2), legend = false)
      Plots.png("parameters")"""
end


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
        #write_to_file(radiation)
        println(@sprintf("Radiation loaded with %.2e packets.", sum(radiation.packets)))

    elseif mode == "atom"

        # ==================================================================
        # LOAD ATOM
        # ==================================================================
        print("--Loading atom.............................")
        atom_parameters = collect_atom_data()
        atom = AtomicLine(collect_atom_data()...)
        nλ_bb, nλ_bf = get_nλ()
        nλ = nλ_bf*2 + nλ_bb
        nλ += 1-nλ%2
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
            println("\n  ITERATION ", n)
            println("="^91)
            populations = new_populations
            # ==================================================================
            # LOAD RADIATION DATA WITH CURRENT POPULATIONS
            # ==================================================================
            print("--Loading radiation data...................")
            radiation_parameters = collect_radiation_data(atmosphere, atom, populations)
            radiation = Radiation(radiation_parameters...)
            #write_to_file(radiation)
            println(@sprintf("Radiation loaded with %.2e packets.", sum(radiation.packets[1,:,:,:])))
        end
    end

    check_parameters(atmosphere, radiation)
end

function average_column(array)
      Statistics.mean(array, dims=[2,3])[:,1,1]
end

run()
