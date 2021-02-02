include("../../src/radiation.jl")
import Plots
import Statistics

function check_parameters(atmosphere::Atmosphere, radiation::Radiation)

      # ==================================================================
      # LOAD PARAMETERS
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
