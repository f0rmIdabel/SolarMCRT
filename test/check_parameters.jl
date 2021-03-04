include("../../src/radiation.jl")
import Plots
import Statistics
using Test

function full_check()

    # =============================================================================
    # ATMOSPHERE
    # =============================================================================
    atmosphere_parameters = collect_atmosphere_data()
    atmosphere = Atmosphere(atmosphere_parameters...)

    check_atmosphere(atmosphere)
    plot_atmosphere(atmosphere)

    # =============================================================================
    # BACKGROUND RADIATION
    # =============================================================================
    λ = get_background_λ()
    radiation_parameters = collect_radiation_data(atmosphere, λ)
    radiationBackground = RadiationBackground(radiation_parameters...)

    check_radiationBackground(radiationBackground)
    plot_radiationBackground(radiationBackground, atmosphere.z)

    # =============================================================================
    # ATOM
    # =============================================================================
    atom_parameters = collect_atom_data(atmosphere)
    atom = Atom(atom_parameters...)

    check_atom(atom)

    # =============================================================================
    # INITIAL POPULATIONS
    # =============================================================================
    populations = collect_initial_populations()

    check_populations(populations)
    plot_populations(populations, atmosphere.z)

    # =============================================================================
    # INITIAL TRANSITION RATES
    # =============================================================================
    Bλ = blackbody_lambda(atom.λ, atmosphere.temperature)
    rate_parameters = calculate_transition_rates(atom, atmosphere, populations, Bλ)
    rates = TransitionRates(rate_parameters...)

    check_rates(rates)
    plot_rates(rates, atmosphere.z)

    # =============================================================================
    # RADIATION
    # =============================================================================
    radiation_parameters = collect_radiation_data(atmosphere, atom, rates, populations)
    radiation = Radiation(radiation_parameters...)

    check_radiation(radiation)
    plot_radiation(radiation, atmosphere.z, atom.λ)
end


function check_atmosphere(atmosphere::Atmosphere)
    # ===========================================================
    # LOAD DATA
    # ===========================================================
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
    @test all( dimension.(z) .== Unitful.𝐋 )
    @test all( dimension.(x) .== Unitful.𝐋 )
    @test all( dimension.(y) .== Unitful.𝐋 )
    @test all( dimension.(T) .== Unitful.𝚯 )
    @test all( dimension.(v) .== Unitful.𝐋 * Untitful.𝐓^-1 )
    @test all( dimension.(vz) .== Unitful.𝐋 * Untitful.𝐓^-1 )
    @test all( dimension.(electron_density) .== Unitful.𝐋^-3 )
    @test all( dimension.(hydrogen_populations) .== Unitful.𝐋^-3 )

    # ===========================================================
    # NO NEGAITVE VALUES
    # ===========================================================
    @test all( ustrip.(T) .>= 0.0 )
    @test all( ustrip.(electron_density) .>= 0.0 )
    @test all( ustrip.(hydrogen_populations) .>= 0.0 )

    # ===========================================================
    # DECREASING Z, INCREASING X AND Y
    # ===========================================================
    dz = z[2:end] .- z[1:end-1]
    dx = x[2:end] .- x[1:end-1]
    dy = y[2:end] .- y[1:end-1]

    @test all( ustrip.(dz) .<= 0.0 )
    @test all( ustrip.(dx) .<= 0.0 )
    @test all( ustrip.(dy) .<= 0.0 )
end

function plot_atmosphere(atmosphere::Atmosphere)
end

function check_radiationBackground(radiationBackground)
    # ===========================================================
    # LOAD DATA
    # ===========================================================
    λ = radiationBackground.λ
    α_continuum = radiationBackground.α_continuum
    ε_continuum = radiationBackground.ε_continuum
    boundary = radiationBackground.boundary
    packets = radiationBackground.packets
    intensity_per_packet = radiationBackground.intensity_per_packet

    # ===========================================================
    # CHECK DIMENSIONS
    # ===========================================================
    size(α_continuum) = nλ, nz, nx, ny

    @assert size(α_continuum) == size(α_continuum) == size(packets)
    @assert size(boundary) == (nλ, nx, ny)
    @assert length(λ) == length(intensity_per_packet) == nλ

    # ===========================================================
    # CHECK UNITS
    # ===========================================================
    @test all( dimension.(λ) .==  Unitful.𝐋 )
    @test all( dimension.(α_continuum) .== Unitful.𝐋^-1 )
    @test all( dimension.(ε_continuum) .== NoDims
    @test all( dimension.(boundary) .== NoDims
    @test all( dimension.(packets) .== NoDims
    @test all( dimension.(intensity_per_packet) .== Unitful.𝐋^-1 * Unitful.𝐌 * Unitful.𝐓^-3 )

    # ===========================================================
    # NO NEGAITVE VALUES
    # ===========================================================
    @test all( ustrip.(λ) .>= 0.0 )
    @test all( ustrip.(α_continuum) .>= 0.0 )
    @test all( ustrip.(intensity_per_packet) .>= 0.0 )
    @test all( ε_continuum .>= 0.0 )
    @test all( boundary .>= 0 )
    @test all( packets .>= 0 )
end

function plot_radiationBackground(radiationBackground, atmosphere.z, λ)
end

function check_atom(atom, atmosphere_size)
    # ===========================================================
    # LOAD DATA
    # ===========================================================
    line = atom.line
    Aul = line.Aul
    Bul = line.Bul
    Blu = line.Bli
    λ0 = line.λ0
    doppler_width = atom.doppler_width
    damping_constant = atom.damping_constant
    λ = atom.λ
    nλ_bb = atom.nλ_bb
    nλ_bf = atom.nλ_bf

    # ===========================================================
    # CHECK DIMENSIONS
    # ===========================================================
    @assert size(doppler_width) == size(damping_constant) == size(Bul)
                                == size(Blu) == size(Aul) == atmosphere_size
    @assert length(nλ) == 2nλ_bf + nλ_bb

    # ===========================================================
    # CHECK UNITS
    # ===========================================================
    @test all( dimension.(Aul) .==  Unitful.𝐓^-1 )
    @test all( dimension.(Bul) .==  Unitful.𝐓^2 * Unitful.𝐌^-1 * Unitful.𝐋)
    @test all( dimension.(Blu) .==  Unitful.𝐓^2 * Unitful.𝐌^-1 * Unitful.𝐋)
    @test all( dimension.(λ)   .==  Unitful.𝐋)
    @test all( dimension.(damping_constant)  .==  Unitful.𝐋^-2)
    @test all( dimension.(doppler_width)  .==  Unitful.𝐋)
    @test dimension(λ0)  ==  Unitful.𝐋

    # ===========================================================
    # NO NEGAITVE VALUES
    # ===========================================================
    @test all( ustrip.(Aul) .>= 0.0 )
    @test all( ustrip.(Bul) .>= 0.0 )
    @test all( ustrip.(Blu) .>= 0.0 )
    @test all( ustrip.(damping_constant) .>= 0.0 )
    @test all( ustrip.(doppler_width) .>= 0.0 )
    @test all( ustrip.(λ) .>= 0.0 )
    @test ustrip.(λ0) .>= 0.0
    @test nλ_bb >= 0
    @test nλ_bf >= 0
end

function check_populations(populations, atmosphere_size)

    # ===========================================================
    # CHECK DIMENSIONS
    # ===========================================================
    @assert size(populations) == atmosphere_size

    # ===========================================================
    # CHECK UNITS
    # ===========================================================
    @test all( dimension.(populations) .==  Unitful.𝐋^-3)

    # ===========================================================
    # NO NEGAITVE VALUES
    # ===========================================================
    @test all( ustrip.(populations) .>= 0.0 )
end

function plot_populations(populations, atmosphere.z)
end

function check_rates(rates, atmosphere_size)
    # ===========================================================
    # LOAD DATA
    # ===========================================================
    R12 = rates.R12
    R13 = rates.R13
    R23 = rates.R23
    R21 = rates.R21
    R31 = rates.R31
    R32 = rates.R32
    C12 = rates.C12
    C13 = rates.C13
    C23 = rates.C23
    C21 = rates.C21
    C31 = rates.C31
    C32 = rates.C32

    # ===========================================================
    # CHECK DIMENSIONS
    # ===========================================================
    @assert size(R12) == atmosphere_size
    @assert size(R13) == atmosphere_size
    @assert size(R23) == atmosphere_size
    @assert size(R21) == atmosphere_size
    @assert size(R31) == atmosphere_size
    @assert size(R32) == atmosphere_size
    @assert size(C12) == atmosphere_size
    @assert size(C13) == atmosphere_size
    @assert size(C23) == atmosphere_size
    @assert size(C21) == atmosphere_size
    @assert size(C31) == atmosphere_size
    @assert size(C32) == atmosphere_size

    # ===========================================================
    # CHECK UNITS
    # ===========================================================
    @test all( dimension.(R12) .==  Unitful.𝐓^-1)

    # ===========================================================
    # NO NEGAITVE VALUES
    # ===========================================================
    @test all( ustrip.(populations) .>= 0.0 )

end

function plot_rates(rates, atmosphere.z)
end

function check_radiation(radiation)
end

function plot_radiation(radiation, atmosphere.z, atom.λ)
end

#### USEFUL JUNK

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


function average_column(array)
      Statistics.mean(array, dims=[2,3])[:,1,1]
end

run()
