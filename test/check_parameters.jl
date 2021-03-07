include("../src/mcrt.jl")
include("../src/populations.jl")
using Test

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
    nz, nx, ny = size(T)

    @assert size(T) == size(v) == size(electron_density) == size(vz) == size(v)
    @assert length(v[1,1,1]) == 3
    @assert size(hydrogen_populations) == (nz, nx, ny, 3)
    @assert length(z) == nz+1
    @assert length(x) == nx+1
    @assert length(y) == ny+1

    # ===========================================================
    # CHECK UNITS
    # ===========================================================
    @test dimension(z[1]) == Unitful.𝐋
    @test dimension(x[1]) == Unitful.𝐋
    @test dimension(y[1]) == Unitful.𝐋
    @test dimension(T[1]) == Unitful.𝚯
    @test dimension(v[1][1]) == Unitful.𝐋 * Unitful.𝐓^-1
    @test dimension(vz[1]) == Unitful.𝐋 * Unitful.𝐓^-1
    @test dimension(electron_density[1]) == Unitful.𝐋^-3
    @test dimension(hydrogen_populations[1]) == Unitful.𝐋^-3

    # ===========================================================
    # NO NEGAITVE OR INFINITE VALUES
    # ===========================================================
    @test all( Inf .> ustrip.(T) .>= 0.0 )
    @test all( Inf .> ustrip.(electron_density) .>= 0.0 )
    @test all( Inf .> ustrip.(hydrogen_populations) .>= 0.0 )

    # ===========================================================
    # DECREASING Z, INCREASING X AND Y
    # ===========================================================
    dz = z[2:end] .- z[1:end-1]
    dx = x[2:end] .- x[1:end-1]
    dy = y[2:end] .- y[1:end-1]

    @test all( ustrip.(dz) .<= 0.0 )
    @test all( ustrip.(dx) .>= 0.0 )
    @test all( ustrip.(dy) .>= 0.0 )
end

function check_radiationBackground(radiationBackground, atmosphere_size)
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
    nλ = length(λ)
    nz, nx, ny = atmosphere_size

    @assert size(α_continuum) == size(α_continuum) == size(packets)
    @assert size(boundary) == (nλ, nx, ny)
    @assert length(λ) == length(intensity_per_packet) == nλ

    # ===========================================================
    # CHECK UNITS
    # ===========================================================
    @test dimension(λ[1]) ==  Unitful.𝐋
    @test dimension(α_continuum[1]) == Unitful.𝐋^-1
    @test dimension(ε_continuum[1]) == NoDims
    @test dimension(boundary[1]) == NoDims
    @test dimension(packets[1]) == NoDims
    @test dimension(intensity_per_packet[1]) == Unitful.𝐋^-1 * Unitful.𝐌 * Unitful.𝐓^-3

    # ===========================================================
    # NO NEGAITVE OR INFINITE VALUES
    # ===========================================================
    @test all(  Inf .> ustrip.(λ) .>= 0.0 )
    @test all(  Inf .> ustrip.(α_continuum) .>= 0.0 )
    @test all(  Inf .> ustrip.(intensity_per_packet) .>= 0.0 )
    @test all(  Inf .> ε_continuum .>= 0.0 )
    @test all(  Inf .> boundary .>= 0 )
    @test all(  Inf .> packets .>= 0 )
end

function check_atom(atom, atmosphere_size)
    # ===========================================================
    # LOAD DATA
    # ===========================================================
    line = atom.line
    Aul = line.Aji
    Bul = line.Bji
    Blu = line.Bij
    λ0 = line.λ0
    χl = atom.χl
    χu = atom.χu
    χ∞ = atom.χ∞
    doppler_width = atom.doppler_width
    damping_constant = atom.damping_constant
    λ = atom.λ
    nλ_bb = atom.nλ_bb
    nλ_bf = atom.nλ_bf

    # ===========================================================
    # CHECK DIMENSIONS
    # ===========================================================
    @assert size(doppler_width) == atmosphere_size
    @assert size(damping_constant) == atmosphere_size
    @assert length(λ) == 2nλ_bf + nλ_bb

    # ===========================================================
    # CHECK UNITS
    # ===========================================================
    @test dimension(Aul) ==  Unitful.𝐓^-1
    @test dimension(Bul) ==  Unitful.𝐓^2 * Unitful.𝐌^-1 * Unitful.𝐋
    @test dimension(Blu) ==  Unitful.𝐓^2 * Unitful.𝐌^-1 * Unitful.𝐋
    @test dimension(λ[1])   ==  Unitful.𝐋
    @test dimension(damping_constant[1])  ==  Unitful.𝐋^-2
    @test dimension(doppler_width[1])  ==  Unitful.𝐋
    @test dimension(λ0)  ==  Unitful.𝐋
    @test dimension(χl) ==  Unitful.𝐌 * Unitful.𝐋^2 * Unitful.𝐓^-2
    @test dimension(χu) ==  Unitful.𝐌 * Unitful.𝐋^2 * Unitful.𝐓^-2
    @test dimension(χ∞) ==  Unitful.𝐌 * Unitful.𝐋^2 * Unitful.𝐓^-2

    # ===========================================================
    # NO NEGAITVE OR INFINITE VALUES
    # ===========================================================
    @test all( Inf > ustrip(Aul) >= 0.0 )
    @test all( Inf > ustrip(Bul) >= 0.0 )
    @test all( Inf > ustrip(Blu) >= 0.0 )
    @test all( Inf .> ustrip.(damping_constant) .>= 0.0 )
    @test all( Inf .> ustrip.(doppler_width) .>= 0.0 )
    @test all( Inf .> ustrip.(λ) .>= 0.0 )
    @test Inf > ustrip(λ0) >= 0.0
    @test Inf > nλ_bb >= 0
    @test Inf > nλ_bf >= 0

    @test ustrip(χl) == 0.0
    @test ustrip(χu) > 0.0
    @test χ∞ > χu
end

function check_populations(populations, atmosphere_size)
    # ===========================================================
    # CHECK DIMENSIONS
    # ===========================================================
    nz, nx, ny = atmosphere_size
    @assert size(populations) == (nz, nx, ny, 3)

    # ===========================================================
    # CHECK UNITS
    # ===========================================================
    @test dimension(populations[1]) ==  Unitful.𝐋^-3

    # ===========================================================
    # NO NEGAITVE VALUES
    # ===========================================================
    @test all( Inf .> ustrip.(populations) .> 0.0 ) # divide by zero problem
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
    @test dimension(R12[1]) ==  Unitful.𝐓^-1
    @test dimension(R13[1]) ==  Unitful.𝐓^-1
    @test dimension(R23[1]) ==  Unitful.𝐓^-1
    @test dimension(R21[1]) ==  Unitful.𝐓^-1
    @test dimension(R31[1]) ==  Unitful.𝐓^-1
    @test dimension(R32[1]) ==  Unitful.𝐓^-1
    @test dimension(C12[1]) ==  Unitful.𝐓^-1
    @test dimension(C13[1]) ==  Unitful.𝐓^-1
    @test dimension(C23[1]) ==  Unitful.𝐓^-1
    @test dimension(C21[1]) ==  Unitful.𝐓^-1
    @test dimension(C31[1]) ==  Unitful.𝐓^-1
    @test dimension(C32[1]) ==  Unitful.𝐓^-1

    # ===========================================================
    # NO NEGAITVE OR INFINITE VALUES
    # ===========================================================
    @test all( Inf .> ustrip.(R12) .>= 0.0 )
    @test all( Inf .> ustrip.(R13) .>= 0.0 )
    @test all( Inf .> ustrip.(R23) .>= 0.0 )
    @test all( Inf .> ustrip.(R21) .>= 0.0 )
    @test all( Inf .> ustrip.(R31) .>= 0.0 )
    @test all( Inf .> ustrip.(R32) .>= 0.0 )
    @test all( Inf .> ustrip.(C12) .>= 0.0 )
    @test all( Inf .> ustrip.(C13) .>= 0.0 )
    @test all( Inf .> ustrip.(C23) .>= 0.0 )
    @test all( Inf .> ustrip.(C21) .>= 0.0 )
    @test all( Inf .> ustrip.(C31) .>= 0.0 )
    @test all( Inf .> ustrip.(C32) .>= 0.0 )
end

function check_radiation(radiation, atom, atmosphere_size)
    # ===========================================================
    # LOAD DATA
    # ===========================================================
    α_continuum = radiation.α_continuum
    ε_continuum = radiation.ε_continuum
    α_line_constant = radiation.α_line_constant
    ε_line = radiation.ε_line

    boundary = radiation.boundary
    packets = radiation.packets
    intensity_per_packet = radiation.intensity_per_packet

    nz, nx, ny = atmosphere_size
    # ===========================================================
    # LOAD DATA
    # ===========================================================
    λ = atom.λ
    nλ = length(λ)

    nλ_bb = atom.nλ_bb
    nλ_bf = atom.nλ_bf
    α_line = Array{PerLength, 4}(undef,nλ_bb,nz,nx,ny)

    for l=1:nλ_bb
        α_line[l,:,:,:] = line_extinction.(λ[2nλ_bf + l], atom.line.λ0, atom.doppler_width, atom.damping_constant, α_line_constant)
    end

    # ===========================================================
    # CHECK DIMENSIONS
    # ===========================================================
    @assert size(α_continuum) == (nλ, nz, nx, ny)
    @assert size(α_continuum) ==  (nλ, nz, nx, ny)
    @assert size(ε_continuum) == (nλ, nz, nx, ny)
    @assert size(packets) == (nλ, nz, nx, ny)
    @assert size(α_line_constant) == (nz, nx, ny)
    @assert size(ε_line) == (nz, nx, ny)
    @assert size(boundary) == (nλ, nx, ny)
    @assert length(intensity_per_packet) == nλ

    # ===========================================================
    # CHECK UNITS
    # ===========================================================
    @test dimension(α_continuum[1]) == Unitful.𝐋^-1
    @test dimension(ε_continuum[1]) == NoDims
    @test dimension(α_line[1]) == Unitful.𝐋^-1
    @test dimension(α_line_constant[1]) == NoDims
    @test dimension(ε_line[1]) == NoDims
    @test dimension(boundary[1]) == NoDims
    @test dimension(packets[1]) == NoDims
    @test dimension(intensity_per_packet[1]) == Unitful.𝐋^-1 * Unitful.𝐌 * Unitful.𝐓^-3

    # ===========================================================
    # NO NEGAITVE VALUES
    # ===========================================================
    @test all( Inf .> ustrip.(α_continuum) .>= 0.0 )
    @test all( Inf .> ε_continuum .>= 0.0 )
    @test all( Inf .> ε_line .>= 0.0 )
    @test all( Inf .> boundary .>= 0 )
    @test all( Inf .> ustrip.(α_line_constant) .>= 0.0 )
    @test all( Inf .> ustrip.(α_line) .>= 0.0 )
    @test all( Inf .> packets .>= 0 )
    @test all( Inf .> ustrip.(intensity_per_packet) .>= 0.0 )
end
