include("../src/mcrt.jl")
include("../src/populations.jl")
using Test

"""
    check_atmosphere(atmosphere::Atmosphere)

Check that atmosphere data has valid units, dimensions and values.
"""
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

"""
    check_radiationBackground(radiationBackground::RadiationBackground,
                              atmosphere_size::Tuple)

Check that background radiation data has valid units, dimensions and values.
"""
function check_radiationContinuum(radiation::Radiation,
                                   atmosphere_size::Tuple, nλ::Int64)
    # ===========================================================
    # LOAD DATA
    # ===========================================================
    α_continuum = radiation.α_continuum
    ε_continuum = radiation.ε_continuum
    boundary = radiation.boundary
    packets = radiation.packets
    intensity_per_packet = radiation.intensity_per_packet

    # ===========================================================
    # CHECK DIMENSIONS
    # ===========================================================
    nz, nx, ny = atmosphere_size

    @assert size(α_continuum) == size(ε_continuum) == size(packets) == (nλ, nz, nx, ny)
    @assert size(boundary) == (nλ, nx, ny)
    @assert length(intensity_per_packet) == nλ

    # ===========================================================
    # CHECK UNITS
    # ===========================================================
    @test dimension(α_continuum[1]) == Unitful.𝐋^-1
    @test dimension(ε_continuum[1]) == NoDims
    @test dimension(boundary[1]) == NoDims
    @test dimension(packets[1]) == NoDims
    @test dimension(intensity_per_packet[1]) == Unitful.𝐋^-1 * Unitful.𝐌 * Unitful.𝐓^-3

    # ===========================================================
    # NO NEGAITVE OR INFINITE VALUES
    # ===========================================================
    @test all(  Inf .> ustrip.(α_continuum) .>= 0.0 )
    @test all(  Inf .> ustrip.(intensity_per_packet) .>= 0.0 )
    @test all(  1.0 .>= ε_continuum .>= 0.0 )
    @test all(  Inf .> boundary .>= 0 )
    @test all(  Inf .> packets .>= 0 )
end

"""
    check_atom(atom::Atom,
               atmosphere_size::Tuple)

Check that atom data has valid units, dimensions and values.
"""
function check_atom(atom::Atom,
                    atmosphere_size::Tuple)
    # ===========================================================
    # LOAD DATA
    # ===========================================================

    density = atom.density
    n_levels = atom.n_levels
    n_lines = atom.n_lines
    χ = atom.χ
    g = atom.g
    Z = atom.Z
    f_value = atom.f_value
    λ = atom.λ
    nλ = atom.nλ

    # ===========================================================
    # CHECK DIMENSIONS
    # ===========================================================
    @assert size(density) == atmosphere_size
    @assert length(χ) == n_levels + 1
    @assert length(g) == n_levels + 1
    @assert length(f_value) == n_lines

    # ===========================================================
    # CHECK UNITS
    # ===========================================================
    @test dimension(λ[1][1])    ==  Unitful.𝐋
    @test dimension(χ[1])       ==  Unitful.𝐌 * Unitful.𝐋^2 * Unitful.𝐓^-2
    @test dimension(density[1]) ==  Unitful.𝐋^-3
    # ===========================================================
    # NO NEGAITVE OR INFINITE VALUES
    # ===========================================================
    @test all( Inf .> ustrip.(density) .>= 0.0)
    @test all( Inf .> f_value .>= 0.0)
    @test all( Inf .> g .>= 0.0)

    @test ustrip(χ[1]) == 0.0
    @test all( ustrip.(χ[2:end] .- χ[1:end-1]) .> 0 )

    for l=1:(n_levels+n_lines)
        @test all(Inf .> ustrip.(λ[l]) .>= 0.0 )
    end

end

"""
    check_line(line::Line,
               atmosphere_size::Tuple)

Check that atom data has valid units, dimensions and values.
"""
function check_line(line::Line,
                    atmosphere_size::Tuple)
    # ===========================================================
    # LOAD DATA
    # ===========================================================
    lineData = line.lineData
    Aul = lineData.Aji
    Bul = lineData.Bji
    Blu = lineData.Bij
    λ0 = lineData.λ0
    doppler_width = line.doppler_width
    damping_constant = line.damping_constant

    # ===========================================================
    # CHECK DIMENSIONS
    # ===========================================================
    @assert size(doppler_width) == atmosphere_size
    @assert size(damping_constant) == atmosphere_size

    # ===========================================================
    # CHECK UNITS
    # ===========================================================
    @test dimension(Aul) ==  Unitful.𝐓^-1
    @test dimension(Bul) ==  Unitful.𝐓^2 * Unitful.𝐌^-1 * Unitful.𝐋
    @test dimension(Blu) ==  Unitful.𝐓^2 * Unitful.𝐌^-1 * Unitful.𝐋
    @test dimension(λ0)   ==  Unitful.𝐋
    @test dimension(damping_constant[1])  ==  Unitful.𝐋^-2
    @test dimension(doppler_width[1])  ==  Unitful.𝐋

    # ===========================================================
    # NO NEGAITVE OR INFINITE VALUES
    # ===========================================================
    @test all( Inf > ustrip(Aul) >= 0.0 )
    @test all( Inf > ustrip(Bul) >= 0.0 )
    @test all( Inf > ustrip(Blu) >= 0.0 )
    @test all( Inf .> ustrip.(damping_constant) .>= 0.0 )
    @test all( Inf .> ustrip.(doppler_width) .>= 0.0 )
    @test Inf > ustrip(λ0) >= 0.0
end

"""
    check_populations(populations::Array{<:NumberDensity,4},
                      atmosphere_size::Tuple)

Check that populations have valid units, dimensions and values.
"""
function check_populations(populations::Array{<:NumberDensity,4},
                           atmosphere_size::Tuple)
    # ===========================================================
    # CHECK DIMENSIONS
    # ===========================================================
    nz, nx, ny = atmosphere_size
    @assert size(populations[:,:,:,1]) == (nz, nx, ny)

    # ===========================================================
    # CHECK UNITS
    # ===========================================================
    @test dimension(populations[1]) ==  Unitful.𝐋^-3

    # ===========================================================
    # NO NEGAITVE VALUES
    # ===========================================================
    @test all( Inf .> ustrip.(populations) .> 0.0 )
end

"""
    check_rates(rates::TransitionRates,
                atmosphere_size::Tuple)

Check that transition rates have valid units, dimensions and values.
"""
function check_rates(rates::TransitionRates,
                     atmosphere_size::Tuple,
                     n_levels::Int64)
    # ===========================================================
    # LOAD DATA
    # ===========================================================
    R = rates.R
    C = rates.C

    # ===========================================================
    # CHECK DIMENSIONS
    # ===========================================================
    for l=1:n_levels
        for u=l+1:n_levels+1
            @assert size(R[l,u,:,:,:]) == atmosphere_size
            @assert size(R[u,l,:,:,:]) == atmosphere_size
            @assert size(C[l,u,:,:,:]) == atmosphere_size
            @assert size(C[u,l,:,:,:]) == atmosphere_size
        end
    end

    # ===========================================================
    # CHECK UNITS
    # ===========================================================
    for l=1:n_levels
        for u=l+1:n_levels+1
            @test dimension(R[l,u,1,1,1]) == Unitful.𝐓^-1
            @test dimension(R[u,l,1,1,1]) == Unitful.𝐓^-1
            @test dimension(C[l,u,1,1,1]) == Unitful.𝐓^-1
            @test dimension(C[u,l,1,1,1]) == Unitful.𝐓^-1
        end
    end

    # ===========================================================
    # NO NEGAITVE OR INFINITE VALUES
    # ===========================================================
    for l=1:n_levels
        for u=l+1:n_levels+1
            @test all( Inf .> ustrip.(R[l,u,:,:,:]) .>= 0.0 )
            @test all( Inf .> ustrip.(R[u,l,:,:,:]) .>= 0.0 )
            @test all( Inf .> ustrip.(C[l,u,:,:,:]) .>= 0.0 )
            @test all( Inf .> ustrip.(C[u,l,:,:,:]) .>= 0.0 )
        end
    end
end

"""
check_radiation(radiation::Radiation,
                atom::Atom,
                atmosphere_size::Tuple)

Check that radiation data has valid units, dimensions and values.
"""
function check_radiationLine(radiation::Radiation,
                             lineRadiation::LineRadiation,
                             λ::Array{<:Unitful.Length, 1},
                             line::Line,
                             atmosphere_size::Tuple)
    # ===========================================================
    # LOAD DATA
    # ===========================================================
    α_continuum = radiation.α_continuum
    ε_continuum = radiation.ε_continuum
    boundary = radiation.boundary
    packets = radiation.packets
    intensity_per_packet = radiation.intensity_per_packet

    α_line_constant = lineRadiation.α_line_constant
    ε_line = lineRadiation.ε_line

    nz, nx, ny = atmosphere_size
    # ===========================================================
    # LOAD DATA
    # ===========================================================
    nλ = length(λ)
    lineData = line.lineData

    α = Array{PerLength, 4}(undef,nλ,nz,nx,ny)
    ε = Array{Float64, 4}(undef,nλ,nz,nx,ny)

    for l=1:nλ
        α_line = line_extinction.(λ[l], lineData.λ0, line.doppler_width,
                                      line.damping_constant, α_line_constant)

        α[l,:,:,:] = α_line .+ α_continuum
        ε[l,:,:,:] = (ε_line .* α_line .+ ε_continuum .* α_continuum) ./ α[l,:,:,:]
    end

    # ===========================================================
    # CHECK DIMENSIONS
    # ===========================================================
    @assert size(α_continuum) == (nz, nx, ny)
    @assert size(ε_continuum) == (nz, nx, ny)
    @assert size(α_line_constant) == (nz, nx, ny)
    @assert size(ε_line) == (nz, nx, ny)

    @assert size(packets) == (nλ, nz, nx, ny)
    @assert size(boundary) == (nλ, nx, ny)
    @assert length(intensity_per_packet) == nλ

    # ===========================================================
    # CHECK UNITS
    # ===========================================================
    @test dimension(α_continuum[1]) == Unitful.𝐋^-1
    @test dimension(ε_continuum[1]) == NoDims
    @test dimension(α[1]) == Unitful.𝐋^-1
    @test dimension(α_line_constant[1]) == NoDims
    @test dimension(ε_line[1]) == NoDims
    @test dimension(boundary[1]) == NoDims
    @test dimension(packets[1]) == NoDims
    @test dimension(intensity_per_packet[1]) == Unitful.𝐋^-1 * Unitful.𝐌 * Unitful.𝐓^-3

    # ===========================================================
    # NO NEGAITVE VALUES
    # ===========================================================
    @test all( Inf .> ustrip.(α_continuum) .>= 0.0 )
    @test all( 1.0 .>= ε_continuum .>= 0.0 )
    @test all( 1.0 .>= ε_line .>= 0.0 )
    @test all( Inf .> ustrip.(α_line_constant) .>= 0.0 )
    @test all( Inf .> ustrip.(α) .>= 0.0 )
    @test all( 1.0 .>= ε .>= 0.0 )
    @test all( Inf .> boundary .>= 0 )
    @test all( Inf .> packets .>= 0 )
    @test all( Inf .> ustrip.(intensity_per_packet) .>= 0.0 )
end
