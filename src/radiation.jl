include("atmosphere.jl")
include("rates.jl")
include("atom.jl")

struct Radiation
    α_continuum::Array{<:PerLength, 4}                     # (nλ, nz, nx, ny)
    ε_continuum::Array{Float64,4}                        # (nλ, nz, nx, ny)
    α_line_constant::Array{Float64, 3}
    ε_line::Array{Float64,3}
    boundary::Array{Int32,3}                             # (nλ, nx, ny)
    packets::Array{Int32,4}                              # (nλ, nz, nx, ny)
    intensity_per_packet::Array{<:UnitsIntensity_λ, 1}     # (nλ)
end

struct RadiationBackground
    λ::Array{<:Unitful.Length, 1}                          # (nλ)
    α_continuum::Array{<:PerLength, 4}                     # (nλ, nz, nx, ny)
    ε_continuum::Array{Float64,4}                        # (nλ, nz, nx, ny)
    boundary::Array{Int32,3}                             # (nλ, nx, ny)
    packets::Array{Int32,4}                              # (nλ, nz, nx, ny)
    intensity_per_packet::Array{<:UnitsIntensity_λ, 1}     # (nλ)
end

"""
    collect_radiation_data(atmosphere::Atmosphere,
                           λ::Unitful.Length,
                           τ_max::Float64,
                           target_packets::Float64)

Collects radition data for background processes at a single wavelength
Returns data to go into structure.
"""
function collect_radiation_data(atmosphere::Atmosphere,
                                λ::Unitful.Length,
                                τ_max::Float64,
                                target_packets::Float64)

    # ==================================================================
    # GET ATMOSPHERE DATA AND WAVELENGTH
    # ==================================================================
    x = atmosphere.x
    y = atmosphere.y
    z = atmosphere.z
    temperature = atmosphere.temperature
    electron_density = atmosphere.electron_density
    hydrogen_populations = atmosphere.hydrogen_populations
    nz, nx, ny = size(temperature)

    λ = [λ]
    # ==================================================================
    # INITIALISE VARIABLES
    # ==================================================================
    α = Array{PerLength,4}(undef, 1, nz, nx, ny)
    ε = Array{Float64,4}(undef, 1, nz, nx, ny)
    boundary = Array{Int32,3}(undef, 1, nx, ny)
    packets = Array{Int32,4}(undef, 1, nz, nx, ny)
    intensity_per_packet = Array{UnitsIntensity_λ,1}(undef, 1)

    # ==================================================================
    # EXTINCTION AND DESTRUCTION PROBABILITY FOR BACKGROUND PROCESSES
    # ==================================================================
    proton_density = hydrogen_populations[:,:,:,end]
    hydrogen_ground_density = hydrogen_populations[:,:,:,1]
    hydrogen_neutral_density = hydrogen_populations[:,:,:,1] .+ hydrogen_populations[:,:,:,2]

    α_continuum_abs = α_cont_abs.(λ, temperature, electron_density, hydrogen_neutral_density, proton_density)
    α_continuum_scat = α_cont_scatt.(λ, electron_density, hydrogen_ground_density)

    α[1,:,:,:] = α_continuum_abs .+ α_continuum_scat
    ε[1,:,:,:] = α_continuum_abs ./ α[1,:,:,:]

    # ==================================================================
    # FIND OPTICAL DEPTH BOUNDARY
    # ==================================================================
    boundary[1,:,:] = optical_depth_boundary(α[1,:,:,:], z, τ_max)

    # ==================================================================
    # FIND DISTRIBUTION OF PACKETS
    # ==================================================================
    packets[1,:,:,:], intensity_per_packet[1] = distribute_packets(λ[1], target_packets, x, y, z,
                                                                   temperature, α_continuum_abs, boundary[1,:,:])

    # ==================================================================
    # CHECK FOR UNVALID VALUES
    # ==================================================================
    @test all(  Inf .> ustrip.(λ) .>= 0.0 )
    @test all(  Inf .> ustrip.(α) .>= 0.0 )
    @test all(  Inf .> ustrip.(intensity_per_packet) .>= 0.0 )
    @test all(  Inf .> ε .>= 0.0 )
    @test all(  Inf .> boundary .>= 0 )
    @test all(  Inf .> packets .>= 0 )

    return λ, α, ε, boundary, packets, intensity_per_packet
end

"""
    collect_radiation_data(atmosphere::Atmosphere,
                           atom::Atom,
                           rates::TransitionRates,
                           populations::Array{<:NumberDensity,4},
                           τ_max::Float64,
                           target_packets::Float64)

Collects radition data wavelength associated with bound-bound and
bound-free processes. Returns data to go into structure.
"""
function collect_radiation_data(atmosphere::Atmosphere,
                                atom::Atom,
                                rates::TransitionRates,
                                populations::Array{<:NumberDensity,4},
                                τ_max::Float64,
                                target_packets::Float64)

    # ==================================================================
    # GET ATMOSPHERE DATA
    # ==================================================================
    x = atmosphere.x
    y = atmosphere.y
    z = atmosphere.z
    temperature = atmosphere.temperature
    electron_density = atmosphere.electron_density
    hydrogen_populations = atmosphere.hydrogen_populations
    nz, nx, ny = size(temperature)
    velocity_z = atmosphere.velocity_z

    # ==================================================================
    # GET ATOM DATA
    # ==================================================================
    line = atom.line
    λ0 = line.λ0
    ΔλD = atom.doppler_width
    damping_constant = atom.damping_constant
    nλ_bf = atom.nλ_bf
    λ = atom.λ
    nλ = length(λ)

    # ==================================================================
    # INITIALISE VARIABLES
    # ==================================================================
    boundary = Array{Int32,3}(undef, nλ, nx, ny)
    packets = Array{Int32,4}(undef, nλ, nz, nx, ny)
    intensity_per_packet =  Array{UnitsIntensity_λ, 1}(undef, nλ)

    # ==================================================================
    # EXTINCTION AND DESTRUCTION PROBABILITY FOR EACH WAVELENGTH
    # ==================================================================
    α_continuum, ε_continuum = continuum_extinction_destruction(atmosphere, atom, rates, populations, λ)
    ε_line = line_destruction(rates)
    α_line_constant = line_extinction_constant.(Ref(line), populations[:,:,:,1], populations[:,:,:,2])

    # ==================================================================
    # FIND OPTICAL DEPTH BOUNDARY AND PACKET DISTRIBUTION FOR EACH λ
    # ==================================================================

    # BF wavelengths
    for l=1:2nλ_bf
        boundary[l,:,:] = optical_depth_boundary(α_continuum[l,:,:,:], z, τ_max)

        α_continuum_abs = α_continuum[l,:,:,:] .* ε_continuum[l,:,:,:]
        packets[l,:,:,:], intensity_per_packet[l] = distribute_packets(λ[l], target_packets, x, y, z,
                                                                      temperature, α_continuum_abs, boundary[l,:,:])
    end

    # BB wavelengths
    for l=2nλ_bf+1:nλ
        α = α_continuum[l,:,:,:] .+ line_extinction.(λ[l], λ0, ΔλD, damping_constant, α_line_constant, velocity_z)
        boundary[l,:,:] = optical_depth_boundary(α, z, τ_max)

        α_line = line_extinction.(λ[l], λ0, ΔλD, damping_constant, α_line_constant)
        α_abs =  α_continuum[l,:,:,:] .* ε_continuum[l,:,:,:] .+ α_line .* ε_line
        packets[l,:,:,:], intensity_per_packet[l] = distribute_packets(λ[l], target_packets, x, y, z,
                                                                       temperature, α_abs, boundary[l,:,:])
    end

    # ==================================================================
    # CHECK FOR UNVALID VALUES
    # ==================================================================
    @test all( Inf .> ustrip.(α_continuum) .>= 0.0 )
    @test all( Inf .> ε_continuum .>= 0.0 )
    @test all( Inf .> ε_line .>= 0.0 )
    @test all( Inf .> boundary .>= 0 )
    @test all( Inf .> ustrip.(α_line_constant) .>= 0.0 )
    @test all( Inf .> packets .>= 0 )
    @test all( Inf .> ustrip.(intensity_per_packet) .>= 0.0 )

    return α_continuum, ε_continuum, α_line_constant, ε_line, boundary, packets, intensity_per_packet
end



# ==================================================================
# EXTINCTION AND DESTRUCTION
# ==================================================================

"""
    continuum_extinction_destruction(atmosphere::Atmosphere,
                                     atom::Atom,
                                     rates::TransitionRates,
                                     atom_populations::Array{<:NumberDensity,4},
                                     λ::Array{<:Unitful.Length, 1})

Collect non-line extinction and destruction for all wavelengths.
Includes H bf, H- bf and ff, H2+ bf and ff, thomson and rayleigh.
"""

function continuum_extinction_destruction(atmosphere::Atmosphere,
                                          atom::Atom,
                                          rates::TransitionRates,
                                          atom_populations::Array{<:NumberDensity,4},
                                          λ::Array{<:Unitful.Length, 1})

    # ==================================================================
    # EXTINCTION AND DESTRUCTION PROBABILITY FROM BACKGROUND PROCESSES
    # ==================================================================

    temperature = atmosphere.temperature
    electron_density = atmosphere.electron_density
    hydrogen_populations = atmosphere.hydrogen_populations

    nλ_bf = atom.nλ_bf
    nλ_bb = atom.nλ_bb

    nz, nx, ny = size(temperature)
    nλ = length(λ)
    λ0 = λ[nλ_bf*2 + (nλ_bb÷2) + 1]

    # ==================================================================
    # EXTINCTION AND DESTRUCTION PROBABILITY FROM BACKGROUND PROCESSES
    # ==================================================================
    α_background = Array{PerLength, 4}(undef, 2nλ_bf, nz, nx, ny)
    ε_background = Array{Float64,4}(undef, 2nλ_bf, nz, nx, ny)
    α_continuum = Array{PerLength, 4}(undef, nλ, nz, nx, ny)
    ε_continuum = Array{Float64,4}(undef, nλ, nz, nx, ny)

    proton_density = hydrogen_populations[:,:,:,end]
    hydrogen_ground_density = hydrogen_populations[:,:,:,1]
    hydrogen_neutral_density = hydrogen_populations[:,:,:,1] .+ hydrogen_populations[:,:,:,2]

    # Background at bound-free wavelengths
    @Threads.threads for l=1:2*nλ_bf
        α_abs = α_cont_abs.(λ[l], temperature, electron_density, hydrogen_neutral_density, proton_density)
        α_scatt = α_cont_scatt.(λ[l], electron_density, hydrogen_ground_density)

        α_background[l,:,:,:] = α_scatt .+ α_abs
        ε_background[l,:,:,:] = α_abs ./ α_background[l,:,:,:]
    end

    # ==================================================================
    # BACKGROUND EXTINCTION AND DESTRUCTION FOR LINE
    # ==================================================================

    # Assume constant background over line profile wavelengths
    α_abs = α_cont_abs.(λ0, temperature, electron_density, hydrogen_neutral_density, proton_density)
    α_scatt =  α_cont_scatt.(λ0, electron_density, hydrogen_ground_density)
    α_background_line = α_abs .+ α_scatt
    ε_background_line = α_abs ./ α_background_line

    @Threads.threads for l=2*nλ_bf+1:nλ
        α_continuum[l,:,:,:] = α_background_line
        ε_continuum[l,:,:,:] = ε_background_line
    end

    # ==================================================================
    # EXTINCTION AND DESTRUCTION FROM ATOM BOUND-FREE
    # ==================================================================
    ν = c_0 ./ λ
    n_eff = sqrt(E_∞ / (atom.χu - atom.χl)) |> u"J/J"

    C31 = rates.C31
    R31 = rates.R31
    C32 = rates.C32
    R32 = rates.R32

    ε_bf_l = C31 ./ (R31 .+ C31)
    ε_bf_u = C32 ./ (R32 .+ C32)

    @Threads.threads for l=1:nλ_bf
        α_bf_l = hydrogenic_bf.(ν[l], ν[nλ_bf],
                               temperature,  atom_populations[:,:,:,1],
                               1.0, n_eff)

        α_bf_u = hydrogenic_bf.(ν[l+nλ_bf], ν[2*nλ_bf],
                               temperature, atom_populations[:,:,:,2],
                               1.0, n_eff)

        α_continuum[l,:,:,:] = α_background[l,:,:,:] .+ α_bf_l
        ε_continuum[l,:,:,:] = ( ε_background[l,:,:,:] .* α_background[l,:,:,:] .+ ε_bf_l .* α_bf_l ) ./  (α_background[l,:,:,:] .+ α_bf_l)

        α_continuum[l+nλ_bf,:,:,:] = α_background[l+nλ_bf,:,:,:] .+ α_bf_u
        ε_continuum[l+nλ_bf,:,:,:] = ( ε_background[l+nλ_bf,:,:,:] .* α_background[l+nλ_bf,:,:,:] .+ ε_bf_u .* α_bf_u ) ./ (α_background[l+nλ_bf,:,:,:] .+ α_bf_u)
    end

    return α_continuum, ε_continuum
end

"""
    line_extinction(λ::Unitful.Length,
                    λ0::Unitful.Length,
                    ΔλD::Unitful.Length,
                    damping_constant::PerArea,
                    α_line_constant::Float64,
                    v_los::Unitful.Velocity=0u"m/s")

Calculate line profile and return bound-bound
extinction contribution for a line wavelength.
"""
function line_extinction(λ::Unitful.Length,
                         λ0::Unitful.Length,
                         ΔλD::Unitful.Length,
                         damping_constant::PerArea,
                         α_line_constant::Float64,
                         v_los::Unitful.Velocity=0u"m/s")

    damping = damping_constant*λ^2 |> u"m/m"
    v = (λ - λ0 .+ λ0 .* v_los ./ c_0) ./ ΔλD
    profile = voigt_profile.(damping, ustrip(v), ΔλD)
    α = α_line_constant * profile

    return α
end

"""
    line_destruction(rates::TransitionRates)

Returns line destruction probability for the two level atom.
"""
function line_destruction(rates::TransitionRates)
    C21 = rates.C21
    R21 = rates.R21
    return C21 ./ (R21 .+ C21)
end

"""
    line_extinction_constant(line::AtomicLine, n_l::NumberDensity, n_u::NumberDensity)

Compute the line extinction constant to be
multiplied by the profile (per length).
"""
function line_extinction_constant(line::AtomicLine, n_l::NumberDensity, n_u::NumberDensity)
    (h * c_0 / (4 * π * line.λ0) * (n_l * line.Bij - n_u * line.Bji)) |> u"m/m"
end

"""
    α_cont_abs(λ::Unitful.Length,
               temperature::Unitful.Temperature,
               electron_density::NumberDensity,
               h_neutral_density::NumberDensity,
               proton_density::NumberDensity)

The extinction from continuum absorption processes for a given λ.
Includes H- ff, H- bf, H ff, H2+ ff and H2+ bf. Credit: Tiago
"""
function α_cont_abs(λ::Unitful.Length,
                    temperature::Unitful.Temperature,
                    electron_density::NumberDensity,
                    h_neutral_density::NumberDensity,
                    proton_density::NumberDensity)

    α = Transparency.hminus_ff_stilley(λ, temperature, h_neutral_density, electron_density)
    α += Transparency.hminus_bf_geltman(λ, temperature, h_neutral_density, electron_density)
    α += hydrogenic_ff(c_0 / λ, temperature, electron_density, proton_density, 1)
    α += h2plus_ff(λ, temperature, h_neutral_density, proton_density)
    α += h2plus_bf(λ, temperature, h_neutral_density, proton_density)
    return α
end

"""
    α_cont_scatt(λ::Unitful.Length,
                 electron_density::NumberDensity,
                 h_ground_density::NumberDensity)

The extincion from Thomson and Rayleigh scattering
for a given λ. Credit: Tiago
"""
function α_cont_scatt(λ::Unitful.Length,
                      electron_density::NumberDensity,
                      h_ground_density::NumberDensity)

    α = thomson(electron_density)
    α += rayleigh_h(λ, h_ground_density)
    return α
end

"""
    optical_depth(α::Array{<:PerLength, 3},
                  z::Array{<:Unitful.Length, 1})

Calculates the vertical optical depth of the atmosphere.
"""
function optical_depth(α::Array{<:PerLength, 3},
                       z::Array{<:Unitful.Length, 1})

    nz, nx, ny = size(α)
    τ = Array{Float64,3}(undef, nz-1, nx, ny)

    # Calculate vertical optical depth for each column
    Threads.@threads for col=1:nx*ny
        j = 1 + (col-1)÷nx
        i = col - (j-1)*nx
        τ[1,i,j] = 0.5(z[1] - z[2]) * (α[1,i,j] + α[2,i,j])

        for k=2:nz-1
            τ[k,i,j] =  τ[k-1,i,j] + 0.5(z[k] - z[k+1]) * (α[k,i,j] + α[k+1,i,j])
        end
    end

    return τ
end

"""
    optical_depth_boundary(α::Array{<:PerLength, 3},
                           z::Array{<:Unitful.Length, 1},
                           τ_max::Real)

Returns 2D array containing the z-indices
where the optical depth reaches τ_max.
"""
function optical_depth_boundary(α::Array{<:PerLength, 3},
                                z::Array{<:Unitful.Length, 1},
                                τ_max::Real)
    nz, nx, ny = size(α)
    columns = nx*ny
    boundary = Array{Int32, 2}(undef, nx, ny)

    # Calculate vertical optical depth for each column
    Threads.@threads for col=1:columns
        j = 1 + (col-1)÷nx
        i = col - (j-1)*nx

        τ = 0
        k = 1

        while τ < τ_max && k < nz
            # Trapezoidal rule
            τ += 0.5(z[k] - z[k+1]) * (α[k,i,j] + α[k+1,i,j])
            k += 1
        end
        boundary[i,j] = k
    end

    return boundary
end

"""
    distribute_packets(λ::Unitful.Length,
                       target_packets::Real,
                       x::Array{<:Unitful.Length, 1},
                       y::Array{<:Unitful.Length, 1},
                       z::Array{<:Unitful.Length, 1},
                       temperature::Array{<:Unitful.Temperature, 3},
                       α::Array{<:PerLength, 3},
                       boundary::Array{Int32,2})

Returns a 3D array of the # of packets to be
generated in each box above the boundary.
As well as the intensity contained in each packet.
"""
function distribute_packets(λ::Unitful.Length,
                            target_packets::Real,
                            x::Array{<:Unitful.Length, 1},
                            y::Array{<:Unitful.Length, 1},
                            z::Array{<:Unitful.Length, 1},
                            temperature::Array{<:Unitful.Temperature, 3},
                            α::Array{<:PerLength, 3},
                            boundary::Array{Int32,2})

    nz, nx, ny = size(α)

    emissivity = blackbody_lambda.(λ, temperature) .* α # u"kW / m^3 / sr / nm"
    box_emission = zeros(Float64,nz,nx,ny)u"kW / sr / nm"
    intensity_per_packet = 0.0u"kW / m^2 / sr / nm"

    Δz = (z[1:end-1] .- z[2:end])
    Δx = (x[2:end] .- x[1:end-1])
    Δy = (y[2:end] .- y[1:end-1])

    @Threads.threads for j=1:ny
        for i=1:nx
            for k=1:boundary[i,j]
                box_emission[k,i,j] = emissivity[k,i,j]*Δz[k]*Δx[i]*Δy[j]
                intensity_per_packet += emissivity[k,i,j]*Δz[k]
            end
        end
    end

    packets_per_box = Int.(round.( (box_emission/sum(box_emission)) * target_packets ))
    intensity_per_packet /= sum(packets_per_box)

    return packets_per_box, intensity_per_packet
end

"""
    blackbody_lambda(λ::Unitful.Length,
                     temperature::Unitful.Temperature)

Calculates the Blackbody (Planck) function per
wavelength, for a given wavelength and temperature.
Returns monochromatic intensity. Credit: Tiago
"""
function blackbody_lambda(λ::Unitful.Length,
                          temperature::Unitful.Temperature)
    B = (2h * c_0^2) / ( λ^5 * (exp((h * c_0 / k_B) / (λ * temperature)) - 1) ) |> u"kW / m^2 / sr / nm"
end

"""
    blackbody_lambda(λ::Array{<:Unitful.Length,1},
                     temperature::Unitful.Temperature)

Calculates the Blackbody (Planck) function per
wavelength, for an array of wavelengths and 3D temperature.
Returns monochromatic intensity.
"""
function blackbody_lambda(λ::Array{<:Unitful.Length,1},
                          temperature::Array{<:Unitful.Temperature,3})
    nλ = length(λ)
    nz, nx, ny = size(temperature)
    B = Array{UnitsIntensity_λ, 4}(undef, nλ, nz, nx, ny)

    for l=1:nλ
        B[l,:,:,:] = (2h * c_0^2) ./ ( λ[l]^5 * (exp.((h * c_0 / k_B) ./ (λ[l] * temperature)) .- 1) ) .|> u"kW / m^2 / sr / nm"
    end

    return B
end


# ==================================================================
#  WRITE TO FILE
# ==================================================================

"""
    write_to_file(radiation::RadiationBackground, output_path)

Write the relevant radition data to the output file.
"""
function write_to_file(radiation::RadiationBackground, output_path::String)
    h5open(output_path, "r+") do file
        file["packets"][:,:,:,:] = ustrip(radiation.packets)
        file["boundary"][:,:,:] = radiation.boundary
        file["intensity_per_packet"][:] = ustrip(radiation.intensity_per_packet)
    end
end

"""
    write_to_file(radiation::Radiation, iteration::Int64, output_path::String)

Write the relevant radition data to the output file.
"""
function write_to_file(radiation::Radiation, iteration::Int64, output_path::String)
    h5open(output_path, "r+") do file
        file["packets"][iteration,:,:,:,:] = ustrip(radiation.packets)
        file["boundary"][iteration,:,:,:] = radiation.boundary
        file["intensity_per_packet"][iteration,:] = ustrip(radiation.intensity_per_packet)
    end
end
