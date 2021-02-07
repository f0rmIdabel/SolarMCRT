include("atmosphere.jl")
include("atom.jl")

struct Radiation
    λ::Array{<:Unitful.Length, 1}                        # (nλ)
    α::Array{<:PerLength, 4}                             # (nλ, nz, nx, ny)
    ε::Array{Float64,4}                                  # (nλ, nz, nx, ny)
    boundary::Array{Int32,3}                             # (nλ, nx, ny)
    packets::Array{Int32,4}                              # (nλ, nz, nx, ny)
    intensity_per_packet::Array{<:UnitsIntensity_λ, 1}   # (nλ)
    max_scatterings::Int64                               # Int64
end

struct RadiationAtom
    λ::Array{<:Unitful.Length, 1}                        # (nλ)
    α_continuum::Array{<:PerLength, 4}                   # (nλ, nz, nx, ny)
    ε::Array{Float64,4}                                  # (nλ, nz, nx, ny)
    boundary::Array{Int32,3}                             # (nλ, nx, ny)
    packets::Array{Int32,4}                              # (nλ, nz, nx, ny)
    intensity_per_packet::Array{<:UnitsIntensity_λ, 1}   # (nλ)
    max_scatterings::Int64                               # Int64
    a::Array
    ΔλD::Array
end


"""
TEST MODE: BACKGROUND PROCESSES
Collects radition data for background processes at a single wavelength
Returns data to go into structure.
"""
function collect_radiation_data(atmosphere::Atmosphere, λ::Unitful.Length)
    # ==================================================================
    # GET KEYWORD INPUT
    # ==================================================================
    max_scatterings = get_max_scatterings()
    target_packets = get_target_packets()
    τ_max = get_cut_off()

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
    α = Array{<:PerLength,4}(undef, 1, nz, nx, ny)
    ε = Array{Float64,4}(undef, 1, nz, nx, ny)
    boundary = Array{Int32,3}(undef, 1, nx, ny)
    packets = Array{Int32,4}(undef, 1, nz, nx, ny)
    intensity_per_packet = Array{<:UnitsIntensity_λ,1}(undef, 1)

    # ==================================================================
    # EXTINCTION AND DESTRUCTION PROBABILITY FOR BACKGROUND PROCESSES
    # ==================================================================
    proton_density = hydrogen_populations[:,:,:,end]
    hydrogen_ground_popuplation = hydrogen_populations[:,:,:,1]

    α_abs = α_cont_abs.(λ, temperature, electron_density, hydrogen_ground_popuplation, proton_density)
    α_scat = α_cont_scatt.(λ, electron_density, hydrogen_ground_popuplation)

    α[1,:,:,:] = α_abs .+ α_scat
    ε[1,:,:,:] = α_abs ./ α[1,:,:,:]

    # ==================================================================
    # FIND OPTICAL DEPTH BOUNDARY
    # ==================================================================
    boundary[1,:,:] = optical_depth_boundary(α[1,:,:,:], z, τ_max)

    # ==================================================================
    # FIND DISTRIBUTION OF PACKETS
    # ==================================================================
    packets[1,:,:,:], intensity_per_packet[1] = distribute_packets(λ[1], target_packets, x, y, z,
                                                                   temperature, α[1,:,:,:], boundary[1,:,:])

    return λ, α, ε, boundary, packets, intensity_per_packet, max_scatterings
end

"""
FULL MODE: POPULATION ITERATION
Collects radition data wavelength associated with bound-bound and bound-free processes.
Returns data to go into structure.
"""
function collect_radiation_data(atmosphere::Atmosphere, atom::AtomicLine, populations::Array{<:PerLength,4})
   # ==================================================================
   # GET KEYWORD INPUT
   # ==================================================================
    max_scatterings = get_max_scatterings()
    target_packets = get_target_packets()
    τ_max = get_cut_off()

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

    # ==================================================================
    # SAMPLE WAVELENGTHS
    # ==================================================================
    λ = sample_λ(atom)
    nλ = length(λ)

    # ==================================================================
    # INITIALISE VARIABLES
    # ==================================================================
    boundary = Array{Int32,3}(undef, nλ, nx, ny)
    packets = Array{Int32,4}(undef, nλ, nz, nx, ny)
    intensity_per_packet =  Array{<:UnitsIntensity_λ, 1}(undef, nλ)  #u"kW / m^2 / sr / nm"

    # ==================================================================
    # EXTINCTION AND DESTRUCTION PROBABILITY FOR EACH WAVELENGTH
    # ==================================================================
    α, α_line, ε, a, ΔλD  = atom_extinction_data(atmosphere, atom, populations, λ)

    # ==================================================================
    # FIND OPTICAL DEPTH BOUNDARY AND PACKET DISTRIBUTION FOR EACH λ
    # ==================================================================
    for l=1:nλ
        boundary[l,:,:] = optical_depth_boundary(α[l,:,:,:], z, τ_max)
        packets[l,:,:,:], intensity_per_packet[l] = distribute_packets(λ[l], target_packets, x, y, z,
                                                                       temperature, α[l,:,:,:], boundary[l,:,:])
    end

    α[2nλ_bf+1:end,:,:,:] -= α_line

    return λ, α_continuum, ε, boundary, packets, intensity_per_packet, max_scatterings, a, ΔλD
end

function sample_λ(atom::Atom)

    # Get atom data
    χl = atom.χl
    χu = atom.χu
    χ∞ = atom.χ∞
    nλ_bb = atom.nλ_bb
    nλ_bf = atom.nλ_bf
    bfl_min = atom.bfl_min
    bfu_min = atom.bfu_min
    λ_bf_l_min = atom.λ_bf_l_min
    λ_bf_l_max = atom.λ_bf_l_max
    λ_bf_u_min = atom.λ_bf_u_min
    λ_bf_u_max = atom.λ_bf_u_max
    λ_bb_center = atom.λ_bb_center
    qwing = atom.qwing
    qcore = atom.qcore

    # Make sure odd # of bb wavelengths
    if nλ_bb > 0 && nλ_bb%2 == 0
        nλ_bb += 1
    end

    # Initialise wavelength array
    nλ = nλ_bf*2 + nλ_bb
    λ = Array{Float64,1}(undef, nλ)u"nm"

    # =================================================
    # Bound-free transitions
    # Linear spacing
    # =================================================
    if nλ_bf > 0
        Δλ_bf_l = (λ_bf_l_max - λ_bf_l_min)/nλ_bf
        Δλ_bf_u = (λ_bf_u_max - λ_bf_u_min)/nλ_bf

        λ[1] = λ_bf_l_min
        λ[nλ_bf+1] = λ_bf_u_min

        for l=2:nλ_bf
            λ[l] = λ[l-1] + Δλ_bf_l
            λ[l+nλ_bf] = λ[l+nλ_bf - 1] + Δλ_bf_u
        end
    end

    # =================================================
    # Bound-bound transition
    # Follows github.com/ITA-Solar/rh/blob/master/getlambda.c
    # =================================================
    if nλ_bb > 0
        vmicro_char = 2.5u"km/s"

        n = nλ_bb/2 # Questionable
        β = qwing/(2*qcore)
        y = β + sqrt(β*β + (β - 1.0)*n + 2.0 - 3.0*β)
        b = 2.0*log(y) / (n - 1)
        a = qwing / (n - 2.0 + y*y)

        center = nλ_bf*2 + (nλ_bb÷2) + 1
        λ[center] = λ_bb_center
        q_to_λ = λ[center] * vmicro_char / c_0

        for l=1:(nλ_bb÷2)
            Δλ = a*(l + (exp(b*l) - 1.0)) * q_to_λ
            λ[center-l] = λ[center] - Δλ
            λ[center+l] = λ[center] + Δλ
        end
    end

    return λ
end

function atom_extinction_data(atmosphere::Atmosphere, atom::AtomicLine, atom_populations::Array{<:PerLength,4}, λ::Array{<:Unitful.Length, 1})

    temperature = atmosphere.temperature
    electron_density = atmosphere.electron_density
    hydrogen_populations = atmosphere.hydrogen_populations

    nλ_bb, nλ_bf = get_nλ()

    nz, nx, ny = size(temperature)
    nλ = length(λ)

    λ0 = λ[nλ_bf*2 + (nλ_bb÷2) + 1]

    # ==================================================================
    # EXTINCTION AND DESTRUCTION PROBABILITY FROM BACKGROUND PROCESSES
    # ==================================================================
    α_background = Array{<:PerLength, 4}(undef, nλ, nz, nx, ny)
    ε_background = Array{Float64,4}(undef, nλ, nz, nx, ny)
    α = Array{<:PerLength, 4}(undef, nλ, nz, nx, ny)
    ε = Array{Float64,4}(undef, nλ, nz, nx, ny)
    α_line = Array{<:PerLength, 4}(undef, nλ, nz, nx, ny)

    proton_density = hydrogen_populations[:,:,:,end]
    hydrogen_ground_popuplation = hydrogen_populations[:,:,:,1]

    # Background at bound-free wavelengths
    @Threads.threads for l=1:2*nλ_bf
        α_abs = α_cont_abs.(λ[l], temperature, electron_density, hydrogen_ground_popuplation, proton_density)
        α_scatt = α_cont_scatt.(λ[l], electron_density, hydrogen_ground_popuplation)

        α_background[l,:,:,:] = α_scatt + α_abs
        ε_background[l,:,:,:] = α_abs/ α_background
    end

    # Assume constant background over line profile wavelengths
    α_abs = α_cont_abs.(λ0, temperature, electron_density, hydrogen_ground_popuplation, proton_density)
    α_scatt =  α_cont_scatt.(λ0, electron_density, hydrogen_ground_popuplation)

    @Threads.threads for l=2*nλ_bf+1:nλ
        α_background[l,:,:,:] = α_scatt + α_abs
        ε_background[l,:,:,:] = α_abs/ α_background
    end

    # ==================================================================
    # EXTINCTION AND DESTRUCTION FROM BOUND-FREE PROCESSES
    # ==================================================================
    ν = c_0 ./ λ
    n_eff = sqrt(E∞ / (atom.χu - atom.χl))

    @Threads.threads for l=1:nλ_bf
        α_bf_l = hydrogenic_bf(ν[l], ν[nλ_bf],
                               temperature,  hydrogen_populations[:,:,:,1],
                               1.0, n_eff)

        α_bf_u = hydrogenic_bf(ν[l+nλ_bf], ν[2*nλ_bf],
                               temperature, hydrogen_populations[:,:,:,2],
                               1.0, n_eff)

        ## ε_bf_l = Cji_bf ./ (Rji_bf .+ Cji_bf)   Needs to be added to Transparency
        ## ε_bf_u = Cji_bf ./ (Rji_bf .+ Cji_bf)

        α[l,:,:,:] = α_background[l,:,:,:] + α_bf_l
        ε[l,:,:,:] = ( ε_background[l,:,:,:] .* α_background[l,:,:,:] .+ ε_bf_l .* α_bf_l ) ./ α[l,:,:,:]

        α[l+nλ_bf,:,:,:] = α_background[l+nλ_bf,:,:,:] + α_bf_u
        ε[l+nλ_bf,:,:,:] = ( ε_background[l+nλ_bf,:,:,:] .* α_background[l+nλ_bf,:,:,:] .+ ε_bf_u .* α_bf_u ) ./ α[l+nλ_bf,:,:,:]
    end

    # ==================================================================
    # EXTINCTION FROM BOUND-BOUND PROCESSES
    # ==================================================================

    unsold_const = γ_unsold_const(line)
    γ = γ_unsold.(unsold_const, temperature, hydrogen_populations[:,:,:,1]) .+ atom.Aji
    ΔλD = doppler_width.(λ[center], atom.atom_weight, temperature)
    a = Array{Float64, 4}(nλ, nz, nx, ny)

    @Threads.threads for l=(nλ_bf*2+1):nλ
        a[l,:,:,:] = damping.(γ, λ[l], ΔλD)
        v = (λ[l] - λ[center]) ./ ΔλD
        profile = voigt_profile.(a, ustrip(v), ΔλD)

        α_line[l,:,:,:] = αline_λ.(Ref(line), profile, atom_populations[:,:,:,2], atom_populations[:,:,:,1])
        #ε_line = Cji ./ (Rji[l] .+ Cji)

        α[l,:,:,:] = α_background[l,:,:,:] + α_line
        ε[l,:,:,:] = ( ε_background[l,:,:,:] .* α_background[l,:,:,:] .+ ε_line .* α_line) ./ α[l,:,:,:]
    end

    return α, ε, α_line, a, ΔλD
end

function α_line(line::AtomicLine, λ::Unitful.Length, λ0::Unitful.Length, ΔλD::Float64, a::Float64, vlos::Unitful.Velocity)

    v = (λ - λ0 + λ0*v_los/c_0 ) ./ ΔλD
    profile = voigt_profile.(a, ustrip(v), ΔλD)
    α_line = αline_λ.(Ref(line), profile, atom_populations[:,:,:,2], atom_populations[:,:,:,1])

    return α
end


"""
The extinction from continuum absorption processes for a given λ.
Includes H- ff, H- bf, H ff, H2+ ff and H2+ bf.
Credit: Tiago
"""
function α_cont_abs(λ::Unitful.Length,
               temperature::Unitful.Temperature,
               electron_density::NumberDensity,
               h_ground_density::NumberDensity,
               proton_density::NumberDensity)

    α = Transparency.hminus_ff_stilley(λ, temperature, h_ground_density, electron_density)
    α += Transparency.hminus_bf_geltman(λ, temperature, h_ground_density, electron_density)
    α += hydrogenic_ff(c_0 / λ, temperature, electron_density, proton_density, 1)
    α += h2plus_ff(λ, temperature, h_ground_density, proton_density)
    α += h2plus_bf(λ, temperature, h_ground_density, proton_density)
    return α
end

"""
The extincion from Thomson and Rayleigh scattering for a given λ.
Credit: Tiago
"""
function α_cont_scatt(λ::Unitful.Length,
                      electron_density::NumberDensity,
                      h_ground_density::NumberDensity)

    α = thomson(electron_density)
    α += rayleigh_h(λ, h_ground_density)
    return α
end

"""
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
        k = 0

        while τ < τ_max && k < nz
            k += 1
            # Trapezoidal rule
            τ += 0.5(z[k] - z[k+1]) * (α[k,i,j] + α[k+1,i,j])
        end
        boundary[i,j] = k
    end

    return boundary
end

"""
Returns a 3D array of the # of packets to be
generated in each box above the boundary.
As well as the scale
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
Calculates the Blackbody (Planck) function per wavelength,
for given arrays of wavelength and temperature.
Returns monochromatic intensity.
Credit: Tiago
"""
function blackbody_lambda(λ::Unitful.Length,
                          temperature::Unitful.Temperature)
    (2h * c_0^2) / ( λ^5 * (exp((h * c_0 / k_B) / (λ * temperature)) - 1) ) |> u"kW / m^2 / sr / nm"
end

function write_to_file(radiation::Radiation)
    h5open("../out/output.h5", "w") do file
        write(file, "lambda", ustrip(radiation.λ))
        write(file, "chi", ustrip(radiation.α))
        write(file, "epsilon", radiation.ε)
        write(file, "packets", ustrip(radiation.packets))
        write(file, "boundary", radiation.boundary)
        write(file, "intensity_per_packet", ustrip(radiation.intensity_per_packet))
    end
end
