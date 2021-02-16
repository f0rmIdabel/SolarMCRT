include("atmosphere.jl")
include("atom.jl")

struct Radiation
    λ::Array{<:Unitful.Length, 1}                        # (nλ)
    α_continuum::Array{<:PerLength, 4}                   # (nλ, nz, nx, ny)
    ε_continuum::Array{Float64,4}                        # (nλ, nz, nx, ny)
    ε_line::Array{Float64,4}
    boundary::Array{Int32,3}                             # (nλ, nx, ny)
    packets::Array{Int32,4}                              # (nλ, nz, nx, ny)
    intensity_per_packet::Array{<:UnitsIntensity_λ, 1}   # (nλ)
    max_scatterings::Int64                               # Int64
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
    ε_line = ε*0.0

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
function collect_radiation_data(atmosphere::Atmosphere, atom::AtomicLine, rates::TransitionRates, populations::Array{<:PerLength,4})
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
    velocity_z = atmosphere.velocity_z

    nz, nx, ny = size(temperature)

    # ==================================================================
    # SAMPLE WAVELENGTHS
    # ==================================================================
    λ = sample_λ(atom)
    nλ = length(λ)
    nλ_bf = atom.nλ_bf

    # ==================================================================
    # INITIALISE VARIABLES
    # ==================================================================
    boundary = Array{Int32,3}(undef, nλ, nx, ny)
    packets = Array{Int32,4}(undef, nλ, nz, nx, ny)
    intensity_per_packet =  Array{<:UnitsIntensity_λ, 1}(undef, nλ)  #u"kW / m^2 / sr / nm"

    # ==================================================================
    # EXTINCTION AND DESTRUCTION PROBABILITY FOR EACH WAVELENGTH
    # ==================================================================
    α_continuum, ε_continuum = continuum_extinction_destruction(atmosphere, atom, rates, populations, λ)
    ε_line = line_destruction(rates)

    # ==================================================================
    # FIND OPTICAL DEPTH BOUNDARY AND PACKET DISTRIBUTION FOR EACH λ
    # ==================================================================

    for l=1:2nλ_bf
        boundary[l,:,:] = optical_depth_boundary(α_continuum[l,:,:,:], z, τ_max)
        packets[l,:,:,:], intensity_per_packet[l] = distribute_packets(λ[l], target_packets, x, y, z,
                                                                      temperature, α_continuum[l,:,:,:], boundary[l,:,:])
    end

    for l=2nλ_bf+1:nλ
        v_los = velocity_z
        α = α_continuum[l,:,:,:] .+ line_extinction(λ[2nλ_bf+1:end], λ0, atom.ΔλD, atom.dc, atom.αlc, v_los)
        boundary[l,:,:] = optical_depth_boundary(α, z, τ_max)

        fill!(v_los, 0u"m/s")
        α = α_continuum[l,:,:,:] .+ line_extinction(λ[2nλ_bf+1:end], λ0, atom.ΔλD, atom.dc, atom.αlc, v_los)
        packets[l,:,:,:], intensity_per_packet[l] = distribute_packets(λ[l], target_packets, x, y, z,
                                                                      temperature, α, boundary[l,:,:])
    end

    return λ, α_continuum, ε_continuum, ε_line, boundary, packets, intensity_per_packet, max_scatterings
end

function sample_λ(atom)

    nλ_bb, nλ_bf = get_nλ()

    atom_file = h5open(get_atom_path(), "r")
    λ_bf_l_min = read(atom_file, "bfl_min")u"nm"
    λ_bf_u_min = read(atom_file, "bfu_min")u"nm"
    qwing = read(atom_file, "qwing")
    qcore = read(atom_file, "qcore")
    close(atom_file)

    λ_bf_l_max = transition_λ(atom.χl, atom.χ∞)
    λ_bf_u_max = transition_λ(atom.χu, atom.χ∞)
    λ_bb_center = transition_λ(atom.χl, atom.χu)

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

function continuum_extinction_destruction(atmosphere::Atmosphere, atom::AtomicLine, rates::TransitionRates, atom_populations::Array{<:PerLength,4}, λ::Array{<:Unitful.Length, 1})

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
    α_background = Array{<:PerLength, 4}(undef, nλ, nz, nx, ny)
    ε_background = Array{Float64,4}(undef, nλ, nz, nx, ny)
    α_continuum = Array{<:PerLength, 4}(undef, nλ, nz, nx, ny)
    ε_continuum = Array{Float64,4}(undef, nλ, nz, nx, ny)
    ε_line = Array{Float64,4}(undef, nλ_bb, nz, nx, ny)

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
    # EXTINCTION AND DESTRUCTION FROM ATOM BOUND-FREE
    # ==================================================================
    ν = c_0 ./ λ
    n_eff = sqrt(E∞ / (atom.χu - atom.χl))

    C31 = rates.C31
    R31 = rates.R31
    C32 = rates.C32
    R32 = rates.R32

    ε_bf_l = C31 ./ (R31.+ C31)
    ε_bf_u = C32 ./ (R32 .+ C32)

    @Threads.threads for l=1:nλ_bf
        α_bf_l = hydrogenic_bf(ν[l], ν[nλ_bf],
                               temperature,  hydrogen_populations[:,:,:,1],
                               1.0, n_eff)

        α_bf_u = hydrogenic_bf(ν[l+nλ_bf], ν[2*nλ_bf],
                               temperature, hydrogen_populations[:,:,:,2],
                               1.0, n_eff)

        α_continuum[l,:,:,:] = α_background[l,:,:,:] + α_bf_l
        ε_continuum[l,:,:,:] = ( ε_background[l,:,:,:] .* α_background[l,:,:,:] .+ ε_bf_l .* α_bf_l ) ./ α[l,:,:,:]

        α_continuum[l+nλ_bf,:,:,:] = α_background[l+nλ_bf,:,:,:] + α_bf_u
        ε_continuum[l+nλ_bf,:,:,:] = ( ε_background[l+nλ_bf,:,:,:] .* α_background[l+nλ_bf,:,:,:] .+ ε_bf_u .* α_bf_u ) ./ α[l+nλ_bf,:,:,:]
    end

    return α_continuum, ε_continuum
end

function line_extinction(λ::Unitful.Length, λ0::Unitful.Length, ΔλD::Float64, dc::Float64, αlc::Float64, v_los::Unitful.Velocity)
    damp = dc*λ^2 |> u"m/m"
    v = (λ - λ0 .+ λ0 .* v_los ./ c_0) ./ ΔλD
    profile = voigt_profile.(damp, ustrip(v), ΔλD)
    α = αlc * profile
    return α
end

function line_destruction(rates::TransitionRates)
    C21 = rates.C21
    R21 = rates.R21
    return C21 ./ (R21 .+ C21)
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

function transition_λ(χ1::Unitful.Energy, χ2::Unitful.Energy)
    ((h * c_0) / (χ2-χ1)) |> u"nm"
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
