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
    α = Array{Float64,4}(undef, 1, nz, nx, ny)u"m^-1"
    ε = Array{Float64,4}(undef, 1, nz, nx, ny)
    boundary = Array{Int32,3}(undef, 1, nx, ny)
    packets = Array{Int32,4}(undef, 1, nz, nx, ny)
    intensity_per_packet = Array{Float64,1}(undef, 1)u"kW / m^2 / sr / nm"

    # ==================================================================
    # EXTINCTION AND DESTRUCTION PROBABILITY FOR BACKGROUND PROCESSES
    # ==================================================================
    proton_density = hydrogen_populations[:,:,:,end]
    hydrogen_ground_popuplation = hydrogen_populations[:,:,:,1]

    α_cont_a = α_cont_abs.(λ, temperature, electron_density, hydrogen_ground_popuplation, proton_density)
    α_cont_s = α_cont_scatt.(λ, electron_density, hydrogen_ground_popuplation)

    α[1,:,:,:] = α_cont_a .+ α_cont_s
    ε[1,:,:,:] = α_cont_a ./ α[1,:,:,:]

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
function collect_radiation_data(atmosphere::Atmosphere, atom::AtomicLine, populations)
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
    intensity_per_packet =  Array{Float64,1}(undef, nλ)u"kW / m^2 / sr / nm"

    # ==================================================================
    # EXTINCTION AND DESTRUCTION PROBABILITY FOR EACH WAVELENGTH
    # ==================================================================
    α, ε = α_and_ε_atom(atmosphere, atom, populations, λ)


    # ==================================================================
    # FIND OPTICAL DEPTH BOUNDARY AND PACKET DISTRIBUTION FOR EACH λ
    # ==================================================================
    for l=1:nλ
        boundary[l,:,:] = optical_depth_boundary(α[l,:,:,:], z, τ_max)
        packets[l,:,:,:], intensity_per_packet[l] = distribute_packets(λ[l], target_packets, x, y, z,
                                                                       temperature, α[l,:,:,:], boundary[l,:,:])
    end

    return λ, α, ε, boundary, packets, intensity_per_packet, max_scatterings
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



function α_and_ε_atom(atmosphere::Atmosphere, atom::AtomicLine, atom_populations, λ)

    vz = atmosphere,velocity_z
    vx = atmosphere.velocity_x
    vy = atmosphere.velocity_y
    temperature = atmosphere.temperature
    electron_density = atmosphere.electron_density
    hydrogen_populations = atmosphere.hydrogen_populations

    nλ_bb, nλ_bf = get_nλ()

    nz, nx, ny = size(temperature)
    nλ = length(λ)

    # For each wavelength, find χ and ε
    α = Array{Float64,4}(undef, nλ, nz, nx, ny)u"m^-1"
    ε = Array{Float64,4}(undef, nλ, nz, nx, ny)

    # Find bound-free continuum
    @Threads.threads for l=1:nλ_bf
        α[l,:,:,:], ε[l,:,:,:] =  α_and_ε_cont(λ[l], temperature, electron_density, hydrogen_populations)
        α[l+nλ_bf,:,:,:], ε[l+nλ_bf,:,:,:] =  α_and_ε_cont(λ[l+nλ_bf], temperature, electron_density, hydrogen_populations)
    end

    # Find bound-bound continuum
    # assume continuum constant over line
    center = nλ_bf*2 + (nλ_bb÷2)
    α_cont, ε_cont =  α_and_ε_cont(λ[center], temperature, electron_density, hydrogen_populations)

    # Compute line extinction (van der Waals + natural broadening)
    unsold_const = γ_unsold_const(atom)
    γ = γ_unsold.(unsold_const, temperature, hydrogen_populations[:,:,:,1]) .+ atom.Aji
    ΔλD = doppler_width.(λ[center], atom.atom_weight, temperature)

    Cji = Cji_RH()                          # replace with atom.Cij

    @Threads.threads for l=(nλ_bf*2+1):nλ

        a = damping.(γ, λ[l], ΔλD)
        v = (λ[l] - λ[center]) ./ ΔλD
        profile = voigt_profile.(a, ustrip(v), ΔλD)
        α_line = αline_λ.(Ref(atom), profile, atom_populations[:,:,:,2], atom_populations[:,:,:,1])


        B = blackbody_lambda.(λ[l], temperature)
        Rji = atom.Aji .+ atom.Bji.*B #should I use J in later iterations
        ε_line = Cji ./ (Rji .+ Cji)
        α[l,:,:,:] = α_line .+ α_cont
        ε[l,:,:,:] = ε_line .* (α_line ./ α[l,:,:,:])  .+ ε_cont .* (α_cont ./ α[l,:,:,:])
    end

    return α, ε
end



function α_and_ε_atom_test(atmosphere::Atmosphere, atom::AtomicLine, atom_populations, λ)

    vz = atmosphere,velocity_z
    vx = atmosphere.velocity_x
    vy = atmosphere.velocity_y
    temperature = atmosphere.temperature
    electron_density = atmosphere.electron_density
    hydrogen_populations = atmosphere.hydrogen_populations

    nλ_bb, nλ_bf = get_nλ()

    nz, nx, ny = size(temperature)
    nλ = length(λ)

    # For each wavelength, find χ and ε
    α_cont = Array{Float64,4}(undef, nλ, nz, nx, ny)u"m^-1"
    ε_cont = Array{Float64,4}(undef, nλ, nz, nx, ny)
    ε_line = Array{Float64,4}(undef, nλ, nz, nx, ny)

    # Find bound-free continuum
    @Threads.threads for l=1:nλ_bf
        α_cont[l,:,:,:], ε_cont[l,:,:,:] =  α_and_ε_cont(λ[l], temperature, electron_density, hydrogen_populations)
        α_cont[l+nλ_bf,:,:,:], ε_cont[l+nλ_bf,:,:,:] =  α_and_ε_cont(λ[l+nλ_bf], temperature, electron_density, hydrogen_populations)
    end

    # Find bound-bound continuum
    # assume continuum constant over line
    center = nλ_bf*2 + (nλ_bb÷2)
    α_cont_line, ε_cont_line =  α_and_ε_cont(λ[center], temperature, electron_density, hydrogen_populations)

    # Compute line extinction (van der Waals + natural broadening)
    unsold_const = γ_unsold_const(atom)
    γ = γ_unsold.(unsold_const, temperature, hydrogen_populations[:,:,:,1]) .+ atom.Aji
    ΔλD = doppler_width.(λ[center], atom.atom_weight, temperature)

    Cji = Cji_RH()                          # replace with atom.Cij

    @Threads.threads for l=(nλ_bf*2+1):nλ

        α_cont[l,:,:,:] = α_cont_line

        a = damping.(γ, λ[l], ΔλD)

        v = (λ[l] - λ[center]) ./ ΔλD
        profile = voigt_profile.(a, ustrip(v), ΔλD)

        αC = αline_λ.(Ref(atom), profile, atom_populations[:,:,:,2], atom_populations[:,:,:,1]) ./ profile

        B = blackbody_lambda.(λ[l], temperature)
        Rji = atom.Aji .+ atom.Bji.*B
        ε_line = Cji ./ (Rji .+ Cji)

    end

    return α_cont, ε_cont, ε_line, αC, ΔλD, a
end

function α_line(λ, αC, ΔλD, a, vlos)

    v = (λ[l] - λ[center]) ./ ΔλD
    profile = voigt_profile.(a, ustrip(v), ΔλD)

    α = αC*profile

    return α
end


"""
DELETE once Cji in Transparency
"""
function Cji_RH()
    rh_aux = h5open("/mn/stornext/u3/idarhan/MScProject/SolarMCRT/run/atoms/Cij_aux.h5", "r")
    Cji = read(rh_aux, "Cij")[:,:,:,4]u"s^-1"
    close(rh_aux)

    # original dimensions of data
    nz, nx, ny = size(Cji)

    # ===========================================================
    # FLIP AXES
    # ===========================================================

    Cji = Cji[:,:,end:-1:1]
    # ===========================================================
    # CUT AND SLICE ATMOSPHERE BY INDEX
    # ===========================================================

    ze, xe, ye = get_stop()
    zs, xs, ys = get_start()
    dz, dx, dy = get_step()

    # Cut z-direction from below

    if ze != nothing && ze < nz
        nz = ze
        Cji = Cji[1:nz,:,:]
    end

    # Cut  z-direction from up top
    if zs > 1
        nz = zs
        Cji = Cji[nz:end,:,:]
    end

    # Cut x-direction from right
    if xe != nothing && xe < nx
        nx = xe
        Cji = Cji[:,1:nx,:]
    end

    # Cut x-direction from the left
    if xs > 1
        nx = xs
        Cji = Cji[:,nx:end,:]
    end

    # Cut y-direction from right
    if ye != nothing && ye < ny
        ny = ye
        Cij = Cij[:,:,1:ny]
    end

    # Cut y-direction from the left
    if ys > 1
        ny = ys
        Cji = Cji[:,:,ny:end]
    end

    # Only keep every dz-th box in z-direction
    if dz > 1
        Cji = Cji[1:dz:end,:,:]
    end

    # Only keep every dx-th box in x-direction
    if dx > 1
        Cji = Cji[:,1:dx:end,:]
    end

    # Only keep every dy-th box in y-direction
    if dy > 1
        Cji = Cji[:,:,1:dy:end]
    end

    return Cji
end



""" DELETE
function α_and_ε_cont(λ, temperature, electron_density, hydrogen_populations)

    proton_density = hydrogen_populations[:,:,:,end]
    hydrogen_ground_popuplation = hydrogen_populations[:,:,:,1]

    # continuum
    α_cont_a = α_cont_abs.(λ, temperature, electron_density, hydrogen_ground_popuplation, proton_density)
    α_cont_s = α_cont_scatt.(λ, electron_density, hydrogen_ground_popuplation)

    α_cont = α_cont_a .+ α_cont_s
    ε_cont = α_cont_a ./ α_cont

    return α_cont, ε_cont
end"""

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
