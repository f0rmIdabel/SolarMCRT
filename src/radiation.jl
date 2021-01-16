include("atmosphere.jl")

struct Radiation
    λ::Array{<:Unitful.Length, 1}                  # (nλ)
    χ::Array{<:Unitful.Quantity, 4}                # (nλ, nz, nx, ny)
    ε::Array{Float64,4}                            # (nλ, nz, nx, ny)
    boundary::Array{Int32,3}                       # (nλ, nx, ny)
    S::Array{Int32,4}                              # (nλ, nz, nx, ny)
    rad_per_packet::Array{<:Unitful.Quantity, 1}   # (nλ)
    max_scatterings::Real                          # Int64
    escape_bins::Array{Int64,1}       # (nϕ, nθ)   # (2)
end



"""
Collects radition data to go into structure.
"""
function collect_radiation_data(atmosphere::Atmosphere, λ::Unitful.Length)
    # Read from input file
    max_scatterings = get_max_scatterings()
    escape_bins = get_escape_bins()

    target_packets = get_target_packets()
    τ_max = get_cut_off()

    # Get atmosphere data
    x = atmosphere.x
    y = atmosphere.y
    z = atmosphere.z
    temperature = atmosphere.temperature
    electron_density = atmosphere.electron_density
    hydrogen_populations = atmosphere.hydrogen_populations

    # Calculate χ and ε
    χ = Array{Float64,4}(undef, 1, nz-1, nx, ny)u"m^-1"
    ε = Array{Float64,4}(undef, 1, nz-1, nx, ny)
    χ[1,:,:,:], ε[1,:,:,:] =  χ_and_ε_cont(λ, temperature, electron_density, hydrogen_populations)

    # Find opticla depth boundary
    boundary = Array{Int32,3}(undef, 1, nx, ny)
    boundary[1,:,:] = optical_depth_boundary(χ, z, τ_max)

    # Calculate distribuion of packets
    S = Array{Int32,4}(undef, 1, nz, nx, ny)
    intensity_per_packet = Array{Float64,1}(undef, 1)
    S[1,:,:,:], intensity_per_packet[1] = distribute_packets(λ, target_packets, x, y, z,
                                                 temperature, χ, boundary)

    return λ, χ, ε, boundary, S, intensity_per_packet, max_scatterings, escape_bins
end

function collect_radiation_data(atomsphere::Atmosphere, atom::AtomicLine)

    # Read from input file
    max_scatterings = get_max_scatterings()
    escape_bins = get_escape_bins()
    target_packets = get_target_packets()
    τ_max = get_cut_off()
    nλ_bb, nλ_bf = get_nλ()

    # Get atmosphere data
    x = atmosphere.x
    y = atmosphere.y
    z = atmosphere.z
    temperature = atmosphere.temperature
    electron_density = atmosphere.electron_density
    hydrogen_populations = atmosphere.hydrogen_populations

    # Get atom data
    χl = atom.χl
    χu = atom.χu
    λ0 = atom.λ0

    # Sample wavelengths
    λ = get_λ(χl, χu, nλ_bb, nλ_bf)
    nλ = length(λ)

    # Get opacity and destruction probability
    # For each wavelength, find χ and ε
    χ, ε = χ_and_ε_atom(atom, λ, nλ_bb, nλ_bf, temperature, electron_density, hydrogen_populations)

    # Get boundary and packet distribuion

    for l=1:nλ
        # Find opticla depth boundary
        boundary[l,:] = optical_depth_boundary(χ[l,:,:,:], z, τ_max)

        # Calculate distribuion of packets
        S[l,:], intensity_per_packet[l] = distribute_packets(λ[l], target_packets, x, y, z,
                                                     temperature, χ[l,:,:,:], boundary[l,:,:])
    end

    return λ, χ, ε, boundary, S, intensity_per_packet, max_scatterings, escape_bins

end


function get_λ(χl, χu, nλ_bb, nλ_bf)

    λ_bf_edge_l = energytolambda...#convert(Quantity{T, Unitful.𝐋}, ((h * c_0) / (χu - χl)) |> u"nm")
    λ_bf_edge_u =
    λ_bb_center = energytolambda

    Δλ_bf = 1.0
    Δλ_bb = 0.1

    if nλ_bb % 2 == 0
        nλ_bb += 1.0
    end

    nλ = nλ_bf*2 + nλ_bb

    λ = Array{Float64,1}(undef,nλ)u"m"

    for l=1:nλ_bf
        λ[l] = λ_bf_edge_l + Δλ_bf*l
        λ[l+nλ_bf] = λ_bf_edge_u + Δλ_bf*l
    end

    center = nλ_bf*2 + (nλ_bb÷2)
    λ[center] = atom.λ0

    for l=1:(nλ_bb÷2)
        λ[center-l] = λ[center - l + 1] - Δλ_bb
        λ[center+l] = λ[center + l - 1] + Δλ_bb
    end

    return λ
end


function χ_and_ε_atom(atom, λ, nλ_bb, nλ_bf, temperature, electron_density, hydrogen_populations)
    nz, nx, ny = shape(temperature)
    nλ = length(λ)

    # For each wavelength, find χ and ε
    χ = Array{Float64,4}(undef, nλ, nz-1, nx, ny)u"m^-1"
    ε = Array{Float64,4}(undef, nλ, nz-1, nx, ny)

    # Find bound-free continuum
    for l=1:nλ_bf
        χ[l,:,;,:], ε[l,:,:,:] =  χ_and_ε_cont(λ[l], temperature, electron_density, hydrogen_populations)
        χ[l+nλ_bf,:,;,:], ε[l+nλ_bf,:,:,:] =  χ_and_ε_cont(λ[l], temperature, electron_density, hydrogen_populations)
    end

    # Find bound-bound continuum
    # assume continuum constant over line
    center = nλ_bf*2 + (nλ_bb÷2)
    χ_cont, ε_cont =  χ_and_ε_cont(λ[center], temperature, electron_density, hydrogen_populations)

    # Compute line extinction (van der Waals + natural broadening)
    unsold_const = γ_unsold_const(atom)
    γ = γ_unsold.(unsold_const, temperature, hydrogen_populations[:, 1]) .+ atom.Aji
    ΔλD = doppler_width.(λ0, atom.atom_weight, temperature)

    for l=(nλ_bf*2+1):nλ

        a = damping.(γ, λ[l], ΔλD)
        v = (λ[l] - λ0) ./ ΔλD
        profile = voigt_profile.(a, v, ΔλD)
        χ_line = αline_λ.(Ref(atom), profile, hydrogen_populations[:, 2], hydrogen_populations[:, 1])

        B = blackbody_lambda(λ[l], temperature)
        Rji = atom.Aji .+ atom.Bji.*B
        Cji = Cji_rh()                          # replace with atom.Cji
        ε_line = Cji ./ (Rji .+ Cji)
        χ[l,:,:,:] = χ_line .+ χ_cont
        ε[l,:,:,:] = ε_line .* (χ_line ./ χ[l,:,:,:])  .+ ε_cont .* (χ_cont ./ χ[l,:,:,:])
        end
    end

    return χ, ε
end


function χ_and_ε_cont(λ, temperature, electron_density, hydrogen_populations)

    proton_density = hydrogen_populations[:,:,;,end]
    hydrogen_ground_popuplation = hydrogen_populations[:,:,:,1] #unclear if I should use all neutral hydrogen

    # continuum
    χ_cont_a = χ_cont_abs.(λ, temperature, electron_density, hydrogen_ground_popuplation, proton_density)
    χ_cont_s = χ_cont_scatt.(λ, electron_density, hydrogen_ground_popuplation)

    χ_cont = χ_cont_a .+ χ_cont_s
    ε_cont = χ_cont_a ./ χ_cont

    return χ_cont, ε_cont
end


"""
DELETE once Cji in Transparency
"""
function Cji_rh()
    rh_aux = h5open("../../../../basement/MScProject/Atmospheres/output_aux.hdf5", "r")
    Cji = read(rh_aux, "atom_CA/Cji_line")[:,:,:,4]
    close(rh_aux)
    return Cji
end


"""
The extinction from continuum absorption processes for a given λ.
Includes H- ff, H- bf, H ff, H2+ ff and H2+ bf.
Credit: Tiago
"""
function χ_cont_abs(λ::Unitful.Length,
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
function χ_cont_scatt(λ::Unitful.Length,
                 electron_density::NumberDensity,
                 h_ground_density::NumberDensity)

    α = thomson(electron_density)
    α += rayleigh_h(λ, h_ground_density)
    return α
end



#############################################################################33

"""
Calculates the vertical optical depth of the atmosphere.
"""
function optical_depth(χ,
                       z::Array{<:Unitful.Length, 1})
    nz, nx, ny = size(χ)
    columns = nx*ny

    τ = Array{Float64,3}(undef, nz-1, nx, ny)

    # Calculate vertical optical depth for each column
    Threads.@threads for col=1:columns
        j = 1 + (col-1)÷nx
        i = col - (j-1)*nx
        τ[1,i,j] = 0.5(z[1] - z[2]) * (χ[1,i,j] + χ[2,i,j])

        for k=2:nz-1
            τ[k,i,j] =  τ[k-1,i,j] + 0.5(z[k] - z[k+1]) * (χ[k,i,j] + χ[k+1,i,j])
        end
    end

    return τ
end

"""
Returns 2D array containing the z-indices where the optical depth reaches τ_max.
"""
function optical_depth_boundary(χ::Array{<:Unitful.Quantity{<:Real, Unitful.𝐋^(-1)}, 3},
                                z::Array{<:Unitful.Length, 1},
                                τ_max::Real)
    nz, nx, ny = size(χ)
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
            τ += 0.5(z[k] - z[k+1]) * (χ[k,i,j] + χ[k+1,i,j])
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
                            x,
                            y,
                            z,
                            temperature,
                            χ,
                            boundary)

    nz, nx, ny = size(χ)

    j = blackbody_lambda.(λ, temperature) .* χ

    box_emission = zeros(Float64,nz,nx,ny)u"kW / sr / nm"
    intensity_per_packet = zeros(Float64,nz,nx,ny)u"kW / m^2 / sr / nm"

    Δz = (z[1:end-1] .- z[2:end])
    Δx = (x[2:end] .- x[1:end-1])
    Δy = (y[2:end] .- y[1:end-1])

    @Threads.threads for j=1:ny
        for i=1:nx
            for k=1:boundary[i,j]
                box_emission[k,i,j] = j[k,i,j]*Δz[k]*Δx[i]*Δy[j]
                intensity_per_packet = j[k,i,j]/(Δx[i]*Δy[j])
            end
        end
    end

    packets_per_box = Int.(round.( (box_emission/sum(box_emission)) * target_packets ))
    intensity_per_packet = sum(rad_per_packet)/sum(packets)

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
