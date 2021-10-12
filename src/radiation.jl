include("atmosphere.jl")
include("rates.jl")
include("atom.jl")


struct Radiation
    α_continuum::Array{<:PerLength, 4}                     # (nλ, nz, nx, ny)
    ε_continuum::Array{Float64,4}                          # (nλ, nz, nx, ny)
    boundary::Array{Int32,3}                               # (nλ, nx, ny)
    packets::Array{Float64,4}                                # (nλ, nz, nx, ny)
    packets_to_intensity::Array{<:UnitsIntensity_λ, 4}     # (nλ)
end

struct LineRadiation
    α_line_constant::Array{Float64, 3}
    ε_line::Array{Float64,3}
end

function collect_radiation(atmosphere::Atmosphere, atom::Atom, rates::TransitionRates,
                           lines, lineRadiations, populations::Array{<:NumberDensity,4},
                           boundary_config, packet_config)
    x = atmosphere.x
    y = atmosphere.y
    z = atmosphere.z

    temperature = atmosphere.temperature
    nz,nx,ny = size(temperature)
    λ = atom.λ
    nλ = atom.nλ

    α_continuum, ε_continuum = opacity_continuum(atmosphere, atom, rates, populations)
    α_line, ε_line = opacity_line(atom,  lines, lineRadiations)

    α_abs = α_continuum .* ε_continuum .+ α_line .* ε_line

    α = α_continuum .+ α_line
    α_scat = α .- α_abs
    ε = α_abs ./ α


    boundary = Array{Int32,3}(undef, nλ, nx, ny)
    packets = Array{Float64,4}(undef, nλ, nz, nx, ny)
    packets_to_intensity = Array{UnitsIntensity_λ,4}(undef, nλ, nz, nx, ny)

    # ==================================================================
    # FIND OPTICAL DEPTH BOUNDARY AND PACKET DISTRIBUTION FOR EACH λ
    # ==================================================================

    if boundary_config[1] == false
        fill!(boundary, nz)
    else
        for w=1:nλ
            boundary[w,:,:] = depth_boundary(α[w,:,:,:], ε[w,:,:,:], z, boundary_config)
        end
    end

    for w=1:nλ
        packets[w,:,:,:], packets_to_intensity[w,:,:,:] = distribute_packets(λ[w], x, y, z, temperature, α_abs[w,:,:,:],
                                                                             ε[w,:,:,:], boundary[w,:,:], packet_config)
    end

    return α_continuum, ε_continuum, boundary, packets, packets_to_intensity
end


function collect_line_radiation_data(line::Line, rates::TransitionRates, populations::Array{<:NumberDensity,4})

    l = line.l
    u = line.u
    lineData = line.lineData

    α_line_constant = line_extinction_constant.(Ref(lineData), populations[:,:,:,l], populations[:,:,:,u])

    Cul = rates.C[u,l,:,:,:]
    Rul = rates.R[u,l,:,:,:]
    ε_line = Cul ./ (Rul .+ Cul)

    return α_line_constant, ε_line
end

function opacity_line(atom::Atom, lines, lineRadiations)
    # ==================================================================
    # GET ATOM AND LINE DATA
    # ==================================================================
    λ = atom.λ
    iλbb = atom.iλbb
    n_levels = atom.n_levels

    nz, nx, ny, = size(atom.density)
    nλ = atom.nλ

    # ==================================================================
    # INITIALISE VARIABLES
    # ==================================================================
    α = zeros(Float64, nλ, nz, nx, ny)u"m^-1"
    ε = zeros(Float64, nλ, nz, nx, ny)

    # ==================================================================
    # EXTINCTION AND DESTRUCTION PROBABILITY FOR BACKGROUND PROCESSES
    # ==================================================================

    for l=1:n_levels-1
        for u=2:n_levels

            line_number = sum((n_levels-l+1):(n_levels-1)) + (u - l)
            line = lines[line_number]
            lineRadiation = lineRadiations[line_number]
            ε_line = lineRadiation.ε_line
            start, stop = iλbb[line_number]

            for w=start:stop
                α_line = line_extinction.(λ[w], line.lineData.λ0, line.doppler_width, line.damping_constant, lineRadiation.α_line_constant)
                α[w,:,:,:] = α_line
                ε[w,:,:,:] = ε_line
            end
        end
    end

    # ==================================================================
    # CHECK FOR UNVALID VALUES
    # ==================================================================
    @test all(  Inf .> ustrip.(α) .>= 0.0 )
    @test all(  1.0 .>= ε .>= 0.0 )

    return α, ε
end

function opacity_continuum(atmosphere::Atmosphere, atom::Atom, rates::TransitionRates, populations::Array{<:NumberDensity,4})
    # ==================================================================
    # GET ATMOSPHERE DATA
    # ==================================================================
    temperature = atmosphere.temperature
    electron_density = atmosphere.electron_density
    hydrogen_populations = atmosphere.hydrogen_populations
    proton_density = hydrogen_populations[:,:,:,end]
    hydrogen_ground_density = hydrogen_populations[:,:,:,1]
    hydrogen_neutral_density = sum(hydrogen_populations, dims=4)[:,:,:,1] .- proton_density

    nz, nx, ny = size(temperature)
    # ==================================================================
    # GET ATOM AND RATE DATA
    # ==================================================================
    λ = atom.λ
    iλbf = atom.iλbf
    n_levels = atom.n_levels
    χ = atom.χ
    Z = atom.Z
    R = rates.R
    C = rates.C

    nλ = length(λ)
    # ==================================================================
    # INITIALISE VARIABLES
    # ==================================================================
    α = zeros(Float64, nλ, nz, nx, ny)u"m^-1"
    α_abs = zeros(Float64, nλ, nz, nx, ny)u"m^-1"

    # ==================================================================
    # EXTINCTION FROM BACKGROUND PROCESSES
    # ==================================================================

    for w=1:nλ
        α_absorp = α_cont_abs.(λ[w], temperature, electron_density, hydrogen_neutral_density, proton_density)
        α_scatt = α_cont_scatt.(λ[w], electron_density, hydrogen_ground_density)
        α[w,:,:,:] = α_scatt .+ α_absorp
        α_abs[w,:,:,:] = α_absorp
    end

    # ==================================================================
    # EXTINCTION AND DESTRUCTION PROBABILITY FROM BOUND-FREE PROCESSES
    # ==================================================================

    for level=1:n_levels

        n_eff = sqrt(E_∞ / (χ[end] - χ[level])) |> u"J/J"
        ε_bf = C[end,level,:,:,:] ./ (C[end,level,:,:,:] .+ R[end,level,:,:,:])

        start,stop = iλbf[level]
        for w=start:stop
            α_bf = hydrogenic_bf.(c_0/λ[w], c_0/λ[stop],
                                  temperature, populations[:,:,:,level],
                                  1.0, n_eff)

            α[w,:,:,:] .+= α_bf
            α_abs[w,:,:,:] .+= α_bf .* ε_bf
        end
    end

    ε = α_abs ./ α

    # ==================================================================
    # CHECK FOR UNVALID VALUES
    # ==================================================================
    @test all(  Inf .> ustrip.(α) .>= 0.0 )
    @test all(  1.0 .>= ε .>= 0.0 )

    return α, ε
end


"""
    collect_radiation_data(atmosphere::Atmosphere,
                           λ::Unitful.Length,
                           τ_max::Float64,
                           target_packets::Float64)

Collects radition data for background processes at a single wavelength.
This is what is called in the test mode.

Returns data to go into structure.
"""
function collect_background_radiation(atmosphere::Atmosphere,
                                      λ::Array{<:Unitful.Length,1},
                                      boundary_config,
                                      packet_config)

    # ==================================================================
    # GET ATMOSPHERE DATA AND WAVELENGTH
    # ==================================================================
    x = atmosphere.x
    y = atmosphere.y
    z = atmosphere.z
    volume = box_volume(z, x, y)
    temperature = atmosphere.temperature
    electron_density = atmosphere.electron_density
    hydrogen_populations = atmosphere.hydrogen_populations
    nz, nx, ny = size(temperature)

    # ==================================================================
    # INITIALISE VARIABLES
    # ==================================================================
    α = Array{PerLength,4}(undef, 1, nz, nx, ny)
    ε = Array{Float64,4}(undef, 1, nz, nx, ny)
    boundary = Array{Int32,3}(undef, 1, nx, ny)
    packets = Array{Float64,4}(undef, 1, nz, nx, ny)
    packets_to_intensity = Array{UnitsIntensity_λ,4}(undef, 1,nz, nx, ny)

    # ==================================================================
    # EXTINCTION AND DESTRUCTION PROBABILITY FOR BACKGROUND PROCESSES
    # ==================================================================
    proton_density = hydrogen_populations[:,:,:,end]
    hydrogen_ground_density = hydrogen_populations[:,:,:,1]
    hydrogen_neutral_density = sum(hydrogen_populations, dims=4)[:,:,:,1] .- proton_density

    α_continuum_abs = α_cont_abs.(λ, temperature, electron_density, hydrogen_neutral_density, proton_density)
    α_continuum_scat = α_cont_scatt.(λ, electron_density, hydrogen_ground_density)

    α[1,:,:,:] = α_continuum_abs .+ α_continuum_scat
    ε[1,:,:,:] = α_continuum_abs ./ α[1,:,:,:]

    # ==================================================================
    # FIND OPTICAL DEPTH BOUNDARY
    # ==================================================================
    if boundary_config[1] == false
        fill!(boundary, nz)
    else
        boundary[1,:,:] = depth_boundary(α[1,:,:,:], ε[1,:,:,:], z, boundary_config)
    end

    # ==================================================================
    # FIND DISTRIBUTION OF PACKETS
    # ==================================================================
    packets[1,:,:,:], packets_to_intensity[1,:,:,:] = distribute_packets(λ[1], x, y, z,temperature, α_continuum_abs,
                                                                         ε[1,:,:,:], boundary[1,:,:], packet_config)

    # ==================================================================
    # CHECK FOR UNVALID VALUES
    # ==================================================================
    @test all(  Inf .> ustrip.(λ) .>= 0.0 )
    @test all(  Inf .> ustrip.(α) .>= 0.0 )
    @test all(  Inf .> ustrip.(packets_to_intensity) .>= 0.0 )
    @test all(  1.0 .>= ε .>= 0.0 )
    @test all(  Inf .> boundary .>= 0 )
    @test all(  Inf .> packets .>= 0 )

    return α, ε, boundary, packets, packets_to_intensity
end


# ==================================================================
# EXTINCTION AND DESTRUCTION
# ==================================================================

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
    for col=1:nx*ny
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
function depth_boundary(α::Array{<:PerLength, 3},
                        ε::Array{Float64, 3},
                        z::Array{<:Unitful.Length, 1},
                        boundary_config,
                        depth_condition=false)
    nz, nx, ny = size(α)
    columns = nx*ny
    boundary = Array{Int32, 2}(undef, nx, ny)
    criterion, depth_exponent = boundary_config

    Δz = (z[1:end-1] .- z[2:end])*0.5

    # Calculate vertical optical depth for each column
    for col=1:columns
        j = 1 + (col-1)÷nx
        i = col - (j-1)*nx

        if depth_condition
            τ = 0
            k = 1

            α = α .* ε.^depth_exponent
            while τ < criterion && k < nz
                # Trapezoidal rule
                τ += 0.5(z[k] - z[k+1]) * (α[k,i,j] + α[k+1,i,j])
                k += 1
            end

            boundary[i,j] = k
        else
            c = 0
            k = nz

            while c < criterion && k > 1
                c = (1.0 .- ε[k,i,j]).^(Δz[k]*α[k,i,j]).^2
                k -= 1
            end

            boundary[i,j] = k
        end
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
                            x::Array{<:Unitful.Length, 1},
                            y::Array{<:Unitful.Length, 1},
                            z::Array{<:Unitful.Length, 1},
                            temperature::Array{<:Unitful.Temperature, 3},
                            α::Array{<:PerLength, 3},
                            ε::Array{Float64, 3},
                            boundary::Array{Int32,2},
                            packet_config)

    nz, nx, ny = size(α)
    Bλ = blackbody_lambda.(λ, temperature)
    emissivity = Bλ .* α # u"kW / m^3 / sr / nm"
    box_emission = zeros(Float64,nz,nx,ny)u"J / s / sr / nm"

    Δz = (z[1:end-1] .- z[2:end]) .|>u"m"
    Δx = (x[2:end] .- x[1:end-1]) .|>u"m"
    Δy = (y[2:end] .- y[1:end-1]) .|>u"m"

    for j=1:ny
        for i=1:nx
            for k=1:boundary[i,j]
                box_emission[k,i,j] = emissivity[k,i,j]*Δz[k]*Δx[i]*Δy[j]
            end
        end
    end

    target_packets, packet_exponent = packet_config

    total_emission = sum(box_emission)
    packets_per_box = round.(box_emission ./ total_emission * target_packets)
    scale = sum(packets_per_box .* ε.^packet_exponent) / sum(packets_per_box)
    packets_per_box = round.(box_emission ./ total_emission * target_packets * scale)

    total_packets = sum(packets_per_box)
    packets_to_intensity = total_emission / total_packets ./ box_surface_scale(Δz, Δx, Δy)

    # treat sub-boundary separately to get J = B
    for j=1:ny
        for i=1:nx
            for k=(boundary[i,j]+1):nz
                packets_per_box[k,i,j] = Bλ[k,i,j] / packets_to_intensity[k,i,j]
            end
        end
    end

    return packets_per_box, packets_to_intensity
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
    B = (2h * c_0^2) / ( λ^5 * (exp((h * c_0 / k_B) / (λ * temperature)) - 1) ) |> u"J / s / m^2 / sr / nm"
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
        B[l,:,:,:] = (2h * c_0^2) ./ ( λ[l]^5 * (exp.((h * c_0 / k_B) ./ (λ[l] * temperature)) .- 1) ) .|> u"J / s / m^2 / sr / nm"
    end

    return B
end

function box_volume(z::Array{<:Unitful.Length,1},
                    x::Array{<:Unitful.Length,1},
                    y::Array{<:Unitful.Length,1})

    nz = length(z) - 1
    nx = length(x) - 1
    ny = length(y) - 1
    volume = Array{Unitful.Volume, 3}(undef, nz, nx, ny)
    for k = 1:nz
        for i=1:nx
            for j=1:ny
                volume[k,i,j] = (z[k] - z[k+1]) * (x[i+1] - x[i]) * (y[i+1] - y[i])
            end
        end
    end
    return volume
end

function box_surface_scale(Δz, Δx, Δy)
    nz = length(Δz)
    nx = length(Δx)
    ny = length(Δy)
    surface_scale = Array{Unitful.Area, 3}(undef, nz, nx, ny)

    for k = 1:nz
        for i=1:nx
            for j=1:ny
                surface_scale[k,i,j] = Δz[k] * sqrt( Δy[j]^2 +  Δx[i]^2) * sqrt(Δz[i]/Δx[i])
            end
        end
    end
    return surface_scale
end

# ==================================================================
#  WRITE TO FILE
# ==================================================================

"""
    write_to_file(radiation::RadiationBackground, output_path)

Write the relevant radition data to the output file.
"""
function write_to_file(radiation::Radiation, output_path::String)
    h5open(output_path, "r+") do file
        file["packets"][1,:,:,:,:] = ustrip(radiation.packets)
        file["boundary"][1,:,:,:] = radiation.boundary
        file["packets_to_intensity"][1,:,:,:,:] = ustrip(radiation.packets_to_intensity)
    end
end

"""
    write_to_file(radiation::Radiation, iteration::Int64, output_path::String)

Write the relevant radition data to the output file.
"""
function write_to_file(radiation, iteration::Int64, output_path::String)
    h5open(output_path, "r+") do file
        file["packets"][iteration,:,:,:,:] = ustrip(radiation.packets)
        file["boundary"][iteration,:,:,:] = radiation.boundary
        file["packets_to_intensity"][iteration,:,:,:,:] = ustrip(radiation.packets_to_intensity)
    end
end
