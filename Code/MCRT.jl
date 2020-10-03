include("MyLibs/physLib.jl")
include("Atmos.jl")

using Random
using ProgressBars
using Unitful
#using LinearAlgebra # Where do I use this?

"""
    simulate(Atmosphere::Module, max_scatterings::Real,
             scale_emission::Real, chi_max::Real,  rng::MersenneTwister)

Simulates the radiation field in a given atmosphere with
a lower optical depth boundary given by tau_max.
"""
function simulate(atmosphere::Atmosphere, wavelength::Unitful.Length,
                  max_scatterings::Real, τ_max::Real,
                  total_packets::Real, num_bins = [4,2])

    # Atmosphere data
    x = atmosphere.x
    y = atmosphere.y
    z = atmosphere.z

    temperature = atmosphere.temperature
    ε = atmosphere.ε_continuum
    χ = atmosphere.χ_continuum

    λ = wavelength

    # Useful quantities
    nx = length(x)
    ny = length(y)
    nz = length(z)
    total_boxes = nx*ny*nz

    edge = [x[1] x[end]
            y[1] y[end]]

    # Initialise variables
    total_destroyed = Threads.Atomic{Int64}(0) #Remove Base. ?
    total_escaped = Threads.Atomic{Int64}(0)
    total_scatterings = Threads.Atomic{Int64}(0)

    boundary = optical_depth_boundary(χ, z, τ_max)
    total_emission_ = total_emission(χ, temperature, x, y, z, boundary, λ)
    scale_emission = total_packets/total_emission_

    # Bin escape angles
    ϕ_bins, θ_bins = num_bins

    surface = zeros(Int64, nx, ny, ϕ_bins, θ_bins)
    J = zeros(Int64, nz, nx, ny)

    nynz = ny*nz

    # Go through all boxes
    Threads.@threads for box in ProgressBar(1:total_boxes)

        # Find (x,y,z) indices of box
        i = 1 + box ÷ (nynz + 1)
        j = 1 + (box - (i-1)*nynz) ÷ (nz + 1)
        k = 1 + (box - (i-1)*nynz - 1) % (nz)

        # Skip boxes beneath boundary
        if k > boundary[i,j]
            continue
        end

        box_id = [k, i, j]

        # Find dimensions of box
        corner = [x[i], y[j], z[k]]
        box_dim = [x[i+1], y[j+1], z[k+1]] .- corner
        box_volume = box_dim[1]*box_dim[2]*(-box_dim[3])


        # Based on condition in box,
        # create certain number of photons packets
        B = blackbody_lambda(λ, temperature[box_id...])
        box_emission = B*χ[box_id...]*box_volume
        packets = Int(round(box_emission*scale_emission))

        # Add to local field
        J[k,i,j] += packets

        for packet=1:packets

            # Initial position uniformely drawn from box
            r = corner .+ (box_dim .* rand(3))

            # Initial box
            box_id = [k, i, j]

            # Scatter each packet until destroyed,
            # escape or reach max_scatterings
            for s=1:Int(max_scatterings)

                Threads.atomic_add!(total_scatterings, 1)

                # Scatter packet once
                r, box_id, escape, escape_angle, destroyed, J = scatter_packet!(x, y, z, χ, edge, box_id, r, J, boundary)

                # Check if escaped
                if escape
                    ϕ, θ  = escape_angle
                    ϕ_bin = 1 + Int(ϕ÷(2π/ϕ_bins))
                    θ_bin = 1 + sum(θ .> ((π/2)/(2 .^(2:θ_bins))))
                    surface[box_id[1], box_id[2], ϕ_bin, θ_bin] += 1
                    Threads.atomic_add!(total_escaped, 1)
                    break
                # Check if destroyed in bottom
                elseif destroyed
                    Threads.atomic_add!(total_destroyed, 1)
                    break
                # Check if destroyed in next particle interaction
                elseif rand() < ε[box_id...]
                    Threads.atomic_add!(total_destroyed, 1)
                    break
                end
            end
        end
    end

    return surface, total_destroyed, total_escaped, total_scatterings, J
end



"""
    scatter(Atmos::Module, box_id::Array{Int,1}, r::Array{<:Unitful.Length, 1}, J::Array{Int, 3})

Scatters photon packet once. Returns new position, box_id,
escape/destroyed-status and an updated mean radiation field J.
"""
function scatter_packet!(x::Array{<:Unitful.Length, 1}, y::Array{<:Unitful.Length, 1}, z::Array{<:Unitful.Length, 1},
                         χ::Array{<:Unitful.Quantity, 3}, edge::Array{<:Unitful.Length, 2}, box_id::Array{Int,1},
                         r::Array{<:Unitful.Length, 1}, J::Array{Int, 3}, boundary::Array{Int, 2})

    # Import relevant atmosphere data
    dim = size(χ) #Atmosphere.dim

    # Keep track of status
    escape = false
    escape_angle = nothing
    destroyed = false

    # Draw scattering depth and direction
    τ = -log(rand())
    ϕ = 2π * rand()
    θ =  π * rand()

    unit_vector = [cos(θ), sin(θ)*cos(ϕ), sin(θ)*sin(ϕ)]

    # Find distance to closest face
    ds, face = next_edge(x, y, z, r, unit_vector, box_id)

    direction = sign.(unit_vector)

    # Add optical depth and update position
    τ_cum = ds * χ[box_id...]
    r += ds*unit_vector

    # If depth target not reached in current box,
    # traverse boxes until target is reached
    while τ > τ_cum

        # Switch to new box                    # Find bins
        if face == 1
            box_id[1] -= direction[1] # Consequence of height array up->down

            # Top escape
            if box_id[1] == 0
                escape = true
                escape_angle = [ϕ, θ]
                break

            # Bottom destruction, Here they will very likely get destroyed anyway, so might not need this
        elseif box_id[1] == boundary[box_id[2], box_id[3]] + 1
                destroyed = true
                break
            end

        else
            box_id[face] += direction[face]

            # Handle side escapes with periodic boundary
            if box_id[face] == 0
                box_id[face] = dim[face]
                r[face] = edge[face,2]

            elseif  box_id[face] == dim[face] + 1
                box_id[face] = 1
                r[face] = edge[face,1]
            end

        end

        # Add to radiation field
        J[box_id...]+= 1

        # Find distance to closest face
        ds, face = next_edge(x, y, z, r, unit_vector, box_id)

        # Update optical depth and position
        τ_cum += ds*χ[box_id...]
        r += ds*unit_vector
    end

    if escape || destroyed
        r = nothing
    else
        # Correct for overshoot in final box
        r -= unit_vector * (τ_cum - τ)/χ[box_id...]
    end

    return r, box_id, escape, escape_angle, destroyed, J
end



"""
    next_edge(Atmosphere::Module, r::Array{Unitful.Length, 1},
              vector::Array{Float64, 1}, box_id::Array{Int,1})

Calculates the distance to the next box and the face that it crosses,
given an initial position and direction of travel.
"""
function next_edge(x::Array{<:Unitful.Length, 1}, y::Array{<:Unitful.Length, 1}, z::Array{<:Unitful.Length, 1},
                   r::Array{<:Unitful.Length, 1}, unit_vector::Array{Float64, 1}, box_id::Array{Int,1})

    direction = unit_vector .> 0.0

    # Find distance to box crossings in all dimensions
    distance = ([x[box_id[1] + direction[1]],
                y[box_id[2] + direction[2]],
                z[box_id[3] + !direction[3]]] .- r) ./unit_vector

    # Closest face
    face = argmin(distance)
    ds = distance[face]

    return ds, face
end
