include("atmos.jl")
include("MyLibs/physLib.jl")
using Random
using ProgressBars

"""
    simulate(Atmosphere::Module, max_scatterings::Real,
             scale_emission::Real, chi_max::Real,  rng::MersenneTwister)

Simulates the radiation field in a given atmosphere with
a lower optical depth boundary given by tau_max.
"""
function simulate(atmosphere::Atmosphere, wavelength::Unitful.Length,
                  max_scatterings::Real, τ_max::Real,
                  target_packets::Real, num_bins = [4,2])

    println("Loading atmosphere...")
    # Atmosphere data
    x = atmosphere.x
    y = atmosphere.y
    z = atmosphere.z
    ε = atmosphere.ε_continuum
    χ = atmosphere.χ_continuum
    temperature = atmosphere.temperature

    # Chosen wavelengths
    λ = wavelength

    println("Doing pre-calculations...")
    # Pre-calculations
    boundary = optical_depth_boundary(χ, z, τ_max)
    box_packets = packets_per_box(x,y,z,χ,temperature,
                                  λ,target_packets,boundary)

    # Number of boxes
    nx, ny = size(boundary)
    nz = maximum(boundary)
    total_boxes = nx*ny*nz

    # Number of escape bins
    ϕ_bins, θ_bins = num_bins

    # Initialise variables
    surface = zeros(Int64, nx, ny, ϕ_bins, θ_bins)
    J = zeros(Int64, nx, ny, nz)

    total_packets = Threads.Atomic{Int64}(0)
    total_destroyed = Threads.Atomic{Int64}(0)
    total_escaped = Threads.Atomic{Int64}(0)
    total_scatterings = Threads.Atomic{Int64}(0)

    # Go through all boxes
    Threads.@threads for box in ProgressBar(1:total_boxes)

        # Find (x,y,z) indices of box
        i = 1 + box ÷ (ny*nz + 1)
        j = 1 + (box - (i-1)*ny*nz) ÷ (nz + 1)
        k = 1 + (box - (i-1)*ny*nz - 1) % (nz)

        # Skip boxes beneath boundary
        if k > boundary[i,j]
            continue
        end

        box_id = [i,j,k]

        # Find dimensions of box
        corner = [x[i], y[j], z[k]]
        box_dim = [x[i+1], y[j+1], z[k+1]] .- corner

        # Add to local field
        packets = box_packets[box_id...]
        J[box_id...] += packets
        Threads.atomic_add!(total_packets, packets)

        for packet=1:packets

            # Initial position uniformely drawn from box
            r = corner .+ (box_dim .* rand(3))

            # Initial box
            box_id = [i,j,k]

            # Scatter each packet until destroyed,
            # escape or reach max_scatterings
            for s=1:Int(max_scatterings)

                Threads.atomic_add!(total_scatterings, 1)

                # Scatter packet once
                box_id, r, escaped, destroyed = scatter_packet(x, y, z, χ, boundary,
                                                               box_id, r, J)

                # Check if escaped
                if escaped[1]
                    ϕ, θ  = escaped[2]
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
function scatter_packet(x::Array{<:Unitful.Length, 1}, y::Array{<:Unitful.Length, 1}, z::Array{<:Unitful.Length, 1},
                        χ::Array{<:Unitful.Quantity{<:Real, Unitful.𝐋^(-1)}, 3}, boundary::Array{Int, 2},
                        box_id::Array{Int,1}, r::Array{<:Unitful.Length, 1}, J::Array{Int, 3})

    # Usefule quantities
    side_dim = size(boundary)
    edge = [x[1] x[end]
            y[1] y[end]]

    # Keep track of status
    escaped = false
    escape_angle = nothing
    destroyed = false

    # Draw scattering depth and direction
    τ = -log(rand())
    ϕ = 2π * rand()
    θ =  π * rand()

    unit_vector = [cos(θ), sin(θ)*cos(ϕ), sin(θ)*sin(ϕ)]

    # Find distance to closest face
    face, ds = next_edge(x, y, z, box_id, r, unit_vector)

    direction = sign.(unit_vector)

    # Add optical depth and update position
    τ_cum = ds * χ[box_id...]
    r += ds * unit_vector

    # If depth target not reached in current box,
    # traverse boxes until target is reached
    while τ > τ_cum

        # Switch to new box
        if face == 3
            box_id[3] -= direction[3] # Consequence of height array up->down

            # Top escape
            if box_id[3] == 0
                escaped = [true, [ϕ, θ]]
                break

            # Bottom destruction
            elseif box_id[3] == boundary[box_id[1], box_id[2]] + 1
                destroyed = true
                break
            end

        else
            box_id[face] += direction[face]

            # Handle side escapes with periodic boundary
            if box_id[face] == 0
                box_id[face] = side_dim[face]
                r[face] = edge[face,2]

            elseif  box_id[face] == side_dim[face] + 1
                box_id[face] = 1
                r[face] = edge[face,1]
            end

        end

        # Add to radiation field
        J[box_id...] += 1

        # Find distance to closest face
        face, ds = next_edge(x, y, z, box_id, r, unit_vector)

        # Update optical depth and position
        τ_cum += ds*χ[box_id...]
        r += ds*unit_vector
    end

    if escaped[1] || destroyed
        r = nothing
    else
        # Correct for overshoot in final box
        r -= unit_vector*(τ_cum - τ)/χ[box_id...]
    end

    return box_id, r, escaped, destroyed
end



"""
    next_edge(Atmosphere::Module, r::Array{Unitful.Length, 1},
              vector::Array{Float64, 1}, box_id::Array{Int,1})

Calculates the distance to the next box and the face that it crosses,
given an initial position and direction of travel.
"""
function next_edge(x::Array{<:Unitful.Length, 1}, y::Array{<:Unitful.Length, 1}, z::Array{<:Unitful.Length, 1},
                   box_id::Array{Int,1}, r::Array{<:Unitful.Length, 1}, unit_vector::Array{Float64, 1})

    direction = unit_vector .> 0.0

    # Find distance to box crossings in all dimensions
    distance = ([x[box_id[1] +  direction[1]],
                 y[box_id[2] +  direction[2]],
                 z[box_id[3] + !direction[3]]] .- r) ./unit_vector

    # Closest face
    face = argmin(distance)
    ds = distance[face]

    return face, ds
end
