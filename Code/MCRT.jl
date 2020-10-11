include("atmos.jl")
include("MyLibs/physLib.jl")
using Random
using Printf
using ProgressMeter

"""
    function simulate(atmosphere::Atmosphere,
                      wavelengths::Unitful.Length,
                      max_scatterings::Real,
                      τ_max::Real,
                      target_packets::Real,
                      num_bins = [4,2])

Simulates the radiation field in a given atmosphere with
a lower optical depth boundary given by τ_max.
"""
function simulate(atmosphere::Atmosphere,
                  wavelengths::Unitful.Length,
                  max_scatterings::Real,
                  τ_max::Real,
                  target_packets::Real,
                  num_bins = [4,2])

    # ==================================================================
    # ATMOSPHERE DATA
    # ==================================================================
    x = atmosphere.x
    y = atmosphere.y
    z = atmosphere.z
    ε = atmosphere.ε_continuum
    χ = atmosphere.χ_continuum
    temperature = atmosphere.temperature

    # ===================================================================
    # CHOSEN WAVELENGTHS AND ESCAPE BINS
    # ===================================================================
    λ = wavelengths
    ϕ_bins, θ_bins = num_bins

    # ===================================================================
    # PRE-CALCULATIONS
    # ===================================================================
    println("--Performing pre-calculations...")

    # Find boundary for given τ_max
    boundary = optical_depth_boundary(χ, z, τ_max)

    # Number of boxes
    nx, ny = size(boundary)
    nz = maximum(boundary)
    total_boxes = nx*ny*nz

    # Initialise variables
    surface_intensity = Tuple([zeros(Int, nx, ny, ϕ_bins, θ_bins) for t in 1:Threads.nthreads()])
    J = Tuple([zeros(Int, nx, ny, nz) for t in 1:Threads.nthreads() ])

    total_destroyed = Threads.Atomic{Int}(0)
    total_scatterings = Threads.Atomic{Int}(0)

    # Find number of packets per box and add to source function
    S = packets_per_box(x,y,z,χ,temperature,
                        λ,target_packets,boundary)

    # ===================================================================
    # SIMULATION
    # ===================================================================
    println(@sprintf("--Starting simulation, using %d thread(s)...\n",
            Threads.nthreads()))

    # Create ProgressMeter working with threads
    p = Progress(total_boxes)
    update!(p,0)
    jj = Threads.Atomic{Int}(0)
    l = Threads.SpinLock()

    # Go through all boxes
    Threads.@threads for box in 1:total_boxes

        # Advance ProgressMeter
        Threads.atomic_add!(jj, 1)
        Threads.lock(l)
        update!(p, jj[])
        Threads.unlock(l)

        # Find (x,y,z) indices of box
        i = 1 + (box-1) ÷ (ny*nz)
        j = 1 + (box - (i-1)*ny*nz - 1) ÷ nz
        k = 1 + (box - (i-1)*ny*nz - 1) % nz

        # Skip boxes beneath boundary
        if k > boundary[i,j]
            continue
        end

        # Packets in box
        packets = S[i,j,k]

        if packets < 1
            continue
        end

        # Dimensions of box
        corner = [x[i], y[j], z[k]]
        box_dim = [x[i+1], y[j+1], z[k+1]] .- corner

        for packet=1:packets

            # Initial box
            box_id = [i,j,k]

            # Initial position uniformely drawn from box
            r = corner .+ (box_dim .* rand(3))

            # Scatter each packet until destroyed,
            # escape or reach max_scatterings
            for s=1:Int(max_scatterings)

                Threads.atomic_add!(total_scatterings, 1)

                # Scatter packet once
                box_id, r, escaped, destroyed = scatter_packet(x, y, z, χ, boundary,
                                                               box_id, r,
                                                               J[Threads.threadid()])
                # Check if escaped
                if escaped[1]
                    ϕ, θ  = escaped[2]
                    ϕ_bin = 1 + Int(ϕ÷(2π/ϕ_bins))
                    θ_bin = 1 + sum(θ .> ((π/2)/(2 .^(2:θ_bins))))

                    surface_intensity[Threads.threadid()][box_id[1], box_id[2], ϕ_bin, θ_bin] += 1
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

    # ===================================================================
    # POST-CALCULATIONS
    # ===================================================================
    println("\n--Performing post-calculations...")

    # Collect results from all threads
    surface_intensity = reduce(+, surface_intensity)
    J = reduce(+, J)

    # Evaluate field
    J = J .+ S
    mean_J, min_J, max_J = field_above_boundary(z, χ, J, τ_max)
    J_data = [J, S, mean_J, min_J, max_J]

    # Actual number of packets generated ≈ target_packets
    total_packets = sum(S)
    total_escaped = sum(surface_intensity)
    packet_data = [total_packets, total_destroyed.value,
                   total_escaped, total_scatterings.value]

    return packet_data, J_data, surface_intensity
end


"""
    function scatter_packet(x::Array{<:Unitful.Length, 1},
                            y::Array{<:Unitful.Length, 1},
                            z::Array{<:Unitful.Length, 1},
                            χ::Array{<:Unitful.Quantity{<:Real, Unitful.𝐋^(-1)}, 3},
                            boundary::Array{Int, 2},
                            box_id::Array{Int,1},
                            r::Array{<:Unitful.Length, 1},
                            J::Array{Int, 3})

Scatters photon packet once. Returns new position, box_id and escape/destroyed-status.
"""
function scatter_packet(x::Array{<:Unitful.Length, 1},
                        y::Array{<:Unitful.Length, 1},
                        z::Array{<:Unitful.Length, 1},
                        χ::Array{<:Unitful.Quantity{<:Real, Unitful.𝐋^(-1)}, 3},
                        boundary::Array{Int, 2},
                        box_id::Array{Int,1},
                        r::Array{<:Unitful.Length, 1},
                        J::Array{Int, 3})

    # Keep track of status
    escaped = [false, nothing]
    destroyed = false

    # Useful quantities
    side_dim = size(boundary)
    side_edge = [x[1] x[end]
                 y[1] y[end]]

    # ===================================================================
    # DRAW DEPTH AND DIRECTION
    # ===================================================================

    # Draw scattering depth and direction
    τ = -log(rand())
    ϕ = 2π * rand()
    θ =  π * rand()

    # Find direction
    unit_vector = [sin(θ)*cos(ϕ), sin(θ)*sin(ϕ), cos(θ)]
    direction = sign.(unit_vector)
    direction[3] = -direction[3] # height array up->down

    # ===================================================================
    # MOVE PACKET TO FIRST BOX INTERSECTION
    # ===================================================================

    # Next face cross in all dimensions
    next_edge = (direction .> 0) .+ box_id

    # Distance to next face cross in all dimensions
    distance = ([x[next_edge[1]],
                 y[next_edge[2]],
                 z[next_edge[3]]] .- r) ./unit_vector

    # Closest face cross
    face = argmin(distance)
    ds = distance[face]

    # Update optical depth and position
    τ_cum = ds * χ[box_id...]
    r += ds * unit_vector

    # ===================================================================
    # TRAVERSE BOXES UNTIL DEPTH TARGET REACHED
    # ===================================================================
    while τ > τ_cum

        # Switch to new box
        box_id[face] += direction[face]
        next_edge[face] += direction[face]

        # Check that within bounds of atmosphere
        if face == 3
            # Top escape
            if box_id[3] < 1 # == 0
                escaped = [true, [ϕ, θ]]
                break

            # Bottom destruction
            elseif box_id[3] > boundary[box_id[1], box_id[2]] #+ 1 wronggggg
                destroyed = true
                break
            end

        # Handle side escapes with periodic boundary
        else
            # Left-going packets
            if box_id[face] == 0
                box_id[face] = side_dim[face]
                next_edge[face] = side_dim[face]
                r[face] = side_edge[face,2]
            # Right-going packets
            elseif box_id[face] == side_dim[face] + 1
                box_id[face] = 1
                next_edge[face] = 2
                r[face] = side_edge[face,1]
            end
        end

        # Add to radiation field
        J[box_id...] += 1

        # Distance to next face cross in all dimensions
        distance = ([x[next_edge[1]],
                     y[next_edge[2]],
                     z[next_edge[3]]] .- r) ./unit_vector

        # Closest face cross
        face = argmin(distance)
        ds = distance[face]

        # Update optical depth and position
        τ_cum += ds*χ[box_id...]
        r += ds*unit_vector
    end

    # ===================================================================
    # CHECK PACKET STATUS BEFORE RETURNING
    # ===================================================================
    if escaped[1] || destroyed
        r = nothing
    else
        # Correct for overshoot in final box
        r -= unit_vector*(τ_cum - τ)/χ[box_id...]
    end

    return box_id, r, escaped, destroyed
end
