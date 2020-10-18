include("atmosphere.jl")
include("radiation.jl")
using Random
using Future # for randjump in rng when using threads
using Printf
using ProgressMeter
using HDF5
using Unitful

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
function mcrt(atmosphere::Atmosphere,
              wavelengths::Unitful.Length,
              S::Array{Int64, 3},
              max_scatterings = 1e9,
              num_bins = [4,2])

    # ==================================================================
    # ATMOSPHERE DATA
    # ==================================================================
    x = atmosphere.x
    y = atmosphere.y
    z = atmosphere.z
    χ = atmosphere.χ
    ε = atmosphere.ε
    boundary = atmosphere.boundary

    # ===================================================================
    # RADIATION DATA
    # ===================================================================
    #
    #
    #
    #
    λ = wavelengths
    ϕ_bins, θ_bins = num_bins

    # ===================================================================
    # SET UP VARIABLES
    # ===================================================================
    nx, ny, nz = size(χ)
    total_boxes = nx*ny*nz

    # Initialise variables
    surface_intensity = Tuple([zeros(Int64, nx, ny, ϕ_bins, θ_bins) for t in 1:Threads.nthreads()])
    J = Tuple([zeros(Int64, nx, ny, nz) for t in 1:Threads.nthreads() ])
    total_destroyed = Threads.Atomic{Int64}(0)
    total_scatterings = Threads.Atomic{Int64}(0)

    # Give each thread seeded rng
    rng = Tuple([Future.randjump(MersenneTwister(1), t*big(10)^20) for t in 1:Threads.nthreads()])

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
            r = corner .+ (box_dim .* rand(rng[Threads.threadid()],3))

            # Scatter each packet until destroyed,
            # escape or reach max_scatterings
            for s=1:Int(max_scatterings)

                Threads.atomic_add!(total_scatterings, 1)

                # Scatter packet once
                box_id, r, lost = scatter_packet(x, y, z, χ, boundary,
                                                 box_id, r,
                                                 J[Threads.threadid()],
                                                 surface_intensity[Threads.threadid()],
                                                 rng[Threads.threadid()],
                                                 num_bins)
                # Check if escaped or lost in bottom
                if lost
                    break
                # Check if destroyed in next particle interaction
                elseif rand(rng[Threads.threadid()]) < ε[box_id...]
                    Threads.atomic_add!(total_destroyed, 1)
                    break
                end
            end
        end
    end

    # ===================================================================
    # WRITE TO FILE
    # ===================================================================
    println("\n--Writing results to file...\n")

    # Collect results from all threads
    surface_intensity = reduce(+, surface_intensity)
    J = reduce(+, J)
    J = J .+ S

    out = h5open("../out/output.hdf5", "w")
    write(out, "total_packets", sum(S))
    write(out, "total_destroyed", total_destroyed.value)
    write(out, "total_escaped", sum(surface_intensity))
    write(out, "SNR", sqrt(maximum(surface_intensity)))
    write(out, "total_scatterings", total_scatterings.value)
    write(out, "S", S)
    write(out, "J", J)
    write(out, "surface_intensity", surface_intensity)
    close(out)
end


"""
    function scatter_packet(x::Array{<:Unitful.Length, 1},
                            y::Array{<:Unitful.Length, 1},
                            z::Array{<:Unitful.Length, 1},
                            χ::Array{<:Unitful.Quantity{<:Real, Unitful.𝐋^(-1)}, 3},
                            boundary::Array{Int, 2},
                            box_id::Array{Int,1},
                            r::Array{<:Unitful.Length, 1},
                            J::Array{Int, 3},
                            rng::MersenneTwister)

Scatters photon packet once. Returns new position, box_id and escape/destroyed-status.
"""
function scatter_packet(x::Array{<:Unitful.Length, 1},
                        y::Array{<:Unitful.Length, 1},
                        z::Array{<:Unitful.Length, 1},
                        χ::Array{PerLength, 3},
                        boundary::Array{Int64, 2},
                        box_id::Array{Int64,1},
                        r::Array{<:Unitful.Length, 1},
                        J::Array{Int64, 3},
                        surface_intensity::Array{Int64,4},
                        rng::MersenneTwister,
                        num_bins::Array{Int64, 1})

    # Keep track of status
    lost = false

    # Useful quantities
    side_dim = size(boundary)
    side_edge = [x[1] x[end]
                 y[1] y[end]]

    # ===================================================================
    # DRAW DEPTH AND DIRECTION
    # ===================================================================

    # Draw scattering depth and direction
    τ = -log(rand(rng))
    ϕ = 2π * rand(rng)
    θ =  π * rand(rng)

    # Find direction
    unit_vector = [sin(θ)*cos(ϕ), sin(θ)*sin(ϕ), cos(θ)]
    direction = Int.(sign.(unit_vector))
    direction[3] = -direction[3] # because height array up->down

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

        # Check if escaped
        if face == 3
            if box_id[3] == 0
                lost = true
                ϕ_bin = 1 + Int(ϕ÷(2π/num_bins[1]))
                θ_bin = 1 + sum(θ .> ((π/2)/(2 .^(2:num_bins[2]))))
                surface_intensity[box_id[1], box_id[2], ϕ_bin, θ_bin] += 1
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

        # Check that above boundary
        if box_id[3] > boundary[box_id[1], box_id[2]]
            lost = true
            break
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
    # CORRECT FOR OVERSHOOT IN FINAL BOX
    # ===================================================================
    if !lost
        r -= unit_vector*(τ_cum - τ)/χ[box_id...]
    end

    return box_id, r, lost
end
