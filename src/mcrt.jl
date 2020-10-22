include("radiation.jl")
using Random
using Future # for randjump in rng when using threads
using Printf
using ProgressMeter
using HDF5
using Unitful
#using StaticArrays

"""
Simulates the radiation field in a given atmosphere with
a lower optical depth boundary given by τ_max.
"""
function mcrt(atmosphere::Atmosphere,
              radiation::Radiation)

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
    λ = radiation.λ
    S = radiation.S
    max_scatterings = radiation.max_scatterings
    num_bins = radiation.escape_bins

    # ===================================================================
    # SET UP VARIABLES
    # ===================================================================
    nλ, nx, ny, nz = size(χ)
    total_boxes = nx*ny*nz

    # Initialise variables
    surface_intensity = Tuple([Array{UInt32,4}(undef, nx, ny, num_bins...) for t in 1:Threads.nthreads()])
    J = Tuple([Array{UInt32,3}(undef, nx, ny, nz) for t in 1:Threads.nthreads() ])
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
    p = Progress(total_boxes*nλ)
    update!(p,0)
    jj = Threads.Atomic{Int}(0)
    ll = Threads.SpinLock()

    for l=1:nλ

        # Put field to zero for each λ
        for t=1:Threads.nthreads()
            fill!(J[t],0)
            fill!(surface_intensity[t], 0)
        end
        total_destroyed = Threads.Atomic{Int64}(0)
        total_scatterings = Threads.Atomic{Int64}(0)

        # Go through all boxes
        Threads.@threads for box=1:total_boxes
            println(box)

            # Advance ProgressMeter
            Threads.atomic_add!(jj, 1)
            Threads.lock(ll)
            update!(p, jj[])
            Threads.unlock(ll)

            # Find (x,y,z) indices of box
            i = 1 + (box-1) ÷ (ny*nz)
            j = 1 + (box - (i-1)*ny*nz - 1) ÷ nz
            k = 1 + (box - (i-1)*ny*nz - 1) % nz

            # Packets in box
            packets = S[l,i,j,k]

            if packets < 1
                continue
            end

            # Dimensions of box
            corner = [x[i], y[j], z[k]]
            box_dim = [x[i+1], y[j+1], z[k+1]] .- corner
            #corner = SA[x[i], y[j], z[k]]
            #box_dim = SA[x[i+1], y[j+1], z[k+1]] .- corner

            for packet=1:packets

                # Initial box
                #box_id = SVector{3, UInt16}([i,j,k])
                box_id = UInt16.([i,j,k])

                # Initial position uniformely drawn from box
                r = corner .+ (box_dim .* rand(rng[Threads.threadid()],3))
                #r = SVector{3, Float64}(corner .+ (box_dim .* rand(rng[Threads.threadid()],3)))

                # Scatter each packet until destroyed,
                # escape or reach max_scatterings
                for s=1:Int(max_scatterings)

                    Threads.atomic_add!(total_scatterings, 1)

                    # Scatter packet once
                    box_id, r, lost = scatter_packet(x, y, z, χ[l,:,:,:], boundary[l,:,:],
                                                     box_id, r,
                                                     J[Threads.threadid()],
                                                     surface_intensity[Threads.threadid()],
                                                     rng[Threads.threadid()],
                                                     num_bins)
                    # Check if escaped or lost in bottom
                    if lost
                        break
                    # Check if destroyed in next particle interaction
                    elseif rand(rng[Threads.threadid()]) < ε[l, box_id...]
                        Threads.atomic_add!(total_destroyed, 1)
                        break
                    end


                end
            end
        end

        # Collect results from all threads
        surface_intensity = reduce(+, surface_intensity)
        J = reduce(+, J)
        J = J .+ S[l,:,:,:]

        # ===================================================================
        # WRITE TO FILE
        # ===================================================================
        println("\n--Writing results to file...\n")
        output(S[l,:,:,:], J, surface_intensity, total_destroyed.value, total_scatterings.value)

    end

end


"""
Scatters photon packet once. Returns new position, box_id and escape/destroyed-status.
"""
function scatter_packet(x::Array{<:Unitful.Length, 1},
                        y::Array{<:Unitful.Length, 1},
                        z::Array{<:Unitful.Length, 1},
                        χ::Array{PerLength, 3},
                        boundary::Array{UInt16, 2},
                        box_id::Array{UInt16,1}, #SArray{3,UInt16}
                        r::Array{<:Unitful.Length, 1}, #SArray{3,Float64}
                        J::Array{UInt32, 3},
                        surface_intensity::Array{UInt32,4},
                        rng::MersenneTwister,
                        num_bins::Array{UInt16, 1})

    # Keep track of status
    lost::Bool = false

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
    unit_vector = [sin(θ)*cos(ϕ), sin(θ)*sin(ϕ), cos(θ)] #SA
    direction = Int.(sign.(unit_vector)) #SArray{3}(Int.(sign.(unit_vector)))
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
