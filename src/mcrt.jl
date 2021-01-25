include("radiation.jl")

"""
Simulates the radiation field in a given atmosphere with
a lower optical depth boundary given by a maximum τ.
"""
function mcrt(atmosphere::Atmosphere,
              radiation::Radiation)

    # ==================================================================
    # ATMOSPHERE DATA
    # ==================================================================
    x = atmosphere.x
    y = atmosphere.y
    z = atmosphere.z
    vx = atmosphere.velocity_x
    vy = atmosphere.velocity_y
    vz = atmosphere.velocity_z

    # ===================================================================
    # RADIATION DATA
    # ===================================================================
    λ = radiation.λ
    α = radiation.α
    ε = radiation.ε
    boundary = radiation.boundary
    packets = radiation.packets
    rad_per_packet = radiation.intensity_per_packet
    max_scatterings = radiation.max_scatterings

    # ===================================================================
    # SET UP VARIABLES
    # ===================================================================
    nλ, nz, nx, ny = size(α)

    # Open output file and initialise variables
    file = h5open("../out/output.h5", "w")
    J = create_dataset(file, "J", datatype(Int32), dataspace(nλ,nz,nx,ny), chunk=(1,nz,nx,ny))
    write(file, "total_destroyed", Array{Int64,1}(undef,nλ))
    write(file, "total_scatterings", Array{Int64,1}(undef,nλ))

    # Initialise placeholder variables
    J_λ = zeros(Int32, nz, nx, ny)
    total_destroyed = Threads.Atomic{Int64}(0)
    total_scatterings = Threads.Atomic{Int64}(0)

    # ===================================================================
    # SIMULATION
    # ===================================================================
    println(@sprintf("--Starting simulation, using %d thread(s)...",
            Threads.nthreads()))

    for λi=1:nλ
        # fill!(surface_intensity, 0.0)
        fill!(J_λ, 0.0)
        total_destroyed = Threads.Atomic{Int64}(0)
        total_scatterings = Threads.Atomic{Int64}(0)

        packets_λ = view(packets, λi,:,:,:)
        α_λ = view(α, λi,:,:,:)
        boundary_λ = view(boundary, λi, :,:)
        ε_λ = view(ε, λi,:,:,:)

        println("\n--[",λi,"/",nλ, "]        ", @sprintf("λ = %.3f nm...............", ustrip(λ[λi])*1e9))

        # Create ProgressMeter working with threads
        p = Progress(ny)
        update!(p,0)
        jj = Threads.Atomic{Int}(0)
        l = Threads.SpinLock()

        # Go through all boxes
        Threads.@threads for j=1:ny

            # Advance ProgressMeter
            Threads.atomic_add!(jj, 1)
            Threads.lock(l)
            update!(p, jj[])
            Threads.unlock(l)

            for i=1:nx
                for k=1:nz

                    # Packets in box
                    pcts = packets_λ[k,i,j]

                    if pcts == 0
                        continue
                    end

                    # Dimensions of box
                    corner = SA[z[k], x[i], y[j]]
                    box_dim = SA[z[k+1], x[i+1], y[j+1]] .- corner

                    for pct=1:pcts

                        # Initial box
                        box_id = [k,i,j]

                        # Initial position uniformely drawn from box
                        r = corner .+ (box_dim .* rand(3))

                        # Scatter each packet until destroyed,
                        # escape or reach max_scatterings
                        for s=1:Int(max_scatterings)

                            Threads.atomic_add!(total_scatterings, 1)

                            # Scatter packet once
                            box_id, r, lost = scatter_packet(x, y, z,
                                                             α_λ,
                                                             boundary_λ,
                                                             box_id, r,
                                                             J_λ,
                                                             vx, vy, vz)

                            # Check if escaped or lost in bottom
                            if lost
                                break
                            # Check if destroyed in next particle interaction
                            elseif rand() < ε_λ[box_id...]
                                Threads.atomic_add!(total_destroyed, 1)
                                break
                            end
                        end
                    end
                end
            end
        end

        # ===================================================================
        # WRITE TO FILE
        # ===================================================================
        J[λi,:,:,:] = J_λ
        file["total_destroyed"][λi] = total_destroyed.value
        file["total_scatterings"][λi] = total_scatterings.value
    end
    close(file)
    write_to_file(radiation)
end

"""
Scatters photon packet once.
Returns new position, box_id and lost-status.
"""
function scatter_packet(x::Array{<:Unitful.Length, 1},
                        y::Array{<:Unitful.Length, 1},
                        z::Array{<:Unitful.Length, 1},
                        α,
                        boundary,
                        box_id::Array{Int64,1},
                        r::Array{<:Unitful.Length, 1},
                        J::Array{Int32, 3},
                        vx, vy, vz)

    # Keep track of status
    lost = false

    # Useful quantities
    side_dim = SVector(size(boundary))
    side_edge = SA[x[1] x[end]
                   y[1] y[end]]

    # ===================================================================
    # DRAW DEPTH AND DIRECTION
    # ===================================================================

    # Draw scattering depth and direction
    τ = -log(rand())
    ϕ = 2π * rand()
    θ =  π * rand()

    # Find direction
    unit_vector = [cos(θ), sin(θ)*cos(ϕ), sin(θ)*sin(ϕ)]
    direction = Int.(sign.(unit_vector))
    direction[1] = -direction[1] # because height array up->down

    # ===================================================================
    # MOVE PACKET TO FIRST BOX INTERSECTION
    # ===================================================================

    # Next face cross in all dimensions
    next_edge = (direction .> 0) .+ box_id

    # Closest face and distance to it
    face, ds = closest_edge([z[next_edge[1]], x[next_edge[2]], y[next_edge[3]]],
                             r, unit_vector)

    τ_cum = ds * α[box_id...]
    r += ds * unit_vector

    # ===================================================================
    # TRAVERSE BOXES UNTIL DEPTH TARGET REACHED
    # ===================================================================
    while τ > τ_cum
        # Switch to new box
        box_id[face] += direction[face]
        next_edge[face] += direction[face]

        # Check if escaped
        if face == 1
            if box_id[1] == 0
                lost = true
                break
            end
        # Handle side escapes with periodic boundary
        else
            # Left-going packets
            if box_id[face] == 0
                box_id[face] = side_dim[face-1]
                next_edge[face] = side_dim[face-1]
                r[face] = side_edge[face-1,2]
            # Right-going packets
            elseif box_id[face] == side_dim[face-1] + 1
                box_id[face] = 1
                next_edge[face] = 2
                r[face] = side_edge[face-1,1]
            end
        end

        # Check that above boundary
        if box_id[1] > boundary[box_id[2], box_id[3]]
            lost = true
            break
        end

        # Add to radiation field
        J[box_id...] += 1

        # Closest face and distance to it
        face, ds = closest_edge([z[next_edge[1]], x[next_edge[2]], y[next_edge[3]]],
                                 r, unit_vector)

        τ_cum += ds * α[box_id...]
        r += ds * unit_vector
    end

    # ===================================================================
    # CORRECT FOR OVERSHOOT IN FINAL BOX
    # ===================================================================
    if !lost
        r -= unit_vector*(τ_cum - τ)/α[box_id...]
    end

    return box_id, r, lost
end


"""
Returns the face (1=z,2=x or 3=y) that
the packet will cross next and the distance to it.
"""
function closest_edge(next_edges::Array{<:Unitful.Length, 1},
                      r::Array{<:Unitful.Length, 1},
                      unit_vector::Array{Float64, 1})

    # Distance to next face cross in all dimensions
    distance = (next_edges .- r) ./unit_vector

    # Closest face cross
    face = argmin(distance)
    ds = distance[face]

    return face, ds
end
