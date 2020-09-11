module MCRT

import PhysLib
using LinearAlgebra
using ProgressBars
using Markdown
using Unitful



"""
    MC(Atmosphere::Module, max_scatterings::Int, scale_emission::Real, chi_max::Real)


Parameters
----------
    - Atmosphere::Module
    - max_scatterings::Int
    - scale_emission::Real
    - tau_max::Real
"""
function simulate(Atmosphere::Module, max_scatterings::Int, scale_emission::Real, tau_max::Real)

    # Import relevant atmosphere data
    x = Atmosphere.x
    y = Atmosphere.y
    z = Atmosphere.height

    temperature = Atmosphere.temperature
    epsilon = Atmosphere.epsilon_continuum
    chi = Atmosphere.chi_continuum
    wavelength = Atmosphere.wavelength

    nx, ny, nz = Atmosphere.dim
    total_boxes = nx*ny*nz

    # Initialise variables
    total_packets = 0
    total_destroyed = 0
    total_escaped = 0
    avg_scatterings = 0
        
    boundary = PhysLib.optical_depth_boundary(chi, z, tau_max)
    
    hits = zeros(Int, nx, ny)
    J = zeros(Int, nx, ny, nz)

    # Go through all boxes
    for box in ProgressBar(1:total_boxes)

        # Find (x,y,z) indices of box
        i = 1 + box ÷ (ny*nz + 1)
        j = 1 + (box - (i-1)*ny*nz) ÷ (nz + 1)
        k = 1 + (box - (i-1)*ny*nz - 1) % (nz)

        # Skip boxes beneath boundary 
        if k > boundary[i,j]
            continue
        end
            
        box_id = [i, j, k]

        # Based on condition in box,
        # create certain numbers of photons packets
        B = PhysLib.blackbody_lambda(wavelength, temperature[box_id...])
        emission = B*chi[box_id...]
                
        photon_packets = Int(floor(ustrip(emission)*scale_emission))
        J[i,j,k] += photon_packets
        total_packets += photon_packets
        
        box_dim = ([x[i+1], y[j+1], z[k+1]] .- [x[i], y[j], z[k]])
        
        for packet=1:photon_packets

            # Initial position uniformely drawn from box
            r =  [x[i], y[j], z[k]] .+ (box_dim .* rand(3))

            # Initial box
            box_id = [i, j, k]

            # Scatter each packet until destroyed,
            # escape or reach max_scatterings
            for s=1:Int(max_scatterings)
                
                avg_scatterings += 1
                
                # Scatter packet
                r, box_id, escape, destroyed, J = scatter_packet!(Atmosphere, box_id, r, J, boundary)

                # Check if escaped
                if escape
                    hits[box_id[1], box_id[2]] += 1
                    total_escaped +=1
                    break
                # Destroyed in bottom
                elseif destroyed
                    total_destroyed += 1
                    break
                # Check if destroyed in new cell
                elseif rand() < epsilon[box_id...]
                    total_destroyed += 1
                    break

                end
            end
        end
    end
    
    avg_scatterings /= total_packets

    return hits, total_destroyed, total_escaped, total_packets, avg_scatterings
end



"""
    scatter(Atmos::Module, box_id::Array{Int,1}, r::Array{<:Unitful.Length, 1}, J::Array{Int, 3})

Parameters
----------
    - Atmosphere::Module
    - r::Array{<:Unitful.Length, 1}
    - box_id::Array{Int, 1}


"""
function scatter_packet!(Atmosphere::Module, box_id::Array{Int,1}, r::Array{<:Unitful.Length, 1}, J::Array{Int, 3}, boundary::Array{Int, 2})

    # Import relevant atmosphere data
    chi = Atmosphere.chi_continuum
    dim = Atmosphere.dim
    edge = Atmosphere.edge

    # Keep track of status
    escape = false
    destroyed = false

    # Draw scattering depth and direction
    τ = -log(rand()) # Check this
    ϕ =  π * rand()
    θ = 2π * rand()

    unit_vector = [sin(θ)*cos(ϕ), sin(θ)*sin(ϕ), cos(θ)]

    # Find distance to closest face
    ds, face = next_edge(Atmosphere, r, unit_vector, box_id)
    direction = sign.(unit_vector)

    # Add optical depth and update position
    τ_cum = ds * chi[box_id...]
    r += ds*unit_vector

    # If depth target not reached in current box,
    # traverse boxes until target is reached
    while τ > τ_cum

        # Switch to new box
        if face == 3
            box_id[3] -= direction[3] # Consequence of height array up->down
            
            # Top escape 
            if box_id[3] == 0
                escape = true 
                box_id[face] = 1
                break
            
            # Bottom destruction, Here they will very likely get destroyed anyway, so might not need this 
            elseif box_id[3] == boundary[box_id[1], box_id[2]] + 1
                destroyed = true
                box_id[3] -= 1
                break
            end
            
        else
            box_id[face] += direction[face]
                
            # Handle side escapes with periodic boundary
            if box_id[face] == 0
                box_id[face] = dim[face]
                r[face] = edge[face][2]

            elseif  box_id[face] == dim[face] + 1
                box_id[face] = 1
                r[face] = edge[face][1]
            end
            
        end
            
        # Add to radiation field 
        J[box_id...] += 1

        # Find distance to closest face
        ds, face = next_edge(Atmosphere, r, unit_vector, box_id)

        # Update optical depth and position
        τ_cum += ds*chi[box_id...]
        r += ds*unit_vector
    end

    if escape || destroyed
        r = nothing
    else
        # Correct for overshoot
        r -= unit_vector * (τ_cum - τ)/chi[box_id...]
    end

    return r, box_id, escape, destroyed, J
end


"""
    next_edge(Atmosphere::Module, r::Array{Unitful.Length, 1},
              vector::Array{Float64, 1}, box_id::Array{Int,1})

Calculates the distance to the next box and the face that it crosses,
given an initial position and direction of travel.

Parameters
----------
    - Atmosphere::Module
      Collection of atmosphere parameters.

    - r::Array{Float64, 1}


    - unit_vector::Array{Float64, 1}


    - box_id


"""
function next_edge(Atmosphere::Module, r::Array{<:Unitful.Length, 1},
                   unit_vector::Array{Float64, 1}, box_id::Array{Int,1})

    x, y, z = Atmosphere.x, Atmosphere.y, Atmosphere.height

    direction = unit_vector .> 0.0

    # Find distance to box crossings
    # in all dimensions
    distance = ([x[box_id[1] + direction[1]],
                y[box_id[2] + direction[2]],
                z[box_id[3] + !direction[3]]] .- r) ./unit_vector

    # Closest face
    face = argmin(distance)
    ds = distance[face]

    return ds, face
end


# End module
end
