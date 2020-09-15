module MCRT

import PhysLib
using Random
using LinearAlgebra
using ProgressBars
using Markdown
using Unitful

"""
    simulate(Atmosphere::Module, max_scatterings::Real, 
             scale_emission::Real, chi_max::Real,  rng::MersenneTwister)

Simulates the radiation field in a given atmosphere with 
a lower optical depth boundary given by tau_max.
"""
function simulate(Atmosphere::Module, max_scatterings::Real, tau_max::Real, total_packets::Real, rng::MersenneTwister)

    # Import relevant atmosphere data
    x = Atmosphere.x
    y = Atmosphere.y
    z = Atmosphere.height

    temperature = Atmosphere.temperature
    epsilon = Atmosphere.epsilon_continuum
    chi = Atmosphere.chi_continuum
    wavelength = Atmosphere.wavelength
    
    edge = Atmosphere.edge
    nx, ny, nz = Atmosphere.dim
    total_boxes = nx*ny*nz

    # Initialise variables
    total_destroyed = 0
    total_escaped = 0
    total_scatterings = 0
        
    boundary = PhysLib.optical_depth_boundary(chi, z, tau_max)
    total_emission = PhysLib.total_emission(chi, temperature, x, y, z, boundary, wavelength)
    scale_emission = total_packets/total_emission
    
    surface = zeros(Int, nx, ny)
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
        
        # Find dimensions of box
        corner = [x[i], y[j], z[k]]
        box_dim = ([x[i+1], y[j+1], z[k+1]] .- corner)
        
        box_volume = box_dim[1]*box_dim[2]*(-box_dim[3])

        # Based on condition in box,
        # create certain number of photons packets
        B = PhysLib.blackbody_lambda(wavelength, temperature[box_id...])
        box_emission = B*chi[box_id...]*box_volume
        packets = Int(round(box_emission*scale_emission))
        
        # Add to local field  
        J[i,j,k] += packets
        
        # Draw random numbers for initial position
        if packets > 0
            rn = rand(rng, packets, 3)
        end
        
        for packet=1:packets

            # Initial position uniformely drawn from box
            r = corner .+ (box_dim .* rn[packet,:])

            # Initial box
            box_id = [i, j, k]

            # Scatter each packet until destroyed,
            # escape or reach max_scatterings
            for s=1:Int(max_scatterings)
                
                total_scatterings += 1
                
                # Scatter packet once 
                r, box_id, escape, destroyed, J = scatter_packet!(x, y, z, chi, edge, box_id, r, J, boundary, rng)

                # Check if escaped
                if escape
                    surface[box_id[1], box_id[2]] += 1
                    total_escaped +=1
                    break
                # Check if destroyed in bottom
                elseif destroyed
                    total_destroyed += 1
                    break
                # Check if destroyed in next particle interaction
                elseif rand() < epsilon[box_id...]
                    total_destroyed += 1
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
                         chi::Array{<:Unitful.Quantity, 3}, edge::Array{<:Unitful.Length, 2}, box_id::Array{Int,1}, 
                         r::Array{<:Unitful.Length, 1}, J::Array{Int, 3}, boundary::Array{Int, 2}, rng::MersenneTwister)

    # Import relevant atmosphere data
    dim = size(chi) #Atmosphere.dim

    # Keep track of status
    escape = false
    destroyed = false

    # Draw scattering depth and direction
    τ = -log(rand(rng))
    ϕ =  π * rand(rng)
    θ = 2π * rand(rng)

    unit_vector = [sin(θ)*cos(ϕ), sin(θ)*sin(ϕ), cos(θ)]
    
    # Find distance to closest face
    ds, face = next_edge(x, y, z, r, unit_vector, box_id)
    
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
                break
            
            # Bottom destruction, Here they will very likely get destroyed anyway, so might not need this 
            elseif box_id[3] == boundary[box_id[1], box_id[2]] + 1
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
        J[box_id...] += 1

        # Find distance to closest face
        ds, face = next_edge(x, y, z, r, unit_vector, box_id)

        # Update optical depth and position
        τ_cum += ds*chi[box_id...]
        r += ds*unit_vector
    end

    if escape || destroyed
        r = nothing
    else
        # Correct for overshoot in final box
        r -= unit_vector * (τ_cum - τ)/chi[box_id...]
    end

    return r, box_id, escape, destroyed, J
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


# End module
end
