include("../src/radiation.jl")
using InteractiveUtils
 print("\n--Loading wavelengths......................")
 λ = get_λ()[1]
 println(@sprintf(" %d wavelength(s) loaded.", length(λ)))

 # ==================================================================
 # LOAD ATMOSPHERE DATA AND CALCULATE BOUNDARY
 # ==================================================================
 print("--Loading atmosphere data..................")
 atmosphere_parameters = collect_atmosphere_data(λ)
 atmosphere = Atmosphere(atmosphere_parameters...)
 println(@sprintf(" Atmosphere loaded with tau_max = %.1f.", get_τ_max()))
 # ==================================================================
 # LOAD RADIATION DATA AND CALCULATE # PACKETS
 # ==================================================================
 print("--Loading radiation data...................")
 radiation_parameters = collect_radiation_data(atmosphere, λ)
 radiation = Radiation(λ, radiation_parameters...)
 println(@sprintf(" Radiation loaded with %.2e packets.", sum(radiation.S)))

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
 nz, nx, ny = size(χ)
 total_boxes = nx*ny*nz

 # Initialise variables
 surface_intensity = zeros(Int64, nx, ny, num_bins...)
 J = zeros(Int64, nz, nx, ny)
 total_destroyed = Threads.Atomic{Int64}(0)
 total_scatterings = Threads.Atomic{Int64}(0)

 # Give each thread seeded rng
 rng = MersenneTwister(1)

 i = j = k = 1

 # Dimensions of box
 corner = SA[z[k], x[i], y[j]]
 box_dim = SA[z[k+1], x[i+1], y[j+1]] .- corner

 # Initial box
 box_id = [k,i,j]

 # Initial position uniformely drawn from box
 r = corner .+ (box_dim .* rand(rng,3))


 """
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
     lost::Bool = false

     # Useful quantities
     side_dim = SVector(size(boundary))
     side_edge = SA[x[1] x[end]
                    y[1] y[end]]

     # ===================================================================
     # DRAW DEPTH AND DIRECTION
     # ===================================================================

     # Draw scattering depth and direction
     τ = -log(rand(rng))
     ϕ = 2π * rand(rng)
     θ =  π * rand(rng)

     # Find direction
     unit_vector = [cos(θ), sin(θ)*cos(ϕ), sin(θ)*sin(ϕ)]
     direction = Int.(sign.(unit_vector))
     direction[1] = -direction[1] # because height array up->down

     # ===================================================================
     # MOVE PACKET TO FIRST BOX INTERSECTION
     # ===================================================================

     # Next face cross in all dimensions
     next_edge = (direction .> 0) .+ box_id

     # Distance to next face cross in all dimensions
     distance = ([z[next_edge[1]],
                  x[next_edge[2]],
                  y[next_edge[3]]] .- r) ./unit_vector

     # Closest face cross
     face = argmin(distance)
     ds = distance[face]

     # Update optical depth and position
     τ_cum::Float64 = ds * χ[box_id...]
     r::Quantity{Float64,𝐋,Unitful.FreeUnits{(m,),𝐋,nothing}} = r + ds * unit_vector

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
                 ϕ_bin = 1 + Int(ϕ÷(2π/num_bins[1]))
                 θ_bin = 1 + sum(θ .> ((π/2)/(2 .^(2:num_bins[2]))))
                 surface_intensity[box_id[2], box_id[3], ϕ_bin, θ_bin] += 1
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

         # Distance to next face cross in all dimensions
         distance = ([z[next_edge[1]],
                      x[next_edge[2]],
                      y[next_edge[3]]] .- r) ./unit_vector

         # Closest face cross
         face = argmin(distance)
         ds = distance[face]

         # Update optical depth and position
         τ_cum += ds*χ[box_id...]
         r::Quantity{Float64,𝐋,Unitful.FreeUnits{(m,),𝐋,nothing}} = r + ds*unit_vector
     end

     # ===================================================================
     # CORRECT FOR OVERSHOOT IN FINAL BOX
     # ===================================================================
     if !lost
         r::Quantity{Float64,𝐋,Unitful.FreeUnits{(m,),𝐋,nothing}} = r - unit_vector*(τ_cum - τ)/χ[box_id...]
     end

     return box_id, r, lost
 end

using Profile
 # Scatter packet once
@code_warntype scatter_packet(x, y, z, χ,
                              boundary,
                              box_id, r,
                              J, surface_intensity,
                              rng, num_bins)
