

boxes = 29719872

tau = [5,10,15,20]
packets = [1e13, 1e14, 1e15]


it_per_sec = [23087.2  2423.6  222.5 
              32648.3  3271.3  330.8 
              36336.1  3714.6  384.7
              36817.2  3813.6  401.7]


destroyed = [151_800_578  1_519_176_904  15_192_919_612
             155_522_073  1_556_446_148  15_565_601_916
             155_869_453  1_560_031_402  15_601_457_687
             156_247_741  1_563_834_092  15_639_477_891]

escaped = [5_269_575  53_708_942  538_073_012
           1_801_327  18_957_539  190_651_365
             894_509   9_763_126   98_696_835 
             553_554   6_265_110   63_710_930]

scatterings = [157_168_769  1_573_873_117  15_740_864_554
               157_400_059  1_576_169_883  15_763_931_160
               156_832_426  1_570_484_434  15_707_054_985
               156_865_591  1_570_750_029  15_709_697_269]

time = boxes ./ it_per_sec

total_photons = destroyed .+ escaped
avg_scatterings = scatterings ./total_photons

escape_rate = escaped ./total_photons

ENV["GKSwstype"]="nul"

# escape rate VS tau

plot(tau, escape_rate[:,1], yscale=:log10, label = "1e13", xlabel="Optical depth boundary", ylabel = "Escape rate", color = "red", linewidth = 4.0, opacity = 0.5, linestyle =:dot)
plot!(tau, escape_rate[:,2], label = "1e14", color = "blue", opacity = 0.5, linestyle =:dot, linewidth = 4.0)
plot!(tau, escape_rate[:,3], label = "1e15", color = "green", opacity = 0.5, linestyle =:dot, linewidth = 4.0)


png("/mn/stornext/u3/idarhan/SolarMCRT/Results/Plots/escapesVStau")


# avg scat VS tau
plot(tau, avg_scatterings[:,1], yscale=:log10, label = "1e13", xlabel="Optical depth boundary", ylabel = "Avg scatterings", color = "red", opacity = 0.5, linestyle =:dot, linewidth = 4.0)
plot!(tau, avg_scatterings[:,2], label = "1e14", color = "blue", opacity = 0.5, linestyle =:dot, linewidth = 4.0)
plot!(tau, avg_scatterings[:,3], label = "1e15", color = "green", markersize = 5, opacity = 0.5, linestyle =:dot, linewidth = 4.0)
png("/mn/stornext/u3/idarhan/SolarMCRT/Results/Plots/scatVStau")


# time VS escaped 
# avg scat VS tau

plot(time[:,1], escaped[:,1], yscale=:log10, label = "1e13", xlabel="Time", ylabel = "Escapes", color = "red", markersize = 5, opacity = 0.5, linestyle =:dot, linewidth = 4.0)
plot!(time[:,2], escaped[:,2], label = "1e14", color = "blue", opacity = 0.5, linestyle =:dot, linewidth = 4.0)
plot!(time[:,3], escaped[:,3], label = "1e15", color = "green", opacity = 0.5, linestyle =:dot, linewidth = 4.0)
png("/mn/stornext/u3/idarhan/SolarMCRT/Results/Plots/timeVsEscapes")




"""
Same as scatter!, but without periodic boundaries
This is slower and gives fewer photons in the field. 
So probably not useful.
"""
function scatter_packet_no_periodic!(Atmosphere::Module, box_id::Array{Int,1}, r::Array{<:Unitful.Length, 1}, J::Array{Int, 3}, boundary::Array{Int, 2})

    # Import relevant atmosphere data
    chi = Atmosphere.chi_continuum
    dim = Atmosphere.dim
    edge = Atmosphere.edge

    # Keep track of status
    escape = false
    destroyed = false
    side_escape = false

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
                break
            
            # Bottom destruction, Here they will very likely get destroyed anyway, so might not need this 
            elseif box_id[3] == boundary[box_id[1], box_id[2]] + 1
                destroyed = true
                break
            end
            
        else
            box_id[face] += direction[face]
                
            # Handle side escapes with periodic boundary
            if box_id[face] == 0 || box_id[face] == dim[face] + 1
                side_escape = true
                break
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

    if escape || destroyed || side_escape
        r = nothing
    else
        # Correct for overshoot
        r -= unit_vector * (τ_cum - τ)/chi[box_id...]
    end

    return r, box_id, escape, destroyed, J, side_escape
end



"""
    optical_depth(chi::Array{<:Unitful.Quantity, 3}, height::Array{<:Unitful.Length, 1})

Calculate optical depth given

Future improvement: remove ustrip

Not used in current code 
"""
function optical_depth(chi::Array{<:Unitful.Quantity, 3}, height::Array{<:Unitful.Length, 1})

    dim = size(chi)
    tau = Array{Float64, 3}(undef, dim...) 
    columns = dim[1]*dim[2]

    # Calculate vertical optical depth for each column
    for col=1:columns
        i = 1 + (col - 1)÷dim[2]
        j = col - (i - 1)*dim[2]
        tau[i,j,:] = -cumul_integrate(ustrip(height[1:end-1]), ustrip(chi[i,j,:]))
    end
    return tau
end
