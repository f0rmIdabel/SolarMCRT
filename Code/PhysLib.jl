module PhysLib

using Unitful
using Markdown
import PhysicalConstants.CODATA2018: c_0, h, k_B


"""
Calculates the Blackbody (Planck) function per wavelength, for given
arrays of wavelength and temperature.

Returns monochromatic intensity (energy per time per area per steradian per wavelength).

(Copied from Tiago, SSB)
"""
function blackbody_lambda(λ::Unitful.Length, temperature::Unitful.Temperature)
    (2h * c_0^2) / (λ^5 * (exp((h * c_0 / k_B) / (λ * temperature)) - 1)) |> u"kW / m^2 / sr / nm"
end

"""
2D array  of indices of
"""
function optical_depth_boundary(χ::Array{<:Unitful.Quantity, 3}, height::Array{<:Unitful.Length, 1}, τ_max::Real)

    dim = size(χ)
    columns = dim[1]*dim[2]
    boundary = Array{Int, 2}(undef, dim[1], dim[2])

    # Calculate vertical optical depth for each column
    for col=1:columns
        i = 1 + (col - 1)÷dim[2]
        j = col - (i - 1)*dim[2]

        τ = 0
        k = 0

        while τ < τ_max
            k += 1
            # Trapezoidal rule
            τ += 0.5(height[k] - height[k+1]) * (χ[i,j,k] + χ[i,j,k+1])
        end
        boundary[i,j] = k

    end
    return boundary
end

"""
The total emission above the optical depth boundary
"""

function total_emission(chi::Array{<:Unitful.Quantity, 3}, temperature::Array{<:Unitful.Temperature, 3},
                        x::Array{<:Unitful.Length, 1}, y::Array{<:Unitful.Length, 1}, z::Array{<:Unitful.Length, 1},
                        boundary::Array{Int,2}, λ::Unitful.Length)

    total_emission = 0.0u"kW / sr / nm"

    dim = size(chi)
    columns = dim[1]*dim[2]

    # Calculate vertical optical depth for each column
    for col=1:columns

        i = 1 + (col - 1)÷dim[2]
        j = col - (i - 1)*dim[2]

        for k=1:boundary[i,j]

            box_volume = (x[i+1] - x[1])*(y[j+1] - y[1])*(z[k] - z[k+1])
            B = blackbody_lambda(λ, temperature[i,j,k])
            total_emission += B*chi[i,j,k]*box_volume
        end
    end

    return total_emission
end




function extract_surface_bin(surface::Array{Int, 4}, bin = :[:,:])

    surface = surface[:,:, eval(bin)[1], eval(bin)[2]]

    dim = length(size(surface))

    if dim == 4
        surface = sum(surface, dims = 4)
        surface = surface[:,:,:,1]
    end

    if dim == 4 || dim == 3
        surface = sum(surface, dims=3)
        surface = surface[:,:,1]
    end
    
    return surface
end


# End of module
end
