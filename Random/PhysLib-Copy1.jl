module PhysLib

using Unitful
using NumericalIntegration
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


function optical_depth_boundary(chi::Array{<:Unitful.Quantity, 3}, height::Array{<:Unitful.Length, 1}, tau_max::Real)

    dim = size(chi)
    columns = dim[1]*dim[2]
    boundary = Array{Int, 2}(undef, dim[1], dim[2])

    # Calculate vertical optical depth for each column
    for col=1:columns
        i = 1 + (col - 1)÷dim[2]
        j = col - (i - 1)*dim[2]
        
        tau = 0 
        k = 0
        
        while tau < tau_max
            k += 1
            
            # Trapezoidal rule 
            tau += 0.5(height[k] - height[k+1]) * (chi[i,j,k] + chi[i,j,k+1])
        end
            
        boundary[i,j] = k
        
    end
    return boundary
end


end
