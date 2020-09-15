module PhysLib2

using Unitful
using Markdown
import PhysicalConstants.CODATA2018: c_0, h, k_B


"""
Calculates the Blackbody (Planck) function per wavelength, for given
arrays of wavelength and temperature.

Returns monochromatic intensity (energy per time per area per steradian per wavelength).

(Copied from Tiago, SSB)
"""
function blackbody_lambda(λ::Unitful.Length, temperature::Unitful.Quantity)
    (2h * c_0^2) / (λ^5 * (exp((h * c_0 / k_B) / (λ * temperature)) - 1)) |> u"kW / m^2 / sr / nm"
end


function optical_depth_boundary(Atmosphere::Module, tau_max::Real)
    
    chi = Atmosphere.chi_continuum
    height = Atmosphere.height
    
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


function total_emission(Atmosphere::Module, boundary::Array{Int,2}, wavelength::Unitful.Length)
    
    x = Atmosphere.x
    y = Atmosphere.y
    z = Atmosphere.height
    
    chi = Atmosphere.chi_continuum
    temperature = Atmosphere.temperature
    
    
    total_emission = 0.0u"kW / sr / nm"
    
    dim = size(chi)
    columns = dim[1]*dim[2]
    
    # Calculate vertical optical depth for each column
    for col=1:columns
        
        i = 1 + (col - 1)÷dim[2]
        j = col - (i - 1)*dim[2]
                
        for k=1:boundary[i,j]
            
            box_volume = (x[i+1] - x[1])*(y[j+1] - y[1])*(z[k] - z[k+1])
            B = blackbody_lambda(wavelength, temperature[i,j,k])
            total_emission += B*chi[i,j,k]*box_volume
        end
    end
                
    return total_emission
end

# End of module 
end
