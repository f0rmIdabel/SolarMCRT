include("atmosphere.jl")
import PhysicalConstants.CODATA2018: c_0, h, k_B

struct Radiation
    λ::Array{<:Unitful.Length}
    S::Array{Int64,3}
    max_scatterings::Real
    num_bins::Array{Int64,2}
end

"""
    function blackbody_lambda(λ::Unitful.Length,
                              temperature::Unitful.Temperature)

Calculates the Blackbody (Planck) function per wavelength, for given
arrays of wavelength and temperature. Returns monochromatic intensity.
Copied from Tiago, SSB
"""
function blackbody_lambda(λ::Unitful.Length,
                          temperature::Unitful.Temperature)
    (2h * c_0^2) / ( λ^5 * (exp((h * c_0 / k_B) / (λ * temperature)) - 1) ) |> u"kW / m^2 / sr / nm"
end


"""


Returns a 3D array of the # of packets to be generated in each box.
"""
function packets_per_box(atmosphere::Atmosphere,
                         λ::Unitful.Length,
                         target_packets::Real)

    x = atmosphere.x
    y = atmosphere.y
    z = atmosphere.z
    χ = atmosphere.χ
    temperature = atmosphere.temperature
    boundary = atmosphere.boundary

    nx, ny = size(boundary)
    nz = maximum(boundary)
    total_boxes = nx*ny*nz

    total_emission = 0.0u"kW / sr / nm"
    #total_emission_ = total_emission(atmosphere, λ)
    box_emission = blackbody_lambda.(λ, temperature) .* χ
    packets = zeros(Int64,nx,ny,nz)

    for box=1:total_boxes
        i = 1 + (box-1) ÷ (ny*nz)
        j = 1 + (box - (i-1)*ny*nz - 1) ÷ nz
        k = 1 + (box - (i-1)*ny*nz - 1) % nz

        # Skip boxes beneath boundary
        if k > boundary[i,j]
            continue
        end

        box_volume = (x[i+1] - x[i])*(y[j+1] - y[j])*(z[k] - z[k+1])
        box_emission[i,j,k] *= box_volume
        total_emission += box_emission[i,j,k]
    end

    scale_emission = target_packets/total_emission_
    packets = Int.(round.(box_emission*scale_emission))

    return packets
end


function get_radiation_data(atmosphere::Atmos)
