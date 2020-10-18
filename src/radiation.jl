include("atmosphere.jl")
import PhysicalConstants.CODATA2018: c_0, h, k_B

struct Radiation
    λ::Unitful.Length
    S::Array{Int64,3}
    max_scatterings::Real
    escape_bins::Array{Int64,1}
end

"""
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

    box_emissivity = blackbody_lambda.(λ, temperature) .* χ
    box_emission = zeros(Float64,nx,ny,nz)u"kW / sr / nm"
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
        box_emission[i,j,k] = box_emissivity[i,j,k]*box_volume
    end

    total_emission = sum(box_emission)
    scale_emission = target_packets/total_emission
    packets = Int.(round.(box_emission*scale_emission))

    return packets
end

"""
Collects radition data that will go into structure
"""
function collect_radiation_data(atmosphere::Atmosphere, λ::Unitful.Length)
    target_packets = get_target_packets()
    max_scatterings = get_max_scatterings()
    escape_bins = get_escape_bins()
    S = packets_per_box(atmosphere, λ, target_packets)
    return S, max_scatterings, escape_bins
end
