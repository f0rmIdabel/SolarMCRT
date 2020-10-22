include("atmosphere.jl")
import PhysicalConstants.CODATA2018: c_0, h, k_B

struct Radiation
    # (nλ)
    λ::Array{<:Unitful.Length,1}
    # (nλ, nx, ny, nz)
    S::Array{UInt32,4}
    max_scatterings::Real
    # (2)
    escape_bins::Array{UInt16,1}  ##SArray{2,UInt16}
end

"""
Collects radition data that will go into structure
"""
function collect_radiation_data(atmosphere::Atmosphere, λ::Array{<:Unitful.Length,1})
    # Read from input file
    target_packets = get_target_packets()
    max_scatterings = get_max_scatterings()
    escape_bins = get_escape_bins()
    # Calculate escapes per box
    S = packets_per_box(atmosphere, λ, target_packets)
    return S, max_scatterings, escape_bins
end

"""
Returns a 3D array of the # of packets to be generated in each box.
"""
function packets_per_box(atmosphere::Atmosphere,
                         λ::Array{<:Unitful.Length,1},
                         target_packets::Real)
    x = atmosphere.x
    y = atmosphere.y
    z = atmosphere.z
    χ = atmosphere.χ
    temperature = atmosphere.temperature
    boundary = atmosphere.boundary

    nλ, nx, ny, nz = size(χ)
    total_boxes = nx*ny*nz

    box_emission = zeros(Float64,nλ,nx,ny,nz)u"kW / sr / nm"

    @Threads.threads for box=1:total_boxes
        k = 1 + (box-1) ÷ (ny*nx)
        j = 1 + (box - (k-1)*ny*nx - 1) ÷ nx
        i = 1 + (box - (k-1)*ny*nx - 1) % nx

        box_volume = (x[i+1] - x[i])*(y[j+1] - y[j])*(z[k] - z[k+1])

        for l=1:nλ
            # Skip boxes beneath boundary
            if k > boundary[l,i,j]
                continue
            end

            box_emission[l,i,j,k] = blackbody_lambda(λ[l], temperature[i,j,k])*χ[l,i,j,k]*box_volume
        end
    end

    total_emission = sum(box_emission)
    scale_emission = target_packets/total_emission
    packets = UInt32.(round.(box_emission*scale_emission))

    return packets
end


"""
Calculates the Blackbody (Planck) function per wavelength, for given
arrays of wavelength and temperature. Returns monochromatic intensity.
Copied from Tiago, SSB# Skip boxes beneath boundary
"""
function blackbody_lambda(λ::Unitful.Length,
                          temperature::Unitful.Temperature)
    (2h * c_0^2) / ( λ^5 * (exp((h * c_0 / k_B) / (λ * temperature)) - 1) ) |> u"kW / m^2 / sr / nm"
end
