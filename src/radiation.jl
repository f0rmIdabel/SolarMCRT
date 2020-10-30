include("atmosphere.jl")

struct Radiation
    λ::Array{Unitful.Length,1}    # (nλ)
    S::Array{Int64,3}             # (nz*, nx, ny)
    max_scatterings::Real
    escape_bins::Array{Int64,1}   # (nϕ, nθ)
end


"""
Collects radition data to go into structure.
"""
function collect_radiation_data(atmosphere::Atmosphere, λ::Unitful.Length)
    # Read from input file
    target_packets = get_target_packets()
    max_scatterings = get_max_scatterings()
    escape_bins = get_escape_bins()

    # Calculate escapes per box
    S = packets_per_box(atmosphere, λ, target_packets)
    return S, max_scatterings, escape_bins
end


"""
Returns a 3D array of the # of packets to be
generated in each box above the boundary.
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

    nz, nx, ny = size(χ)
    total_boxes = nz*nx*ny

    box_emissivity = blackbody_lambda.(λ, temperature) .* χ
    box_emission = zeros(Float64,nz,nx,ny)u"kW / sr / nm"
    packets = zeros(Int64,nz,nx,ny)

    @Threads.threads for box=1:total_boxes
        j = 1 + (box-1) ÷ (nx*nz)
        i = 1 + (box - (j-1)*nx*nz - 1) ÷ nz
        k = 1 + (box - (j-1)*nx*nz - 1) % nz

        # Skip boxes beneath boundary
        if k > boundary[i,j]
            continue
        end

        box_volume = (x[i+1] - x[i])*(y[j+1] - y[j])*(z[k] - z[k+1])
        box_emission[k,i,j] = box_emissivity[k,i,j]*box_volume
    end

    total_emission = sum(box_emission)
    scale_emission = target_packets/total_emission
    packets = Int.(round.(box_emission*scale_emission))

    return packets
end


"""
Calculates the Blackbody (Planck) function per wavelength,
for given arrays of wavelength and temperature.
Returns monochromatic intensity.
Credit: Tiago
"""
function blackbody_lambda(λ::Unitful.Length,
                          temperature::Unitful.Temperature)
    (2h * c_0^2) / ( λ^5 * (exp((h * c_0 / k_B) / (λ * temperature)) - 1) ) |> u"kW / m^2 / sr / nm"
end
