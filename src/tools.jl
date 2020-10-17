using Unitful
import PhysicalConstants.CODATA2018: c_0, h, k_B
@derived_dimension NumberDensity Unitful.𝐋^-3
@derived_dimension PerLength Unitful.𝐋^-1

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
    function optical_depth_boundary(χ::Array{<:Unitful.Quantity{<:Real, Unitful.𝐋^(-1)}, 3},
                                    z::Array{<:Unitful.Length, 1},
                                    τ_max::Real)

Returns 2D array containing the k-indices where the optical depth reaches τ_max.
"""
function optical_depth_boundary(χ::Array{PerLength, 3},
                                z::Array{<:Unitful.Length, 1},
                                τ_max::Real)

    nx, ny, nz = size(χ)
    columns = nx*ny
    boundary = Array{Int, 2}(undef, nx, ny)

    # Calculate vertical optical depth for each column
    Threads.@threads for col=1:columns
        i = 1 + (col-1)÷ny
        j = col - (i-1)*ny

        τ = 0
        k = 0

        while τ < τ_max && k < ny
            k += 1
            # Trapezoidal rule
            τ += 0.5(z[k] - z[k+1]) * (χ[i,j,k] + χ[i,j,k+1])
        end
        boundary[i,j] = k
    end

    return boundary
end

"""
    function total_emission(χ::Array{<:Unitful.Quantity{<:Real, Unitful.𝐋^(-1)}, 3},
                        temperature::Array{<:Unitful.Temperature, 3},
                        x::Array{<:Unitful.Length, 1},
                        y::Threads.@threads Array{<:Unitful.Length, 1},
                        z::Array{<:Unitful.Length, 1},
                        boundary::Array{Int,2},
                        λ::Unitful.Length)

Reurns the total emission above the optical depth boundary in intensity units.
"""
function total_emission(χ::Array{PerLength, 3},
                        temperature::Array{<:Unitful.Temperature, 3},
                        x::Array{<:Unitful.Length, 1},
                        y::Array{<:Unitful.Length, 1},
                        z::Array{<:Unitful.Length, 1},
                        boundary::Array{Int,2},
                        λ::Unitful.Length)

    total_emission = 0.0u"kW / sr / nm"

    nx,ny = size(boundary)
    columns = nx*ny

    # Calculate vertical optical depth for each column
    for col=1:columns

        i = 1 + (col-1)÷ny
        j = col - (i-1)*ny

        for k=1:boundary[i,j]
            box_volume = (x[i+1] - x[i])*(y[j+1] - y[j])*(z[k] - z[k+1])
            B = blackbody_lambda(λ, temperature[i,j,k])
            total_emission += B*χ[i,j,k]*box_volume
        end
    end

    return total_emission
end

"""
    packets_per_box(x::Array{<:Unitful.Length, 1},
                         y::Array{<:Unitful.Length, 1},
                         z::Array{<:Unitful.Length, 1},
                         χ::Array{<:Unitful.Quantity{<:Real, Unitful.𝐋^(-1)}, 3},
                         temperature::Array{<:Unitful.Temperature, 3},
                         λ::Unitful.Length,
                         target_packets::Real,
                         boundary::Array{Int64,2})

Returns a 3D array of the # of packets to be generated in each box.
"""
function packets_per_box(x::Array{<:Unitful.Length, 1},
                         y::Array{<:Unitful.Length, 1},
                         z::Array{<:Unitful.Length, 1},
                         χ::Array{PerLength, 3},
                         temperature::Array{<:Unitful.Temperature, 3},
                         λ::Unitful.Length,
                         target_packets::Real,
                         boundary::Array{Int64,2})

    nx,ny = size(boundary)
    nz = maximum(boundary)
    total_boxes = nx*ny*nz

    total_emission_ = total_emission(χ, temperature, x, y, z, boundary, λ)
    scale_emission = target_packets/total_emission_

    packets = zeros(Int64,nx,ny,nz)

    for box=1:total_boxes
        i = 1 + (box-1) ÷ (ny*nz)
        j = 1 + (box - (i-1)*ny*nz - 1) ÷ nz
        k = 1 + (box - (i-1)*ny*nz - 1) % nz

        # Skip boxes beneath boundary
        if k > boundary[i,j]
            continue
        end

        B = blackbody_lambda(λ, temperature[i,j,k])
        box_volume = (x[i+1] - x[i])*(y[j+1] - y[j])*(z[k] - z[k+1])
        box_emission = B*χ[i,j,k]*box_volume
        packets[i,j,k] = Int(round(box_emission*scale_emission))
    end
    return packets
end
