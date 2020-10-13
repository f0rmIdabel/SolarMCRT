using Unitful
import PhysicalConstants.CODATA2018: c_0, h, k_B

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
function optical_depth_boundary(χ::Array{<:Unitful.Quantity{<:Real, Unitful.𝐋^(-1)}, 3},
                                z::Array{<:Unitful.Length, 1},
                                τ_max::Real)

    nx, ny, nz = size(χ)
    columns = nx*ny
    boundary = Array{Int, 2}(undef, nx, ny)

    # Calculate vertical optical depth for each column
    for col=1:columns
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
                        y::Array{<:Unitful.Length, 1},
                        z::Array{<:Unitful.Length, 1},
                        boundary::Array{Int,2},
                        λ::Unitful.Length)

Reurns the total emission above the optical depth boundary in intensity units.
"""
function total_emission(χ::Array{<:Unitful.Quantity{<:Real, Unitful.𝐋^(-1)}, 3},
                        temperature::Array{<:Unitful.Temperature, 3},
                        x::Array{<:Unitful.Length, 1},
                        y::Array{<:Unitful.Length, 1},
                        z::Array{<:Unitful.Length, 1},
                        boundary::Array{Int,2},
                        λ::Unitful.Length)

    total_emission = 0.0u"kW / sr / nm"

    dim = size(χ)
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
    function extract_surface_bin(surface::Array{Int, 4},
                                 bin = :[:,:])

Slices the surface for the given escape bins of ϕ and θ.
"""
function extract_surface_bin(surface::Array{Int, 4},
                             bin = :[:,:])

    # Pick out bin
    surface = surface[:,:, eval(bin)[1], eval(bin)[2]]
    dim = length(size(surface))

    # If no bin specified
    if dim == 4
        surface = sum(surface, dims = 4)
        surface = surface[:,:,:,1]
    end

    # If less than one bin specified
    if dim == 4 || dim == 3
        surface = sum(surface, dims=3)
        surface = surface[:,:,1]
    end

    return surface
end

"""
    function field_above_boundary(z::Array{<:Unitful.Length, 1},
                                  χ::Array{<:Unitful.Quantity{<:Real, Unitful.𝐋^(-1)}, 3},
                                  J::Array{Int, 3},
                                  τ_max::Real)

Calculates the average, maximum and minimum intensity value of the field above the boundary.
"""
function field_above_boundary(z::Array{<:Unitful.Length, 1},
                              χ::Array{<:Unitful.Quantity{<:Real, Unitful.𝐋^(-1)}, 3},
                              J::Array{Int, 3},
                              τ_max::Real)

    boundary = optical_depth_boundary(χ, z, τ_max)

    dim = size(χ)
    columns = dim[1]*dim[2]

    # Initialize variables
    meanJ = 0.0
    minJ =  J[1,1,1]
    maxJ =  J[1,1,1]

    total_boxes = 0

    # Calculate mean, min and max J above boundary
    for col=1:columns

        i = 1 + (col - 1)÷dim[2]
        j = col - (i - 1)*dim[2]

        for k=1:boundary[i,j]
            meanJ += J[i,j,k]
            total_boxes += 1

            if J[i,j,k] > maxJ
                maxJ = J[i,j,k]
            elseif J[i,j,k] < minJ
                minJ = J[i,j,k]
            end
        end
    end

    meanJ = meanJ/total_boxes
    return meanJ, minJ, maxJ
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
                         χ::Array{<:Unitful.Quantity{<:Real, Unitful.𝐋^(-1)}, 3},
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

        box_dim = [x[i+1], y[j+1], z[k+1]] .- [x[i], y[j], z[k]]
        box_volume = box_dim[1]*box_dim[2]*(-box_dim[3])

        # Based on condition in box,
        # create certain number of photons packets
        B = blackbody_lambda(λ, temperature[i,j,k])
        box_emission = B*χ[i,j,k]*box_volume
        packets[i,j,k] = Int(round(box_emission*scale_emission))
    end
    return packets
end
