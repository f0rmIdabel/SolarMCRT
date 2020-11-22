"""
File from Tiago,

changed from (x,y,z) to (z,x,y) to fit with convention used in rest of project
"""

using Unitful
using Images: warp, indices_spatial
using CoordinateTransformations, Interpolations


"""
    shift_image(image::Array, shift_x::Real, shift_y::Real)

Shift a 2D array by an amount of `shift_x` and `shift_y` pixels,
in the first and seecond dimensions. Assumes the image is horizontally
periodic, so the returned array has the same dimensions as `image`.
"""
function shift_image(image::Array, shift_x::Real, shift_y::Real)
    transl = Translation(shift_x, shift_y)

    u = unit(image[1,1])
    image = ustrip(image)

    warp(image, transl, indices_spatial(image), Linear(), Periodic())
    image *= u
end


"""
    shift_variable!(var::Array, height::Array{<:Unitful.Length, 1},
                         pixel_size::Unitful.Length, μ::Real)

Shift (or translate) a 3D array that is horizontally periodic in the first two dimensions
according to a polar angle θ given by μ = cos(θ).

This function is easier to understand, but about 5x slower than `translate!`.
"""
function shift_variable!(var::Array, height::Array{<:Unitful.Length, 1},
                         pixel_size::Unitful.Length, μ::Real)
    θ = acos(μ)
    shift_pix = uconvert.(Unitful.NoUnits, height .* tan(θ) ./ pixel_size)
    for i=1:length(height)
        var[i,:, :] .= shift_image(var[i, :, :], shift_pix[i], 0.)
    end
end



"""
    translate!(data::Array, height::Array{<:Unitful.Length, 1},
                    pixel_size::Unitful.Length, μ::Real, φ::Real)

Shift (or translate) a 3D array that is horizontally periodic in the first two dimensions
according to a polar angle θ given by μ = cos(θ), and an azimuthal angle φ given
in radians. Uses cubic splines, modifies data in-place.

Adapted from the fortran routine `trtnslt` by Åke Nordlund and Bob Stein.

Valid
μ ∈ (0, 1)
φ ∈
"""
function translate!(data::Array, height::Array{<:Unitful.Length, 1},
                    pixel_size::Unitful.Length, μ::Real, φ::Real)
    tanθ = sin(acos(μ)) / μ
    dxdz = tanθ * cos(φ)
    dydz = tanθ * sin(φ)

    dx = ustrip(pixel_size)
    z = ustrip(height)
    height_pix = uconvert.(Unitful.NoUnits, height / pixel_size)

    tmp = similar(data[1, :, :])
    eps = 1.e-6
    nz, nx, ny = size(data)

    if abs(dxdz) > eps
        for n=1:nz
            xk = dxdz * height_pix[n] / nx
            xk = mod(xk, 1.)
            if xk < 0
                xk += 1.
            end
            xk = nx * xk
            k = floor(Int, xk)
            p = xk - k
            k += nx
            q = 1. - p
            af = q + p * q * (q-p)
            bf = p - p * q * (q-p)
            ad = p * q * q * 0.5
            bd = -p * q * p * 0.5
            ac = af - bd
            bc = bf + ad

            for m=1:ny, l=1:nx
                tmp[l, m] = data[n, l, m]
            end

            for l=1:nx
                lm1 = mod(l + k - 2, nx) + 1
                lp0 = mod(l + k - 1, nx) + 1
                lp1 = mod(l + k, nx) + 1
                lp2 = mod(l + k + 1, nx) + 1
                for m=1:ny
                    data[n,l,m] = ac*tmp[lp0,m] + bc*tmp[lp1,m] - ad*tmp[lm1,m] + bd*tmp[lp2,m]
                end
            end
        end
    end

    if abs(dydz) > eps
        for n=1:nz
            yk = dydz * height_pix[n] / ny
            yk = mod(yk, 1.)
            if yk < 0
                yk += 1.
            end
            yk = ny * yk
            k = floor(Int, yk)
            p = yk - k
            k += ny
            q = 1. - p
            af = q + p * q * (q-p)
            bf = p - p * q * (q-p)
            ad = p * q * q * 0.5
            bd = -p * q * p * 0.5
            ac = af - bd
            bc = bf + ad

            for m=1:ny, l=1:nx
                tmp[l, m] = data[n, l, m]
            end

            for m=1:ny
                mm1 = mod(m + k - 2, ny) + 1
                mp0 = mod(m + k - 1, ny) + 1
                mp1 = mod(m + k, ny) + 1
                mp2 = mod(m + k + 1, ny) + 1
                for l=1:nx
                    data[n,l,m] = ac*tmp[l,mp0] + bc*tmp[l,mp1] - ad*tmp[l,mm1] + bd*tmp[l,mp2]
                end
            end
        end
    end
end
