include("../../src/radiation.jl")
include("shift_tools.jl")
using FastGaussQuadrature

"""
1.5D feautrier calculation.
"""
function feautrier(atmosphere::Atmosphere, λ::Unitful.Length, nμ::Int64, nϕ::Int64)

    # ==================================================================
    # ATMOSPHERE DATA
    # ==================================================================
    x = atmosphere.x
    z = atmosphere.z
    χ = atmosphere.χ
    temperature = atmosphere.temperature

    # ===================================================================
    # CALCULATE BB SOURCE FUNCTION
    # ===================================================================
    S = blackbody_lambda.(λ, temperature)

    # ===================================================================
    # SET UP VARIABLES
    # ===================================================================
    μ, w = gausslegendre(nμ)
    μ = μ ./2.0 .+ 0.5
    ϕ = range(0, step=2π/nϕ, length=nϕ)

    pixel_size = abs(x[2] - x[1])
    nz, nx, ny = size(χ)

    D = Array{Float64,3}(undef,nz-1,nx,ny)
    E = Array{Float64,3}(undef,nz-1,nx,ny)u"kW / m^2 / sr / nm"
    P = zeros(nz,nx,ny)u"kW / m^2 / sr / nm"
    J = zeros(nz,nx,ny)u"kW / m^2 / sr / nm"

    # ==================================================================
    # FEAUTRIER
    # ==================================================================
    Threads.@threads for direction=1:nμ*nϕ
        m = 1 + (direction - 1)÷nϕ
        p = direction - (m-1)*nϕ

        S_ = copy(S)
        χ_ = copy(χ)

        println(m, " ", p)

        # Shift atmosphere
        translate!(S_, z[1:end-1], pixel_size, μ[m], ϕ[p])
        translate!(χ_, z[1:end-1], pixel_size, μ[m], ϕ[p])
        τ = optical_depth(χ_, z) / μ[m]

        # Feautrier
        P[end,:,:] = forward(D, E, S_, τ, μ[m])
        backward(P, D, E)

        # Shift back
        translate!(P, z[1:end-1], pixel_size, -μ[m], -ϕ[p])

        # Add contribution to radiation field
        J += w[m]*P/nϕ

        #J = J .+ w[m]*P/nϕ
    end

    # ==================================================================
    # WRITE TO FILE
    # ==================================================================
    out = h5open("../../out/output.hdf5", "cw")
    write(out, "J_feautrier", ustrip(J))
    close(out)
end

"""
Forward-propagation to find the coefficients.
"""
function forward(D::Array{Float64, 3},
                 E::Array{<:Unitful.Quantity,3},
                 S::Array{<:Unitful.Quantity,3},
                 τ::Array{Float64, 3},
                 μ::Float64)

    nz, nx, ny =  size(τ)
    P_end = Array{Float64,2}(undef,nx,ny)u"kW / m^2 / sr / nm"

    for j=1:ny
        for i=1:nx
            Δτ = τ[2:end,i,j] .- τ[1:end-1, i,j]

            # From boundary condition at the top
            E[1,i,j] = 0.0u"kW / m^2 / sr / nm"
            D[1,i,j] = 1.0/(Δτ[1]/μ + 1.0)

            #forward
            for k=2:length(Δτ)
                A = 2μ^2 / (Δτ[k-1]*(Δτ[k-1] + Δτ[k]))
                B = 1.0 + 2μ^2 / (Δτ[k]*Δτ[k-1])            #should use steins trick here
                C = 2μ^2 /(Δτ[k]*(Δτ[k-1] + Δτ[k]))

                D[i] = C / (B - A*D[k-1,i,j])
                E[i] = (S[k,i,j] + A*E[k-1,i,j]) / (B - A*D[k-1,i,j])
            end

            # From boundary condition at the bottom
            P_end[i,j] = ( E[end,i,j] + Δτ[end]/μ * S[end,i,j] + (S[end,i,j] - S[end-1,i,j]) ) / (1.0 - D[end,i,j] + Δτ[end]/μ)
        end
    end

    return P_end
end

"""
Back-propagation to find the P.
"""
function backward(P::Array{<:Unitful.Quantity,3},
                  D::Array{Float64, 3},
                  E::Array{<:Unitful.Quantity,3})

    nz, nx, ny = size(D)

    for j=1:ny
        for i=1:nz
            for k in range(nz, step=-1,stop=1)
                P[k,i,j] = D[k,i,j]*P[k+1,i,j] + E[k,i,j]
            end
        end
    end
end
