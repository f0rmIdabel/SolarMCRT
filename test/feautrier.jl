include("../src/radiation.jl")
include("shift_tools.jl")
"""
1.5D feautrier calculation.
"""
function feautrier(atmosphere::Atmosphere, λ::Unitful.Length)

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

    pixel_size = abs(x[2] - x[1])

    # Using 12 angles
    μ = [-0.7745966692414833770359, 0, 0.7745966692414833770359] ./2.0 .+ 0.5
    w = [0.5555555555555555555556, 0.8888888888888888888889, 0.555555555555555555556]

    ϕ = [0, π/2, π, 3/2π]

    nz, nx, ny = size(χ)
    nμ = length(μ)
    nϕ = length(ϕ)

    P = zeros(nz,nx,ny)u"kW / m^2 / sr / nm"
    D = Array{Float64,3}(undef,nz-1, nx, ny)
    E = Array{Float64,3}(undef,nz-1, nx, ny)u"kW / m^2 / sr / nm"

    J = zeros(nz, nx, ny)u"kW / m^2 / sr / nm"

    # ==================================================================
    # FEAUTRIER, COLUMN BY COLUMN
    # ==================================================================
    for m=1:nμ

        #fill!(P, 0.0u"kW / m^2 / sr / nm") #can possibly be removed

        for p=1:nϕ
            translate!(S, z, pixel_size, μ[m], ϕ[p])
            translate!(χ, z, pixel_size, μ[m], ϕ[p])
            τ = optical_depth(χ, z) / μ[m]

            P[end,i,j] = forward(D, E, S[:,i,j], τ[:,i,j], μ[m])
            backward(P, D, E)

            # Shift back
            translate!(P, z[1:end-1], pixel_size, -μ[m], -ϕ[p])
        end

        J = J .+ w[m]*P/nϕ
    end

    # ==================================================================
    # WRITE TO FILE
    # ==================================================================
    out = h5open("../out/output.hdf5", "cw")
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

    for j=1:ny
        for i=1:nx
            Δτ = τ[2:end,i,j] .- τ[1:end-1, i,j]

            # From boundary condition at the top
            E[1,i,j] = 0.0u"kW / m^2 / sr / nm"
            D[1,i,j] = 1.0/(Δτ[1]/μ + 1.0)

            #forward
            for k=2:length(Δτ)
                A = 2μ^2 / (Δτ[k-1]*(Δτ[k-1] + Δτ[k]))
                B = 1.0 + 2μ^2 / (Δτ[k]*Δτ[k-1]) #use steins trick
                C = 2μ^2 /(Δτ[k]*(Δτ[k-1] + Δτ[k]))

                D[i] = C / (B - A*D[k-1,i,j])
                E[i] = (S[k,i,j] + A*E[k-1,i,j]) / (B - A*D[k-1,i,j])
            end

            # From boundary condition at the bottom
            P_end = ( E[end,i,j] + Δτ[end]/μ * S[end,i,j] + (S[end,i,j] - S[end-1,i,j]) ) / (1.0 - D[end,i,j] + Δτ[end]/μ)
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

"""
Gauss-Legendre with 12 angles.
"""
function gauss_legendre12(P::Array{<:Unitful.Quantity,2},
                          μ::Array{Float64, 1})

    nμ, nz = size(P)
    J = zeros(nz) * u"kW / m^2 / sr / nm"

    w = [0.0471753363865118271946
         0.1069393259953184309603
         0.1600783285433462263347
         0.2031674267230659217491
         0.233492536538354808761
         0.2491470458134027850006
         0.2491470458134027850006
         0.233492536538354808761
         0.203167426723065921749
         0.160078328543346226335
         0.1069393259953184309603
         0.0471753363865118271946] ./2

    for i=1:nμ
        J = J .+ w[i]*P[i,:]
    end

    return J
end
