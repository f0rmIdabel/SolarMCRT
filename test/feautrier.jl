include("../src/atmosphere.jl")
include("../src/radiation.jl")
using Unitful


"""
THIS DOES NOT WORK. YET...
"""

function feautrier()
    λ = 500u"nm"
    μ = [1.0]

    atmosphere_parameters = collect_atmosphere_data(λ)
    atmosphere = Atmosphere(atmosphere_parameters...)

    χ = atmosphere.χ
    z = atmosphere.z
    boundary = atmosphere.boundary
    temperature = atmosphere.temperature

    nx,ny,nz = size(χ)
    nμ = length(μ)

    nz -= 50
    τ = optical_depth(χ, z, boundary)[:,:,1:nz]
    S = blackbody_lambda.(λ, temperature)[:,:,1:nz]
    P = Array{Float64,1}(undef,nz)

    nx = ny = 1

    D = Array{Float64,1}(undef,nz-1)
    E = Array{Float64,1}(undef,nz-1) * u"kW / m^2 / sr / nm"
    P = Array{Float64,1}(undef,nz) * u"kW / m^2 / sr / nm"
    Pμ = Array{Float64,2}(undef,nμ, nz) * u"kW / m^2 / sr / nm"

    for i=1:nx
        for j=1:ny
            for m=1:nμ
                τ = τ[i,j,:]
                S = S[i,j,:]

                forward(P, D, E, S, τ, μ[m])
                backward(P, D, E)

                Pμ[m,:] = P
            end
            #J[i,j,:] = quadrature(Pμ, μ)
        end
    end
    print(Pμ)
end

function forward(P, D, E, S, τ, μ)
    Δτ = τ[2:end] .- τ[1:end-1]

    E[1] = 0.0 * u"kW / m^2 / sr / nm"
    D[1] = 1.0/(Δτ[1]/μ + 1.0)

    #forward
    for i=2:length(Δτ)
        A = 2μ^2 / (Δτ[i-1]*(Δτ[i-1] + Δτ[i]))
        B = 1.0 + 2μ^2 / (Δτ[i]*Δτ[i-1]) #use steins trick
        C = 2μ^2 /(Δτ[i]*(Δτ[i-1] + Δτ[i]))

        D[i] = C / (B - A*D[i-1])
        E[i] = (S[i] + A*E[i-1]) / (B - A*D[i-1])
    end

    P[end] = ( E[end] + Δτ[end]/μ * S[end] + (S[end] - S[end-1]) ) / (1.0 - D[end] + Δτ[end]/μ)

    println(Δτ)
end

function backward(P, D, E)
    nz = length(D)
    for i=nz-1:1
        P[i] = D[i]*P[i+1] + E[i]
    end
end

function quadrature()

end

"""
Returns 2D array containing the k-indices where the optical depth.
"""

feautrier()
