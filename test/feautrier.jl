include("../src/atmosphere.jl")
include("../src/radiation.jl")
using Unitful
using Plots

function feautrier()
    λ = 500u"nm"
    μ = [-0.9815606342467192506906  -0.9041172563704748566785  -0.769902674194304687037
         -0.5873179542866174472967  -0.3678314989981801937527  -0.1252334085114689154724
          0.1252334085114689154724   0.3678314989981801937527   0.5873179542866174472967
          0.7699026741943046870369   0.9041172563704748566785   0.9815606342467192506906] ./2.0 .+ 0.5

    atmosphere_parameters = collect_atmosphere_data(λ)
    atmosphere = Atmosphere(atmosphere_parameters...)

    χ = atmosphere.χ
    z = atmosphere.z
    temperature = atmosphere.temperature

    τ = optical_depth(χ, z)
    S = blackbody_lambda.(λ, temperature)

    nz, nx, ny = size(τ)
    nμ = length(μ)

    P = Array{Float64,1}(undef,nz)
    D = Array{Float64,1}(undef,nz-1)
    E = Array{Float64,1}(undef,nz-1) * u"kW / m^2 / sr / nm"
    P = Array{Float64,1}(undef,nz) * u"kW / m^2 / sr / nm"
    Pμ = Array{Float64,2}(undef,nμ, nz) * u"kW / m^2 / sr / nm"
    J = Array{Float64,3}(undef,nz, nx, ny) * u"kW / m^2 / sr / nm"
    for j=1:ny
        for i=1:nx
            for m=1:nμ

                P[end] = forward(D, E, S[:,i,j], τ[:,i,j], μ[m])
                backward(P, D, E)

                Pμ[m,:] = P
            end
            J[:,i,j] .= quadrature(Pμ, μ)
        end
    end

    surface = ustrip(J[1,:,:])
    # To avoid ssh display problems
    ENV["GKSwstype"]="nul"
    heatmap(1:size(surface,1), 1:size(surface,2), surface, c=:grays, aspect_ratio=:equal)
    plot!(size=(410,400))
    fig = "/mn/stornext/u3/idarhan/MScProject/Analysis/feautrier"
    png(fig)


end

function forward(D, E, S, τ, μ)
    Δτ = τ[2:end] .- τ[1:end-1] #nz - 1

    # From boundary at the top
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
    # From boundary condition at bottom
    P_end = ( E[end] + Δτ[end]/μ * S[end] + (S[end] - S[end-1]) ) / (1.0 - D[end] + Δτ[end]/μ)

    return P_end
end

function backward(P, D, E)
    n = length(D)
    for i in range(n, step=-1,stop=1)
        P[i] = D[i]*P[i+1] + E[i]
    end
end

function quadrature(P, μ)

    nμ, nz = size(P)
    J = zeros(nz) * u"kW / m^2 / sr / nm"

    w = [0.0471753363865118271946  0.1069393259953184309603  0.1600783285433462263347  0.2031674267230659217491
         0.233492536538354808761   0.2491470458134027850006  0.2491470458134027850006  0.233492536538354808761
         0.203167426723065921749   0.160078328543346226335   0.1069393259953184309603  0.0471753363865118271946] ./2

    for i=1:nμ
        J = J .+ w[i]*P[i,:]
    end

    return J

end

feautrier()
