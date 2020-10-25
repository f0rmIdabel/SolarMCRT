include("../src/atmosphere.jl")
include("../src/radiation.jl")
using Unitful
using Plots

function feautrier()
    λ = 500u"nm"
    μ = [1.0]

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
    nx = ny = 50

    D = Array{Float64,1}(undef,nz-1)
    E = Array{Float64,1}(undef,nz-1) * u"kW / m^2 / sr / nm"
    P = Array{Float64,1}(undef,nz) * u"kW / m^2 / sr / nm"
    #Pμ = Array{Float64,2}(undef,nμ, nz) * u"kW / m^2 / sr / nm"
    J = Array{Float64,3}(undef,nz, nx, ny) * u"kW / m^2 / sr / nm"
    for j=1:ny
        for i=1:nx
            for m=1:nμ

                P[end] = forward(D, E, S[:,i,j], τ[:,i,j], μ[m])
                backward(P, D, E)

                #Pμ[m,:] = P
            end
            #J[:,i,j] = quadrature(Pμ, μ)
            J[:,i,j] = P
        end
        println(j)
    end


    surface = ustrip(J[1,:,:])
    # To avoid ssh display problems
    ENV["GKSwstype"]="nul"
    heatmap(1:size(surface,1), 1:size(surface,2), surface, c=:grays, aspect_ratio=:equal)
    plot!(size=(410,400))
    fig = @sprintf("/mn/stornext/u3/idarhan/MScProject/Analysis/bf_tau%.1f_pcts%.0e_bin%s",
                   τ_max, total_packets, string(bin))
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

    for i=1:nμ
        J += w[i]*P[i]
    end

    -0.5773502691896257645092
    0.5773502691896257645092

end


"""
Returns 2D array containing the k-indices where the optical depth.
"""

feautrier()
