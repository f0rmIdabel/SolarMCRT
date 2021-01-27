include("feautrier.jl")
include("../../src/radiation.jl")
using LinearAlgebra

"""
From Tiago, translated from Python to Julia
Computes S and J from a lambda iteration.
"""
function lambda_iteration(atmosphere::Atmosphere, radiation::Radiation, nμ=3, nϕ=4, max_iterations=1000)
    # ==================================================================
    # ATMOSPHERE DATA
    # ==================================================================
    x = atmosphere.x
    z = atmosphere.z
    temperature = atmosphere.temperature
    pixel_size = abs(x[2] - x[1])

    # ==================================================================
    # RADIATION DATA
    # ==================================================================
    λ = radiation.λ[1]
    χ = radiation.χ[1,:,:,:]
    ε = radiation.ε[1,:,:,:]

    nλ, nz, nx, ny = size(χ)

    # ===================================================================
    # CALCULATE BB SOURCE FUNCTION
    # ===================================================================

    J = Array{Float64,3}(undef, nλ, nz, nx, ny)
    B = Array{Float64,3}(undef, nλ, nz, nx, ny)
    S = Array{Float64,3}(undef, nλ, nz, nx, ny)

    Threads.@threads for l=1:nλ
        B[l,:,:,:] = blackbody_lambda.(λ[l], temperature)
    end

    S = copy(B)

    for n=1:max_iterations
        println("\nIteration ", n)

        Threads.@threads for l=1:nλ
            J[l,:,:,:] = feautrier(S[l,:,:,:], χ[l,:,:,:], z, nμ, nϕ, pixel_size)
        end

        S_new = (1.0 - ε) .* J + ε .* B

        if convergence(S, S_new, 1e-3)
            println("Convergence at iteration n = ", n)
            S = S_new
            iterations = n
            break
        else
            S = S_new
        end
    end

    # ==================================================================
    # WRITE TO FILE
    # ==================================================================
    if iterations < max_iterations
        out = h5open("../../out/output_integral.h5", "cw")
        write(out, "lambda", ustrip(λ))
        write(out, "J", ustrip(J))
        write(out, "S", ustrip(S))
        write(out, "B", ustrip(B))
        write(out, "iterations", n)
        close(out)
    end
end

"""
Check if the relative difference between two arrays
S1 and S2 is smaller than a given criterion.
"""
function convergence(S1, S2, criterion = 1e-4)

    converge = false
    c = norm( abs.(S1 .- S2) ./S1 )

    if c < criterion
        converge = true
    else
        println("c = ", c)
    end

    return converge
end
