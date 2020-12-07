include("feautrier.jl")
include("../../src/radiation.jl")
using LinearAlgebra

"""
From Tiago, translated from Python to Julia
Computes S and J from a lambda iteration.
"""
function lambda_iteration(atmosphere::Atmosphere, λ::Unitful.Length, nμ=3, nϕ=4, max_iterations=1000)
    # ==================================================================
    # ATMOSPHERE DATA
    # ==================================================================
    x = atmosphere.x
    z = atmosphere.z
    χ = atmosphere.χ
    ε = atmosphere.ε
    temperature = atmosphere.temperature

    # ===================================================================
    # CALCULATE BB SOURCE FUNCTION
    # ===================================================================
    B = blackbody_lambda.(λ, temperature)
    S = copy(B)
    nz, nx, ny = size(χ)
    pixel_size = abs(x[2] - x[1])

    J = Array{Float64,3}(undef, nz, nx, ny)

    for n=1:max_iterations
        println("\nIteration ", n)
        J = feautrier(S, χ, z, nμ, nϕ, pixel_size)
        Snew = (1.0 .- ε) .* J + ε .* B

        if convergence(S, Snew, 1e-4)
            S = copy(Snew)
            println("Convergence at iteration n = ", n)
            break
        else
            S = copy(Snew)
        end
    end

    # ==================================================================
    # WRITE TO FILE
    # ==================================================================
    out = h5open("../../out/output.hdf5", "cw")
    write(out, "J_integral", ustrip(J))
    write(out, "S_integral", ustrip(S))
    close(out)
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
