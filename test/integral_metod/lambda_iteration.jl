include("feautrier.jl")
include("../../src/radiation.jl")
using LinearAlgebra

"""
From Tiago, translated from Python to Julia
Computes S and J from a lambda iteration.
"""
function lambda_iteration(atmosphere::Atmosphere, radiation::Radiation, nμ=3, nϕ=4, max_iterations=100)
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
    λ = radiation.λ
    α = radiation.α
    ε = radiation.ε

    nλ, nz, nx, ny = size(α)

    # ==================================================================
    # SET UP TIMER
    # ==================================================================

    time = Array{Float64,2}(undef, max_iterations, nλ)
    error = Array{Float64,2}(undef, max_iterations, nλ)

    # ===================================================================
    # CALCULATE BB SOURCE FUNCTION
    # ===================================================================

    J = Array{Float64,4}(undef, nλ, nz, nx, ny)u"kW / m^2 / sr / nm"
    B = Array{Float64,4}(undef, nλ, nz, nx, ny)u"kW / m^2 / sr / nm"
    S = Array{Float64,4}(undef, nλ, nz, nx, ny)u"kW / m^2 / sr / nm"

    Threads.@threads for l=1:nλ
        B[l,:,:,:] = blackbody_lambda.(λ[l], temperature)
    end

    S = copy(B)
    println("--Starting λ-iteration.....................\n")
    for n=1:max_iterations
        print("--Iteration ", n, "..............................")

        Threads.@threads for l=1:nλ
            f = @timed feautrier(S[l,:,:,:], α[l,:,:,:], z, nμ, nϕ, pixel_size)
            J[l,:,:,:] = f.value
            time[n,l] = f.time
        end

        S_new = (1.0 .- ε) .* J + ε .* B

        converged = check_converging(S, S_, error, n)

        if converged
            println("--Convergence at iteration n = ", n, ". λ-iteration finished.")
            S = S_new
            time = time[1:n,:]
            error = error[1:n,:]
            break
        else
            S = S_new
        end
    end

    # ==================================================================
    # WRITE TO FILE
    # ==================================================================
    out = h5open("../../out/output_integral.h5", "w")
    write(out, "lambda", ustrip(λ))
    write(out, "J", ustrip(J))
    write(out, "S", ustrip(S))
    write(out, "B", ustrip(B))
    write(out, "time", time)
    write(out, "error", error)
    close(out)
end

"""
Check if the relative difference between two arrays
S1 and S2 is smaller than a given criterion.
"""
function check_converged(S1, S2, error, n, criterion = 1e-4)
    nλ, nz, nx, ny = size(S1)
    err = Array{Float64, 1}(undef, nλ)
    for l=1:nλ
        err[l] = norm( abs.(S1[l,:,:,:] .- S2[l,:,:,:]) ./S1[l,:,:,:] )
    end

    error[n,:] = err
    println(@sprintf("Relative error = %.2e.", maximum(err)))

    converged = false

    if maximum(err) < criterion
        converged = true
    end

    return converged
end
