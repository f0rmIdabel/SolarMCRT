#include("atom.jl")
#include("rates.jl")

function collect_initial_populations()

    pop = h5open(get_initial_populations_path(), "r")
    initial_populations = read(pop, "initial_populations")u"m^-3"
    close(pop)

    return initial_populations
end

function LTE_populations(atom::Atom,
                         temperature::Array{<:Unitful.Temperature, 3},
                         electron_density::Array{<:NumberDensity, 3})

    χl = atom.χl
    χu = atom.χu
    χ∞ = atom.χ∞
    gl = atom.gl
    gu = atom.gu
    g∞ = atom.g∞

    atom_density = sum(atom.populations, dims=4)[:,:,:,1]

    nz, nx, ny = size(temperature)
    populations = Array{Float64, 4}(undef, nz, nx, ny, 3)u"m^-3"

    C = 2π*m_e*k_B/h^2
    U1 = gl * exp.(-χl/k_B./temperature)
    U2 = gu * exp.(-χu/k_B./temperature)
    U3 = g∞ * exp.(-χ∞/k_B./temperature)

    K = 1 ./electron_density .* 2 .* U3 .^2 ./ (U1 .+ U2) / g∞ .* (C*temperature).^(1.5)

    populations[:,:,:,3] = K .* atom_density ./ (1.0 .+ K)
    populations[:,:,:,1] = (atom_density .- populations[:,:,:,3]) ./ (1.0 .+ U2 ./ U1)
    populations[:,:,:,2] = atom_density .- populations[:,:,:,1] .- populations[:,:,:,3]

    return populations
end

function check_population_convergence(populations::Array{<:PerLength, 4},
                                      new_populations::Array{<:PerLength, 4},
                                      criterion::Real = 1e-3)
    N = length(populations)
    error = norm( abs.(populations .- new_population) ./populations ./N)

    println(@sprintf("--Relative error = %.2e.", err))

    converged = false

    if error < criterion
        converged = true
    end

    return converged
end

function get_revised_populations(atom::Atom,
                                 rates::TransitionRates)

    # Transition probabilities
    P12 = rates.R12 .+ rates.C12
    P13 = rates.R13 .+ rates.C13
    P23 = rates.R23 .+ rates.C23
    P21 = rates.R21 .+ rates.C21
    P31 = rates.R31 .+ rates.C31
    P32 = rates.R32 .+ rates.C32

    # Revised populations
    atom_density = sum(atom.populations, dims=4)[:,:,:,1]
    revised_populations[:,:,:,3] = n3(atom_density, P12, P13, P21, P23, P31, P32)
    revised_populations[:,:,:,2] = n2(atom_density, revised_populations[:,:,:,3], P12, P21, P23, P32)
    revised_populations[:,:,:,1] = n1(atom_density, revised_populations[:,:,:,3], revised_populations[:,:,:,2])

    return revised_populations
end

function n3(N::Array{<:NumberDensity,3},
            P12::Array{<:Unitful.Frequency,3},
            P13::Array{<:Unitful.Frequency,3},
            P21::Array{<:Unitful.Frequency,3},
            P23::Array{<:Unitful.Frequency,3},
            P31::Array{<:Unitful.Frequency,3},
            P32::Array{<:Unitful.Frequency,3})
    a = N .* P12 ./ (P21 .+ P23 .+ P12)
    b = (P32 .- P12) ./ (P21 .+ P23 .+ P12)
    c = N .* (P12 .+ P13) .- a .* (P21 .+ P12 .+ P13)
    d = b .* (P21 .+ P12 .+ P13) .+ P31 .+ P12 .+ P13

    return c ./ d
end

function n2(N::Array{<:NumberDensity,3},
            n3::Array{<:NumberDensity,3},
            P12::Array{<:Unitful.Frequency,3},
            P21::Array{<:Unitful.Frequency,3},
            P23::Array{<:Unitful.Frequency,3},
            P32::Array{<:Unitful.Frequency,3})
    return (N .* P12  + n3 .* (P32 .- P12) ) ./ (P21 .+ P23 .+ P12)
end

function n1(N::Array{<:NumberDensity,3},
            n3::Array{<:NumberDensity,3},
            n2::Array{<:NumberDensity,3})
    return N .- n3 .- n2
end

function write_to_file(populations::Array{<:NumberDensity,4})
    h5open("../out/output.h5", "w") do file
        write(file, "populations", ustrip(populations))
    end
end
