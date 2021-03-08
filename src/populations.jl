include("rates.jl")

function collect_initial_populations()

    pop = h5open(get_initial_populations_path(), "r")
    initial_populations = read(pop, "populations")u"m^-3"
    close(pop)

    return initial_populations
end

function check_population_convergence(populations::Array{<:NumberDensity, 4},
                                      new_populations::Array{<:NumberDensity, 4},
                                      n::Integer,
                                      criterion::Real = 1e-3)
    N = length(populations)
    error = sum( abs.(populations .- new_populations) ./populations ) ./N
    write_error(n, error)

    converged = false
    if error < criterion
        converged = true
    end

    return converged
end

function get_revised_populations(atom::Atom,
                                 rates::TransitionRates,
                                 populations::Array{<:NumberDensity, 4})

    # Transition probabilities
    P12 = rates.R12 .+ rates.C12
    P13 = rates.R13 .+ rates.C13
    P23 = rates.R23 .+ rates.C23
    P21 = rates.R21 .+ rates.C21
    P31 = rates.R31 .+ rates.C31
    P32 = rates.R32 .+ rates.C32

    nz, nx, ny, nl = size(populations)
    revised_populations = Array{Float64, 4}(undef, nz, nx, ny, nl)u"m^-3"

    atom_density = sum(populations, dims=4)[:,:,:,1]
    revised_populations[:,:,:,3] = n3(atom_density, P12, P13, P21, P23, P31, P32)
    revised_populations[:,:,:,2] = n2(atom_density, revised_populations[:,:,:,3], P12, P21, P23, P32)
    revised_populations[:,:,:,1] = n1(atom_density, revised_populations[:,:,:,3], revised_populations[:,:,:,2])

    @test all( Inf .> ustrip.(revised_populations) .> 0.0 )

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
    h5open("../out/output.h5", "r+") do file
        write(file, "populations", ustrip(populations))
    end
end

function write_error(n, error)

    if n == 1
        file = h5open("../out/output_error.h5", "w")
        write(file, "error", Array{Float64,1}(undef,get_max_iterations()))
        file["error"][1] = error
        close(file)
    else
        file = h5open("../out/output_error.h5", "r+")
        file["error"][n] = error
        close(file)
    end
end
