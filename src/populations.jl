include("rates.jl")

function collect_initial_populations(atom::Atom,
                                     temperature::Array{<:Unitful.Temperature, 3},
                                     electron_density::Array{<:NumberDensity,3})

    density = atom.density
    population_mode = get_population_distribution()

    if population_mode == "LTE"
        initial_populations = LTE_populations(atom, temperature, electron_density)
    elseif population_mode == "zero_radiation"
        initial_populations = zero_radiation_populations(atom, temperature, electron_density)
    end

    return initial_populations
end


function zero_radiation_populations(atom::Atom,
                                    temperature::Array{<:Unitful.Temperature, 3},
                                    electron_density::Array{<:NumberDensity,3})
    nλ = length(atom.λ)
    nz,nx,ny = size(temperature)

    J0 = zeros(Float64,nλ,nz,nx,ny)u"J/s/nm/m^2/sr"
    rates0 = TransitionRates(calculate_transition_rates(atom, temperature, electron_density, J0)...)
    populations = get_revised_populations(rates0, atom.density)
    return populations
end


function check_population_convergence(populations::Array{<:NumberDensity, 4},
                                      new_populations::Array{<:NumberDensity, 4},
                                      criterion::Real = 1e-3)
    N = length(populations)
    error = sum( abs.(populations .- new_populations) ./populations ) ./N

    converged = false
    if error < criterion
        converged = true
    end

    return converged, error
end


function get_revised_populations(rates::TransitionRates,
                                 atom_density::Array{<:NumberDensity, 3},
                                 iteration::Int64,
                                 output_path::String)

    # Transition probabilities
    P12 = rates.R12 .+ rates.C12
    P13 = rates.R13 .+ rates.C13
    P23 = rates.R23 .+ rates.C23
    P21 = rates.R21 .+ rates.C21
    P31 = rates.R31 .+ rates.C31
    P32 = rates.R32 .+ rates.C32

    nz, nx, ny = size(P12)
    revised_populations = Array{Float64, 4}(undef, nz, nx, ny, 3)u"m^-3"

    revised_populations[:,:,:,3] = n3(atom_density, P12, P13, P21, P23, P31, P32)
    revised_populations[:,:,:,2] = n2(atom_density, revised_populations[:,:,:,3], P12, P21, P23, P32)
    revised_populations[:,:,:,1] = n1(atom_density, revised_populations[:,:,:,3], revised_populations[:,:,:,2])


    @test all( sum(revised_populations, dims=4)[:,:,:,1] .≈ atom_density )
    @test all( Inf .> ustrip.(revised_populations[:,:,:,3]) .>= 0.0 )
    @test all( Inf .> ustrip.(revised_populations[:,:,:,2]) .>= 0.0 )
    @test all( Inf .> ustrip.(revised_populations[:,:,:,1]) .>= 0.0 )

    write_to_file(revised_populations, iteration, output_path)

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


function write_to_file(populations::Array{<:NumberDensity,4},
                       iteration::Int64,
                       output_path::String)
    h5open(output_path, "r+") do file
        file["populations"][iteration+1,:,:,:,:] = ustrip.(populations)
    end
end
