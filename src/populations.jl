include("rates.jl")

"""
    collect_initial_populations(atom::Atom,
                                temperature::Array{<:Unitful.Temperature, 3},
                                electron_density::Array{<:NumberDensity,3},
                                population_mode=nothing)

Collect initial population distribution. Either LTE or zero-radiation.
"""
function collect_initial_populations(atmosphere,
                                     atom::Atom)
                                     #temperature::Array{<:Unitful.Temperature, 3},
                                     #electron_density::Array{<:NumberDensity,3})

    temperature = atmosphere.temperature
    electron_density = atmosphere.electron_density
    density = atom.density
    population_mode = get_population_distribution()

    if population_mode == "LTE"
        initial_populations = LTE_populations(atom, temperature, electron_density)
    elseif population_mode == "zero_radiation"
        initial_populations = zero_radiation_populations(atmosphere, atom)
    end

    return initial_populations
end

"""
    zero_radiation_populations(atom::Atom,
                               temperature::Array{<:Unitful.Temperature, 3},
                               electron_density::Array{<:NumberDensity,3})

For a given atom density, calculate the populations according to zero-radiation.
"""
function zero_radiation_populations(atmosphere::Atmosphere, atom::Atom)
                                    #temperature::Array{<:Unitful.Temperature, 3},
                                    #electron_density::Array{<:NumberDensity,3})
    nλ = atom.nλ
    nz,nx,ny = size(atmosphere.temperature)

    λ = atom.λ
    J = []
    for t=1:length(λ)
        nλ = length(λ[t])
        j = zeros(Float64,nλ,nz,nx,ny)u"J/s/nm/m^2/sr"
        #j = zeros(UnitsIntensity_λ, nλ,nz,nx,ny)#u"J/s/nm/m^2/sr"
        append!(J, [j])
    end

    zero_rates = TransitionRates(calculate_transition_rates(atmosphere, atom, J)...)
    populations = get_revised_populations(zero_rates, atom.density)

    @test all( Inf .> ustrip.(populations) .>= 0.0 )

    return populations
end

"""
    check_population_convergence(populations::Array{<:NumberDensity, 4},
                                 new_populations::Array{<:NumberDensity, 4},
                                 criterion::Real = 1e-3)
Check if the relative difference between two populations satisfy
the convergence criterion. Return convergence status and error.
"""
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

"""
    get_revised_populations(rates::TransitionRates,
                            atom_density::Array{<:NumberDensity, 3},
                            iteration::Int64,
                            output_path::String)

Calculate the population distribution using statistical equlibrium.
"""
function get_revised_populations(rates::TransitionRates,
                                 atom_density::Array{<:NumberDensity, 3})


    # Transition probabilities
    P12 = rates.Rlu[3,:,:,:] + rates.Clu[3,:,:,:]
    P13 = rates.Rlu[1,:,:,:] + rates.Clu[1,:,:,:]
    P23 = rates.Rlu[2,:,:,:] + rates.Clu[2,:,:,:]

    P21 = rates.Rul[3,:,:,:] + rates.Cul[3,:,:,:]
    P31 = rates.Rul[1,:,:,:] + rates.Cul[1,:,:,:]
    P32 = rates.Rul[2,:,:,:] + rates.Cul[2,:,:,:]


    nz, nx, ny = size(P12)
    revised_populations = Array{Float64, 4}(undef, nz, nx, ny, 3)u"m^-3"

    revised_populations[:,:,:,3] = n3(atom_density, P12, P13, P21, P23, P31, P32)
    revised_populations[:,:,:,2] = n2(atom_density, revised_populations[:,:,:,3], P12, P21, P23, P32)
    revised_populations[:,:,:,1] = n1(atom_density, revised_populations[:,:,:,3], revised_populations[:,:,:,2])

    return revised_populations
end

"""
    n3(N::Array{<:NumberDensity,3},
       P12::Array{<:Unitful.Frequency,3},
       P13::Array{<:Unitful.Frequency,3},
       P21::Array{<:Unitful.Frequency,3},
       P23::Array{<:Unitful.Frequency,3},
       P31::Array{<:Unitful.Frequency,3},
       P32::Array{<:Unitful.Frequency,3})

Given the atom density and se    qcore, qwing = get_qvalues()
lected transition rates,
calculate the ionised population density.
"""
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

"""
    n2(N::Array{<:NumberDensity,3},
       n3::Array{<:NumberDensity,3},
       P12::Array{<:Unitful.Frequency,3},
       P21::Array{<:Unitful.Frequency,3},
       P23::Array{<:Unitful.Frequency,3},
       P32::Array{<:Unitful.Frequency,3})

Given the atom density, ionised density and selected rates,
calculate the excited level population.
"""
function n2(N::Array{<:NumberDensity,3},
            n3::Array{<:NumberDensity,3},
            P12::Array{<:Unitful.Frequency,3},
            P21::Array{<:Unitful.Frequency,3},
            P23::Array{<:Unitful.Frequency,3},
            P32::Array{<:Unitful.Frequency,3})
    return (N .* P12  + n3 .* (P32 .- P12) ) ./ (P21 .+ P23 .+ P12)
end

"""
    n1(N::Array{<:NumberDensity,3},
       n3::Array{<:NumberDensity,3},
       n2::Array{<:NumberDensity,3})

Given the atom density, ionised density and excited density,
get the remaining ground density.
"""
function n1(N::Array{<:NumberDensity,3},
            n3::Array{<:NumberDensity,3},
            n2::Array{<:NumberDensity,3})
    return N .- n3 .- n2
end

"""
    write_to_file(populations::Array{<:NumberDensity,4},
                  iteration::Int64,
                  output_path::String)

Write the populations for a given iteration to the output file.
"""
function write_to_file(populations::Array{<:NumberDensity,4},
                       iteration::Int64,
                       output_path::String)
    h5open(output_path, "r+") do file
        file["populations"][iteration+1,:,:,:,:] = ustrip.(populations)
    end
end
