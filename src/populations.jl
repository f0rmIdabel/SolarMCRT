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
    nλ = atom.nλ
    nz,nx,ny = size(atmosphere.temperature)

    λ = atom.λ
    J = []
    for t=1:length(λ)
        nλ = length(λ[t])
        j = zeros(Float64,nλ,nz,nx,ny)u"J/s/nm/m^2/sr"
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


function get_revised_populations(rates::TransitionRates, atom_density::Array{<:NumberDensity,3})

    P = rates.R .+ rates.C
    n_levels = size(P)[1] - 1
    nz,nx,ny = size(atom_density)

    A = Array{Float64, 5}(undef, n_levels, n_levels, nz, nx, ny)u"s^-1"
    b = Array{Float64, 4}(undef, n_levels, nz, nx, ny)u"s^-1 * m^-3"
    populations = Array{Float64, 4}(undef, nz, nx, ny, n_levels+1)u"m^-3"

    for r=1:n_levels
        A[r,r,:,:,:] = P[1,r+1,:,:,:]
        for c=setdiff(1:n_levels, r)
            A[c,r,:,:,:] = P[1,r+1,:,:,:] .- P[c+1,r+1,:,:,:]
            A[r,r,:,:,:] .+= P[r+1,c+1,:,:,:]
        end

        b[r,:,:,:] = atom_density .* P[1,r+1,:,:,:]
    end

    for k=1:nz
        for i=1:nx
            for j=1:ny
                populations[k,i,j,2:end] .= inv(A[:,:,k,i,j]) * b[:,k,i,j]
            end
        end
    end

    populations[:,:,:,1] = atom_density .- sum(populations[:,:,:,2:end],dims=4)[:,:,:,1]

    return populations

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


"""
    write_to_file(populations::Array{<:NumberDensity,4},
                  iteration::Int64,
                  output_path::String)

Write the populations for a given iteration to the output file.
"""
function write_to_file(error::Float64,
                       iteration::Int64,
                       output_path::String)
    h5open(output_path, "r+") do file
        file["populations"][iteration+1] = ustrip.(populations)
    end
end



"""
    get_revised_populations(rates::TransitionRates,
                            atom_density::Array{<:NumberDensity, 3},
                            iteration::Int64,
                            output_path::String)

Calculate the population distribution using statistical equlibrium.
"""
function get_revised_populations2(rates::TransitionRates,
                                 atom_density::Array{<:NumberDensity, 3})


    # Transition probabilities
    P = rates.R+ rates.C

    nz, nx, ny = size(atom_density)
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
