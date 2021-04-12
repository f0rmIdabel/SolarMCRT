include("rates.jl")

"""
    collect_initial_populations(atom::Atom,
                                temperature::Array{<:Unitful.Temperature, 3},
                                electron_density::Array{<:NumberDensity,3},
                                population_mode=nothing)

Collect initial population distribution. Either LTE or zero-radiation.
"""
function collect_initial_populations(atmosphere::Atmosphere,
                                     atom::Atom,
                                     distribution::String)

    temperature = atmosphere.temperature
    electron_density = atmosphere.electron_density
    density = atom.density

    if distribution == "LTE"
        initial_populations = LTE_populations(atom, temperature, electron_density)
    elseif distribution == "zero_radiation"
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


function get_revised_populations(rates::TransitionRates, atom_density::Array{<:NumberDensity,3})

    P = rates.R .+ rates.C
    n_levels = size(P)[1] - 1
    nz,nx,ny = size(atom_density)

    A = Array{Float64, 5}(undef, n_levels, n_levels, nz, nx, ny)u"s^-1"
    b = Array{Float64, 4}(undef, n_levels, nz, nx, ny)u"s^-1 * m^-3"
    populations = Array{Float64, 4}(undef, nz, nx, ny, n_levels+1)u"m^-3"

    for r=1:n_levels
        A[r,r,:,:,:] = P[1,r+1,:,:,:] .+ P[r+1,1,:,:,:]
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
