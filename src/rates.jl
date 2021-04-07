include("atmosphere.jl")
include("atom.jl")

struct TransitionRates
    Rlu::Array{<:Unitful.Frequency, 4}
    Rul::Array{<:Unitful.Frequency, 4}
    Clu::Array{<:Unitful.Frequency, 4}
    Cul::Array{<:Unitful.Frequency, 4}
end

"""
    calculate_transition_rates(atom::Atom,
                               temperature::Array{<:Unitful.Temperature, 3},
                               electron_density::Array{<:NumberDensity,3},
                               J::Array{<:UnitsIntensity_λ, 4})

Given the radiation field, calculate all transition rates for
the excitation, de-excitation, ionisations and re-combiantions
of a two level atom. Level 1,2 and 3 represent the ground level,
excited level and ionised level respectively.
"""
function calculate_transition_rates(atmosphere::Atmosphere,
                                    atom::Atom,
                                    J::Array{Any,1})

    # ==================================================================
    # LOAD ATMOSPHERE PARAMETERS
    # ==================================================================
    temperature = atmosphere.temperature
    electron_density = atmosphere.electron_density
    hydrogen_populations = atmosphere.hydrogen_populations

    LTE_pops = LTE_populations(atom, temperature, electron_density)

    nz,nx,ny,nl = size(LTE_pops)
    n_levels = nl - 1
    # ==================================================================
    # LOAD WAVELENGTHS
    # ==================================================================
    λ = atom.λ

    n_levels = atom.n_levels
    n_lines = atom.n_lines

    n_transitions = n_levels + n_lines
    # ==================================================================
    # CALCULATE RADIATIVE RATES
    # ==================================================================

    Rlu = Array{Unitful.Frequency,4}(undef,n_transitions,nz,nx,ny)
    Rul = Array{Unitful.Frequency,4}(undef,n_transitions,nz,nx,ny)
    Clu = Array{Unitful.Frequency,4}(undef,n_transitions,nz,nx,ny)
    Cul = Array{Unitful.Frequency,4}(undef,n_transitions,nz,nx,ny)

    # BF

    for level = 1:n_levels
        σ = σic(level, atom, λ[level])
        G = Gij(level, n_levels+1, λ[level], temperature, LTE_pops)

        Rlu[level,:,:,:] = Rij(J[level], σ, λ[level])
        Rul[level,:,:,:] = Rji(J[level], σ, G, λ[level])

        Clu[level,:,:,:] = Cij(level, n_levels+1, electron_density, temperature, LTE_pops)
        Cul[level,:,:,:] = Cij(n_levels+1, level, electron_density, temperature, LTE_pops)
    end


    line_number = 0
    for l=1:n_levels-1
        for u=(l+1):n_levels

            line_parameters = collect_line_data(atmosphere, atom, u, l)
            line = Line(line_parameters...)

            line_number += 1

            σ = σij(l, u, line, λ[n_levels + line_number])
            G = Gij(l, u, λ[n_levels+line_number], temperature, LTE_pops)

            Rlu[n_levels + line_number,:,:,:] = Rij(J[n_levels+line_number], σ, λ[n_levels+line_number])
            Rul[n_levels + line_number,:,:,:] = Rji(J[n_levels+line_number], σ, G, λ[n_levels+line_number])

            Clu[n_levels + line_number,:,:,:] = Cij(l, u, electron_density, temperature, LTE_pops)
            Cul[n_levels + line_number,:,:,:] = Cij(u, l, electron_density, temperature, LTE_pops)
        end
    end


    return Rlu, Rul, Clu, Cul
end

"""
    Rij(J::Array{<:UnitsIntensity_λ, 4},
        σij::Array{<:Unitful.Area, 4},
        λ::Array{<:Unitful.Length, 1})

Raditive rate for excitation transitions.
"""
function Rij(J::Array{<:UnitsIntensity_λ, 4},
             σij::Array{<:Unitful.Area, 4},
             λ::Array{<:Unitful.Length, 1})

    nλ, nz, nx, ny = size(J)
    R = Array{Unitful.Frequency,3}(undef,nz,nx,ny)
    fill!(R,0.0u"s^-1")

    for l=1:(nλ-1)
        R += 2π/hc * (( λ[l]   * σij[l,:,:,:]   .* J[l,:,:,:] .+
                        λ[l+1] * σij[l+1,:,:,:] .* J[l+1,:,:,:]) .* (λ[l+1] - λ[l]))
    end

    return R
end

"""
    Rij(J::Array{<:UnitsIntensity_λ, 4},
        σij::Array{<:Unitful.Area, 1},
        λ::Array{<:Unitful.Length, 1})

Radiative rate for ionisation transitions.
"""
function Rij(J::Array{<:UnitsIntensity_λ, 4},
             σij::Array{<:Unitful.Area, 1},
             λ::Array{<:Unitful.Length, 1})

    nλ, nz, nx, ny = size(J)
    R = Array{Unitful.Frequency,3}(undef,nz,nx,ny)
    fill!(R,0.0u"s^-1")

    for l=1:(nλ-1)
        R +=  2π/hc * ( λ[l]   * σij[l]   .* J[l,:,:,:]   .+
                        λ[l+1] * σij[l+1] .* J[l+1,:,:,:]  ) .* (λ[l+1] - λ[l])
    end

    return R
end

"""
    Rji(J::Array{<:UnitsIntensity_λ, 4},
        σij::Array{<:Unitful.Area, 4},
        Gij::Array{Float64, 4},
        λ::Array{<:Unitful.Length, 1})

Radiative rate for de-excitation transitions.
"""
function Rji(J::Array{<:UnitsIntensity_λ, 4},
             σij::Array{<:Unitful.Area, 4},
             Gij::Array{Float64, 4},
             λ::Array{<:Unitful.Length, 1})

    nλ, nz, nx, ny = size(J)
    R = Array{Unitful.Frequency,3}(undef,nz,nx,ny)
    fill!(R,0.0u"s^-1")

    # Trapezoid rule
    for l=1:(nλ-1)
        R += 4π/hc * (σij[l,:,:,:]   .* Gij[l,:,:,:]   .* λ[l]   .* (2*h*c_0^2 ./ λ[l]^5   .+ J[l,:,:,:] ) .+
                      σij[l+1,:,:,:] .* Gij[l+1,:,:,:] .* λ[l+1] .* (2*h*c_0^2 ./ λ[l+1]^5 .+ J[l+1,:,:,:] )) .* (λ[l+1] - λ[l])
    end

    return R
end

"""
    Rji(J::Array{<:UnitsIntensity_λ, 4},
        σij::Array{<:Unitful.Area, 1},
        Gij::Array{Float64, 4},
        λ::Array{<:Unitful.Length, 1})

Radiative rate for recombination transitions.
"""
function Rji(J::Array{<:UnitsIntensity_λ, 4},
             σij::Array{<:Unitful.Area, 1},
             Gij::Array{Float64, 4},
             λ::Array{<:Unitful.Length, 1})

    nλ, nz, nx, ny = size(J)
    R = Array{Unitful.Frequency,3}(undef,nz,nx,ny)
    fill!(R,0.0u"s^-1")

    # Trapezoid rule
    for l=1:(nλ-1)
        R += 4π/hc * (σij[l]   .* Gij[l,:,:,:]   .* λ[l]   .* (2*hc*c_0 ./ λ[l]^5   .+ J[l,:,:,:] ) .+
                      σij[l+1] .* Gij[l+1,:,:,:] .* λ[l+1] .* (2*hc*c_0 ./ λ[l+1]^5 .+ J[l+1,:,:,:] )) .* (λ[l+1] - λ[l])
    end

    return R
end

"""
    σij(i::Integer,
        j::Integer,
        atom::Atom,
        λ::Array{<:Unitful.Length, 1})

Calculates the bound-bound cross section.
"""
function σij(i::Integer,
             j::Integer,
             line::Line,
             λ::Array{<:Unitful.Length, 1})

    damping_constant = line.damping_constant
    ΔλD = line.doppler_width
    λ0 = line.lineData.λ0
    Bij = line.lineData.Bij
    σ_constant = h*c_0/(4π*λ0) * Bij
    nλ = length(λ)
    nz, nx, ny = size(ΔλD)
    σ = Array{Unitful.Area, 4}(undef, nλ, nz, nx, ny)

    @Threads.threads for l=1:nλ
        damping = (λ[l]^2 * damping_constant) .|> u"m/m"
        v = (λ[l] - λ0) ./ ΔλD
        profile = voigt_profile.(damping, ustrip(v), ΔλD)
        σ[l,:,:,:] = σ_constant .* profile
    end

    return σ
end

"""
    σic(i::Integer,
        atom::Atom,
        λ::Array{<:Unitful.Length, 1})

Calculates the bound-free cross-section.
"""
function σic(i::Integer,
             atom::Atom,
             λ::Array{<:Unitful.Length, 1})

    λ_edge = λ[end]
    λ3_ratio = (λ ./ λ_edge).^3
    n_eff = sqrt(E_∞ / (atom.χ[end] - atom.χ[i])) |>u"J/J" # should be χu - χl
    charge = atom.Z
    σ_constant = (4 * e^2 / (3 * π * sqrt(3) * ε_0 * m_e * c_0^2 * R_∞)) |> u"m^2"

    σ = (σ_constant * charge^4 * n_eff * λ3_ratio .* gaunt_bf.(λ, charge, n_eff))

    return σ
end

"""
    Gij(i::Integer,
        j::Integer,
        λ::Array{<:Unitful.Length, 1},
        temperature::Array{<:Unitful.Temperature, 3},
        LTE_populations::Array{<:NumberDensity, 4})

Factor to go into the de-excitation and recombination expressions.
"""
function Gij(i::Integer,
             j::Integer,
             λ::Array{<:Unitful.Length, 1},
             temperature::Array{<:Unitful.Temperature, 3},
             LTE_populations::Array{<:NumberDensity, 4})

    nλ = length(λ)
    nz, nx, ny = size(temperature)
    G = Array{Float64, 4}(undef, nλ, nz, nx, ny)

    n_ratio = LTE_populations[:,:,:,i] ./LTE_populations[:,:,:,j]

    @Threads.threads for l=1:nλ
        G[l,:,:,:] =  n_ratio .* exp.(- hc ./ (k_B*λ[l]*temperature))
    end

    return G
end

"""
    Cij(i::Integer,
        j::Integer,
        electron_density::Array{<:NumberDensity,3},
        temperature::Array{<:Unitful.Temperature,3},
        LTE_populations::Array{<:NumberDensity,4})

Calculates the collisional rates for all possible
two-level atom transitions.
"""
function Cij(i::Integer,
             j::Integer,
             electron_density::Array{<:NumberDensity,3},
             temperature::Array{<:Unitful.Temperature,3},
             LTE_populations::Array{<:NumberDensity,4})

    ionisation_level = size(LTE_populations)[end]

    # If UP
    if i < j
        if j < ionisation_level
            C = coll_exc_hydrogen_johnson.(i, j, electron_density, temperature)
        elseif j == ionisation_level
            C = coll_ion_hydrogen_johnson.(i, electron_density, temperature)
        end

    # If DOWN
    elseif i > j
        if i < ionisation_level
            C = coll_exc_hydrogen_johnson.(j, i, electron_density, temperature)
        elseif i == ionisation_level
            C = coll_ion_hydrogen_johnson.(j, electron_density, temperature)
        end
        C = C .* ( LTE_populations[:,:,:,j] ./ LTE_populations[:,:,:,i] )
    end

    return C
end

"""
    gaunt_bf(charge::Int, n_eff::Number, λ::Unitful.Length)::Float64

Compute bound-free Gaunt factor for a given charge, effective principal
quantum number and wavelength λ. Taken from RH. Formula from
[Seaton (1960), Rep. Prog. Phys. 23, 313](https://ui.adsabs.harvard.edu/abs/1960RPPh...23..313S/abstract),
page 316. Recipes from Seaton. (This is copied from github.com/tiagopereira/Transparency.jl)
"""
function gaunt_bf(λ::Unitful.Length,
                  charge::Real,
                  n_eff::Real)::Float64

    x = ustrip(1 / (λ * R_∞ * charge^2) |> u"m/m")
    x3 = x^(1/3)
    nsqx = 1 / (n_eff^2 * x)
    g_bf = 1 + 0.1728 * x3 * (1 - 2 * nsqx) - 0.0496 * x3^2 * (1 - (1 - nsqx) * 0.66666667 * nsqx)
    @assert g_bf >= 0 "gaunt_bf negative, calculation will not be reliable"
    return g_bf
end

"""
    LTE_populations(atom::Atom,
                    temperature::Array{<:Unitful.Temperature, 3},
                    electron_density::Array{<:NumberDensity, 3})
Given the atom density, calculate the atom populations according to LTE.
"""
function LTE_populations2(atom::Atom,
                          temperature::Array{<:Unitful.Temperature, 3},
                          electron_density::Array{<:NumberDensity, 3})

    atom_density = atom.density
    χ = atom.χ
    g = atom.g
    U = atom.U
    n_levels = length(g)

    nz, nx, ny = size(temperature)
    nl = length(χ)
    populations = Array{Float64, 4}(undef, nz, nx, ny, nl)u"m^-3"

    z = Array{Float64, 4}(undef, nl, nz, nx, ny)

    for l=1:n_levels
        z[l,:,:,:] = g[l] * exp.(-χ[l]/k_B./temperature)
    end

    c = ( 2π*m_e*k_B/h^2 .* temperature ).^(1.5) .* 2.0 ./ electron_density * U[2] ./ U[1] .* exp.(-χ[end]/k_B./temperature)

    n3  = (c .* atom_density ./ (1.0 .+ c)) .|> u"m^-3"
    n2 = ((atom_density .- n3) ./ (z[1] ./ z[2] .+ 1.0) ) .|> u"m^-3"
    n1 = (atom_density .- n2 .- n3) .|> u"m^-3"

    populations[:,:,:,1] = n1
    populations[:,:,:,2] = n2
    populations[:,:,:,3] = n3

    @test all( Inf .> ustrip.(populations) .>= 0.0 )

    return populations
end

function LTE_populations(atom::Atom,
                         temperature::Array{<:Unitful.Temperature, 3},
                         electron_density::Array{<:NumberDensity, 3})
    χ = atom.χ
    g = atom.g
    atom_density = atom.density
    nz,nx,ny = size(atom_density)

    n_levels = length(χ)
    n_relative = ones(Float64, nz,nx,ny, n_levels)

    saha_const = (k_B / h) * (2π * m_e) / h
    saha_factor = 2 * ((saha_const * temperature).^(3/2) ./ electron_density) .|> u"m/m"

    for i=2:n_levels
        ΔE = χ[i] - χ[1]
        n_relative[:,:,:,i] = g[i] / g[1] * exp.(-ΔE ./ (k_B * temperature))
    end

    # Last level is ionised stage (H II)
    n_relative[:,:,:,n_levels] .*= saha_factor
    n_relative[:,:,:,1] = 1 ./ sum(n_relative, dims=4)[:,:,:,1]
    n_relative[:,:,:,2:end] .*= n_relative[:,:,:,1]

    return n_relative .* atom_density
end




"""
    write_to_file(populations::Array{<:NumberDensity,4},
                  iteration::Int64,
                  output_path::String)

Write the populations for a given iteration to the output file.
"""
function write_to_file(rates::TransitionRates,
                       iteration::Int64,
                       output_path::String)
    h5open(output_path, "r+") do file
        file["Rlu"][iteration+1,:,:,:,:] = ustrip.(rates.Rlu)
        file["Rul"][iteration+1,:,:,:,:] = ustrip.(rates.Rul)
        file["Clu"][iteration+1,:,:,:,:] = ustrip.(rates.Clu)
        file["Cul"][iteration+1,:,:,:,:] = ustrip.(rates.Cul)
    end
end
