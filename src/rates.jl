include("atmosphere.jl")
include("atom.jl")

struct TransitionRates
    R12::Array{<:Unitful.Frequency,3}
    R13::Array{<:Unitful.Frequency,3}
    R23::Array{<:Unitful.Frequency,3}
    R21::Array{<:Unitful.Frequency,3}
    R31::Array{<:Unitful.Frequency,3}
    R32::Array{<:Unitful.Frequency,3}
    C12::Array{<:Unitful.Frequency,3}
    C13::Array{<:Unitful.Frequency,3}
    C23::Array{<:Unitful.Frequency,3}
    C21::Array{<:Unitful.Frequency,3}
    C31::Array{<:Unitful.Frequency,3}
    C32::Array{<:Unitful.Frequency,3}
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
function calculate_transition_rates(atom::Atom,
                                    temperature::Array{<:Unitful.Temperature, 3},
                                    electron_density::Array{<:NumberDensity,3},
                                    J::Array{<:UnitsIntensity_λ, 4})

    # ==================================================================
    # LOAD ATMOSPHERE PARAMETERS
    # ==================================================================
    LTE_pops = LTE_populations(atom, temperature, electron_density)

    # ==================================================================
    # LOAD WAVELENGTHS AND SEPERATE BB AND BF
    # ==================================================================
    λ = atom.λ
    nλ_bb = atom.nλ_bb
    nλ_bf = atom.nλ_bf

    λ_bf_l = λ[1:nλ_bf]
    λ_bf_u = λ[nλ_bf+1:2nλ_bf]
    λ_bb = λ[2nλ_bf+1:end]

    J_bf_l = J[1:nλ_bf,:,:,:]
    J_bf_u = J[nλ_bf+1:2nλ_bf,:,:,:]
    J_bb = J[2nλ_bf+1:end,:,:,:]

    # ==================================================================
    # CALCULATE RADIATIVE RATES (nz, nx, ny)
    # ==================================================================
    σ12 = σij(1, 2, atom, λ_bb)
    σ13 = σic(1, atom, λ_bf_l)
    σ23 = σic(2, atom, λ_bf_u)

    G12 = Gij(1, 2, λ_bb, temperature, LTE_pops)
    G13 = Gij(1, 3, λ_bf_l, temperature, LTE_pops)
    G23 = Gij(2, 3, λ_bf_u, temperature, LTE_pops)

    R12 = Rij(J_bb, σ12, λ_bb)
    R13 = Rij(J_bf_l, σ13, λ_bf_l)
    R23 = Rij(J_bf_u, σ23, λ_bf_u)
    R21 = Rji(J_bb, σ12, G12, λ_bb)
    R31 = Rji(J_bf_l, σ13, G13, λ_bf_l)
    R32 = Rji(J_bf_u, σ23, G23, λ_bf_u)

    # ==================================================================
    # CALCULATE COLLISIONAL RATES (nz, nx, ny)
    # ==================================================================
    C12 = Cij(1, 2, electron_density, temperature, LTE_pops)
    C13 = Cij(1, 3, electron_density, temperature, LTE_pops)
    C23 = Cij(2, 3, electron_density, temperature, LTE_pops)
    C21 = Cij(2, 1, electron_density, temperature, LTE_pops)
    C31 = Cij(3, 1, electron_density, temperature, LTE_pops)
    C32 = Cij(3, 2, electron_density, temperature, LTE_pops)

    # ==================================================================
    # CHECK FOR UNVALID VALUES
    # ==================================================================
    @test all( Inf .> ustrip.(R12) .>= 0.0 )
    @test all( Inf .> ustrip.(R13) .>= 0.0 )
    @test all( Inf .> ustrip.(R23) .>= 0.0 )
    @test all( Inf .> ustrip.(R21) .>= 0.0 )
    @test all( Inf .> ustrip.(R31) .>= 0.0 )
    @test all( Inf .> ustrip.(R32) .>= 0.0 )
    @test all( Inf .> ustrip.(C12) .>= 0.0 )
    @test all( Inf .> ustrip.(C13) .>= 0.0 )
    @test all( Inf .> ustrip.(C23) .>= 0.0 )
    @test all( Inf .> ustrip.(C21) .>= 0.0 )
    @test all( Inf .> ustrip.(C31) .>= 0.0 )
    @test all( Inf .> ustrip.(C32) .>= 0.0 )

    return R12, R13, R23, R21, R31, R32,
           C12, C13, C23, C21, C31, C32
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
             atom::Atom,
             λ::Array{<:Unitful.Length, 1})

    damping_constant = atom.damping_constant
    ΔλD = atom.doppler_width
    λ0 = atom.line.λ0
    Bij = atom.line.Bij
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
    n_eff = sqrt(E_∞ / (atom.χu - atom.χl)) |>u"J/J" # should be χu - χl
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
        C = C .* ( LTE_populations[:,:,:,i] ./ LTE_populations[:,:,:,j] )
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
function LTE_populations(atom::Atom,
                         temperature::Array{<:Unitful.Temperature, 3},
                         electron_density::Array{<:NumberDensity, 3})

    χl = atom.χl
    χu = atom.χu
    χ∞ = atom.χ∞
    gl = atom.gl
    gu = atom.gu
    g∞ = atom.g∞
    U0 = atom.U0
    U1 = atom.U1

    atom_density = atom.density
    nz, nx, ny = size(temperature)
    populations = Array{Float64, 4}(undef, nz, nx, ny, 3)u"m^-3"

    z1 = gl * exp.( -χl/k_B./temperature)
    z2 = gu * exp.(- χu/k_B./temperature)

    c = ( 2π*m_e*k_B/h^2 .* temperature ).^(1.5) .* 2.0 ./ electron_density * U1 ./ U0 .* exp.(-χ∞/k_B./temperature)

    n3  = (c .* atom_density ./ (1.0 .+ c)) .|> u"m^-3"
    n2 = ((atom_density .- n3) ./ (z1 ./ z2 .+ 1.0) ) .|> u"m^-3"
    n1 = (atom_density .- n2 .- n3) .|> u"m^-3"

    populations[:,:,:,1] = n1
    populations[:,:,:,2] = n2
    populations[:,:,:,3] = n3

    @test all( Inf .> ustrip.(populations) .>= 0.0 )

    return populations
end
