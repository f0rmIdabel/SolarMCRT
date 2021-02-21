include("io.jl")

struct Atom
    line::AtomicLine
    χl::Unitful.Energy
    χu::Unitful.Energy
    χ∞::Unitful.Energy
    gl::Int64
    gu::Int64
    g∞::Int64
    Z::Int64

    λ::Array{<:Unitful.Length, 1}                    # (nλ)
    nλ_bb::Int64
    nλ_bf::Int64

    doppler_width::Array{<:Unitful.Length, 3}      # (nz, nx, ny)
    damping_constant::Array{<:PerArea,3}             # (nz, nx, ny)
    populations::Array{<:NumberDensity, 4}           # (nz, nx, ny, nl=3)
end

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
Collection of 2 level atoms.
"""
function collect_atom_data(atmosphere::Atmosphere)

    atom = h5open(get_atom_path(), "r")
    χl = read(atom, "chi_l")u"cm^-1"
    χu = read(atom, "chi_u")u"cm^-1"
    χ∞ = read(atom, "chi_inf")u"cm^-1"
    gl = read(atom, "gl")
    gu = read(atom, "gu")
    g∞ = read(atom, "ginf")
    f_value = read(atom, "f_value")
    atom_weight = read(atom, "atom_weight")u"kg"
    Z = read(atom, "Z")
    close(atom)

    # Collect initial populations
    pop = h5open(get_initial_populations_path(), "r")
    initial_populations = read(pop, "initial_populations")u"m^-3"
    close(atom)

    nλ_bb, nλ_bf = get_nλ()

    line = AtomicLine(χu, χl, χ∞, gu, gl, f_value, atom_weight, Z)

    χl = line.χi
    χu = line.χj
    χ∞ = line.χ∞

    # Make sure odd number
    nλ_bb += 1-nλ_bb%2
    nλ = 2nλ_bf + nλ_bb

    # ==================================================================
    # SAMPLE WAVELENGTHS
    # ==================================================================
    λ = sample_λ(nλ_bb, nλ_bf, χl, χu, χ∞)

    # Calculate atmosphere quantities for line
    unsold_const = γ_unsold_const(line)
    γ = γ_unsold.(unsold_const, atmosphere.temperature, atmosphere.hydrogen_populations[:,:,:,1]) .+ line.Aji
    ΔλD = doppler_width.(line.λ0, atom_weight, atmosphere.temperature)
    damping_const = damping_constant.(γ, ΔλD)

    return line,
           χu, χl, χ∞,
           gu, gl, g∞,
           Z,
           λ, nλ_bb, nλ_bf,
           ΔλD, damping_const,
           initial_populations
end

function sample_λ(nλ_bb, nλ_bf, χl, χu, χ∞)

    atom_file = h5open(get_atom_path(), "r")
    λ_bf_l_min = read(atom_file, "bfl_min")u"nm"
    λ_bf_u_min = read(atom_file, "bfu_min")u"nm"
    qwing = read(atom_file, "qwing")
    qcore = read(atom_file, "qcore")
    close(atom_file)

    λ_bf_l_max  = transition_λ(χl, χ∞)
    λ_bf_u_max  = transition_λ(χu, χ∞)
    λ_bb_center = transition_λ(χl, χu)

    # Make sure odd # of bb wavelengths
    if nλ_bb > 0 && nλ_bb%2 == 0
        nλ_bb += 1
    end

    # Initialise wavelength array
    nλ = nλ_bf*2 + nλ_bb
    λ = Array{Float64,1}(undef, nλ)u"nm"

    # =================================================
    # Bound-free transitions
    # Linear spacing
    # =================================================
    if nλ_bf > 0
        Δλ_bf_l = (λ_bf_l_max - λ_bf_l_min)/nλ_bf
        Δλ_bf_u = (λ_bf_u_max - λ_bf_u_min)/nλ_bf

        λ[1] = λ_bf_l_min
        λ[nλ_bf+1] = λ_bf_u_min

        for l=2:nλ_bf
            λ[l] = λ[l-1] + Δλ_bf_l
            λ[l+nλ_bf] = λ[l+nλ_bf - 1] + Δλ_bf_u
        end
    end

    # =================================================
    # Bound-bound transition
    # Follows github.com/ITA-Solar/rh/blob/master/getlambda.c
    # =================================================
    if nλ_bb > 0
        vmicro_char = 2.5u"km/s"

        n = nλ_bb/2 # Questionable
        β = qwing/(2*qcore)
        y = β + sqrt(β*β + (β - 1.0)*n + 2.0 - 3.0*β)
        b = 2.0*log(y) / (n - 1)
        a = qwing / (n - 2.0 + y*y)

        center = nλ_bf*2 + (nλ_bb÷2) + 1
        λ[center] = λ_bb_center
        q_to_λ = λ[center] * vmicro_char / c_0

        for l=1:(nλ_bb÷2)
            Δλ = a*(l + (exp(b*l) - 1.0)) * q_to_λ
            λ[center-l] = λ[center] - Δλ
            λ[center+l] = λ[center] + Δλ
        end
    end

    return λ
end

"""
Compute damping parameter constant .
"""
function damping_constant(γ::Unitful.Frequency,
                          ΔλD::Unitful.Length)
    (γ / (4 * π * c_0 * ΔλD))
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

function calculate_transition_rates(atom::Atom,
                                    J::Array{<:UnitsIntensity_λ, 4},
                                    temperature::Array{<:Unitful.Temperature, 3},
                                    electron_density::Array{<:NumberDensity, 3})

    # ==================================================================
    # LOAD λ AND LTE POPULATIONS
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
    LTE_pops = LTE_populations(atom, temperature, electron_density)

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

    return R12, R13, R23, R21, R31, R32,
           C12, C13, C23, C21, C31, C32
end

function Rij(J::Array{<:UnitsIntensity_λ, 4},
             σij::Array{<:Unitful.Area, 4},
             λ::Array{<:Unitful.Length, 1})
    # Rij  = ∫4πλ/(hc) σij J dλ

    nz, nx, ny = size(J)
    nλ = size(λ)
    R = Array{Unitful.Frequency,3}(undef,nz,nx,ny)

    for l=1:(nλ-1)
        R +=  ( λ[l]   * σij[l,:,:,:]   .* J[l,:,:,:]   .+
                λ[l+1] * σij[l+1,:,:,:] .* J[l+1,:,:,:]  ) .* (λ[l+1] - λ[l])
    end

    R *= 2π/hc

    return R
end

function Rji(J::Array{<:UnitsIntensity_λ, 4},
             σij::Array{<:Unitful.Area, 4},
             Gij::Array{Float64, 4},
             λ::Array{<:Unitful.Length, 1})
    # Rji  = ∫4π/(hν) σij Gij (2hν^3/c_0^2 + J)dν
    # Rji = ∫4π/λ σij Gij (2c/λ^3 + Jλ^2/(hc)) dλ

    nz, nx, ny = size(J)
    R = Array{Unitful.Frequency,3}(undef,nz,nx,ny)

    # Trapezoid rule
    for l=1:(nλ-1)
        R += (σij[l,:,:,:]   .* Gij[l,:,:,:]   ./ λ[l]   .* (2c_0 ./ λ[l]^3   .+ J[l,:,:,:]   .* λ[l]^2   ./ hc ) .+
              σij[l+1,:,:,:] .* Gij[l+1,:,:,:] ./ λ[l+1] .* (2c_0 ./ λ[l+1]^3 .+ J[l+1,:,:,:] .* λ[l+1]^2 ./ hc )  ) .* (λ[l+1] - λ[l])
    end

    R *= 2π

    return R
end

function Rij(J::Array{<:UnitsIntensity_λ, 4},
             σij::Array{<:Unitful.Area, 1},
             λ::Array{<:Unitful.Length, 1})
    # Rij  = ∫4πλ/(hc) σij J dλ

    nλ, nz, nx, ny = size(J)
    R = Array{Unitful.Frequency,3}(undef,nz,nx,ny)

    for l=1:(nλ-1)
        R +=  ( λ[l]   * σij[l]   .* J[l,:,:,:]   .+
                λ[l+1] * σij[l+1] .* J[l+1,:,:,:]  ) .* (λ[l+1] - λ[l])
    end

    R *= 2π/hc

    return R
end

function Rji(J::Array{<:UnitsIntensity_λ, 4},
             σij::Array{<:Unitful.Area, 1},
             Gij::Array{Float64, 4},
             λ::Array{<:Unitful.Length, 1})
    # Rji  = ∫4π/(hν) σij Gij (2hν^3/c_0^2 + J)dν
    # Rji = ∫4π/λ σij Gij (2c/λ^3 + Jλ^2/(hc)) dλ

    nλ, nz, nx, ny = size(J)
    R = Array{Unitful.Frequency,3}(undef,nz,nx,ny)

    # Trapezoid rule
    for l=1:(nλ-1)
        R += (σij[l]   .* Gij[l,:,:,:]   ./ λ[l]   .* (2c_0 ./ λ[l]^3   .+ J[l,:,:,:]   .* λ[l]^2   / hc ) .+
              σij[l+1] .* Gij[l+1,:,:,:] ./ λ[l+1] .* (2c_0 ./ λ[l+1]^3 .+ J[l+1,:,:,:] .* λ[l+1]^2 ./ hc )  ) .* (λ[l+1] - λ[l])
    end

    R *= 2π

    return R
end

function σij(i::Integer,
             j::Integer,
             atom::Atom,
             λ::Array{<:Unitful.Length, 1})

    # σij = hc/(4π*λij)*Bij ϕλ
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

function σic(i::Integer,
             atom::Atom,
             λ::Array{<:Unitful.Length, 1})

    # σij = σic(ν)
    n_eff = sqrt(E_∞ / (atom.χl - atom.χu)) |>u"J/J" # should be χu - χl
    charge = atom.Z
    sigma_constant = 2.815e29 * charge^4/ i^5 * c_0^3
    nλ = length(λ)
    σ = Array{Unitful.Area, 1}(undef, nλ)

    gbf = gaunt_bf.(λ, charge, n_eff)
    σ = sigma_constant * gbf ./ λ.^3

    return σ
end

function Gij(i::Integer,
             j::Integer,
             λ::Array{<:Unitful.Length, 1},
             temperature::Array{<:Unitful.Temperature, 3},
             LTE_populations::Array{<:NumberDensity, 4})

    # Gij = [ni/nj]LTE exp(-hc/λkT)
    nλ = length(λ)
    nz, nx, ny = size(temperature)
    G = Array{Float64, 4}(undef, nλ, nz, nx, ny)

    n_ratio = LTE_populations[:,:,:,i] ./LTE_populations[:,:,:,j]

    @Threads.threads for l=1:nλ
        G[l,:,:,:] =  n_ratio .* exp.(- hc ./ (k_B*λ[l]*temperature))
    end

    return G
end

function Cij(n_i::Integer,
             n_j::Integer,
             electron_density::Array{<:NumberDensity,3},
             temperature::Array{<:Unitful.Temperature,3},
             LTE_populations::Array{<:NumberDensity,4})

    ionisation_level = size(LTE_populations)[end]

    # If UP
    if n_i < n_j
        if n_j < ionisation_level
            C = coll_exc_hydrogen_johnson(n_i, n_j, electron_density, temperature)
        elseif n_j == ionisation_level
            C = coll_ion_hydrogen_johnson(n_i, electron_density, temperature)
        end

    # If DOWN
    elseif n_i > n_j
        if n_i < ionisation_level
            C = coll_exc_hydrogen_johnson(n_j, n_i, electron_density, temperature)
        elseif n_i == ionisation_level
            C = coll_ion_hydrogen_johnson(n_j, electron_density, temperature)
        end

        C *= LTE_populations[:,:,:,n_i] ./ LTE_populations[:,:,:,n_j]
    end

    return C
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

"""
Stolen from Transparency repo
Recipes from Seaton
    gaunt_bf(charge::Int, n_eff::Number, λ::Unitful.Length)::Float64
Compute bound-free Gaunt factor for a given charge, effective principal
quantum number and wavelength λ. Taken from RH. Formula from
[Seaton (1960), Rep. Prog. Phys. 23, 313](https://ui.adsabs.harvard.edu/abs/1960RPPh...23..313S/abstract),
page 316.
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

function write_to_file(λ::Array{<:Unitful.Length,1})
    h5open("../out/output.h5", "w") do file
        write(file, "wavelength", ustrip(λ))
    end
end
