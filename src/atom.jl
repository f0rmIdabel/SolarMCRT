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

    λ::Array{Unitful.Length, 1}                    # (nλ)
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
function damping_constant(γ::Unitful.Frequency, ΔλD::Unitful.Length)
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

    return ustrip(populations)
end

function check_population_convergence(populations::Array{<:PerLength, 4}, new_populations::Array{<:PerLength, 4}, criterion = 1e-3)
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
    LTE_pops = LTE_populations(atom, temperature, electron_density)

    # ==================================================================
    # CALCULATE RADIATIVE RATES (nz, nx, ny)
    # ==================================================================
    σ12 = σij(1, 2, atom)
    σ13 = σic(1, atom)
    σ23 = σic(2, atom)

    G12 = Gij(1, 2, λ, temperature, LTE_pops)
    G13 = Gij(1, 3, λ, temperature, LTE_pops)
    G23 = Gij(2, 3, λ, temperature, LTE_pops)

    R12 = Rij(J, σ12, λ)
    R13 = Rij(J, σ13, λ)
    R23 = Rij(J, σ23, λ)
    R21 = Rji(J, σ12, G12, λ)
    R31 = Rji(J, σ13, G13, λ)
    R32 = Rji(J, σ23, G23, λ)

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

function Rij(Jν, σij, ν)
    # Rij  = ∫4π/(hν) σij J dν
    R = zeros(Float64, nz,nx,ny)

    for i=nν
        R +=  0.5(σij[i,:,:,:] / ν[i] .* Jν[i,:,:,:] * + σij[i+1,:,:,:] / ν[i+1] .* Jν[i+1,:,:,:]) * (ν[i+1] - ν[i])
    end

    R *= 4π/h

    return R
end

function Rji(Jν, σij, Gij, ν)
    # Rji  = ∫4π/(hν) σij Gij (2hν^3/c_0^2 + J)dν

    # m^2 sr
    # [1 / s / m^2 / sr]

    R = zeros(Float64, nz,nx,ny)

    for i=nν
        R +=  0.5(σij[i,:,:,:]   .* Gij[i,:,:,:]   / ν[i]   .* (2h*ν[i]^3/c_0^2   + Jν[i,:,:,:]) .+
                  σij[i+1,:,:,:] .* Gij[i+1,:,:,:] / ν[i+1] .* (2h*ν[i+1]^3/c_0^2 + Jν[i+1,:,:,:])) * (ν[i+1] - ν[i])
    end

    R *= 4π/h

    return R
end

function σij(i::Integer, j::Integer, atom::Atom)
    # σij = hc/(4π*λij)*Bij ϕλ

    damping_constant = atom.damping_constant
    ΔλD = atom.ΔλD
    λ0 = atom.line.λ0
    nλ_bb = atom.nλ_bb
    nλ_bf = atom.nλ_bf
    Bij = atom.line.Bij
    λ = atom.λ[2nλ_bf+1:end]

    σ_constant = h*c_0/(4π*λ0) .* Bij

    for l=1:nλ_bb
        damping = damping_constant .* λ[l]^2 |> u"m/m"
        v = (λ[l] - λ0) ./ ΔλD
        profile = voigt_profile.(damping, ustrip(v), ΔλD)
        σ[l,:,:,:] = σ_constant .* profile
    end

    return σ
end

function σic(i::Integer, atom::Atom)
    # σij = σic(ν)

    n_eff = sqrt(E_∞ / (atom.χu - atom.χl))
    charge = atom.Z
    ν = c_0/atom.λ[1+(i-1)*nλ_bf:i*nλ_bf]
    nλ_bf = atom.nλ_bf
    sigma_constant = 2.815e29 * charge^4/ i^5

    for l=1:nλ_bf
        gbf = gaunt_bf(ν[l], charge, n_eff)
        σ[l,:,:,:] = sigma_constant * gbf / ν[l]^3
    end

    return σ
end

function Gij(i::Integer, j::Integer, λ, temperature, LTE_populations)
    # Gij = [ni/nj]LTE exp(-hc/λkT)

    nλ = length(λ)
    n_ratio = LTE_populations[:,:,:,i] ./LTE_populations[:,:,:,j]
    C = h*c_0/k

    for l=1:nλ
        Gij[l,:,:,:] =  n_ratio .* exp.(-C/(λ[l]*temperature))
    end

    return Gij
end

"""
General for all ij
Excitation/Ionisation
De-excitation/Re-combination
"""
function Cij(n_i::Integer,
             n_j::Integer,
             electron_density::NumberDensity,
             temperature::Unitful.Temperature,
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


function n3(N, P12, P13, P21, P23, P31, P32)
    a = N .* P12 ./ (P21 .+ P23 .+ P12)
    b = (P32 .- P12) ./ (P21 .+ P23 .+ P12)

    c = N .* (P12 .+ P13) .- a .* (P21 .+ P12 .+ P13)
    d = b .* (P21 .+ P12 .+ P13) .+ P31 .+ P12 .+ P13

    return c ./ d
end

function n2(N, n3, P12, P21, P23, P32)
    return (N .* P12  + n3 .* (P32 .- P12) ) ./ (P21 .+ P23 .+ P12)
end

function n1(N, n3, n2)
    return N .- n3 .- n2
end

function write_to_file(λ::Array{<:Unitful.Length,1})
    h5open("../out/output.h5", "w") do file
        write(file, "wavelength", ustrip(λ))
    end
end
