include("io.jl")

struct Atom
    line::AtomicLine
    χl::Unitful.Energy
    χu::Unitful.Energy
    χ∞::Unitful.Energy
    gl::Int64
    gu::Int64
    g∞::Int64

    nλ_bb::Int64
    nλ_bf::Int64

    ΔλD::Array{<:Unitful.Length, 3}                   # (nz, nx, ny)
    dc::Array{Float64,3}                              # (nz, nx, ny)
    αlc::Array{Float64,3}                             # (nz, nx, ny)
    populations::Array{<:Unitful.Length, 4}           # (nz, nx, ny, nl)
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
function collect_atom_data(atmosphere::Atmosphere, populations)

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
    close(atmos)

    nλ_bb, nλ_bf = get_nλ()

    line = AtomicLine(χu, χl, χ∞, gu, gl, f_value, atom_weight, Z)

    χl = line.χi
    χu = line.χu
    χ∞ = line.χ∞

    # Make sure odd number
    nλ_bb += 1-nλ_bb%2
    nλ = 2nλ_bf + nλ_bb

    # Calculate line-atmosphere quantities
    unsold_const = γ_unsold_const(line)
    γ = γ_unsold.(unsold_const, atmosphere.temperature, atmosphere.hydrogen_populations[:,:,:,1]) .+ line.Aji
    ΔλD = doppler_width.(line.λ0, atom_weight, atmosphere.temperature)
    dc = damping_const.(γ, ΔλD)
    αlc = α_line_const(line, initial_populations[:,:,:,1], initial_populations[:,:,:,2])

    return line,
           χu, χl, χ∞,
           gu, gl, g∞,
           nλ_bb, nλ_bf,
           ΔλD, dc, αlc,
           initial_populations
end

"""
Compute damping parameter constant .
"""
function damping_const(γ::Unitful.Frequency, ΔλD::Unitful.Length)
    (γ / (4 * π * c_0 * ΔλD)) |> u"m/m"
end

"""
Compute line extinction given an `AtomicLine` struct, `profile` defined per wavelength,
and upper and lower population densities `n_u` and `n_l`.
"""
function α_line_const(line::AtomicLine, n_u::NumberDensity, n_l::NumberDensity)
    (h * c_0 / (4 * π * line.λ0) * (n_l * line.Bij - n_u * line.Bji)) |> u"m/m"
end

function LTE_populations(atom::Atom, atom_density::NumberDensity, temperature::Unitful.Temperature, electron_density::NumberDensity)

    χl = atom.χi
    χu = atom.χj
    χ∞ = atom.χ∞

    gu = atom.gu
    gl = atom.gl
    g∞ = atom.g∞

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
                                 rates::TransitionRates,
                                 LTE_populations::Array{<:NumberDensity,4},
                                 λ::Array{<:Untiful.Length,1},
                                 temperature::Array{<:Unitful.Temperature,3},
                                 electron_density::Array{<:NumberDensity,3})

    # Transition probabilities
    P12 = rates.R12 .+ rates.C12
    P13 = rates.R13 .+ rates.C13
    P23 = rates.R23 .+ rates.C23
    P21 = rates.R21 .+ rates.C21
    P31 = rates.R31 .+ rates.C31
    P32 = rates.R32 .+ rates.C32

    # Revised populations
    atom_density = sum(LTE_populations, dims=4)[:,:,:,1]
    revised_populations[:,:,:,3] = n3(atom_density, P12, P13, P21, P23, P31, P32)
    revised_populations[:,:,:,2] = n2(atom_density, revised_populations[:,:,:,3], P12, P21, P23, P32)
    revised_populations[:,:,:,1] = n1(atom_density, revised_populations[:,:,:,3], revised_populations[:,:,:,2])

    return revised_populations
end

function calculate_transition_rates(atom, LTE_populations, temperature, electron_density)

    # Read output from simulation
    out = h5open("../out/output.h5", "r")
    J = read(out, "J")
    intensity_per_packet = read(out, "intensity_per_packet")u"kW / m^2 / sr / nm"
    close(out)

    # Convert to Frequency
    n_eff = sqrt(E∞ / (atom.χu - atom.χl))
    ν = c_0/λ
    nλ = length(intensity_per_packet)
    frequency_per_packet = ...
    for l=1:nλ
        J[l,:,:,:] *= frequency_per_packet[l]
    end

    # BB
    σ12 = σij(atom, B12, ν)
    G12 = Gij(1, 2, ν, temperature, LTE_populations)

    # BF
    σ13 = σic(1, ν, charge, n_eff)
    σ23 = σic(2, ν, charge, n_eff)
    G13 = Gij(1, 3, ν, temperature, LTE_populations)
    G23 = Gij(2, 3, ν, temperature, LTE_populations)

    # Radiative rates
    R12 = Rij(J, σ12, ν)
    R13 = Rij(J, σ13, ν)
    R23 = Rij(J, σ23, ν)
    R21 = Rji(J, σ12, G12, ν)
    R31 = Rji(J, σ13, G13, ν)
    R32 = Rji(J, σ23, G23, ν)

    # Collisional rates (nz, nx, ny)
    C12 = Cij(1, 2, electron_density, temperature, LTE_populations)
    C13 = Cij(1, 3, electron_density, temperature, LTE_populations)
    C23 = Cij(2, 3, electron_density, temperature, LTE_populations)
    C21 = Cij(2, 1, electron_density, temperature, LTE_populations)
    C31 = Cij(3, 1, electron_density, temperature, LTE_populations)
    C32 = Cij(3, 2, electron_density, temperature, LTE_populations)

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
    [1 / s / m^2 / sr]

    R = zeros(Float64, nz,nx,ny)

    for i=nν
        R +=  0.5(σij[i,:,:,:]   .* Gij[i,:,:,:]   / ν[i]   .* (2h*ν[i]^3/c_0^2   + Jν[i,:,:,:]) .+
                  σij[i+1,:,:,:] .* Gij[i+1,:,:,:] / ν[i+1] .* (2h*ν[i+1]^3/c_0^2 + Jν[i+1,:,:,:])) * (ν[i+1] - ν[i])
    end

    R *= 4π/h

    reurn R
end

function σij(atom, λ, Bij, ν)
    # σij = h*νij/(4π)*Bij ϕνμ
    dc = atom.dc
    ΔλD = atom.ΔλD
    λ0 = atom.line.λ0
    nλ_bb = atom.nλ_bb

    for l=1:nν_bb
        damp = dc*λ[l]^2 |> u"m/m"
        v = (λ[l] - λ0) ./ ΔλD
        ϕ = voigt_profile.(damp, ustrip(v), ΔλD)
        σ[l,:,:,:] = h/4π * ν[l]*Bij*ϕ
    end

    return σ
end

function σic(i, ν, charge, n_eff)
    # σij = σic(ν)
    c = 2.815*10e29 * charge^4/ i^5

    for l=1:nν_bf
        gbf = gaunt_bf(ν[l], charge, n_eff)
        σ[l,:,:,:] = c * gbf / ν[l]^3
    end

    return
end

function Gij(i, j, ν, temperature, LTE_populations)
    # Gij = [ni/nj]LTE exp(-hν/kT)

    for l=1:nν
        Gij[l,:,:,:] = LTE_populations[:,:,:,i] ./LTE_populations[:,:,:,j]  * exp.(-h*ν[l]/k/temperature)
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

    c = N *. (P12 .+ P13) - a .* (P21 .+ P12 .+ P13)
    d = b .* (P21 .+ P12 .+ P13) .+ P31 .+ P12 .+ P13

    return c ./ d
end

function n2(N, n3, P12, P21, P23, P32)
    return (N .* P12  + n3 .* (P32 .- P12) ) ./ (P21 .+ P23 .+ P12)
end

function n1(N, n3, n2)
    return N .- n3 .- n2
end
