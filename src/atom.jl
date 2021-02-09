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

function check_converge(populations::Array{<:PerLength, 4}, new_populations::Array{<:PerLength, 4}, criterion = 1e-3)
    N = length(populations)
    error = norm( abs.(populations .- new_population) ./populations ./N)

    println(@sprintf("--Relative error = %.2e.", err))

    converged = false

    if error < criterion
        converged = true
    end

    return converged
end

function get_revised_populations(atom::Atom, LTE_populations::Array{<:NumberDensity,4}, λ::Array{<:Untiful.Length,1}, temperature::Array{<:Unitful.Temperature, 3})

    # Read output from simulation
    out = h5open("../out/output.h5", "r")
    J = read(out, "J")
    intensity_per_packet = read(out, "intensity_per_packet")u"kW / m^2 / sr / nm"
    close(out)

    # Convert to Frequency
    ν = c_0/λ
    frequency_per_packet = ...
    nλ = length(intensity_per_packet)
    for l=1:nλ
        J[l,:,:,:] *= frequency_per_packet[l]
    end

    σ12 = ...
    σ13 =

    G12 = ...


    P12 = Rij(J, σ12, ν) + C12
    P21 = Rji(J, σ12, G12, ν) + C21
    ...

    # BB
    # Pij = Rij + Cij
    # Pij = Aij + Bij*J + Cij

    # Rad rates

    # BB
    # σij = h*νij/(4π)*Bij ϕνμ
    # Gij = [ni/nj]LTE exp(-hν/kT)

    #BF
    # σij = σic(ν)
    # Gij = [ni/nc]LTE exp(-hν/kT)

    revised_populations[:,:,:,3] = n3(atom_density, n3, P12, P13, P21, P23, P31, P32)
    revised_populations[:,:,:,2] = n2(atom_density, n3, P12, P21, P23, P32)
    revised_populations[:,:,:,1] = n1(atom_density, n3, n2)
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

    R = zeros(Float64, nz,nx,ny)

    for i=nν
        R +=  0.5(σij[i,:,:,:]   .* Gij[i,:,:,:]   / ν[i]   .* (2h*ν[i]^3/c_0^2   + Jν[i,:,:,:]) .+
                  σij[i+1,:,:,:] .* Gij[i+1,:,:,:] / ν[i+1] .* (2h*ν[i+1]^3/c_0^2 + Jν[i+1,:,:,:])) * (ν[i+1] - ν[i])
    end

    R *= 4π/h

    reurn R
end

function σij(Bij, νij, ϕ)
    # σij = h*νij/(4π)*Bij ϕνμ

    for l=1:nν_bb
        σ[l,:,:,:] = ν[l]*Bij*ϕ[l]
    end
    σ *= h/4π
    return σ
end

function Gij(ni_LTE, nj_LTE, temperature, ν)
    # Gij = [ni/nj]LTE exp(-hν/kT)

    for n=1:nν
        Gij[l,:,:,:] = ni_LTE ./nj_LTE * exp.(-h*ν[l]/k/temperature)
    end
    return Gij
end

function n3(N, n3, P12, P13, P21, P23, P31, P32)
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
