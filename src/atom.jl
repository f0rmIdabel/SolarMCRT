include("io.jl")

struct Atom
    line::AtomicLine
    χu::Unitful.Energy
    χl::Unitful.Energy
    χ∞::Unitful.Energy
    gu::Int64
    gl::Int64
    g∞::Int64

    f_value::Float64           # remove?
    atom_weight::Unitful.Mass  # remove?
    Z::Unitful.Int64           # remove?

    qwing::Float64
    qcore::Float64
    λ_bf_l_min::Unitful.Length
    λ_bf_l_max::Unitful.Length
    λ_bf_u_min::Unitful.Length
    λ_bf_u_max::Unitful.Length
    nλ_bb::Int64
    nλ_bf::Int64
end

"""
Collection of 2 level atoms.
"""
function collect_atom_data()

    atom = h5open(get_atom_path(), "r")
    χu = read(atom, "chi_u")u"cm^-1"
    χl = read(atom, "chi_l")u"cm^-1"
    χ∞ = read(atom, "chi_inf")u"cm^-1"
    gu = read(atom, "gu")
    gl = read(atom, "gl")
    f_value = read(atom, "f_value")
    atom_weight = read(atom, "atom_weight")u"kg"
    Z = read(atom, "Z")
    g∞ = read(atom, "ginf")
    λ_bf_l_min = read(atom, "bfl_min")u"nm"
    λ_bf_u_min = read(atom, "bfu_min")u"nm"
    qwing = read(atom, "qwing")
    qcore = read(atom, "qcore")
    close(atom)

    line = AtomicLine(χu, χl, χ∞, gu, gl, f_value, atom_weight, Z)

    #Converted to energy
    χl = line.χi
    χu = line.χj
    χ∞ = line.χ∞

    nλ_bb, nλ_bf = get_nλ()
    λ_bf_l_max = transition_λ(χl, χ∞)
    λ_bf_u_max = transition_λ(χu, χ∞)
    λ_bb_center = transition_λ(χl, χu)

    # Make sure odd number
    nλ_bf += 1-nλ_bf%2

    return line,
           χu, χl, χ∞,
           gu, gl, g∞,
           f_value, atom_weight, Z,
           qwing, qcore,
           λ_bf_l_min, λ_bf_l_max,
           λ_bf_u_min, λ_bf_u_max,
           nλ_bb, nλ_bf
end

function transition_λ(χ1::Unitful.Energy, χ2::Unitful.Energy)
    ((h * c_0) / (χ2-χ1)) |> u"nm"
end

function collect_initial_populations(atmosphere::Atmosphere, atom::Atom)

    pop = h5open(get_initial_populations_path(), "r")
    atom_populations = read(pop, "initial_populations")u"m^-3"
    close(atmos)

    nz, nx, ny, nl = size(atom_populations)

    # If less than 3 levels given, use Saha-Boltzmann
    if nl < 3
        temperature = atmosphere.temperature
        electron_density = atmosphere.electron_density
        atom_populations = sum(a, dims=4)[:,:,:,1]
        populations = LTE_populations(atom, temperature, atom_populations[:,:,:,1], electron_density)

    # If more than three levels given, pick out two first and last entry
    elseif nl > 3
        populations = Array{Float64, 4}(undef, nz, nx, ny, 3)u"m^-3"
        populations[:,:,:,1:2] = hydrogen_populations[:,:,:,1:2]
        populations[:,:,:,3] = hydrogen_populations[:,:,:,end]
    end

    return populations
end

function LTE_populations(atom::Atom, temperature::Unitful.Temperature, atom_density::NumberDensity, electron_density::NumberDensity)

    χl = atom.χl
    χu = atom.χu
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


function get_revised_populations(atomsphere::Atmosphere, atom::Atom, atom_density::Array{<:NumberDensity,3})

    out = h5open("../out/output.h5", "r")
    J = read(out, "J")
    intensity_per_packet = read(out, "intensity_per_packet")u"kW / m^2 / sr / nm"
    close(out)

    frequency_per_packet = ...

    nλ = length(intensity_per_packet)

    for l=1:nλ
        J[l,:,:,:] *= frequency_per_packet[l]
    end

    temperature = atmosphere.temperature
    electron_density = atmosphere.electron_density
    LTE_populations =  LTE_populations(atom, temperature, atom_density, electron_density)

    σ12 =


    P12 = Rij(J, )

    # BB
    Pij = Rij + Cij
    Pij = Aij + Bij*J + Cij

    # Rad rates

    # BB
    σij = h*νij/(4π)*Bij ϕνμ
    Gij = [ni/nj]LTE exp(-hν/kT)

    #BF
    σij = σic(ν)
    Gij = [ni/nc]LTE exp(-hν/kT)

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

function σij(B12, ν12, profile)
    # σij = h*νij/(4π)*Bij ϕνμ


    σ *= h/4π

end

function Gij(i,j)

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
