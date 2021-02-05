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



function check_converge(populations, new_populations, error, n, criterion = 1e-3)

    err = norm( abs.(populations .- new_population) ./populations)

    println(@sprintf("--Relative error = %.2e.", err))

    converged = false

    if err < criterion
        converged = true
    end

    return converged
end

"""
function get_revised_populations(atom, atomsphere, J)

    n3 = atmosphere.electron_density
    T = atmosphere.temperature

    n1 P12 + n1 P13 = n2 P21 + n3 P31
    n2 P21 + n2 P23 = n1 P12 + n3 P32
    n3 P31 + n3 P32 = n1 P13 + n2 P23


    Pij = Rij + Cij
    Pij = Aij + Bij*J + Cij


    # Rad rates
    Rij  = ∫4π/(hν) σij J dν
    Rji  = ∫4π/(hν) σij Gij (2hν^3/c^2 + J)dν


    # BB
    σij = h*νij/(4π)*Bij ϕνμ
    Gij = [ni/nj]LTE exp(-hν/kT)

    #BF
    σij = σic(ν)
    Gij = [ni/nc]LTE exp(-hν/kT)

end"""
