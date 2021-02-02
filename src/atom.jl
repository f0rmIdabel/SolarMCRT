include("io.jl")

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

    close(atom)
    return χu, χl, χ∞, gu, gl, f_value, atom_weight, Z
end

"""
Only valid for two-level
This needs to be generalised to read a file for other atoms
"""
function collect_initial_populations(hydrogen_populations)

    nz, nx, ny, nl = size(hydrogen_populations)
    populations = Array{Float64, 4}(undef, nz, nx, ny, 3)u"m^-3"
    populations[:,:,:,1:2] = hydrogen_populations[:,:,:,1:2]
    populations[:,:,:,3] = hydrogen_populations[:,:,:,end]

    return populations
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
