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

# Should be generalised
function get_initial_populations()
end


function update_populations(λ, J)
end


"""
χl = 82259.158u"cm^-1"
χ∞ = 109677.617u"cm^-1"
gu = 18
gl = 8
f_value = 6.411e-01
atom_weight = 1.008*m_u
Z = 1
"""
