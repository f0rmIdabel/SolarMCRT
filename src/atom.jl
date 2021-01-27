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


function update_populations(λ, J)
end
