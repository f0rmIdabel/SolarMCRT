include("../../src/atom.jl")

"""
Reworks a typical atmosphere
Should be made more flexible
"""
function rework_atmosphere_data()

    # ===========================================================
    # READ ATMOSPHERE FILE AND SCLICE SNAPSHOT
    # ===========================================================
    atmos = h5open(get_atmosphere_path(), "r+")

    x = read(atmos, "x") #u"m"
    y = read(atmos, "y") #u"m"
    z = read(atmos, "z")[:,1] #u"m"  # slice due to snapshot

    # Change to single variable
    velocity_x = read(atmos, "velocity_x")[:,:,:,1] #u"m/s"
    velocity_y = read(atmos, "velocity_y")[:,:,:,1] #u"m/s"
    velocity_z = read(atmos, "velocity_z")[:,:,:,1] #u"m/s"

    temperature = read(atmos, "temperature")[:,:,:,1] #u"K"
    electron_density = read(atmos, "electron_density")[:,:,:,1] #u"m^-3"
    hydrogen_populations = read(atmos, "hydrogen_populations")[:,:,:,:,1] #u"m^-3"
    close(atmos)

    # original dimensions of data
    nz, nx, ny = size(temperature)

    # ===========================================================
    # FLIP AXES, MOVE TO TEST
    # ===========================================================

    # Make sure x and y are increasing and z decreasing
    if z[1] < z[end]
        z = reverse(z)
        velocity_z = velocity_z[end:-1:1,:,:]
        velocity_x = velocity_x[end:-1:1,:,:]
        velocity_y = velocity_y[end:-1:1,:,:]
        temperature = temperature[end:-1:1,:,:]
        electron_density = electron_density[end:-1:1,:,:]
        hydrogen_populations = hydrogen_populations[end:-1:1,:,:,:]
    end

    if x[1] > x[end]
        x = reverse(x)
        velocity_z = velocity_z[:,end:-1:1,:]
        velocity_x = velocity_x[:,end:-1:1,:]
        velocity_y = velocity_y[:,end:-1:1,:]
        temperature = temperature[:,end:-1:1,:]
        electron_density = electron_density[:,end:-1:1,:]
        hydrogen_populations = hydrogen_populations[:,end:-1:1,:,:]
    end

    if y[1] > y[end]
        y = reverse(y)
        velocity_z = velocity_z[:,:,end:-1:1]
        velocity_x = velocity_x[:,:,end:-1:1]
        velocity_y = velocity_y[:,:,end:-1:1]
        temperature = temperature[:,:,end:-1:1]
        electron_density = electron_density[:,:,end:-1:1]
        hydrogen_populations = hydrogen_populations[:,:,end:-1:1,:]
    end

    # ===========================================================
    # ADD BOX END POINTS
    # ===========================================================

    if length(z) == nz
        z = push!(z, 2*z[end] - z[end-1])
    end
    if length(x) == nx
        x = push!(x, 2*x[end] - x[end-1])
    end
    if length(y) == ny
        y = push!(y, 2*y[end] - y[end-1])
    end

    # ===========================================================
    # LTE POPULATIONS
    # ===========================================================
    n_levels = size(hydrogen_populations)[4]

    if n_levels == 1
        atom = Atom(collect_atom_data()...)
        hydrogen_populations = LTE_populations(atom, temperature*u"K", hydrogen_populations*u"m^-3", electron_density*u"m^-3")
    elseif n_levels > 3
        populations = Array{Float64, 4}(undef, nz, nx, ny, 3)
        populations[:,:,:,1:2] = hydrogen_populations[:,:,:,1:2]
        populations[:,:,:,3] = hydrogen_populations[:,:,:,end]
        hydrogen_populations = populations
    end


    # ===========================================================
    # WRITE TO FILE
    # ===========================================================
    new_atmos = h5open("/mn/stornext/u3/idarhan/basement/MScProject/Atmospheres/bifrost_qs006023_s525_quarter2.hdf5", "w")
    write(new_atmos, "x", x)
    write(new_atmos, "y", y)
    write(new_atmos, "z", z)
    write(new_atmos, "velocity_x", velocity_x)
    write(new_atmos, "velocity_y", velocity_y)
    write(new_atmos, "velocity_z", velocity_z)
    write(new_atmos, "temperature", temperature)
    write(new_atmos, "electron_density", electron_density)
    write(new_atmos, "hydrogen_populations", hydrogen_populations)
    close(new_atmos)

    # ===========================================================
    # LTE POPULATIONS
    # ===========================================================
    pop = h5open(get_initial_populations_path(), "w")
    write(pop, "initial_populations", hydrogen_populations)
    close(pop)
end

function LTE_populations(atom, temperature, hydrogen_density, electron_density)

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

    populations[:,:,:,3] = K .* hydrogen_density ./ (1.0 .+ K)
    populations[:,:,:,1] = (hydrogen_density .- populations[:,:,:,3]) ./ (1.0 .+ U2 ./ U1)
    populations[:,:,:,2] = hydrogen_density .- populations[:,:,:,1] .- populations[:,:,:,3]

    return ustrip(populations)
end

"""
REMOVE
Only valid for two-level
This needs to be generalised to read a file for other atoms
"""
function collect_initial_populations(atmosphere, atom)
    hydrogen_populations = atmosphere.hydrogen_populations

    nz, nx, ny, nl = size(hydrogen_populations)

    # If only number of H-atoms given, use Saha-Boltzmann
    if nl == 1
        temperature = atmosphere.temperature
        electron_density = atmosphere.electron_density
        populations = LTE_populations(atom, temperature, hydrogen_populations[:,:,:,1], electron_density)

    # If several levels and ionisation stage given, pick out two first and last entry
    elseif nl > 1
        populations = Array{Float64, 4}(undef, nz, nx, ny, 3)u"m^-3"
        populations[:,:,:,1:2] = hydrogen_populations[:,:,:,1:2]
        populations[:,:,:,3] = hydrogen_populations[:,:,:,end]
    end

    return populations
end

#rework_atmosphere_data()
