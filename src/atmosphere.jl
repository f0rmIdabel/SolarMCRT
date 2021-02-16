include("io.jl")

struct Atmosphere
    z::Array{<:Unitful.Length, 1}                       # (nz + 1)
    x::Array{<:Unitful.Length, 1}                       # (nx + 1)
    y::Array{<:Unitful.Length, 1}                       # (ny + 1)
    velocity::Array{Array{<:Unitful.Velocity, 1}, 3}    # (nx, ny, nz)
    velocity_z::Array{<:Unitful.Velocity, 3}            # (nx, ny, nz)
    temperature::Array{<:Unitful.Temperature, 3}        # (nx, ny, nz)
    electron_density::Array{<:NumberDensity, 3}      # (nx, ny, nz)
    hydrogen_populations::Array{<:NumberDensity, 4}  # (nx, ny, nz, nl)
end

"""
Reads and slices atmosphere parameters accoring to inputs.
"""

function collect_atmosphere_data()

    # ===========================================================
    # READ ATMOSPHERE FILE
    # ===========================================================
    atmos = h5open(get_atmosphere_path(), "r")
    x = read(atmos, "x")u"m"
    y = read(atmos, "y")u"m"
    z = read(atmos, "z")u"m"

    # Change to single variable
    velocity_x = read(atmos, "velocity_x")u"m/s"
    velocity_y = read(atmos, "velocity_y")u"m/s"
    velocity_z = read(atmos, "velocity_z")u"m/s"

    temperature = read(atmos, "temperature")u"K"
    electron_density = read(atmos, "electron_density")u"m^-3"
    hydrogen_populations = read(atmos, "hydrogen_populations")u"m^-3"
    close(atmos)

    # original dimensions of data
    nz, nx, ny = size(temperature)

    # ===========================================================
    # CUT AND SLICE ATMOSPHERE BY INDEX
    # ===========================================================

    ze, xe, ye = get_stop()
    zs, xs, ys = get_start()
    dz, dx, dy = get_step()

    # Cut z-direction from below
    if ze != nothing && ze < nz
        nz = ze
        z = z[1:nz+1]
        velocity_x = velocity_x[1:nz,:,:]
        velocity_y = velocity_y[1:nz,:,:]
        velocity_z = velocity_z[1:nz,:,:]
        temperature = temperature[1:nz,:,:]
        electron_density = electron_density[1:nz,:,:]
        hydrogen_populations = hydrogen_populations[1:nz,:,:,:]
    end

    # Cut  z-direction from up top
    if zs > 1
        nz = zs
        z = z[nz:end]
        velocity_x = velocity_x[nz:end,:,:]
        velocity_y = velocity_y[nz:end,:,:]
        velocity_z = velocity_z[nz:end,:,:]
        temperature = temperature[nz:end,:,:]
        electron_density = electron_density[nz:end,:,:]
        hydrogen_populations = hydrogen_populations[nz:end,:,:,:]
    end

    # Cut x-direction from right
    if xe != nothing && xe < nx
        nx = xe
        x = x[1:nx+1]
        velocity_x = velocity_x[:,1:nx,:]
        velocity_y = velocity_y[:,1:nx,:]
        velocity_z = velocity_z[:,1:nx,:]
        temperature = temperature[:,1:nx,:]
        electron_density = electron_density[:,1:nx,:]
        hydrogen_populations = hydrogen_populations[:,1:nx,:,:]
    end

    # Cut x-direction from the left
    if xs > 1
        nx = xs
        x = x[nx:end]
        velocity_x = velocity_x[:,nx:end,:]
        velocity_y = velocity_y[:,nx:end,:]
        velocity_z = velocity_z[:,nx:end,:]
        temperature = temperature[:,nx:end,:]
        electron_density = electron_density[:,nx:end,:]
        hydrogen_populations = hydrogen_populations[:,nx:end,:,:]
    end

    # Cut y-direction from right
    if ye != nothing && ye < ny
        ny = ye
        y = y[1:ny+1]
        velocity_x = velocity_x[:,:,1:ny]
        velocity_y = velocity_y[:,:,1:ny]
        velocity_z = velocity_z[:,:,1:ny]
        temperature = temperature[:,:,1:ny]
        electron_density = electron_density[:,:,1:ny]
        hydrogen_populations = hydrogen_populations[:,:,1:ny,:]
    end

    # Cut y-direction from the left
    if ys > 1
        ny = ys
        y = y[ny:end]
        velocity_x = velocity_x[:,:,ny:end]
        velocity_y = velocity_y[:,:,ny:end]
        velocity_z = velocity_z[:,:,ny:end]
        temperature = temperature[:,:,ny:end]
        electron_density = electron_density[:,:,ny:end]
        hydrogen_populations = hydrogen_populations[:,:,ny:end,:]
    end

    # Only keep every dz-th box in z-direction
    if dz > 1
        z = z[1:dz:end]
        velocity_x = velocity_x[1:dz:end,:,:]
        velocity_y = velocity_y[1:dz:end,:,:]
        velocity_z = velocity_z[1:dz:end,:,:]
        temperature = temperature[1:dz:end,:,:]
        electron_density = electron_density[1:dz:end,:,:]
        hydrogen_populations = hydrogen_populations[1:dz:end,:,:,:]
    end

    # Only keep every dx-th box in x-direction
    if dx > 1
        x = x[1:dx:end]
        velocity_x = velocity_x[:,1:dx:end,:]
        velocity_y = velocity_y[:,1:dx:end,:]
        velocity_z = velocity_z[:,1:dx:end,:]
        temperature = temperature[:,1:dx:end,:]
        electron_density = electron_density[:,1:dx:end,:]
        hydrogen_populations = hydrogen_populations[:,1:dx:end,:,:]
    end

    # Only keep every dy-th box in y-direction
    if dy > 1
        y = y[1:dy:end]
        velocity_x = velocity_x[:,:,1:dy:end]
        velocity_y = velocity_y[:,:,1:dy:end]
        velocity_z = velocity_z[:,:,1:dy:end]
        temperature = temperature[:,:,1:dy:end]
        electron_density = electron_density[:,:,1:dy:end]
        hydrogen_populations = hydrogen_populations[:,:,1:dy:end,:]
    end

    # new size of data
    nz, nx, ny = size(temperature)

    # ===========================================================
    # ADD BOX END POINTS IN CASE LOST AFTER CUTTING
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
    # COLLECT VELOCITY IN ONE VARIABLE
    # ===========================================================

    velocity = Array{Array{<:Unitful.Velocity, 1}, 3}(undef,nz,nx,ny)

    for k=1:nz
        for i=1:nx
            for j=1:ny
                velocity[k,i,j] = [velocity_z[k,i,j], velocity_x[k,i,j], velocity_y[k,i,j]]
            end
        end
    end

    return z, x, y, velocity, velocity_z, temperature, electron_density, hydrogen_populations
end
