include("io.jl")

struct Atmosphere
    z::Array{<:Unitful.Length, 1}                  # (nz + 1)
    x::Array{<:Unitful.Length, 1}                  # (nx + 1)
    y::Array{<:Unitful.Length, 1}                  # (ny + 1)
    velocity_z::Array{<:Unitful.Velocity, 3}       # (nx, ny, nz)
    velocity_x::Array{<:Unitful.Velocity, 3}       # (nx, ny, nz)
    velocity_y::Array{<:Unitful.Velocity, 3}       # (nx, ny, nz)
    temperature::Array{<:Unitful.Temperature, 3}   # (nx, ny, nz)
    χ::Array{PerLength, 3}                         # (nx, ny, nz)
    ε::Array{Real, 3}                              # (nx, ny, nz)
    boundary::Array{Int32, 2}                      # (nx, ny)
end

"""
Reads atmosphere parameters and reworks them to fit simulation.
"""
function collect_atmosphere_data(λ, line=false)

    # ===========================================================
    # READ ATMOSPHERE FILE
    # ===========================================================
    atmosphere_path = get_atmosphere_path()
    atmos = h5open(atmosphere_path, "r")

    x = read(atmos, "x")u"m"
    y = read(atmos, "y")u"m"
    z = read(atmos, "z")[:,1]u"m"  # slice due to snapshot

    velocity_x = read(atmos, "velocity_x")[:,:,:,1]u"m/s"
    velocity_y = read(atmos, "velocity_y")[:,:,:,1]u"m/s"
    velocity_z = read(atmos, "velocity_z")[:,:,:,1]u"m/s"

    temperature = read(atmos, "temperature")[:,:,:,1]u"K"
    electron_density = read(atmos, "electron_density")[:,:,:,1]u"m^-3"
    hydrogen_populations = read(atmos, "hydrogen_populations")[:,:,:,:,1]u"m^-3"
    close(atmos)

    # ===========================================================
    # CALCULATE χ AND ϵ
    # ===========================================================

    if line
        χ, ε = χ_ε_line()
    else
        χ, ε = χ_ε_continuum(λ, temperature, electron_density, hydrogen_populations)
    end

    """
    χ_c, ε_c = χ_ε_continuum(λ, temperature, electron_density, hydrogen_populations)

    if line
        χ_l, ε_l = χ_ε_line()

        χ = χ_l .+ χ_c
        ε = ε_l .* (χ_l ./ χ)  .+ ε_c .* (χ_c ./ χ)
    else
        χ = χ_c
        ε = ε_c
    end
    """


    # ===========================================================
    # FLIP AXES
    # ===========================================================

    # Make sure x and y are increasing and z decreasing
    if z[1] < z[end]
        z = reverse(z)
        velocity_z = velocity_z[end:-1:1,:,:]
        velocity_x = velocity_x[end:-1:1,:,:]
        velocity_y = velocity_y[end:-1:1,:,:]
        temperature = temperature[end:-1:1,:,:]
        χ = χ[end:-1:1,:,:]
        ε = ε[end:-1:1,:,:]
    end

    if x[1] > x[end]
        x = reverse(x)
        velocity_z = velocity_z[:,end:-1:1,:]
        velocity_x = velocity_x[:,end:-1:1,:]
        velocity_y = velocity_y[:,end:-1:1,:]
        temperature = temperature[:,end:-1:1,:]
        χ = χ[:,end:-1:1,:]
        ε = ε[:,end:-1:1,:]
    end


    if y[1] > y[end]
        y = reverse(y)
        velocity_z = velocity_z[:,:,end:-1:1]
        velocity_x = velocity_x[:,:,end:-1:1]
        velocity_y = velocity_y[:,:,end:-1:1]
        temperature = temperature[:,:,end:-1:1]
        χ = χ[:,:,end:-1:1]
        ε = ε[:,:,end:-1:1]
    end

    # ===========================================================
    # CALCULATE OPTICAL DEPTH BOUNDARY AND CUT OFF DATA
    # ===========================================================

    # dimensions of data
    nz, nx, ny = size(temperature)

    cut_off = get_cut_off()

    if cut_off == false
        boundary = Array{Int32,2}(undef,nx,ny)
        fill!(boundary, nz)
    else
        boundary = optical_depth_boundary(χ, z, cut_off)
        nz = maximum(boundary)
        z = z[1:nz+1]
        velocity_x = velocity_x[1:nz,:,:]
        velocity_y = velocity_y[1:nz,:,:]
        velocity_z = velocity_z[1:nz,:,:]
        temperature = temperature[1:nz,:,:]
        χ = χ[1:nz,:,:]
        ε = ε[1:nz,:,:]
    end

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
        χ = χ[1:nz,:,:]
        ε = ε[1:nz,:,:]
    end

    # Cut  z-direction from up top
    if zs > 1
        nz = zs
        z = z[nz:end]
        velocity_x = velocity_x[nz:end,:,:]
        velocity_y = velocity_y[nz:end,:,:]
        velocity_z = velocity_z[nz:end,:,:]
        temperature = temperature[nz:end,:,:]
        χ = χ[nz:end,:,:]
        ε = ε[nz:end,:,:]
        boundary = boundary .- nz
    end

    # Cut x-direction from right
    if xe != nothing && xe < nx
        nx = xe
        x = x[1:nx+1]
        velocity_x = velocity_x[:,1:nx,:]
        velocity_y = velocity_y[:,1:nx,:]
        velocity_z = velocity_z[:,1:nx,:]
        temperature = temperature[:,1:nx,:]
        χ = χ[:,1:nx,:]
        ε = ε[:,1:nx,:]
        boundary = boundary[1:nx,:]
    end

    # Cut x-direction from the left
    if xs > 1
        nx = xs
        x = x[nx:end]
        velocity_x = velocity_x[:,nx:end,:]
        velocity_y = velocity_y[:,nx:end,:]
        velocity_z = velocity_z[:,nx:end,:]
        temperature = temperature[:,nx:end,:]
        χ = χ[:,nx:end,:]
        ε = ε[:,nx:end,:]
        boundary = boundary[nx:end,:]
    end

    # Cut y-direction from right
    if ye != nothing && ye < ny
        ny = ye
        y = y[1:ny+1]
        velocity_x = velocity_x[:,:,1:ny]
        velocity_y = velocity_y[:,:,1:ny]
        velocity_z = velocity_z[:,:,1:ny]
        temperature = temperature[:,:,1:ny]
        χ = χ[:,:,1:ny]
        ε = ε[:,:,1:ny]
        boundary = boundary[:,1:ny]
    end

    # Cut y-direction from the left
    if ys > 1
        ny = ys
        y = y[ny:end]
        velocity_x = velocity_x[:,:,ny:end]
        velocity_y = velocity_y[:,:,ny:end]
        velocity_z = velocity_z[:,:,ny:end]
        temperature = temperature[:,:,ny:end]
        χ = χ[:,:,ny:end]
        ε = ε[:,:,ny:end]
        boundary[:,ny:end]
    end

    # Only keep every dz-th box in z-direction
    if dz > 1
        z = z[1:dz:end]
        velocity_x = velocity_x[1:dz:end,:,:]
        velocity_y = velocity_y[1:dz:end,:,:]
        velocity_z = velocity_z[1:dz:end,:,:]
        temperature = temperature[1:dz:end,:,:]
        χ = χ[1:dz:end,:,:]
        ε = ε[1:dz:end,:,:]
        boundary = boundary ./ dz  # fix this
    end

    # Only keep every dx-th box in x-direction
    if dx > 1
        x = x[1:dx:end]
        velocity_x = velocity_x[:,1:dx:end,:]
        velocity_y = velocity_y[:,1:dx:end,:]
        velocity_z = velocity_z[:,1:dx:end,:]
        temperature = temperature[:,1:dx:end,:]
        χ = χ[:,1:dx:end,:]
        ε = ε[:,1:dx:end,:]
        boundary = boundary[1:dx:end,:]
    end

    # Only keep every dy-th box in y-direction
    if dy > 1
        y = y[1:dy:end]
        velocity_x = velocity_x[:,:,1:dy:end]
        velocity_y = velocity_y[:,:,1:dy:end]
        velocity_z = velocity_z[:,:,1:dy:end]
        temperature = temperature[:,:,1:dy:end]
        χ = χ[:,:,1:dy:end]
        ε = ε[:,:,1:dy:end]
        boundary = boundary[:,1:dy:end]
    end

    nz, nx, ny = size(χ)

    # Add endpoints for box calculations
    if length(z) == nz
        z = push!(z, 2*z[end] - z[end-1])
    end
    if length(x) == nx
        x = push!(x, 2*x[end] - x[end-1])
    end
    if length(y) == ny
        y = push!(y, 2*y[end] - y[end-1])
    end

    return z, x, y, velocity_z, velocity_x, velocity_y,
           temperature, χ, ε, boundary
end


"""
To be replaced by Tiago's Transparency functions
"""
function χ_ε_line_rh()
    rh_ray = h5open("../../../../basement/MScProject/Atmospheres/output_ray.hdf5", "r")
    χ_l = read(rh_ray, "chi_line")[2,:,:,:]u"m^-1"
    χ_c = read(rh_ray, "chi_continuum")[2,:,:,:]u"m^-1"
    ε_c = read(rh_ray, "epsilon_continuum")[2,:,:,:]
    close(rh_ray)

    rh_aux = h5open("../../../../basement/MScProject/Atmospheres/output_aux.hdf5", "r")
    Rji = read(rh_aux, "atom_CA/Rji_line")[:,:,:,4]
    Cji = read(rh_aux, "atom_CA/Cji_line")[:,:,:,4]
    close(rh_aux)

    ε_l = Cji ./ (Rji .+ Cji)
    χ = χ_l .+ χ_c

    ε = ε_l .* (χ_l ./ χ)  .+ ε_c .* (χ_c ./ χ)

    return χ, ε
end


function χ_ε_line(temperature, z, electron_density, hydrogen_populations)

    wavelength =

    # Line data
    Hα = AtomicLine(97492.304u"cm^-1", 82259.158u"cm^-1", 109677.617u"cm^-1",
                    18, 8, 6.411e-01, 1.008 * m_u, 1)

    unsold_const = γ_unsold_const(Hα)

    # Compute continuum extinction
    # α_cont = αcont.(Hα.λ0, temp, electron_density, hpops[:, 1], hpops[:, end])
    # j_cont = α_cont .* (blackbody_λ.(Hα.λ0, temp) .|> u"kW / (m^2 * nm)")

    # Compute line extinction (van der Waals + natural broadening)
    γ = γ_unsold.(unsold_const, temperature, hydrogen_populations[:, 1]) .+ Hα.Aji
    ΔλD = doppler_width.(Hα.λ0, Hα.atom_weight, temperature)
    intensity = zeros(length(wavelength))u"kW / (m^2 * nm)"

    for (i, λ) in enumerate(wavelength)
        a = damping.(γ, λ, ΔλD)
        v = (λ - Hα.λ0) ./ ΔλD
        profile = voigt_profile.(a, v, ΔλD)
        χ = αline_λ.(Ref(Hα), profile, hydrogen_populations[:, 3], hydrogen_populations[:, 2])
        #ε = Hα.Aji .... lacking data

        #α_total = α_cont .+ αline_λ.(Ref(Hα), profile, hydrogen_populations[:, 3], hydrogen_populations[:, 2])
        #j_total = j_cont .+ jline_λ.(Ref(Hα), profile, hpops[:, 3])
        #source_function = j_total ./ α_total
        #intensity[i] = calc_intensity(-z, α_total, source_function)
    end

    return χ, ε
end


function χ_ε_continuum(λ, temperature, electron_density, hydrogen_populations)

    ionised_hydrogen_density = hydrogen_populations[:,:,:,end]
    neutral_hydrogen_density = sum(hydrogen_populations, dims = 4) .- ionised_hydrogen_density

    # slice for λ, change later
    χ_a = χ_abs.(λ, temperature, electron_density, neutral_hydrogen_density, ionised_hydrogen_density)[:,:,:,1]
    χ_s = χ_scatt.(λ, electron_density, neutral_hydrogen_density)[:,:,:,1]
    χ = χ_a .+ χ_s
    ε = χ_a ./ χ

    return χ, ε
end

"""
The extinction from continuum absorption processes for a given λ.
Includes H- ff, H- bf, H ff, H2+ ff and H2+ bf.
Credit: Tiago
"""
function χ_abs(λ::Unitful.Length,
               temperature::Unitful.Temperature,
               electron_density::NumberDensity,
               h_ground_density::NumberDensity,
               proton_density::NumberDensity)

    α = Transparency.hminus_ff_stilley(λ, temperature, h_ground_density, electron_density)
    α += Transparency.hminus_bf_geltman(λ, temperature, h_ground_density, electron_density)
    α += hydrogenic_ff(c_0 / λ, temperature, electron_density, proton_density, 1)
    α += h2plus_ff(λ, temperature, h_ground_density, proton_density)
    α += h2plus_bf(λ, temperature, h_ground_density, proton_density)
    return α
end

"""
The extincion from Thomson and Rayleigh scattering for a given λ.
Credit: Tiago
"""
function χ_scatt(λ::Unitful.Length,
                 electron_density::NumberDensity,
                 h_ground_density::NumberDensity)

    α = thomson(electron_density)
    α += rayleigh_h(λ, h_ground_density)
    return α
end

"""
Calculates the vertical optical depth of the atmosphere.
"""
function optical_depth(χ,
                       z::Array{<:Unitful.Length, 1})
    nz, nx, ny = size(χ)
    columns = nx*ny

    τ = Array{Float64,3}(undef, nz-1, nx, ny)

    # Calculate vertical optical depth for each column
    Threads.@threads for col=1:columns
        j = 1 + (col-1)÷nx
        i = col - (j-1)*nx
        τ[1,i,j] = 0.5(z[1] - z[2]) * (χ[1,i,j] + χ[2,i,j])

        for k=2:nz-1
            τ[k,i,j] =  τ[k-1,i,j] + 0.5(z[k] - z[k+1]) * (χ[k,i,j] + χ[k+1,i,j])
        end
    end

    return τ
end

"""
Returns 2D array containing the z-indices where the optical depth reaches τ_max.
"""
function optical_depth_boundary(χ::Array{<:Unitful.Quantity{<:Real, Unitful.𝐋^(-1)}, 3},
                                z::Array{<:Unitful.Length, 1},
                                τ_max::Real)
    nz, nx, ny = size(χ)
    columns = nx*ny
    boundary = Array{Int32, 2}(undef, nx, ny)

    # Calculate vertical optical depth for each column
    Threads.@threads for col=1:columns
        j = 1 + (col-1)÷nx
        i = col - (j-1)*nx

        τ = 0
        k = 0

        while τ < τ_max && k < nz
            k += 1
            # Trapezoidal rule
            τ += 0.5(z[k] - z[k+1]) * (χ[k,i,j] + χ[k+1,i,j])
        end
        boundary[i,j] = k
    end

    return boundary
end
