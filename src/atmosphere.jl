include("io.jl")

struct Atmosphere
    # Dimensions - position of box edges
    x::Array{<:Unitful.Length, 1} # Increasing
    y::Array{<:Unitful.Length, 1} # Increasing
    z::Array{<:Unitful.Length, 1} # Decreasing

    # Local box properties
    # (nx, ny, nz)
    velocity_x::Array{<:Unitful.Velocity, 3}
    velocity_y::Array{<:Unitful.Velocity, 3}
    velocity_z::Array{<:Unitful.Velocity, 3}
    temperature::Array{<:Unitful.Temperature, 3}

    # (λ, nx, ny, nz)
    χ::Array{PerLength, 3}
    ε::Array{Real, 3}

    # (λ, nx, ny)
    boundary::Array{Int64, 2}
end


"""
   function get_atmosphere_data(atmos_data,
                                rh_output)
Reads atmosphere parameters and reworks them to fit simulation.
"""
function collect_atmosphere_data(λ)

    # ===========================================================
    # READ INPUT FILE
    # ===========================================================
    atmosphere_path = get_atmosphere_path()
    τ_max = get_τ_max()

    # ===========================================================
    # READ ATMOSPHERE FILE
    # ===========================================================
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
    # CALCULATE ϵ and χ
    # ===========================================================
    ionised_hydrogen_density = hydrogen_populations[:,:,:,end]
    neutral_hydrogen_density = sum(hydrogen_populations, dims = 4) .- ionised_hydrogen_density

    # slice for λ, change later
    χ_a = χ_abs.(λ, temperature, electron_density, neutral_hydrogen_density, ionised_hydrogen_density)[:,:,:,1]
    χ_s = χ_scatt.(λ, electron_density, neutral_hydrogen_density)[:,:,:,1]
    χ = χ_a .+ χ_s
    ε = χ_a ./ χ

    # ===========================================================
    # RE-WORK DIMENSIONS TO FIT SIMULATION
    # ===========================================================

    # Make sure x and y are increasing and z decreasing
    if z[1] < z[end]
        z = reverse(z)
        velocity_x = velocity_x[end:-1:1,:,:]
        velocity_y = velocity_y[end:-1:1,:,:]
        velocity_z = velocity_z[end:-1:1,:,:]
        temperature = temperature[end:-1:1,:,:]
        χ = χ[end:-1:1,:,:]
        ε = ε[end:-1:1,:,:]
    end

    if x[1] > x[end]
        x = reverse(x)
        velocity_x = velocity_x[:,end:-1:1,:]
        velocity_y = velocity_y[:,end:-1:1,:]
        velocity_z = velocity_z[:,end:-1:1,:]
        temperature = temperature[:,end:-1:1,:]
        χ = χ[:,end:-1:1,:]
        ε = ε[:,end:-1:1,:]
    end

    if y[1] > y[end]
        y = reverse(y)
        velocity_x = velocity_x[:,:,end:-1:1]
        velocity_y = velocity_y[:,:,end:-1:1]
        velocity_z = velocity_z[:,:,end:-1:1]
        temperature = temperature[:,:,end:-1:1]
        χ = χ[:,:,end:-1:1]
        ε = ε[:,:,end:-1:1]
    end

    # Add endpoints for box calculations
    z = push!(z, 2*z[end] - z[end-1])
    x = push!(x, 2*x[end] - x[end-1])
    y = push!(y, 2*y[end] - y[end-1])

    # ===========================================================
    # CALCULATE OPTICAL DEPTH BOUNDARY AND CUT OFF DATA
    # ===========================================================
    boundary = optical_depth_boundary(χ, z, τ_max)

    nz = maximum(boundary)
    z = z[1:nz+1]
    velocity_x = velocity_x[1:nz,:,:]
    velocity_y = velocity_y[1:nz,:,:]
    velocity_z = velocity_z[1:nz,:,:]
    temperature = temperature[1:nz,:,:]
    χ = χ[1:nz,:,:]
    ε = ε[1:nz,:,:]

    return x, y, z, velocity_x, velocity_y, velocity_z,
           temperature, χ, ε, boundary
end

"""
From Tiago
"""
function χ_abs(λ::Unitful.Length,
               temperature::Unitful.Temperature,
               electron_density::NumberDensity,
               h_ground_density::NumberDensity,
               proton_density::NumberDensity)

    α = Transparency.hminus_ff_stilley(λ, temperature, h_ground_density, electron_density)
    α += Transparency.hminus_bf_geltman(λ, temperature, h_ground_density, electron_density)
    α += hydrogenic_ff(c_0 / λ, temperature, electron_density, proton_density, 1)
    #α += h2plus_ff(λ, temperature, h_ground_density, proton_density)
    #α += h2plus_bf(λ, temperature, h_ground_density, proton_density)
    return α
end

"""
From Tiago
"""
function χ_scatt(λ::Unitful.Length,
                 electron_density::NumberDensity,
                 h_ground_density::NumberDensity)

    α = thomson(electron_density)
    α += rayleigh_h(λ, h_ground_density)
    return α
end

function optical_depth(χ, z)
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
    function optical_depth_boundary(χ::Array{<:Unitful.Quantity{<:Real, Unitful.𝐋^(-1)}, 3},
                                    z::Array{<:Unitful.Length, 1},
                                    τ_max::Real)
Returns 2D array containing the k-indices where the optical depth reaches τ_max.
"""
function optical_depth_boundary(χ::Array{<:Unitful.Quantity{<:Real, Unitful.𝐋^(-1)}, 3},
                                z::Array{<:Unitful.Length, 1},
                                τ_max::Real)
    nz, nx, ny = size(χ)
    columns = nx*ny
    boundary = Array{Int64, 2}(undef, nx, ny)

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
