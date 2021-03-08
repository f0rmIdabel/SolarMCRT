import Plots
using UnitfulRecipes
import Statistics

function plot_atmosphere(atmosphere::Atmosphere)
    # ===========================================================
    # LOAD DATA
    # ===========================================================
    z = atmosphere.z[1:end-1]
    T = atmosphere.temperature
    electron_density = atmosphere.electron_density
    hydrogen_populations = atmosphere.hydrogen_populations
    v = atmosphere.velocity
    speed = velocity_to_speed(v)

    # ===========================================================
    # GET AVERAGE COLUMN
    # ===========================================================

    mean_T = average_column(T)
    mean_electron_density = average_column(electron_density)
    mean_speed = average_column(speed)
    mean_h1 = average_column(hydrogen_populations[:,:,:,1])
    mean_h2 = average_column(hydrogen_populations[:,:,:,2])
    mean_h3 = average_column(hydrogen_populations[:,:,:,3])
    total =  mean_h1 .+ mean_h2 .+ mean_h3

    # ===========================================================
    # PLOT
    # ===========================================================

    ENV["GKSwstype"]="nul"
    p1 = Plots.plot(ustrip.(mean_T), ustrip.(z),
                    xlabel = "temperature [K]", ylabel = "z [m]",
                    xscale=:log10, legend = false)
    p2 = Plots.plot(ustrip.(mean_speed), ustrip.(z),
                    xlabel = "speed [m/s]", ylabel = "z [m]",
                    legend = false)
    p3 = Plots.plot(ustrip.(mean_electron_density), ustrip.(z),
                    xlabel = "electron density [m^-3]", ylabel = "z [m]",
                    xscale=:log10, legend = false)
    p4 = Plots.plot([mean_h1./total, mean_h2./total, mean_h3./total], ustrip.(z),
                     xlabel = "population density [m^-3]", ylabel = "z [m]",
                     xscale=:log10, label=permutedims(["ground","excited","ionised"]),
                     legendfontsize=6)
    Plots.plot(p1, p2, p3, p4, layout = (2, 2))
    Plots.png("plots/atmosphere")
end

function plot_populations(populations, z)
    # ===========================================================
    # GET AVERAGE COLUMN
    # ===========================================================
    mean_p1 = average_column(populations[:,:,:,1])
    mean_p2 = average_column(populations[:,:,:,2])
    mean_p3 = average_column(populations[:,:,:,3])
    total =  mean_p1 .+ mean_p2 .+ mean_p3

    # ===========================================================
    # PLOT
    # ===========================================================

    Plots.plot([mean_p1./total, mean_p2./total, mean_p3./total], ustrip.(z),
                xlabel = "population density [m^-3]", ylabel = "z [m]",
                xscale=:log10, label=permutedims(["ground","excited","ionised"]),
                legendfontsize=6)
    Plots.png("plots/initial_populations")
end

function plot_radiationBackground(radiationBackground, z)
    #λ = radiationBackground.λ
    α_continuum = radiationBackground.α_continuum[1,:,:,:]
    ε_continuum = radiationBackground.ε_continuum[1,:,:,:]
    boundary = radiationBackground.boundary[1,:,:]
    packets = radiationBackground.packets[1,:,:,:]
    #intensity_per_packet = radiationBackground.intensity_per_packet

    mean_α = average_column(ustrip.(α_continuum))u"m^-1"
    mean_ε = average_column(ε_continuum)
    mean_packets = average_column(packets)

    nx, ny = size(boundary)
    x = 1:nx
    y = 1:ny
    f(x,y) = ustrip(z[boundary[x,y]])

    ENV["GKSwstype"]="nul"
    p1 = Plots.plot(ustrip.(mean_α), ustrip.(z),
                    xlabel = "Extinction [m^-1]", ylabel = "z [m]",
                    xscale=:log10, legend = false)
    p2 = Plots.plot(mean_ε, ustrip.(z), xlabel = "Destruction probability", ylabel = "z [m]",
                    xscale=:log10, legend = false)
    p3 = Plots.plot(mean_packets, ustrip.(z), xlabel = "Packets", ylabel = "z [m]",
                    legend = false)
    p4 = Plots.surface(x, y, f,zlim = [ustrip(z[end]), ustrip(z[1])], legend = false, camera=(45,10))
    Plots.plot(p1, p2, p3, p4)
    Plots.png("plots/radiation_background")
end

function plot_rates(rates, z)

    z = ustrip.(z)
    R12 = average_column(ustrip.(rates.R12))
    R13 = average_column(ustrip.(rates.R13))
    R23 = average_column(ustrip.(rates.R23))
    R21 = average_column(ustrip.(rates.R21))
    R31 = average_column(ustrip.(rates.R31))
    R32 = average_column(ustrip.(rates.R32))
    C12 = average_column(ustrip.(rates.C12))
    C13 = average_column(ustrip.(rates.C13))
    C23 = average_column(ustrip.(rates.C23))
    C21 = average_column(ustrip.(rates.C21))
    C31 = average_column(ustrip.(rates.C31))
    C32 = average_column(ustrip.(rates.C32))

    p1 = Plots.plot([R12, C12], z,
                     xlabel = "rates [s^-1]", ylabel = "z [m]",
                     xscale=:log10, label=permutedims(["R12","C12"]),
                     legendfontsize=6)

    p2 = Plots.plot([R13, C13], z,
                     xlabel = "rates [s^-1]", ylabel = "z [m]",
                     xscale=:log10, label=permutedims(["R13","C13"]),
                     legendfontsize=6)

    p3 = Plots.plot([R23, C23], z,
                     xlabel = "rates [s^-1]", ylabel = "z [m]",
                     xscale=:log10, label=permutedims(["R23","C23"]),
                     legendfontsize=6)

    p4 = Plots.plot([R21, C21], z,
                      xlabel = "rates [s^-1]", ylabel = "z [m]",
                      xscale=:log10, label=permutedims(["R21","C21"]),
                      legendfontsize=6)

    p5 = Plots.plot([R31, C31], z,
                      xlabel = "rates [s^-1]", ylabel = "z [m]",
                      xscale=:log10, label=permutedims(["R31","C31"]),
                      legendfontsize=6)

    p6 = Plots.plot([R32, C32], z,
                      xlabel = "rates [s^-1]", ylabel = "z [m]",
                      xscale=:log10, label=permutedims(["R32","C32"]),
                      legendfontsize=6)

    Plots.plot(p1, p2, p3, p4, p5, p6, tickfontsize=6)
    Plots.png("plots/transition_rates")
end

function plot_radiation(radiation, atom, z)

    z = ustrip.(z[1:end-1])

    # ===========================================================
    # LOAD RADIATION DATA
    # ===========================================================
    α_continuum = ustrip.(radiation.α_continuum)
    ε_continuum = radiation.ε_continuum
    α_line_constant = radiation.α_line_constant
    ε_line = radiation.ε_line

    boundary = radiation.boundary
    packets = radiation.packets
    nλ, nz, nx, ny = size(packets)
    # ===========================================================
    # LOAD ATOM DATA AND GET LINE OPACITY/DESTRUCTION
    # ===========================================================
    λ = atom.λ
    nλ_bb = atom.nλ_bb
    nλ_bf = atom.nλ_bf
    α_line = Array{Float64, 4}(undef,nλ_bb,nz,nx,ny)

    for l=1:nλ_bb
        α_line[l,:,:,:] = ustrip.(line_extinction.(λ[2nλ_bf + l], atom.line.λ0, atom.doppler_width, atom.damping_constant, α_line_constant))
    end

    α_total = copy(α_continuum)
    ε_total = copy(ε_continuum)

    for l=1:nλ_bb
        α_total[(2nλ_bf + l),:,:,:] += α_line[l,:,:,:]
        ε_total[(2nλ_bf + l),:,:,:] = (ε_continuum[(2nλ_bf + l),:,:,:].*α_continuum[(2nλ_bf + l),:,:,:]  .+  ε_line.*α_line[l,:,:,:]) ./ α_total[(2nλ_bf + l),:,:,:]
    end

    mean_α = Array{Float64, 2}(undef,nλ, nz)
    #mean_α_line = Array{Float64, 2}(nλ_bb, nz)
    #mean_α_continuum = Array{Float64, 2}(nλ, nz)
    mean_ε = Array{Float64, 2}(undef,nλ, nz)

    mean_packets_bf_l = zeros(nz)
    mean_packets_bf_u = zeros(nz)
    mean_packets_bb = zeros(nz)

    mean_boundary = Array{Float64, 1}(undef,nλ)
    max_boundary = Array{Float64, 1}(undef,nλ)
    min_boundary = Array{Float64, 1}(undef,nλ)

    #mean_α[l,:] = average_column(α_total[l,:,:,:])
    #mean_ε[l,:] = average_column(ε_total[l,:,:,:])

    for l=1:nλ_bf
        mean_packets_bf_l += average_column(packets[l,:,:,:])
        mean_packets_bf_u += average_column(packets[nλ_bf+l,:,:,:])
    end

    for l=(2nλ_bf+1):nλ
        mean_packets_bb += average_column(packets[l,:,:,:])
    end

    for l=1:nλ
        mean_boundary[l] = mean(boundary[l,:,:])
        max_boundary[l] = maximum(boundary[l,:,:])
        min_boundary[l] = minimum(boundary[l,:,:])
    end

    λ = ustrip.(λ)
    #box = [70,]

    ENV["GKSwstype"]="nul"
    p1 = Plots.plot([ mean_α[1,:], mean_α[(nλ_bf+1),:], mean_α[(2nλ_bf+nλ_bb÷2+1),:] ], z,
                    xlabel = "Extinction [m^-1]", ylabel = "z [m]",
                    xscale=:log10,
                    label=permutedims(["bf lower", "bf upper", "bb center"]))
    p2 = Plots.plot([ mean_ε[1,:], mean_ε[(nλ_bf+1),:], mean_ε[(2nλ_bf+nλ_bb÷2+1),:] ], z,
                     xlabel = "Destruction probability", ylabel = "z [m]",
                     xscale=:log10,
                     label=permutedims(["bf lower", "bf upper", "bb center"]))
    p3 = Plots.plot([ mean_packets_bf_l, mean_packets_bf_u, mean_packets_bb], z,
                     xlabel = "Packets", ylabel = "z [m]",
                     label=permutedims(["bf lower", "bf upper", "bb center"]))

    p4 = Plots.plot(λ[1:nλ_bf], [mean_boundary[1:nλ_bf], max_boundary[1:nλ_bf], min_boundary[1:nλ_bf]],
                        ylabel = "Optical depth boundary", xlabel = "wavelength [nm]", yflip = true,
                        label=permutedims(["mean", "minimum", "maximum"]))
    p5 = Plots.plot(λ[(nλ_bf+1):2nλ_bf], [mean_boundary[(nλ_bf+1):2nλ_bf], max_boundary[(nλ_bf+1):2nλ_bf], min_boundary[(nλ_bf+1):2nλ_bf]],
                    #ylabel = "Optical depth boundary", xlabel = "wavelength [nm]", yflip = true,
                    yaxis = ("Lol", :flip),
                    label=permutedims(["mean", "minimum", "maximum"]))
    p6 = Plots.plot(λ[(2nλ_bf+1):nλ], [mean_boundary[(2nλ_bf+1):nλ], max_boundary[(2nλ_bf+1):nλ], min_boundary[(2nλ_bf+1):nλ]],
                    ylabel = "Optical depth boundary", xlabel = "wavelength [nm]", yflip = true,
                    label=permutedims(["mean", "minimum", "maximum"]))

    #p5 = Plots.plot([α_total[:, box...], α_continuum[:,box...]], λ, )
    Plots.plot(p1, p2, p3, tickfontsize=6, legendfontsize=6)
    Plots.png("plots/radiation_atmosphere")
    Plots.plot(p4, p5, p6, tickfontsize=6, legendfontsize=6)
    Plots.png("plots/radiation_boundary", layout=(3,1))
end

function average_column(array)
      Statistics.mean(array, dims=[2,3])[:,1,1]
end

function velocity_to_speed(velocity::Array{Array{<:Unitful.Velocity, 1}, 3})
    nz, nx, ny = size(velocity)
    speed = Array{Float64, 3}(undef, nz,nx,ny)

    for j=1:ny
        for i=1:nx
            for k=1:nz
                speed[k,i,j] = ustrip(sqrt(velocity[k,i,j][1]^2 + velocity[k,i,j][2]^2 + velocity[k,i,j][3]^2))
            end
        end
    end

    return speed*u"m/s"
end
