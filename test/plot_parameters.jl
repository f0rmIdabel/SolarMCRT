import Plots
using UnitfulRecipes
import Statistics

function plot_atmosphere(atmosphere::Atmosphere)
    z = atmosphere.z[1:end-1]
    T = atmosphere.temperature
    electron_density = atmosphere.electron_density
    hydrogen_populations = atmosphere.hydrogen_populations
    v = atmosphere.velocity
    speed = velocity_to_speed(v)

    mean_T = average_column(T)
    mean_electron_density = average_column(electron_density)
    mean_speed = average_column(speed)

    mean_h1 = average_column(hydrogen_populations[:,:,:,1])
    mean_h2 = average_column(hydrogen_populations[:,:,:,2])
    mean_h3 = average_column(hydrogen_populations[:,:,:,3])
    total =  mean_h1 .+ mean_h2 .+ mean_h3

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
                     legendfontsize=5)
    Plots.plot(p1, p2, p3, p4, layout = (2, 2))
    Plots.png("plots/atmosphere")
end

function plot_radiationBackground(radiationBackground, z)
    #λ = radiationBackground.λ
    z = z[1:end-1]
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
                    xlabel = "α [m^-1]", ylabel = "z [m]",
                    xscale=:log10, legend = false)
    p2 = Plots.plot(mean_ε, ustrip.(z), xlabel = "ε", ylabel = "z [m]",
                    xscale=:log10, legend = false)
    p3 = Plots.plot(mean_packets, ustrip.(z), xlabel = "Packets", ylabel = "z [m]",
                    xscale=:log10, legend = false)
    p4 = Plots.surface(x, y, f,
                      zlim = [ustrip(z[end]), ustrip(z[1])],
                      camera=(-45,camera_tilt), legend=false)

    Plots.plot(p1, p2, p3, p4)
    Plots.png("plots/radiationBackground")
end

function plot_populations(populations)
    mean_p1 = average_column(populations[:,:,:,1])
    mean_p2 = average_column(populations[:,:,:,2])
    mean_p3 = average_column(populations[:,:,:,3])
    total =  mean_p1 .+ mean_p2 .+ mean_p3

    Plots.plot([mean_p1./total, mean_p2./total, mean_p3./total],
                xlabel = "population density [m^-3]", ylabel = "z [m]",
                yscale=:log10)
    Plots.png("plots/initial_populations")
end


function plot_rates(rates, z)

    R12 = average_column(rates.R12)
    R13 = average_column(rates.R13)
    R23 = average_column(rates.R23)
    R21 = average_column(rates.R21)
    R31 = average_column(rates.R31)
    R32 = average_column(rates.R32)
    C12 = average_column(rates.C12)
    C13 = average_column(rates.C13)
    C23 = average_column(rates.C23)
    C21 = average_column(rates.C21)
    C31 = average_column(rates.C31)
    C32 = average_column(rates.C32)

    p1 = Plot


end


function plot_radiation(radiation, z, λ)
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
