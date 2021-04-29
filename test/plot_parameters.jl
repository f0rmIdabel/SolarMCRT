using Plots
using Unitful
using UnitfulRecipes
import Statistics


"""
    plot_atmosphere(atmosphere::Atmosphere)

Plot column averaged temperature, electron density,
hydrogen populations and speed for an atmosphere.
"""
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

    mean_T = average_column(ustrip.(T))
    mean_electron_density = average_column(ustrip.(electron_density))
    mean_speed = average_column(ustrip.(speed))
    hydrogen_populations = ustrip.(hydrogen_populations)
    mean_h1 = average_column(hydrogen_populations[:,:,:,1])
    mean_h2 = average_column(hydrogen_populations[:,:,:,2])
    mean_h3 = average_column(hydrogen_populations[:,:,:,3])
    total =  mean_h1 .+ mean_h2 .+ mean_h3

    # ===========================================================
    # PLOT
    # ===========================================================

    ENV["GKSwstype"]="nul"
    z = ustrip.(z .|>u"Mm")

    p1 = Plots.plot(z, ustrip.(mean_T),
                    xlabel = "z (Mm)", ylabel = "temperature (K)",
                    yscale=:log10, legend = false)
    p2 = Plots.plot(z, ustrip.(mean_speed),
                    xlabel = "z (Mm)", ylabel = "speed (m/s)",
                    legend = false)
    p3 = Plots.plot(z, ustrip.(mean_electron_density),
                    xlabel = "z (Mm)", ylabel = "electron density (m^-3)",
                    yscale=:log10, legend = false)
    p4 = Plots.plot(z, [mean_h1./total, mean_h2./total, mean_h3./total],
                    xlabel = "z (Mm)", ylabel = "relativ HI populations",
                    yscale=:log10,
                    label=permutedims(["ground","excited","ionised"]),
                    legendfontsize=6)

    Plots.plot(p1, p2, p3, p4, layout = (2, 2))
    Plots.png("plots/atmosphere")
end

"""
    plot_populations(populations::Array{<:NumberDensity, 3},
                     z::Array{<:Unitful.Length,1})

Plot the population distribution for a 2-level atom with continuum.
"""
function plot_populations(populations_LTE::Array{<:NumberDensity, 4},
                          populations_ZR::Array{<:NumberDensity, 4},
                          z::Array{<:Unitful.Length,1})
    # ===========================================================
    # GET AVERAGE COLUMN
    # ===========================================================
    np = size(populations_LTE)[end]
    mean_LTE = []
    mean_ZR = []
    lbl = []

    Np = average_column(sum(populations_LTE, dims=4)[:,:,:,1])

    for i=1:np
        append!(mean_LTE, [average_column(populations_LTE[:,:,:,i]) ./ Np])
        append!(mean_ZR, [average_column(populations_ZR[:,:,:,i]) ./ Np])
        append!(lbl, string(i))
    end

    # ===========================================================
    # PLOT
    # ===========================================================
    ENV["GKSwstype"]="nul"
    z = ustrip.(z .|>u"Mm")

    p1 = Plots.plot(z, mean_LTE,
                    xlabel = "z (Mm)", ylabel = "Relative populations",
                    yscale=:log10, label=permutedims(lbl))
    p2 = Plots.plot(z, mean_ZR,
                    xlabel = "z (Mm)", ylabel = "Relative populations",
                    yscale=:log10, label=permutedims(lbl))
    Plots.plot(p1, p2, title=permutedims(["LTE", "Zero-radiation"]), layout=(2,1))
    Plots.png("plots/initial_populations")
end

"""
    plot_radiationBackground(radiationBackground::RadiationBackground,
                             z::Array{<:Unitful.Length, 1})

Plots column averaged extincion, destruction probability,
package creation and depth boundary for a background wavelength.
"""
function plot_radiationBackground(radiation::Radiation,
                                  z::Array{<:Unitful.Length, 1})

    α_continuum = radiation.α_continuum[1,:,:,:]
    ε_continuum = radiation.ε_continuum[1,:,:,:]
    boundary = radiation.boundary[1,:,:]
    packets = radiation.packets[1,:,:,:]
    #intensity_per_packet = radiationBackground.intensity_per_packet

    nz, nx, ny = size(packets)

    for i=1:nx
        for j=1:ny
            for k=(boundary[i,j]+1):nz
                packets[k,i,j] = 0
            end
        end
    end

    mean_α = average_column(ustrip.(α_continuum))u"m^-1"
    mean_ε = average_column(ε_continuum)
    mean_packets = average_column(packets)

    x = 1:nx
    y = 1:ny
    f(x,y) = ustrip(z[boundary[x,y]])

    ENV["GKSwstype"]="nul"
    z = ustrip.(z .|>u"Mm")
    mean_α = ustrip.(mean_α)

    p1 = Plots.plot(z, mean_α,
                    xlabel = "z (Mm)", ylabel = "Extinction (m^-1)",
                    yscale=:log10, legend = false)
    p2 = Plots.plot(z, mean_ε,
                    xlabel = "z (Mm)", ylabel = "Destruction",
                    yscale=:log10)
    p3 = Plots.plot(z, mean_packets,
                    xlabel = "z (Mm)", ylabel = "Packets")
    p4 = Plots.surface(x, y, f, zlim = [z[end], z[1]], camera=(45,10))
    Plots.plot(p1, p2, p3, p4, legend=false)
    Plots.png("plots/radiation_background")
end


"""
    plot_rates(rates::TransitionRates,
               z::Array{<:Unitful.Length, 1})

Plot the column averaged transition rates.
"""
function plot_rates(rates::TransitionRates,
                    z::Array{<:Unitful.Length, 1})

    # ===========================================================
    # LOAD DATA
    # ===========================================================
    R = rates.R
    C = rates.C
    n_levels = size(R)[1] - 1

    # ===========================================================
    # CHECK DIMENSIONS
    # ===========================================================

    ENV["GKSwstype"]="nul"
    z = ustrip.(z .|>u"Mm")
    plots = []

    for l=1:n_levels
        for u=l+1:n_levels+1
            Rlu = average_column(ustrip.(R[l,u,:,:,:]))
            Rul = average_column(ustrip.(R[u,l,:,:,:]))
            Clu = average_column(ustrip.(C[l,u,:,:,:]))
            Cul = average_column(ustrip.(C[u,l,:,:,:]))

            up = Plots.plot(z, [Rlu, Clu],
                            ylabel = "rates (s^-1)", xlabel = "z (Mm)", yscale=:log10,
                            label=permutedims(["R"*string(l)*string(u), "C"*string(l)*string(u)]),
                            legendfontsize=6)

            down = Plots.plot(z, [Rul, Cul],
                            ylabel = "rates (s^-1)", xlabel = "z (Mm)", yscale=:log10,
                            label=permutedims(["R"*string(u)*string(l), "C"*string(u)*string(l)]),
                            legendfontsize=6)

            append!(plots, [up])
            append!(plots, [down])
        end
    end

    Plots.plot(plots..., tickfontsize=6)
    Plots.png("transition_rates")

end


"""
    plot_radiation(radiation::Radiation,
                   atom::Atom,
                   z::Array{<:Unitful.Length, 1})

For bb-center and bf-edge wavelengths, plot average column extinction,
average column destruction probability, mean boundary depth and average
number of packets created at each height.
"""
function plot_radiation(radiation::Radiation,
                        λ::Array{<:Unitful.Length, 1},
                        level::Int64,
                        z::Array{<:Unitful.Length, 1})

    # ===========================================================
    # LOAD RADIATION DATA
    # ===========================================================
    α = radiation.α_continuum
    ε = radiation.ε_continuum
    boundary = radiation.boundary
    packets = radiation.packets

    nλ, nz, nx, ny = size(packets)

    for l =1:nλ
        for i=1:nx
            for j=1:ny
                for k=(boundary[l,i,j]+1):nz
                    packets[l,k,i,j] = 0
                end
            end
        end
    end

    # ===========================================================
    # LOAD ATOM DATA AND GET LINE OPACITY/DESTRUCTION
    # ===========================================================

    mean_boundary = Array{Float64, 1}(undef,nλ)
    max_boundary = Array{Float64, 1}(undef,nλ)
    min_boundary = Array{Float64, 1}(undef,nλ)

    for l=1:nλ
        mean_boundary[l] = mean(boundary[l,:,:])
        max_boundary[l] = maximum(boundary[l,:,:])
        min_boundary[l] = minimum(boundary[l,:,:])
    end

    α = ustrip.(α)
    mean_α = mean(α, dims=(3,4))[:,:,1,1]
    mean_ε = mean(ε, dims=(3,4))[:,:,1,1]
    mean_packets = mean(packets, dims=(3,4))[:,:,1,1]

    z = ustrip.(z .|>u"Mm")
    λ = ustrip.(λ .|>u"nm")

    c = nλ ÷ 2 + 1

    ENV["GKSwstype"]="nul"
    p1 = Plots.plot(z, [mean_α[1,:], mean_α[c,:], mean_α[end,:]],
                    ylabel = "Extinction (m^-1)", xlabel = "z (Mm)",
                    yscale=:log10,
                    label=permutedims(["low", "mid", "edge"]))

    p2 = Plots.plot(z, [mean_ε[1,:], mean_ε[c,:], mean_ε[end,:] ],
                    ylabel = "Destruction", xlabel = "z (Mm)",
                    yscale=:log10,
                    label=permutedims(["low", "mid", "edge"]))

    p3 = Plots.plot(z, [mean_packets[1,:], mean_packets[c,:], mean_packets[end,:]],
                    ylabel = "Mean packets", xlabel = "z (Mm)",
                    label=permutedims(["low", "mid", "edge"]))

    p4 = Plots.plot(λ, [mean_boundary, max_boundary, min_boundary],
                    xlabel = "wavelength (nm)", yflip = true,
                    label=permutedims(["mean", "minimum", "maximum"]))

    Plots.plot(p1, p2, p3, p4, tickfontsize=6, legendfontsize=6, layout=(2,2))
    Plots.png("plots/radiation_"*string(level)*"c")
end


"""
    plot_radiation(radiation::Radiation,
                   atom::Atom,
                   z::Array{<:Unitful.Length, 1})

For bb-center and bf-edge wavelengths, plot average column extinction,
average column destruction probability, mean boundary depth and average
number of packets created at each height.
"""
function plot_radiation(radiation::Radiation,
                        atom::Atom,
                        lines,
                        lineRadiations,
                        z::Array{<:Unitful.Length, 1})

    # ===========================================================
    # LOAD RADIATION DATA
    # ===========================================================
    α_continuum = radiation.α_continuum
    ε_continuum = radiation.ε_continuum
    boundary = radiation.boundary
    packets = radiation.packets
    nλ, nz, nx, ny = size(packets)

    for l =1:nλ
        for i=1:nx
            for j=1:ny
                for k=(boundary[l,i,j]+1):nz
                    packets[l,k,i,j] = 0
                end
            end
        end
    end

    λ = atom.λ
    iλbf = atom.iλbf
    iλbb = atom.iλbb
    n_lines = atom.n_lines
    n_levels = atom.n_levels

    println(iλbf, " ", iλbb)

    # ===========================================================
    # LOAD ATOM DATA AND GET LINE OPACITY/DESTRUCTION
    # ===========================================================

    for ln=1:n_lines

        lineRadiation = lineRadiations[ln]
        line = lines[ln]

        ε_line = lineRadiation.ε_line

        start, stop = iλbb[ln]

        λi = λ[start:stop]
        nλi = length(λi)

        αci = α_continuum[start:stop,:,:,:]
        εci = ε_continuum[start:stop,:,:,:]
        bi = boundary[start:stop,:,:]

        α = Array{PerLength,4}(undef, nλi, nz, nx, ny)
        ε = Array{Float64,4}(undef, nλi, nz, nx, ny)


        for l=1:nλi
            α_line = line_extinction.(λi[l], line.lineData.λ0, line.doppler_width, line.damping_constant, lineRadiation.α_line_constant)
            α[l,:,:,:] = αci[l,:,:,:] .+ α_line
            ε[l,:,:,:] = (εci[l,:,:,:].*αci[l,:,:,:]  .+  ε_line.*α_line) ./ α[l,:,:,:]
        end

        mean_boundary = Array{Float64, 1}(undef,nλi)
        max_boundary = Array{Float64, 1}(undef,nλi)
        min_boundary = Array{Float64, 1}(undef,nλi)

        for l=1:nλi
            mean_boundary[l] = mean(bi[l,:,:])
            max_boundary[l] = maximum(bi[l,:,:])
            min_boundary[l] = minimum(bi[l,:,:])
        end

        α = ustrip.(α)
        mean_α = mean(α, dims=(3,4))[:,:,1,1]
        mean_ε = mean(ε, dims=(3,4))[:,:,1,1]
        mean_packets = mean(packets, dims=(3,4))[:,:,1,1]

        z = ustrip.(z .|>u"Mm")
        λi = ustrip.(λi .|>u"nm")

        c = nλi ÷ 2 + 1

        ENV["GKSwstype"]="nul"
        p1 = Plots.plot(z, [mean_α[1,:], mean_α[c÷2+1,:], mean_α[c,:]],
                        ylabel = "Extinction (m^-1)", xlabel = "z (Mm)",
                        yscale=:log10,
                        label=permutedims(["wing", "mid", "center"]))

        p2 = Plots.plot(z, [mean_ε[1,:], mean_ε[c÷2+1,:], mean_ε[c,:] ],
                        ylabel = "Destruction", xlabel = "z (Mm)",
                        yscale=:log10,
                        label=permutedims(["wing", "mid", "center"]))

        p3 = Plots.plot(z, [mean_packets[1,:], mean_packets[c÷2+1,:], mean_packets[c,:]],
                        ylabel = "Mean packets", xlabel = "z (Mm)",
                        label=permutedims(["wing", "mid", "center"]))

        p4 = Plots.plot(λi, [mean_boundary, max_boundary, min_boundary],
                        xlabel = "wavelength (nm)", yflip = true,
                        label=permutedims(["mean", "maximum", "minimum"]))

        Plots.plot(p1, p2, p3, p4, tickfontsize=6, legendfontsize=6, layout=(2,2))
        Plots.png("plots/radiation_"*string(line.u)*string(line.l))
    end

    for ln=1:n_levels
        start, stop = iλbf[ln]
        λi = λ[start:stop]
        nλi = length(λi)

        α = α_continuum[start:stop,:,:,:]
        ε = ε_continuum[start:stop,:,:,:]
        b = boundary[start:stop,:,:]

        mean_boundary = Array{Float64, 1}(undef,nλi)
        max_boundary = Array{Float64, 1}(undef,nλi)
        min_boundary = Array{Float64, 1}(undef,nλi)

        for l=1:nλi
            mean_boundary[l] = mean(b[l,:,:])
            max_boundary[l] = maximum(b[l,:,:])
            min_boundary[l] = minimum(b[l,:,:])
        end

        α = ustrip.(α)
        mean_α = mean(α, dims=(3,4))[:,:,1,1]
        mean_ε = mean(ε, dims=(3,4))[:,:,1,1]
        mean_packets = mean(packets, dims=(3,4))[:,:,1,1]

        λi = ustrip.(λi .|>u"nm")
        c = nλi ÷ 2 + 1

        println(mean_packets[1,:])
        println(mean_packets[c,:])
        println(mean_packets[end,:])

        ENV["GKSwstype"]="nul"
        p1 = Plots.plot(z, [mean_α[1,:], mean_α[c,:], mean_α[end,:]],
                        ylabel = "Extinction (m^-1)", xlabel = "z (Mm)",
                        yscale=:log10,
                        label=permutedims(["tail", "mid", "edge"]))

        p2 = Plots.plot(z, [mean_ε[1,:], mean_ε[c,:], mean_ε[end,:] ],
                        ylabel = "Destruction", xlabel = "z (Mm)",
                        yscale=:log10,
                        label=permutedims(["tail", "mid", "edge"]))

        p3 = Plots.plot(z, [mean_packets[1,:], mean_packets[c,:], mean_packets[end,:]],
                        ylabel = "Mean packets", xlabel = "z (Mm)",
                        label=permutedims(["tail", "mid", "edge"]))

        p4 = Plots.plot(λi, [mean_boundary, max_boundary, min_boundary],
                        xlabel = "wavelength (nm)", yflip = true,
                        label=permutedims(["mean", "minimum", "maximum"]))

        Plots.plot(p1, p2, p3, p4, tickfontsize=6, legendfontsize=6, layout=(2,2))
        Plots.png("plots/radiation_c"*string(ln))
    end
end

"""
    average_column(array::Array{Real,3})

Given a (nz, nx, ny)-dimensional array,
get the (nz)-dimensional average column.
"""
function average_column(array::Array)
      Statistics.mean(array, dims=[2,3])[:,1,1]
end


"""
    velocity_to_speed(velocity::Array{Array{<:Unitful.Velocity, 1}, 3})

Given a (nz, nx, ny)-dimensional array containing
3-dimensional velocity [vz,vx,vy] arrays,
calculate the (nz, nx,ny)-dimensional speed.
"""
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
