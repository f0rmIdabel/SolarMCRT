include("physLib.jl")

using Printf
using Plots
using UnitfulRecipes

function plot_boundary_height(z::Array{<:Unitful.Length, 1},
                         χ::Array{<:Unitful.Quantity{<:Real, Unitful.𝐋^(-1)}, 3},
                         τ_max::Real, camera_tilt::Real)

    #Plots.pyplot()

    # To avoid ssh display problems
    ENV["GKSwstype"]="nul"

    b = optical_depth_boundary(χ, z, τ_max)
    nx, ny = size(b)

    x = 1:nx
    y = 1:ny
    f(x,y) = ustrip(z[b[x,y]])

    surface(x, y, f, zlim = [ustrip(z[end]),
            ustrip(z[1])], camera=(-45,camera_tilt))

    fig = @sprintf("/mn/stornext/u3/idarhan/SolarMCRT/Results/Plots/Boundary/boundary_%.1f_%g",
                    τ_max, camera_tilt)
    png(fig)
end


function plot_surface_intensity(surface::Array{Int64, 4}, τ_max::Real,
                           total_packets::Real, bin=:[:,:])

    surface = extract_surface_bin(surface, bin)

    # To avoid ssh display problems
    ENV["GKSwstype"]="nul"

    heatmap(1:size(surface,1), 1:size(surface,2), surface, c=:grays, aspect_ratio=:equal)
    plot!(size=(410,400))
    fig = @sprintf("/mn/stornext/u3/idarhan/SolarMCRT/Results/Plots/Surface/bf_tau%.1f_pcts%.0e_bin%s",
                   τ_max, total_packets, string(bin))
    png(fig)
end


function plot_escape_direction(surface::Array{Int64, 4}, τ_max::Real, total_packets::Real)

    ϕ_bins, θ_bins = size(surface)[3:4]

    ϕ = [0, (2π * (1:(ϕ_bins-1)) / ϕ_bins)...]
    θ = [0, (π/2 * (2 .^(1:(θ_bins-1))) / 2^θ_bins)...]

    ϕ_hits = Array{Int, 1}(undef, ϕ_bins)
    θ_hits = Array{Int, 1}(undef, θ_bins)

    # ϕ direction
    for i=1:ϕ_bins
        ϕ_hits[i] = sum(extract_surface_bin(surface, :[$i,:]))
    end

    # θ direction
    for i=1:θ_bins
        θ_hits[i] = sum(extract_surface_bin(surface, :[:,$i]))
    end

    h1 = bar(ϕ, ϕ_hits, xlabel = "ϕ")
    h2 = bar(θ, θ_hits, xlabel = "θ")

    plot(h1, h2, layout = (1, 2), legend = false)

    fig = @sprintf("/mn/stornext/u3/idarhan/SolarMCRT/Results/Plots/EscapeDirections/escapes_tau%.1f_pcts%.0e",
                   τ_max, total_packets)
    png(fig)
end


function plot_thread_performance(threads::Array{Int64, 1}, time::Array{Float64, 1})
    plot(threads, time, legend=false)
    xlabel!("Threads")
    ylabel!("Time [s]")
    fig = "/mn/stornext/u3/idarhan/SolarMCRT/Results/Plots/time_thread"
    png(fig)
end


function plot_slice(slice, y::Int64, save::Bool)

    heatmap(1:size(slice,1), 1:size(slice,2), slice, c=:grays, aspect_ratio=:equal)
    plot!(size=(410,400))
    if save
        fig = @sprintf("/mn/stornext/u3/idarhan/SolarMCRT/Results/Plots/FieldSlice/field_y%d", y)
        png(fig)
    end
end


function traverse_field_gif(field, ny::Int64)

    anim = @animate for y=1:ny
        slice = field[:,y,:]
        plot_slice(slice, y,false)
    end

    gif(anim, "/mn/stornext/u3/idarhan/SolarMCRT/Results/Plots/anim_fps15.gif", fps = 15)
end
