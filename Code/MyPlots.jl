module MyPlots

using Plots
import PhysLib
import Bifrost
using Unitful
using Printf

"""
function surface_plot(Atmosphere::Module, tau_max::Real, camera_tilt::Real)

"""

function boundary_height(Atmosphere::Module, τ_max::Real, camera_tilt::Real)

    #Plots.pyplot()

    # To avoid ssh display problems
    ENV["GKSwstype"]="nul"

    height = Atmosphere.height
    χ = Atmosphere.chi_continuum

    b = PhysLib.optical_depth_boundary(χ, height, τ_max)
    nx, ny = size(b)

    x = 1:nx
    y = 1:ny
    f(x,y) = ustrip(height[b[x,y]])

    surface(x, y, f, zlim = [ustrip(height[end]), ustrip(height[1])], camera=(-45,camera_tilt))

    fig = @sprintf("/mn/stornext/u3/idarhan/SolarMCRT/Results/Plots/Boundary/boundary_%.1f_%g", τ_max, camera_tilt)
    png(fig)

end


"""
"""
function surface_intensity(surface::Array{Int, 4}, τ_max::Real, total_packets::Real, bin=:[:,:])

    surface = PhysLib.extract_surface_bin(surface, bin)

    # To avoid ssh display problems
    ENV["GKSwstype"]="nul"

    heatmap(1:size(surface,1), 1:size(surface,2), surface, c=:grays, aspect_ratio=:equal)
    plot!(size=(410,400))
    fig = @sprintf("/mn/stornext/u3/idarhan/SolarMCRT/Results/Plots/Surface/bf_tau%.1f_pcts%.0e_bin%s",
                   τ_max, total_packets, string(bin))
    png(fig)

end


function escape_direction(surface::Array{Int, 4}, τ_max::Real, total_packets::Real)

    ϕ_bins, θ_bins = size(surface)[3:4]

    ϕ = [0, (2π * (1:(ϕ_bins-1)) / ϕ_bins)...]
    θ = [0, (π/2 * (2 .^(1:(θ_bins-1))) / 2^θ_bins)...]

    ϕ_hits = Array{Int, 1}(undef, ϕ_bins)  # or use zeros maybe
    θ_hits = Array{Int, 1}(undef, θ_bins)

    # ϕ direction
    for i=1:ϕ_bins
        ϕ_hits[i] = sum(PhysLib.extract_surface_bin(surface, :[$i,:]))
    end

    # θ direction
    for i=1:θ_bins
        θ_hits[i] = sum(PhysLib.extract_surface_bin(surface, :[:,$i]))
    end

    h1 = bar(ϕ, ϕ_hits, xlabel = "ϕ")
    h2 = bar(θ, θ_hits, xlabel = "θ")

    plot(h1, h2, layout = (1, 2), legend = false)

    fig = @sprintf("/mn/stornext/u3/idarhan/SolarMCRT/Results/Plots/EscapeDirections/escapes_tau%.1f_pcts%.0e",
                   τ_max, total_packets)
    png(fig)
end

end
