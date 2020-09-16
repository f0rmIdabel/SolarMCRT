module MyPlots

using Plots
import PhysLib
import Bifrost
using Unitful
using Printf

"""
function surface_plot(Atmosphere::Module, tau_max::Real, camera_tilt::Real) 

"""

function boundary_height(Atmosphere::Module, tau_max::Real, camera_tilt::Real) 

    #Plots.pyplot()

    # To avoid ssh display problems 
    ENV["GKSwstype"]="nul"

    height = Atmosphere.height
    chi = Atmosphere.chi_continuum

    b = PhysLib.optical_depth_boundary(chi, height, tau_max)
    nx, ny = size(b)

    x = 1:nx
    y = 1:ny
    f(x,y) = ustrip(height[b[x,y]])

    surface(x, y, f, zlim = [ustrip(height[end]), ustrip(height[1])], camera=(-45,camera_tilt))

    fig = @sprintf("/mn/stornext/u3/idarhan/SolarMCRT/Results/Plots/Boundary/boundary_%.1f_%g", tau_max, camera_tilt)
    png(fig)

end



function surface_intensity(surface::Array{Int, 2}, tau_max::Real, total_packets::Real)
    
    # To avoid ssh display problems 
    ENV["GKSwstype"]="nul"

    heatmap(1:size(surface,1), 1:size(surface,2), surface, c=:grays, title = "", aspect_ratio=:equal)
    fig = @sprintf("/mn/stornext/u3/idarhan/SolarMCRT/Results/Plots/Surface/bf_%.1f_%.0e", tau_max, total_packets)
    png(fig)

end

end
