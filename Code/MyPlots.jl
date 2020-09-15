module MyPlots

using Plots
import PhysLib
import Bifrost
using Unitful
using Printf

"""
function surface_plot(Atmosphere::Module, tau_max::Real, camera_tilt::Real) 

"""

function surface_plot(Atmosphere::Module, tau_max::Real, camera_tilt::Real) 

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


#surface_plot(Bifrost, 5, 30)




end
