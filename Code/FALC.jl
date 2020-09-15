module FALC 

import IOLib
using Unitful

PATH = "/mn/stornext/d18/RoCS/idarhan/MScProject/Renditions/R1/data/FALC/"

dim = [5, 5, 68]
wavelength = 500u"nm"


# Atmosphere dimension
x_ = IOLib.read_parameter(PATH*"x_FALC.txt")u"m" * 100_000
x = push!(x_, 2*x_[dim[1]] - x_[dim[1]-1]) # Add enddpoint

y_ = IOLib.read_parameter(PATH*"y_FALC.txt")u"m" * 100_000
y = push!(y_, 2*y_[dim[2]] - y_[dim[2]-1]) # Add enddpoint

edge = [x[1] x[end]
        y[1] y[end]]

height_ = IOLib.read_parameter(PATH*"height_FALC.txt")u"m"
height = push!(height_, 2*height_[dim[3]] - height_[dim[3]-1]) # Add enddpoint

intensity = IOLib.read_parameter(PATH*"intensity_FALC.txt")u"kW / m^2 / sr / nm"

temperature = IOLib.read_parameter(PATH*"temperature_FALC.txt")u"K"
epsilon_continuum = IOLib.read_parameter(PATH*"epsilon_continuum_FALC.txt")
chi_continuum = IOLib.read_parameter(PATH*"chi_continuum_FALC.txt")u"m^-1"
velocity_z = IOLib.read_parameter(PATH*"velocity_z_FALC.txt")u"ms^-1"
electron_density = IOLib.read_parameter(PATH*"electron_density_FALC.txt")u"m^-3"


end