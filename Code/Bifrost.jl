module Bifrost 

import IOLib
using Unitful

PATH = "/mn/stornext/d18/RoCS/idarhan/MScProject/Renditions/R1/data/Bifrost/"
dim = [504, 504, 117]
wavelength = 500u"nm"

x_ = IOLib.read_parameter(PATH*"x_Bifrost.txt")u"m"
x = push!(x_, 2*x_[dim[1]] - x_[dim[1]-1]) # Add endpoint

y_ = IOLib.read_parameter(PATH*"y_Bifrost.txt")u"m"
y = -push!(y_, 2*y_[dim[2]] - y_[dim[2]-1]) # Add endpoint

edge = [x[1] x[end]
        y[1] y[end]]

height_ = IOLib.read_parameter(PATH*"height_Bifrost.txt")u"m"
height = push!(height_, 2*height_[dim[3]] - height_[dim[3]-1]) # Add endpoint

intensity = IOLib.read_parameter(PATH*"intensity_Bifrost.txt")u"kW / m^2 / sr / nm"

temperature = IOLib.read_parameter(PATH*"temperature_Bifrost.txt")u"K"
epsilon_continuum = IOLib.read_parameter(PATH*"epsilon_continuum_Bifrost.txt")
chi_continuum = IOLib.read_parameter(PATH*"chi_continuum_Bifrost.txt")u"m^-1"

#electron_density = IOLib.read_parameter(PATH*"electron_density_Bifrost.txt")u"m^-3"
#velocity_z = IOLib.read_parameter(PATH*"velocity_z_Bifrost.txt")u"m / s^1"

end