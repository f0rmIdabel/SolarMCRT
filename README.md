# SolarMCRT

**This is not up to date.**

Simple Monte Carlo radiative transfer code to be used on solar atmosphere models.




## Input
The simulation parameters can be changed in the file *keywords.input* in the run/ directory.

The code requires an atmosphere file containing

    * Dimensions
        * x [m] (nx)
        * y [m] (ny)
        * z [m] (nz)
    * Velocities
        * velocity_x [m/s] (nz, nx, ny)
        * velocity_y [m/s] (nz, nx, ny)
        * velocity_z [m/s] (nz, nx, ny)
    * Temperature
        * temperature [K] (nz, nx, ny)
    * Densities
        * electron_density [m^-3] (nz, nx, ny)
        * hydrogen_populations [m^-3] (nz, nx, ny, levels)

When ran in full *atom mode*,  the program needs a two-level atom file containing

	* Physical properties
		* Atomic weight, atom_weight [kg]
		* Ion charge, Z
		* Ground level wavenumber, chi_l [cm^-1]
		* Second level wavenumber, chi_u [cm^-1]
		* Ionisation wavenumber, chi_inf [cm^-1]
		* Ground level statistical weight, gl
		* Second level statistical weight, gu
		* Ionised level statistical weight, ginf
		* Oscillator strength, f_value
	* λ-sampling parameters
		* Doppler width core, qcore (?)
		* Doppler width wing, qwing (?)
		* Minimum λ to be sampled for lower bf-transition, bfl_min [nm]
		* Minimum λ to be sampled for upper bf-transition, bfl_min [nm]

The initial atom populations are required as well.

	* Initial atom populations, initial_populations [m^-3] (nz, nx, ny, levels=3)

## Running the code
With appropriate input parameters, the code can be executed from the run/ directory with *julia run.jl*.


## Output
The code outputs

    * Wavelengths, λ [nm] (nλ)
	* Opacity, chi [m^-1] (nλ, nz, nx, ny)
	* Destruction probability, epsilon (nλ, nz, nx, ny)
    * Packet distribution, packets (nλ, nz, nx, ny)
	* Sourface boundary, boundary (nλ, nx, ny)
	* Intensity per packet, intensity_per_packet [kW / m^2 / sr / nm] (nλ)

    * Mean radiation field, J (nλ, nz, nx, ny)
    * Total destroyed packets, total_destroyed (nλ)
    * Total scatterings, total_scatterings (nλ)

The get the radiation field in units of intensity, you need to multiply it by the intensity_per_packet variable. All output is collected in the file *output.h5* in the out/ directory. For one wavelength in a ~500x500x450 atmosphere with no cut-offs, this will be around 2 GBs of data.
