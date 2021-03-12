# SolarMCRT

**This README is not up to date :)**

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
        * hydrogen_populations [m^-3] (nz, nx, ny, levels=3)

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
		* Partition function neutral atom, U0
		* Partition function ionised atom, U1
		* Oscillator strength, f_value
		* Atom density, density
	* λ-sampling parameters
		* Doppler width core, qcore
		* Doppler width wing, qwing
		* Minimum λ to be sampled for lower bf-transition, bfl_min [nm]
		* Minimum λ to be sampled for upper bf-transition, bfl_min [nm]


## Running the code
With appropriate input parameters, the code can be executed from the run/ directory with *julia run.jl*.


## Output
The code outputs

	* Simulation result from every iteration
    	* Mean radiation field, J (ni, nλ, nz, nx, ny)
    	* Total destroyed packets, total_destroyed (ni, nλ)
    	* Total scatterings, total_scatterings (ni, nλ)
		* Atom populations, populations (ni, nz, nx, ny, nl)
		* Time, time (ni, nλ)
	* Radiation data from last iteration
    	* Wavelengths, wavelengths [nm] (nλ)
		* Sourface boundary, boundary (nλ, nx, ny)
    	* Packet distribution, packets (nλ, nz, nx, ny)
		* Intensity per packet, intensity_per_packet [kW / m^2 / sr / nm] (nλ)

To get the radiation field in units of intensity, you need to multiply it by the intensity_per_packet variable. All output is collected in the file *output.h5* in the out/ directory. For one wavelength in a ~500x500x450 atmosphere with no cut-offs, this will be around X GBs of data.

## Potential problems

	* ...
