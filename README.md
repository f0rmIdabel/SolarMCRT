# SolarMCRT

**This is not up to date.**

Simple Monte Carlo radiative transfer code to be used on solar atmosphere models.




## Input
The simulation parameters can be changed in the file *keywords.input* in the run/ directory.

The code requires an atmosphere file containing

    * Dimensions
        * x (nx)
        * y (ny)
        * z (nz)
    * Velocities
        * velocity_x (nz, nx, ny)
        * velocity_y (nz, nx, ny)
        * velocity_z (nz, nx, ny)
    * Temperature
        * temperature (nz, nx, ny)
    * Densities
        * electron_density (nz, nx, ny)
        * hydrogen_populations (nz, nx, ny, levels)

When ran in full *atom mode*,  the program needs a two-level atom file containing

	* Atomic weight, atom_weight
	* Ion charge, Z
	* Ionisation energy, chi_inf
	* Ground level energy, chi_l
	* Second level energy, chi_u
	* Ground level statistical weight, gl
	* Second level statistical weight, gu
	* Oscillator strength, f_value

For other atoms than hydrogen, the initial atom populations is required as well.

## Running the code
With appropriate input parameters, the code can be executed from the run/ directory with *julia run.jl*.


## Output
The code outputs

    * Wavelengths, λ [nm]       						                (Float64, nλ)
	* Opacity, chi [m^-1]
	* Destruction probability
    * Source function, S                    						    (Int32, nλ, nz, nx, ny)
    * Mean radiation field, J            						        (Int32, nλ, nz, nx, ny)
    * Surface boundary, boundary                                        (Int32, nx, ny)
    * Intensity per packet, intensity_per_packet [kW / m^2 / sr / nm]   (Float64, nλ)
    * Total packets, total_packets              (Int64, nλ)
    * Total destroyed packets, total_destroyed  (Int64)
    * Total escaped packets, total_escaped      (Int64)
    * Total scatterings, total_scatterings      (Int64)

For each wavelength λ, this is collected in the file *output_λ.hdf5* in the out/ directory. For one wavelength in a ~500x500x450 atmosphere with no boundary cut-off, this will be around 2 GBs of data.
