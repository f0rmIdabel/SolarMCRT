# SolarMCRT

Simple Monte Carlo radiative transfer code to be used on solar atmosphere models.


## Input
The simulation parameters can be changed in the file *keywords.input* in the run/ directory, while the wavelengths are contained in *wavelengths.input*.

The code requires an atmosphere file that contains

    * Dimensions
        * *x* (nx)
        * *y* (ny)
        * *z* (nz)
    * Velocities
        * *velocity_x* (nz, nx, ny)
        * *velocity_y* (nz, nx, ny)
        * *velocity_z* (nz, nx, ny)
    * Temperature
        * *temperature* (nz, nx, ny)
    * Densities
        * *electron_density* (nz, nx, ny)
        * *hydrogen_populations* (nz, nx, ny, levels)

## Running the code
With appropriate input parameters, the code can be executed from the run/ directory with a simple *julia run.jl*.


## Output
The code outputs

    * Source function, S (nx, ny, nz*)
    * Mean radiation field, J (nx, ny, nz*)
    * Surface intensity, surface_intensity (nx, ny, nz*, phi_bins, theta_bins)
    * Total packets, total_packets
    * Total destroyed, total_destroyed
    * Total escaped, total_escaped
    * Total scatterings, total_scatterings
    * Signal to noise ratio, SNR

This is collected in the file *output.hdf5* in the out/ directory.
