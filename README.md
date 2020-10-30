# SolarMCRT

Simple Monte Carlo radiative transfer code to be used on solar atmosphere models.


## Input
The simulation parameters can be changed in the file *keywords.input* in the run/ directory, while the wavelengths are contained in *wavelengths.input*.

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

## Running the code
With appropriate input parameters, the code can be executed from the run/ directory with *julia run.jl*.


## Output
The code outputs

    * Source function, S                        (nlambda, nx, ny, nz*)
    * Mean radiation field, J                   (nlambda, nx, ny, nz*)
    * Surface intensity, surface_intensity      (nlambda, nx, ny, nz*, phi_bins, theta_bins)
    * Total packets, total_packets              (nlambda)
    * Total destroyed packets, total_destroyed  (nlambda)
    * Total escaped packets, total_escaped      (nlambda)
    * Total scatterings, total_scatterings      (nlambda)

This is collected in the file *output.hdf5* in the out/ directory. For on wavelength in a ~500x500x450 atmosphere with no boundary cut-off, this will be around 2 GBs of data.
