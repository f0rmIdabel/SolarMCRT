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

    * Wavelengths, λ                            (Float64 [nm], nλ)
    * Source function, S                        (Int32, nλ, nz*, nx, ny)
    * Mean radiation field, J                   (Int32, nλ, nz*, nx, ny)
    * Surface intensity, surface_intensity      (Int32, nx, ny, phi_bins, theta_bins)
    * Surface boundary, boundary                (Int32, nx, ny)
    * Energy per packet, intensity_per_packet   (Float64 [kW / m^2 / sr / nm], nλ)
    * Total packets, total_packets              (Int64, nλ)
    * Total destroyed packets, total_destroyed  (Int64)
    * Total escaped packets, total_escaped      (Int64)
    * Total scatterings, total_scatterings      (Int64)

For each wavelength λ, this is collected in the file *output_λ.hdf5* in the out/ directory. For one wavelength in a ~500x500x450 atmosphere with no boundary cut-off, this will be around 2 GBs of data.
