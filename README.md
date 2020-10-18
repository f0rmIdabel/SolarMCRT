# SolarMCRT

Simple Monte Carlo radiative transfer code to be used on solar atmosphere models.

## Running the code

With an appropriate atmosphere file, the code can be executed from the run/ directory with a simple *julia run.jl*. All simulation parameters can be changed in the file *keywords.input*.

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

This is collected in the file output.hdf5 in the out/ directory.
