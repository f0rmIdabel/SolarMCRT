#include("../src/mcrt.jl")
#include("../src/radiation.jl") # remove
include("../src/atmosphere.jl")
using Printf

function run(τ_max, target_packets)

    println("\n","-"^42,"\n MCRT calculation in the solar atmosphere \n","-"^42)

    # ==================================================================
    # LOAD WAVELENGTHS
    # ==================================================================
    λ = 499.86u"nm"

    # ==================================================================
    # LOAD ATMOSPHERE DATA
    # ==================================================================
    println("\n--Reading atmosphere model...")
    #atmosphere_parameters = get_atmosphere_data(λ)
    atmosphere_parameters = get_atmosphere_data("bifrost_cb24bih_s385_fullv.ncdf", λ, τ_max)
    atmos = Atmosphere(atmosphere_parameters...)

    # ==================================================================
    # LOAD RADIATION DATA
    # ==================================================================
    #radiation_parameters = get_radiation_data(atmosphere, λ)
    #radiation = Radiation(radiation_parameters...)


    # Find number of packets per box
    S = @timed packets_per_box(atmosphere, λ, target_packets)
    println(S.time)
    # ==================================================================
    # SIMULATION
    # ==================================================================
    #mcrt(atmosphere, boundary, wavelength, S)
    #mcrt(atmosphere, radiation)
end


run(30, 1e7)
