include("../src/mcrt.jl")

function run()

    println("\n","-"^42,"\n MCRT calculation in the solar atmosphere \n","-"^42)

    # ==================================================================
    # LOAD WAVELENGTHS
    # ==================================================================
    println("\n--Loading wavelengths...")
    λ = 499.86u"nm"

    # ==================================================================
    # LOAD ATMOSPHERE DATA AND CALCULATE BOUNDARY
    # ==================================================================
    println("--Loading atmosphere data...")
    atmosphere_parameters = collect_atmosphere_data(λ)
    atmosphere = Atmosphere(atmosphere_parameters...)

    # ==================================================================
    # LOAD RADIATION DATA
    # ==================================================================
    println("--Loading radiation data...")
    radiation_parameters = collect_radiation_data(atmosphere, λ)
    radiation = Radiation(λ, radiation_parameters...)

    # ==================================================================
    # SIMULATION
    # ==================================================================
    mcrt(atmosphere, radiation)
end


run()
