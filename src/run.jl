include("mcrt.jl")
include("tools.jl")
include("atmos.jl")

function run(τ_max=nothing, target_packets=nothing)

    println("\n","-"^42,"\n MCRT calculation in the solar atmosphere \n","-"^42)
    println("\n--Reading atmosphere model...")

    # ==================================================================
    # LOAD ATMOSPHERE DATA
    # ==================================================================
    parameters = get_atmosphere_data("bifrost_cb24bih_s385_fullv.ncdf",
                                     "output_ray.hdf5")
    atmosphere = Atmosphere(parameters...)

    # ==================================================================
    # CHOOSE PARAMETERS
    # ==================================================================
    # Choose wavelengths
    wavelength = 499.86u"nm" #################################################################

    if τ_max == nothing
        # Get user input: τ_max, # packets
        print("--Choose maximum optical depth: ")
        τ_max = readline()
        τ_max = parse(Float64, τ_max)
    else
        println(@sprintf("--Maximum optical depth set to %.1f.", τ_max))
    end

    if target_packets == nothing
        print("--Choose number of photon packets: ")
        target_packets = readline()
        target_packets = parse(Float64, target_packets)
    else
        println(@sprintf("--Generating %.1e photon packets.", target_packets))
    end

    # ==================================================================
    # PRE-CALCULATIONS
    # ==================================================================
    println("--Setting up simulation...")

    # Find boundary for given τ_max
    boundary = optical_depth_boundary(atmosphere.χ_continuum, atmosphere.z, τ_max)

    # Find number of packets per box
    S = packets_per_box(atmosphere.x, atmosphere.y, atmosphere.z,
                        atmosphere.χ_continuum, atmosphere.temperature,
                        wavelength, target_packets, boundary)

    # ==================================================================
    # SIMULATION
    # ==================================================================
    mcrt(atmosphere, boundary, wavelength, S)
end


run(30, 1e8)
