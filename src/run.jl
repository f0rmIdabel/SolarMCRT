include("mcrt.jl")
include("tools.jl")

function run(τ_max=nothing, target_packets=nothing)

    println("\n","-"^42,"\n MCRT calculation in the solar atmosphere \n","-"^42)
    # ==================================================================
    # LOAD ATMOSPHERE DATA
    # ================================================================== MOVE THIS TO MCRT file

    println("\n--Reading atmosphere model...")

    parameters = get_atmosphere_data("bifrost_cb24bih_s385_fullv.ncdf", 499.86u"nm")
    atmosphere = Atmosphere(parameters...)

    # ==================================================================
    # CHOOSE PARAMETERS
    # ================================================================== MOVE THIS TO INPUT FILE
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
    # ================================================================== MOVE THIS TO DIFFERENT FILE
    println("--Setting up simulation...")

    # Find boundary for given τ_max
    boundary = optical_depth_boundary(atmosphere.χ, atmosphere.z, τ_max)

    # Find number of packets per box
    S = packets_per_box(atmosphere.x, atmosphere.y, atmosphere.z,
                        atmosphere.χ, atmosphere.temperature,
                        wavelength, target_packets, boundary)

    # ==================================================================
    # SIMULATION
    # ==================================================================
    mcrt(atmosphere, boundary, wavelength, S)
end


run(30, 1e9)
