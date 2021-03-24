using DelimitedFiles
using BenchmarkTools
using ProgressMeter
using Transparency
using StaticArrays
using AtomicData
using Unitful
using Random
using Printf
using Test
using HDF5

import PhysicalConstants.CODATA2018: c_0, h, k_B, m_u, m_e, R_∞, ε_0, e
const E_∞ = R_∞ * c_0 * h
const hc = h * c_0

@derived_dimension NumberDensity Unitful.𝐋^-3
@derived_dimension PerLength Unitful.𝐋^-1
@derived_dimension PerArea Unitful.𝐋^-2
@derived_dimension PerTime Unitful.𝐓^-1
@derived_dimension UnitsIntensity_λ Unitful.𝐋^-1 * Unitful.𝐌 * Unitful.𝐓^-3

"""
    test_mode()

Get the test mode status.
"""
function background_mode()
    input_file = open(f->read(f, String), "../run/keywords.input")
    i = findfirst("background_mode", input_file)[end] + 1
    file = input_file[i:end]
    i = findfirst("=", file)[end]
    j = findfirst("\n", file[i+1:end])[end] + i
    tm = parse(Bool, file[i+1:j-1])
    return tm
end


"""
    get_output_path()

Get the path to the output file.
"""
function get_output_path()

    if background_mode()
        path = "../out/output_" * string(ustrip.(get_background_λ())) * "nm_" * string(get_target_packets()) * "pcs.h5"
    else
        nλ_bb, nλ_bf = get_nλ()
        nλ_bb += 1-nλ_bb%2
        nλ = 2nλ_bf + nλ_bb
        pop_distrib = get_population_distribution()
        if pop_distrib == "LTE"
            d = "_LTE"
        elseif pop_distrib == "zero_radiation"
            d = "_ZR"
        end
        path = "../out/output_nw" * string(nλ) * "_" * string(get_target_packets()) * "pcs" * d * ".h5"
    end

    return path
end

"""
    get_iteration_mode()

Get the iteration mode (populations or lambda)
"""
function get_iteration_mode()
    input_file = open(f->read(f, String), "../run/keywords.input")
    i = findfirst("iteration_mode", input_file)[end] + 1
    file = input_file[i:end]
    i = findfirst("=", file)[end]
    j = findfirst("\n", file[i+1:end])[end] + i
    im = strip(file[i+1:j-1])
    return im
end

"""
    get_max_iterations()

Get the maximum number of iterations.
"""
function get_max_iterations()
    input_file = open(f->read(f, String), "../run/keywords.input")
    i = findfirst("max_iterations", input_file)[end] + 1
    file = input_file[i:end]
    i = findfirst("=", file)[end] + 1
    j = findfirst("\n", file)[end] - 1
    max_iterations = parse(Int64, file[i:j])
    return max_iterations
end

# =============================================================================
# RADIATION
# =============================================================================

"""
    get_target_packets()

Get the number of packets to be created for each wavelength.
"""
function get_target_packets()
    input_file = open(f->read(f, String), "../run/keywords.input")
    i = findfirst("target_packets", input_file)[end] + 1
    file = input_file[i:end]
    i = findfirst("=", file)[end] + 1
    j = findfirst("\n", file)[end] - 1
    target_packets = parse(Float64, file[i:j])
    return target_packets # move
end

"""
    get_max_scatterings()

Get the maximum number of scatterings
before moving on to the next photon.
"""
function get_max_scatterings()
    input_file = open(f->read(f, String), "../run/keywords.input")
    i = findfirst("max_scatterings", input_file)[end] + 1
    file = input_file[i:end]
    i = findfirst("=", file)[end] + 1
    j = findfirst("\n", file)[end] - 1
    max_scatterings = parse(Float64, file[i:j])
    return max_scatterings # move
end

"""
    get_cut_off()

Get optical depth where the atmosphere will be cut.
"""
function get_cut_off()
    input_file = open(f->read(f, String), "../run/keywords.input")
    i = findfirst("cut_off", input_file)[end] + 1
    file = input_file[i:end]
    i = findfirst("=", file)[end] + 1
    j = findfirst("\n", file)[end] - 1

    cut_off = nothing

    try
        cut_off = parse(Float64, file[i:j])
    catch
        cut_off = parse(Bool, file[i:j])
    end

    return cut_off # move
end

"""
    get_background_λ()

Get the wavelength where to perform the background radiation MC test.
"""
function get_background_λ()
    input_file = open(f->read(f, String), "../run/keywords.input")
    i = findfirst("background_wavelength", input_file)[end] + 1
    file = input_file[i:end]
    i = findfirst("=", file)[end] + 1
    j = findfirst("\n", file)[end] - 1

    λ = parse(Float64, file[i:j])u"nm"

    return λ
end

"""
    get_Jλ(output_path::String, iteration::Int64, intensity_per_packet::Array{<:UnitsIntensity_λ,1})

Get the radiation field from the output file.
"""
function get_Jλ(output_path::String, iteration::Int64, intensity_per_packet::Array{<:UnitsIntensity_λ,1})
    J = nothing
    h5open(output_path, "r") do file
        J = read(file, "J")[iteration,:,:,:,:]
    end

    nλ, nz, nx, ny = size(J)
    Jλ = Array{UnitsIntensity_λ,4}(undef, nλ, nz, nx, ny)

    for l=1:nλ
        Jλ[l,:,:,:] = intensity_per_packet[l] .* J[l,:,:,:]
    end

    return Jλ
end

"""
    get_nλ()

Get the number of bound-bound and bound-free wavelengths to sample.
"""
function get_nλ()
    input_file = open(f->read(f, String), "../run/keywords.input")
    i = findfirst("nλ_bb", input_file)[end] + 1
    file = input_file[i:end]
    i = findfirst("=", file)[end] + 1
    j = findfirst("\n", file)[end] - 1
    nλ_bb = parse(Int64, file[i:j])

    i = findfirst("nλ_bf", input_file)[end] + 1
    file = input_file[i:end]
    i = findfirst("=", file)[end] + 1
    j = findfirst("\n", file)[end] - 1
    nλ_bf = parse(Int64, file[i:j])
    return nλ_bb, nλ_bf
end

# =============================================================================
# ATOM
# =============================================================================

"""
    get_atom_path()

Get path to atom file.
"""
function get_atom_path()
    input_file = open(f->read(f, String), "../run/keywords.input")
    i = findfirst("atom_path", input_file)[end] + 1
    file = input_file[i:end]
    i = findfirst("\"", file)[end]
    j = findfirst("\"", file[i+1:end])[end] + i
    atmosphere_path = string(file[i+1:j-1])
    return atmosphere_path
end

"""
    get_population_distribution()

Get the initial population configuration (LTE or zero-radiation)
"""
function get_population_distribution()
    input_file = open(f->read(f, String), "../run/keywords.input")
    i = findfirst("population_distribution", input_file)[end] + 1
    file = input_file[i:end]
    i = findfirst("\"", file)[end]
    j = findfirst("\"", file[i+1:end])[end] + i
    distribution = string(file[i+1:j-1])
    return distribution
end


"""
    test_mode()

Check whether to write rates or not.
"""
function get_write_rates()
    input_file = open(f->read(f, String), "../run/keywords.input")
    i = findfirst("write_rates", input_file)[end] + 1
    file = input_file[i:end]
    i = findfirst("=", file)[end]
    j = findfirst("\n", file[i+1:end])[end] + i
    tm = parse(Bool, file[i+1:j-1])
    return tm
end

# =============================================================================
# ATMOSPHERE
# =============================================================================
"""
    get_atmosphere_path()

Fetch the atmosphere path.
"""
function get_atmosphere_path()
    input_file = open(f->read(f, String), "../run/keywords.input")
    i = findfirst("atmosphere_path", input_file)[end] + 1
    file = input_file[i:end]
    i = findfirst("\"", file)[end]
    j = findfirst("\"", file[i+1:end])[end] + i
    atmosphere_path = string(file[i+1:j-1])
    return atmosphere_path
end

"""
    get_step()

Fetch the step-size in the z, x and y-direction.
"""
function get_step()
    input_file = open(f->read(f, String), "../run/keywords.input")
    i = findfirst("step", input_file)[end] + 1
    file = input_file[i:end]
    i = findfirst("[", file)[end] + 1
    j = findfirst(",", file)[end] - 1
    dz = parse(UInt16, file[i:j])

    i = j + 2
    file = file[i:end]
    j = findfirst(",", file)[end] - 1
    dx = parse(UInt16, file[1:j])

    i = j + 2
    file = file[i:end]
    j = findfirst("]", file)[end] - 1
    dy = parse(UInt16, file[1:j])
    steps = [dz, dx, dy]
    return steps
end

"""
    get_stop()

Fetch the stop indices in the z, x and y-direction.
"""
function get_stop()

    nz = nothing
    nx = nothing
    ny = nothing

    input_file = open(f->read(f, String), "../run/keywords.input")
    i = findfirst("stop", input_file)[end] + 1
    file = input_file[i:end]
    i = findfirst("[", file)[end] + 1
    j = findfirst(",", file)[end] - 1

    try
        nz = parse(UInt16, file[i:j])
    catch
        nz = nothing
    end

    i = j + 2
    file = file[i:end]
    j = findfirst(",", file)[end] - 1

    try
        nx = parse(UInt16, file[1:j])
    catch
        nx = nothing
    end

    i = j + 2
    j = findfirst("]", file)[end] - 1

    try
        ny = parse(UInt16, file[1:j])
    catch
        ny = nothing
    end

    stop = [nz, nx, ny]

    return stop
end

"""
    get_start()

Fetch the start indices in the z, x and y-direction.
"""
function get_start()
    input_file = open(f->read(f, String), "../run/keywords.input")
    i = findfirst("start", input_file)[end] + 1
    file = input_file[i:end]
    i = findfirst("[", file)[end] + 1
    j = findfirst(",", file)[end] - 1

    nz = parse(UInt16, file[i:j])

    i = j + 2
    file = file[i:end]
    j = findfirst(",", file)[end] - 1

    nx = parse(UInt16, file[1:j])

    i = j + 2
    file = file[i:end]
    j = findfirst("]", file)[end] - 1

    ny = parse(UInt16, file[1:j])

    start = [nz, nx, ny]

    return start
end

# =============================================================================
# OUTPUT FILE
# =============================================================================

"""
    create_output_file(output_path::String, max_iterations::Int64, nλ::Int64, atmosphere_size::Tuple, write_rates::Bool)
Initialise all output variables for the full atom mode.
"""
function create_output_file(output_path::String, max_iterations::Int64, nλ::Int64, atmosphere_size::Tuple, write_rates::Bool)

    nz, nx, ny = atmosphere_size

    h5open(output_path, "w") do file
        write(file, "J", Array{Int32,5}(undef, max_iterations, nλ, nz, nx,ny))
        write(file, "total_destroyed", Array{Int64,2}(undef, max_iterations, nλ))
        write(file, "total_scatterings", Array{Int64,2}(undef,max_iterations, nλ))
        write(file, "time", Array{Float64,2}(undef,max_iterations, nλ))

        write(file, "packets", Array{Int32,5}(undef, max_iterations, nλ, nz, nx, ny))
        write(file, "boundary", Array{Int32,4}(undef,max_iterations, nλ, nx, ny))
        write(file, "intensity_per_packet", Array{Float64,2}(undef,max_iterations, nλ))

        write(file, "populations", Array{Float64,5}(undef, max_iterations+1, nz, nx, ny, 3))

        if write_rates
            write(file, "R12", Array{Float64}(undef, max_iterations+1, nz,nx,ny))
            write(file, "R13", Array{Float64}(undef, max_iterations+1, nz,nx,ny))
            write(file, "R23", Array{Float64}(undef, max_iterations+1, nz,nx,ny))
            write(file, "R21", Array{Float64}(undef, max_iterations+1, nz,nx,ny))
            write(file, "R31", Array{Float64}(undef, max_iterations+1, nz,nx,ny))
            write(file, "R32", Array{Float64}(undef, max_iterations+1, nz,nx,ny))
            write(file, "C12", Array{Float64}(undef, max_iterations+1, nz,nx,ny))
            write(file, "C13", Array{Float64}(undef, max_iterations+1, nz,nx,ny))
            write(file, "C23", Array{Float64}(undef, max_iterations+1, nz,nx,ny))
            write(file, "C21", Array{Float64}(undef, max_iterations+1, nz,nx,ny))
            write(file, "C31", Array{Float64}(undef, max_iterations+1, nz,nx,ny))
            write(file, "C32", Array{Float64}(undef, max_iterations+1, nz,nx,ny))
        end
    end
end

"""
    create_output_file(output_path::String, nλ::Int64, atmosphere_size::Tuple)

Initialise all output variables for the test mode.
"""
function create_output_file(output_path::String, nλ::Int64, atmosphere_size::Tuple)

    nz, nx, ny = atmosphere_size

    h5open(output_path, "w") do file
        write(file, "J", Array{Int32,4}(undef, nλ, nz, nx,ny))
        write(file, "total_destroyed", Array{Int64,1}(undef, nλ))
        write(file, "total_scatterings", Array{Int64,1}(undef, nλ))
        write(file, "time", Array{Float64,1}(undef, nλ))

        write(file, "packets", Array{Int32,4}(undef, nλ, nz, nx, ny))
        write(file, "boundary", Array{Int32,3}(undef, nλ, nx, ny))
        write(file, "intensity_per_packet", Array{Float64,1}(undef, nλ))
    end
end

"""
    cut_output_file(output_path::String, final_iteration::Int64, write_rates::Bool)

Cut output data at a given iteration.
"""
function cut_output_file(output_path::String, final_iteration::Int64, write_rates::Bool)
    h5open(output_path, "r+") do file
        # Slice
        J_new = read(file, "J")[1:final_iteration,:,:,:,:]
        total_destroyed_new = read(file, "total_destroyed")[1:final_iteration,:]
        total_scatterings_new = read(file, "total_scatterings")[1:final_iteration,:]
        time_new = read(file, "time")[1:final_iteration,:]
        packets_new = read(file, "packets")[1:final_iteration,:,:,:,:]
        boundary_new = read(file, "boundary")[1:final_iteration,:,:,:]
        intensity_per_packet_new = read(file, "intensity_per_packet")[1:final_iteration,:]
        populations_new = read(file, "populations")[1:final_iteration+1,:,:,:,:]

        # Delete
        delete_object(file, "J")
        delete_object(file, "total_destroyed")
        delete_object(file, "total_scatterings")
        delete_object(file, "time")
        delete_object(file, "packets")
        delete_object(file, "boundary")
        delete_object(file, "intensity_per_packet")
        delete_object(file, "populations")

        # Write
        write(file, "J", J_new)
        write(file, "total_destroyed", total_destroyed_new)
        write(file, "total_scatterings", total_scatterings_new)
        write(file, "time", time_new)
        write(file, "packets", packets_new)
        write(file, "boundary", boundary_new)
        write(file, "intensity_per_packet", intensity_per_packet_new)
        write(file, "populations", populations_new)

        if write_rates
            R12_new = read(file, "R12")[1:final_iteration,:,:,:]
            R13_new = read(file, "R13")[1:final_iteration,:,:,:]
            R23_new = read(file, "R23")[1:final_iteration,:,:,:]
            R21_new = read(file, "R21")[1:final_iteration,:,:,:]
            R31_new = read(file, "R31")[1:final_iteration,:,:,:]
            R32_new = read(file, "R32")[1:final_iteration,:,:,:]
            C12_new = read(file, "C12")[1:final_iteration,:,:,:]
            C13_new = read(file, "C13")[1:final_iteration,:,:,:]
            C23_new = read(file, "C23")[1:final_iteration,:,:,:]
            C21_new = read(file, "C21")[1:final_iteration,:,:,:]
            C31_new = read(file, "C31")[1:final_iteration,:,:,:]
            C32_new = read(file, "C32")[1:final_iteration,:,:,:]

            delete_object(file, "R12")
            delete_object(file, "R13")
            delete_object(file, "R23")
            delete_object(file, "R21")
            delete_object(file, "R31")
            delete_object(file, "R32")
            delete_object(file, "C12")
            delete_object(file, "C13")
            delete_object(file, "C23")
            delete_object(file, "C21")
            delete_object(file, "C31")
            delete_object(file, "C32")

            write(file, "R12", R12_new)
            write(file, "R13", R13_new)
            write(file, "R23", R23_new)
            write(file, "R21", R21_new)
            write(file, "R31", R31_new)
            write(file, "R32", R32_new)
            write(file, "C12", C12_new)
            write(file, "C13", C13_new)
            write(file, "C23", C23_new)
            write(file, "C21", C21_new)
            write(file, "C31", C31_new)
            write(file, "C32", C32_new)
        end
    end
end

"""
    how_much_data(max_iterations::Int64, nλ::Int64, atmosphere_size::Tuple, write_rates::Bool)

Returns the maximum amount of GBs written to file if the
simulation runs for max_iterations.
"""
function how_much_data(nλ::Int64, atmosphere_size::Tuple, max_iterations::Int64, write_rates::Bool)

    nz, nx, ny = atmosphere_size
    boxes = nz*nx*ny
    slice = nx*ny

    λ_data = 8*nλ + 8*2

    # Iteration data
    J_data   = 4boxes*nλ
    sim_data = 4nλ + 2 * 8nλ
    rad_data = 4boxes*nλ + 4slice*nλ + 8nλ
    pop_data = 8boxes*3

    # Rates
    rate_data = 8*12boxes

    max_data = ( λ_data +
               ( J_data + sim_data + rad_data) * max_iterations +
                                      pop_data * (max_iterations + 1) ) / 1e9

    if write_rates
        max_data += rate_data * (max_iterations+1)/1e9
    end

    return max_data
end

"""
    how_much_data(max_iterations::Int64, nλ::Int64, atmosphere_size::Tuple)

Returns the maximum amount of GBs written
to file for the background mode.
"""
function how_much_data(nλ::Int64, atmosphere_size::Tuple)

    nz, nx, ny = atmosphere_size
    boxes = nz*nx*ny
    slice = nx*ny

    λ_data = 8*nλ + 8*2

    # Iteration data
    J_data   = 4boxes*nλ
    sim_data = 4nλ + 2 * 8nλ
    rad_data = 4boxes*nλ + 4slice*nλ + 8nλ

    max_data = ( λ_data +  J_data + sim_data + rad_data ) / 1e9

    return max_data
end
