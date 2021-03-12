"""
Collection of all imports and functions to read simulation input.
"""

using DelimitedFiles
using BenchmarkTools
using ProgressMeter
using Transparency
using StaticArrays
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

function test_mode()
    input_file = open(f->read(f, String), "../run/keywords.input")
    i = findfirst("test_mode", input_file)[end] + 1
    file = input_file[i:end]
    i = findfirst("=", file)[end]
    j = findfirst("\n", file[i+1:end])[end] + i
    tm = parse(Bool, file[i+1:j-1])
    return tm
end

function get_output_path()

    if test_mode()
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

function get_max_iterations()
    input_file = open(f->read(f, String), "../run/keywords.input")
    i = findfirst("max_iterations", input_file)[end] + 1
    file = input_file[i:end]
    i = findfirst("=", file)[end] + 1
    j = findfirst("\n", file)[end] - 1
    max_iterations = parse(Int64, file[i:j])
    return max_iterations
end

function get_error(output_path, n)
    error = nothing
    h5open(output_path, "r") do file
        error = read(file, "error")
    end
    return error[n] # DELETE
end

# =============================================================================
# RADIATION
# =============================================================================

function get_target_packets()
    input_file = open(f->read(f, String), "../run/keywords.input")
    i = findfirst("target_packets", input_file)[end] + 1
    file = input_file[i:end]
    i = findfirst("=", file)[end] + 1
    j = findfirst("\n", file)[end] - 1
    target_packets = parse(Float64, file[i:j])
    return target_packets # move
end

function get_max_scatterings()
    input_file = open(f->read(f, String), "../run/keywords.input")
    i = findfirst("max_scatterings", input_file)[end] + 1
    file = input_file[i:end]
    i = findfirst("=", file)[end] + 1
    j = findfirst("\n", file)[end] - 1
    max_scatterings = parse(Float64, file[i:j])
    return max_scatterings # move
end

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

function get_background_λ()
    input_file = open(f->read(f, String), "../run/keywords.input")
    i = findfirst("background_wavelength", input_file)[end] + 1
    file = input_file[i:end]
    i = findfirst("=", file)[end] + 1
    j = findfirst("\n", file)[end] - 1

    λ = parse(Float64, file[i:j])u"nm"

    return λ
end

function get_Jλ(output_path, iteration, intensity_per_packet)
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

function get_atom_path()
    input_file = open(f->read(f, String), "../run/keywords.input")
    i = findfirst("atom_path", input_file)[end] + 1
    file = input_file[i:end]
    i = findfirst("\"", file)[end]
    j = findfirst("\"", file[i+1:end])[end] + i
    atmosphere_path = string(file[i+1:j-1])
    return atmosphere_path
end

function get_population_distribution()
    input_file = open(f->read(f, String), "../run/keywords.input")
    i = findfirst("population_distribution", input_file)[end] + 1
    file = input_file[i:end]
    i = findfirst("\"", file)[end]
    j = findfirst("\"", file[i+1:end])[end] + i
    distribution = string(file[i+1:j-1])
    return distribution
end

function get_initial_populations_path()
    input_file = open(f->read(f, String), "../run/keywords.input")
    i = findfirst("initial_populations_path", input_file)[end] + 1
    file = input_file[i:end]
    i = findfirst("\"", file)[end]
    j = findfirst("\"", file[i+1:end])[end] + i
    initial_populations_path = string(file[i+1:j-1])
    return initial_populations_path  # DELETE
end

# =============================================================================
# ATMOSPHERE
# =============================================================================

function get_atmosphere_path()
    input_file = open(f->read(f, String), "../run/keywords.input")
    i = findfirst("atmosphere_path", input_file)[end] + 1
    file = input_file[i:end]
    i = findfirst("\"", file)[end]
    j = findfirst("\"", file[i+1:end])[end] + i
    atmosphere_path = string(file[i+1:j-1])
    return atmosphere_path
end

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

function create_output_file(output_path, max_iterations, nλ, atmosphere_size)

    nz, nx, ny = atmosphere_size

    h5open(output_path, "w") do file
        #J = create_dataset(file, "J", datatype(Int32), dataspace(nλ,nz,nx,ny), chunk=(1,nz,nx,ny))
        write(file, "J", Array{Int32,5}(undef, max_iterations, nλ, nz, nx,ny))
        write(file, "total_destroyed", Array{Int32,2}(undef, max_iterations, nλ))
        write(file, "total_scatterings", Array{Int32,2}(undef,max_iterations, nλ))
        write(file, "time", Array{Float64,2}(undef,max_iterations, nλ))

        write(file, "packets", Array{Int32,5}(undef, max_iterations, nλ, nz, nx, ny))
        write(file, "boundary", Array{Int32,4}(undef,max_iterations, nλ, nx, ny))
        write(file, "intensity_per_packet", Array{Float64,2}(undef,max_iterations, nλ))

        write(file, "populations", Array{Float64,5}(undef, max_iterations+1, nz, nx, ny, 3))
        write(file, "error", Array{Float64,1}(undef, max_iterations))
    end
end


function create_output_file(output_path, nλ, atmosphere_size)

    nz, nx, ny = atmosphere_size

    h5open(output_path, "w") do file
        #J = create_dataset(file, "J", datatype(Int32), dataspace(nλ,nz,nx,ny), chunk=(1,nz,nx,ny))
        write(file, "J", Array{Int32,5}(undef, nλ, nz, nx,ny))
        write(file, "total_destroyed", Array{Int32,2}(undef, nλ))
        write(file, "total_scatterings", Array{Int64,2}(undef, nλ))
        write(file, "time", Array{Float64,2}(undef, nλ))

        write(file, "packets", Array{Int32,5}(undef, nλ, nz, nx, ny))
        write(file, "boundary", Array{Int32,4}(undef, nλ, nx, ny))
        write(file, "intensity_per_packet", Array{Float64,2}(undef, nλ))
    end
end


function cut_output_file(output_path, final_iteration)
    h5open(output_path, "w") do file
        #J = create_dataset(file, "J", datatype(Int32), dataspace(nλ,nz,nx,ny), chunk=(1,nz,nx,ny))
        J_new = read(file, "J")[1:final_iteration,:,:,:,:]
        delete_object(file, "J")
        write(file, "J", J_new)
        # repeat....
        write(file, "total_destroyed", Array{Int32,2}(undef, max_iterations, nλ))
        write(file, "total_scatterings", Array{Int32,2}(undef,max_iterations, nλ))
        write(file, "time", Array{Float64,2}(undef,max_iterations, nλ))

        write(file, "packets", Array{Int32,5}(undef, max_iterations, nλ, nz, nx, ny))
        write(file, "boundary", Array{Int32,4}(undef,max_iterations, nλ, nx, ny))
        write(file, "intensity_per_packet", Array{Float64,2}(undef,max_iterations, nλ))

        write(file, "populations", Array{Float64,5}(undef, max_iterations, nz, nx, ny, 3))
    end
end

function how_much_data(max_iterations, nλ, atmosphere_size)
    λ_data = 8*nλ + 8*2

    nz, nx, ny = atmosphere_size
    boxes = nz*nx*ny
    slice = nx*ny

    # Iteration data
    J_data = 4boxes*nλ
    sim_data = 2 * 4nλ + 8nλ
    rad_data = 4boxes*nλ + 4slice*nλ + 8nλ
    pop_data = 8boxes*3

    min_data = ( λ_data +  J_data + sim_data + rad_data + pop_data ) /1e9
    max_data = ( λ_data + (J_data + sim_data + rad_data + pop_data) * max_iterations ) / 1e9

    return max_data
end
