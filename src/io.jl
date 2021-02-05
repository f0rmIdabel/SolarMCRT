"""
Collection of all imports and funcitons to read simulation input.
"""

using DelimitedFiles
using BenchmarkTools
using ProgressMeter
using Transparency
using StaticArrays
using Unitful
using Random
using Printf
using HDF5

import PhysicalConstants.CODATA2018: c_0, h, k_B, m_u, m_e
@derived_dimension NumberDensity Unitful.𝐋^-3
@derived_dimension PerLength Unitful.𝐋^-1
@derived_dimension UnitsIntensity_λ Unitful.𝐋^-1 * Unitful.𝐌 * Unitful.𝐓^-3

function test_mode()
    input_file = open(f->read(f, String), "/mn/stornext/u3/idarhan/MScProject/SolarMCRT/run/keywords.input")
    i = findfirst("test_mode", input_file)[end] + 1
    file = input_file[i:end]
    i = findfirst("\"", file)[end]
    j = findfirst("\"", file[i+1:end])[end] + i
    tm = parse(Bool, file[i+1:j-1])
    return tm
end

function get_nλ()
    input_file = open(f->read(f, String), "/mn/stornext/u3/idarhan/MScProject/SolarMCRT/run/keywords.input")
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

function get_background_λ()
    input_file = open(f->read(f, String), "/mn/stornext/u3/idarhan/MScProject/SolarMCRT/run/keywords.input")
    i = findfirst("background_wavelength", input_file)[end] + 1
    file = input_file[i:end]
    i = findfirst("=", file)[end] + 1
    j = findfirst("\n", file)[end] - 1

    λ = parse(Float64, file[i:j])u"nm"

    return λ
end

function get_cut_off()
    input_file = open(f->read(f, String), "/mn/stornext/u3/idarhan/MScProject/SolarMCRT/run/keywords.input")
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

    return cut_off
end

function get_atmosphere_path()
    input_file = open(f->read(f, String), "/mn/stornext/u3/idarhan/MScProject/SolarMCRT/run/keywords.input")
    i = findfirst("atmosphere_path", input_file)[end] + 1
    file = input_file[i:end]
    i = findfirst("\"", file)[end]
    j = findfirst("\"", file[i+1:end])[end] + i
    atmosphere_path = string(file[i+1:j-1])
    return atmosphere_path
end

function get_atom_path()
    input_file = open(f->read(f, String), "/mn/stornext/u3/idarhan/MScProject/SolarMCRT/run/keywords.input")
    i = findfirst("atom_path", input_file)[end] + 1
    file = input_file[i:end]
    i = findfirst("\"", file)[end]
    j = findfirst("\"", file[i+1:end])[end] + i
    atmosphere_path = string(file[i+1:j-1])
    return atmosphere_path
end

function get_initial_populations_path()
    input_file = open(f->read(f, String), "/mn/stornext/u3/idarhan/MScProject/SolarMCRT/run/keywords.input")
    i = findfirst("initial_populations_path", input_file)[end] + 1
    file = input_file[i:end]
    i = findfirst("\"", file)[end]
    j = findfirst("\"", file[i+1:end])[end] + i
    initial_populations_path = string(file[i+1:j-1])
    return initial_populations_path
end

function get_target_packets()
    input_file = open(f->read(f, String), "/mn/stornext/u3/idarhan/MScProject/SolarMCRT/run/keywords.input")
    i = findfirst("target_packets", input_file)[end] + 1
    file = input_file[i:end]
    i = findfirst("=", file)[end] + 1
    j = findfirst("\n", file)[end] - 1
    target_packets = parse(Float64, file[i:j])
    return target_packets
end

function get_max_scatterings()
    input_file = open(f->read(f, String), "/mn/stornext/u3/idarhan/MScProject/SolarMCRT/run/keywords.input")
    i = findfirst("max_scatterings", input_file)[end] + 1
    file = input_file[i:end]
    i = findfirst("=", file)[end] + 1
    j = findfirst("\n", file)[end] - 1
    max_scatterings = parse(Float64, file[i:j])
    return max_scatterings
end

function get_step()
    input_file = open(f->read(f, String), "/mn/stornext/u3/idarhan/MScProject/SolarMCRT/run/keywords.input")
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

    input_file = open(f->read(f, String), "/mn/stornext/u3/idarhan/MScProject/SolarMCRT/run/keywords.input")
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
    input_file = open(f->read(f, String), "/mn/stornext/u3/idarhan/MScProject/SolarMCRT/run/keywords.input")
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

function get_λ_min()
    input_file = open(f->read(f, String), "/mn/stornext/u3/idarhan/MScProject/SolarMCRT/run/keywords.input")
    i = findfirst("λ_min_bf", input_file)[end] + 1
    file = input_file[i:end]
    i = findfirst("[", file)[end] + 1
    j = findfirst(",", file)[end] - 1

    λ_min1 = parse(Float64, file[i:j])u"nm"

    i = j + 2
    file = file[i:end]
    j = findfirst("]", file)[end] - 1

    λ_min2 = parse(Float64, file[1:j])u"nm"

    return [λ_min1, λ_min2]
end
