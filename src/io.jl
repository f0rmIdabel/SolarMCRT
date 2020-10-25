import PhysicalConstants.CODATA2018: c_0, h, k_B
using DelimitedFiles
using ProgressMeter
using Transparency
using StaticArrays
using Unitful
using Future
using Random
using Printf
using HDF5

@derived_dimension NumberDensity Unitful.𝐋^-3
@derived_dimension PerLength Unitful.𝐋^-1

#using Future # for randjump in rng when using threads
#using StaticArrays

function output(S::Array{Int64,3},
                J::Array{Int64,3},
                surface_intensity::Array{Int64,4},
                total_destroyed::Int64,
                total_scatterings::Int64)
    out = h5open("../out/output.hdf5", "w")
    write(out, "S", S)
    write(out, "J", J)
    write(out, "surface_intensity", surface_intensity)
    write(out, "total_packets", sum(S))
    write(out, "total_destroyed", total_destroyed)
    write(out, "total_escaped", sum(surface_intensity))
    write(out, "total_scatterings", total_scatterings)
    close(out)
end

function get_λ()
    λ = readdlm("../run/wavelengths.input", ' ')[1,:]u"nm"
    return λ
end

function get_τ_max()
    input_file = open(f->read(f, String), "../run/keywords.input")
    i = findfirst("tau_max", input_file)[end] + 1
    file = input_file[i:end]
    i = findfirst("=", file)[end] + 1
    j = findfirst("\n", file)[end] - 1
    τ_max = parse(Float64, file[i:j])
    return τ_max
end

function get_atmosphere_path()
    input_file = open(f->read(f, String), "../run/keywords.input")
    i = findfirst("atmosphere_path", input_file)[end] + 1
    file = input_file[i:end]
    i = findfirst("\"", file)[end]
    j = findfirst("\"", file[i+1:end])[end] + i
    atmosphere_path = string(file[i+1:j-1])
    return atmosphere_path
end

function get_target_packets()
    input_file = open(f->read(f, String), "../run/keywords.input")
    i = findfirst("target_packets", input_file)[end] + 1
    file = input_file[i:end]
    i = findfirst("=", file)[end] + 1
    j = findfirst("\n", file)[end] - 1
    target_packets = parse(Float64, file[i:j])
    return target_packets
end

function get_max_scatterings()
    input_file = open(f->read(f, String), "../run/keywords.input")
    i = findfirst("max_scatterings", input_file)[end] + 1
    file = input_file[i:end]
    i = findfirst("=", file)[end] + 1
    j = findfirst("\n", file)[end] - 1
    max_scatterings = parse(Float64, file[i:j])
    return max_scatterings
end

function get_escape_bins()
    input_file = open(f->read(f, String), "../run/keywords.input")
    i = findfirst("escape_bins", input_file)[end] + 1
    file = input_file[i:end]
    i = findfirst("[", file)[end] + 1
    j = findfirst(",", file)[end] - 1
    ϕ_bin = parse(UInt16, file[i:j])
    i = j + 2
    j = findfirst("]", file)[end] - 1
    θ_bin = parse(UInt16, file[i:j])
    escape_bins = [ϕ_bin, θ_bin] #SA[ϕ_bin, θ_bin]
    return escape_bins
end
