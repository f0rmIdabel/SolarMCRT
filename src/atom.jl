include("io.jl")

struct Atom
    density::Array{<:NumberDensity,3}                  # (nz, nx, ny)
    n_levels::Int64
    n_lines::Int64
    χ::Array{<:Unitful.Energy,1}
    g::Array{Int64,1}
    Z::Int64
    λ::Array{Any, 1}
    nλ::Int64
end

struct Line
    u::Int
    l::Int
    lineData::AtomicLine
    doppler_width::Array{<:Unitful.Length, 3}        # (n_lines, nz, nx, ny)
    damping_constant::Array{<:PerArea,3}             # (n_lines, nx, ny)
end

"""
    collect_atom_data(atmosphere::Atmosphere)

Reads a two-level atom file, samples wavelengths from its
bound-bound and bound-free transitions and returns all
relevant data for the simulation.
"""
function collect_atom_data(atmosphere::Atmosphere)

    # ==================================================================
    # READ ATOM FILE
    # ==================================================================
    atom = h5open(get_atom_path(), "r")
    stage = read(atom, "stage")
    element = read(atom, "element")
    Z = read(atom, "Z")
    density = read(atom, "density")u"m^-3"
    χ = read(atom, "chi")u"J"
    χ_ion = read(atom, "chi_ion")u"J"
    g = read(atom, "g")
    close(atom)

    # ==================================================================
    # RELEVANT ATMOSPHERE DATA
    # ==================================================================
    temperature = atmosphere.temperature
    electron_density = atmosphere.electron_density
    neutral_hydrogen_density = sum(atmosphere.hydrogen_populations[:,:,:,1:end-1], dims=4)[:,:,:,1]

    nz, nx, ny = size(temperature)

    # ==================================================================
    # ATOMIC STAGE
    # ==================================================================
    n_levels = get_nlevels()
    n_lines = Int(n_levels*(n_levels-1)/2)

    χ = χ[1:n_levels]
    append!(χ, χ_ion)

    g = g[1:n_levels]
    append!(g, 1)

    # ==================================================================
    # SAMPLE ATOM TRANSITION WAVELENGTHS
    # ==================================================================
    nλ_bf = get_nλ_bf()
    λ_min = get_λ_min() .* u"nm"

    nλ_bb = get_nλ_bb()
    qcore = get_qcore()
    qwing = get_qwing()

    λ = []
    nλ = 0

    # Sample bound free transitions
    for level=1:n_levels
        λ_bf = sample_λ_boundfree(nλ_bf[level], λ_min[level], χ[level], χ[end])
        append!(λ, [λ_bf])
        nλ += length(λ_bf)
    end

    # Sample bound-bound transitions
    line = 0
    for l=1:n_levels-1
        for u=(l+1):n_levels
            line += 1
            λ_bb = sample_λ_line(nλ_bb[line], χ[l], χ[u], qwing[line], qcore[line])
            append!(λ, [λ_bb])
            nλ +=length(λ_bb)
        end
    end

    # ===========================================================
    # NO NEGAITVE OR INFINITE VALUES
    # ===========================================================
    @test all( Inf .> ustrip(density) .>= 0.0 )
    @test all( Inf .> ustrip.(χ) .>= 0.0 )
    @test all( Inf .> g .> 0)

    for l=1:(n_levels+n_lines)
        @test all(Inf .> ustrip.(λ[l]) .>= 0.0 )
    end

    return density, n_levels, n_lines, χ, g, Z, λ, nλ
end

"""
    collect_line_data(atmosphere::Atmosphere)

Reads a two-level atom file, samples wavelengths from its
bound-bound and bound-free transitions and returns all
relevant data for the simulation.
"""
function collect_line_data(atmosphere::Atmosphere, atom::Atom, u::Int, l::Int)

    # ==================================================================
    # READ ATOM FILE
    # ==================================================================
    atom_file = h5open(get_atom_path(), "r")
    atom_weight = read(atom_file, "atom_weight")u"kg"
    f_value = read(atom_file, "f_value")
    Z = read(atom_file, "Z")
    close(atom_file)

    # ==================================================================
    # RELEVANT ATMOSPHERE DATA
    # ==================================================================
    temperature = atmosphere.temperature
    electron_density = atmosphere.electron_density
    neutral_hydrogen_density = sum(atmosphere.hydrogen_populations[:,:,:,1:end-1], dims=4)[:,:,:,1]

    nz, nx, ny = size(temperature)

    # ==================================================================
    # RELEVANT LINE DATA
    # ==================================================================
    χ = atom.χ
    g = atom.g

    # ==================================================================
    # COLLECT LINE DATA
    # ==================================================================

    lineData = AtomicLine(χ[u], χ[l], χ[end], g[u], g[l], f_value, atom_weight, Z)
    unsold_const = const_unsold(lineData)
    quad_stark_const = const_quadratic_stark(lineData)

    γ = γ_unsold.(unsold_const, temperature, neutral_hydrogen_density)
    γ .+= lineData.Aji
    γ .+= γ_linear_stark.(electron_density, u, l)
    γ .+= γ_quadratic_stark.(electron_density, temperature, stark_constant=quad_stark_const)

    ΔλD = doppler_width.(lineData.λ0, atom_weight, temperature)
    damping_const = damping_constant.(γ, ΔλD)

    # ===========================================================
    # NO NEGAITVE OR INFINITE VALUES
    # ===========================================================
    @test all( Inf .> ustrip.(damping_const) .>= 0.0 )
    @test all( Inf .> ustrip.(ΔλD) .>= 0.0 )

    return u, l, lineData, ΔλD, damping_const
end


"""
    sample_λ(nλ_bb::Int64, nλ_bf::Int64,
             χl::Unitful.Energy, χu::Unitful.Energy, χ∞::Unitful.Energy)

Get sampling wavelengths. Bound free wavelengths are
linearly sampled, while th bound-bound follow the
log-sampling from github.com/ITA-Solar/rh.
"""
function sample_λ_line(nλ::Int64,
                       χl::Unitful.Energy,
                       χu::Unitful.Energy,
                       qwing::Float64,
                       qcore::Float64)

    atom_file = h5open(get_atom_path(), "r")

    close(atom_file)

    λ_center = transition_λ(χl, χu)

    # Make sure odd # of bb wavelengths
    if nλ > 0 && nλ%2 == 0
        nλ += 1
    end

    # Either 1 or five or more wavelengths
    if 1 < nλ < 5
        nλ = 5
    end

    # Initialise wavelength array
    λ = Array{Float64,1}(undef, nλ)u"nm"

    # =================================================
    # Bound-bound transition
    # Follows github.com/ITA-Solar/rh/blob/master/getlambda.c
    # =================================================
    if nλ == 1
        λ[1] = λ_center

    elseif nλ >= 5
        vmicro_char = 2.5u"km/s"

        n = nλ/2 # Questionable
        β = qwing/(2*qcore)
        y = β + sqrt(β*β + (β - 1.0)*n + 2.0 - 3.0*β)
        b = 2.0*log(y) / (n - 1)
        a = qwing / (n - 2.0 + y*y)

        center = (nλ÷2) + 1
        λ[center] = λ_center
        q_to_λ = λ[center] * vmicro_char / c_0

        for w=1:(nλ÷2)
            Δλ = a*(w + (exp(b*w) - 1.0)) * q_to_λ
            λ[center-w] = λ[center] - Δλ
            λ[center+w] = λ[center] + Δλ
        end
    end

    return λ
end

"""
    sample_λ(nλ_bb::Int64, nλ_bf::Int64,
             χl::Unitful.Energy, χu::Unitful.Energy, χ∞::Unitful.Energy)

Get sampling wavelengths. Bound free wavelengths are
linearly sampled, while th bound-bound follow the
log-sampling from github.com/ITA-Solar/rh.
"""
function sample_λ_boundfree(nλ::Int64,
                            λ_min::Unitful.Length,
                            χl::Unitful.Energy,
                            χ∞::Unitful.Energy)


    λ_max  = transition_λ(χl, χ∞)


    # Initialise wavelength array
    λ = Array{Float64,1}(undef, nλ)u"nm"

    # =================================================
    # Bound-free transitions
    # Linear spacing
    # =================================================
    if nλ == 1

        λ[1] = λ_max

    elseif nλ > 1
        Δλ = (λ_max - λ_min)/(nλ-1)
        λ[1] = λ_min

        for w=2:nλ
            λ[w] = λ[w-1] + Δλ
        end
    end

    return λ
end

"""
    transition_λ(χ1::Unitful.Energy, χ2::Unitful.Energy)

Get the corresponding wavelength for
the energy difference between two levels.
"""
function transition_λ(χ1::Unitful.Energy, χ2::Unitful.Energy)
    ((h * c_0) / (χ2-χ1)) |> u"nm"
end

"""
    damping_constant(γ::Unitful.Frequency,
                     ΔλD::Unitful.Length)

Get daping constant to be multiplied with λ^2.
"""
function damping_constant(γ::Unitful.Frequency,
                          ΔλD::Unitful.Length)
    (γ / (4 * π * c_0 * ΔλD))
end



# ==================================================================
#  WRITE TO FILE
# ==================================================================


"""
    write_to_file(atom::Atom, output_path::String)

Writes wavelengths to the output file.
"""
function write_to_file(λ::Array{<:Unitful.Length,1}, output_path::String)
    h5open(output_path, "r+") do file
        write(file, "wavelength", ustrip(λ))
    end
end

"""
    write_to_file(atom::Atom, output_path::String)

Writes wavelength and number of bound-bound and
bound-free wavelengths to the output file.
"""
function write_to_file(λ, output_path::String)
    λ_flat = Array{Float64,1}(undef,0)
    for t=1:length(λ)
        append!(λ_flat, ustrip.(λ[t]))
    end
    h5open(output_path, "r+") do file
        write(file, "wavelength", λ_flat)
    end
end
