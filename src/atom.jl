include("io.jl")

struct Atom
    line::AtomicLine
    χl::Unitful.Energy
    χu::Unitful.Energy
    χ∞::Unitful.Energy
    gl::Int64
    gu::Int64
    g∞::Int64
    U0::Float64
    U1::Float64
    Z::Int64

    λ::Array{<:Unitful.Length, 1}                    # (nλ)
    nλ_bb::Int64
    nλ_bf::Int64

    doppler_width::Array{<:Unitful.Length, 3}      # (nz, nx, ny)
    damping_constant::Array{<:PerArea,3}             # (nz, nx, ny)
    density::Array{<:NumberDensity}
end

"""
Collection of 2 level atoms.
"""
function collect_atom_data(atmosphere::Atmosphere)

    # ==================================================================
    # READ ATOM FILE
    # ==================================================================

    atom = h5open(get_atom_path(), "r")
    χl = read(atom, "chi_l")u"cm^-1"
    χu = read(atom, "chi_u")u"cm^-1"
    χ∞ = read(atom, "chi_inf")u"cm^-1"
    gl = read(atom, "gl")
    gu = read(atom, "gu")
    g∞ = read(atom, "ginf")
    f_value = read(atom, "f_value")
    atom_weight = read(atom, "atom_weight")u"kg"
    Z = read(atom, "Z")
    U0 = read(atom, "U_neutral")
    U1 = read(atom, "U_ionised")
    density = read(atom, "density")u"m^-3"
    close(atom)

    line = AtomicLine(χu, χl, χ∞, gu, gl, f_value, atom_weight, Z)

    χl = line.χi
    χu = line.χj
    χ∞ = line.χ∞

    # ==================================================================
    # SAMPLE ATOM TRANSITION WAVELENGTHS
    # ==================================================================
    nλ_bb, nλ_bf = get_nλ()
    nλ_bb += 1-nλ_bb%2  # Make sure odd number
    λ = sample_λ(nλ_bb, nλ_bf, χl, χu, χ∞)

    # ==================================================================
    # CALCULATE DAMPING AND DOPPLER WIDTH
    # ==================================================================
    temperature = atmosphere.temperature
    electron_density = atmosphere.electron_density
    neutral_hydrogen_density = atmosphere.hydrogen_populations[:,:,:,1] .+
                               atmosphere.hydrogen_populations[:,:,:,2]

    #unsold_const = γ_unsold_const(line)
    unsold_const = const_unsold(line)
    quad_stark_const = const_quadratic_stark(line)

    γ = γ_unsold.(unsold_const, temperature, neutral_hydrogen_density)
    γ .+= line.Aji
    γ .+= γ_linear_stark.(electron_density, 3, 2)
    γ .+= γ_quadratic_stark.(electron_density, temperature, stark_constant=quad_stark_const)

    ΔλD = doppler_width.(line.λ0, atom_weight, temperature)
    damping_const = damping_constant.(γ, ΔλD)

    # ===========================================================
    # NO NEGAITVE OR INFINITE VALUES
    # ===========================================================
    @test all( Inf > ustrip(line.Aji) >= 0.0 )
    @test all( Inf > ustrip(line.Bji) >= 0.0 )
    @test all( Inf > ustrip(line.Bij) >= 0.0 )
    @test all( Inf .> ustrip.(damping_const) .>= 0.0 )
    @test all( Inf .> ustrip.(ΔλD) .>= 0.0 )
    @test all( Inf .> ustrip.(λ) .>= 0.0 )
    @test Inf > ustrip(line.λ0) >= 0.0
    @test Inf > nλ_bb >= 0
    @test Inf > nλ_bf >= 0
    @test ustrip(χl) == 0.0
    @test ustrip(χu) > 0.0
    @test χ∞ > χu

    return line,
           χl, χu, χ∞,
           gu, gl, g∞,
           U0, U1, Z,
           λ, nλ_bb, nλ_bf,
           ΔλD, damping_const,
           density
end

function sample_λ(nλ_bb, nλ_bf, χl, χu, χ∞)

    atom_file = h5open(get_atom_path(), "r")
    λ_bf_l_min = read(atom_file, "bfl_min")u"nm"
    λ_bf_u_min = read(atom_file, "bfu_min")u"nm"
    qwing = read(atom_file, "qwing")
    qcore = read(atom_file, "qcore")
    close(atom_file)

    λ_bf_l_max  = transition_λ(χl, χ∞)
    λ_bf_u_max  = transition_λ(χu, χ∞)
    λ_bb_center = transition_λ(χl, χu)

    # Make sure odd # of bb wavelengths
    if nλ_bb > 0 && nλ_bb%2 == 0
        nλ_bb += 1
    end

    # Initialise wavelength array
    nλ = nλ_bf*2 + nλ_bb
    λ = Array{Float64,1}(undef, nλ)u"nm"

    # =================================================
    # Bound-free transitions
    # Linear spacing
    # =================================================
    if nλ_bf > 1
        Δλ_bf_l = (λ_bf_l_max - λ_bf_l_min)/nλ_bf
        Δλ_bf_u = (λ_bf_u_max - λ_bf_u_min)/nλ_bf

        λ[1] = λ_bf_l_min
        λ[nλ_bf+1] = λ_bf_u_min

        for l=2:nλ_bf
            λ[l] = λ[l-1] + Δλ_bf_l
            λ[l+nλ_bf] = λ[l+nλ_bf - 1] + Δλ_bf_u
        end
    elseif nλ_bf == 1
        λ[1] = λ_bf_l_max
        λ[2] = λ_bf_u_max
    end

    # =================================================
    # Bound-bound transition
    # Follows github.com/ITA-Solar/rh/blob/master/getlambda.c
    # =================================================
    if nλ_bb > 1
        vmicro_char = 2.5u"km/s"

        n = nλ_bb#/2 # Questionable
        β = qwing/(2*qcore)
        y = β + sqrt(β*β + (β - 1.0)*n + 2.0 - 3.0*β)
        b = 2.0*log(y) / (n - 1)
        a = qwing / (n - 2.0 + y*y)

        center = nλ_bf*2 + (nλ_bb÷2) + 1
        λ[center] = λ_bb_center
        q_to_λ = λ[center] * vmicro_char / c_0

        for l=1:(nλ_bb÷2)
            Δλ = a*(l + (exp(b*l) - 1.0)) * q_to_λ
            λ[center-l] = λ[center] - Δλ
            λ[center+l] = λ[center] + Δλ
        end
    elseif nλ_bb == 1
        λ[2nλ_bf + 1] = λ_bb_center
    end

    return λ
end

function damping_constant(γ::Unitful.Frequency,
                          ΔλD::Unitful.Length)
    (γ / (4 * π * c_0 * ΔλD))
end

function write_to_file(λ::Array{<:Unitful.Length,1})
    h5open("../out/output.h5", "r+") do file
        write(file, "wavelength", ustrip(λ))
    end
end

function write_to_file(atom::Atom, output_path::String)
    h5open(output_path, "r+") do file
        write(file, "wavelength", ustrip(atom.λ))
        write(file, "nlambda_bb", atom.nλ_bb)
        write(file, "nlambda_bf", atom.nλ_bf)
    end
end
