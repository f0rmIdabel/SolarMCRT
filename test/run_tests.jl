include("check_parameters.jl")
include("plot_parameters.jl")

function run_tests(check=true, plot=false)
    target_packets = get_target_packets()
    target_packets_bf = get_target_packets_bf()
    target_packets_bb = get_target_packets_bb()
    cut_off = get_cut_off()

    # =============================================================================
    # ATMOSPHERE
    # =============================================================================
    atmosphere_parameters = collect_atmosphere()
    atmosphere = Atmosphere(atmosphere_parameters...)
    atmosphere_size = size(atmosphere.temperature)
    z = atmosphere.z[1:end-1]

    # =============================================================================
    # BACKGROUND RADIATION
    # =============================================================================
    λ = [get_background_λ()]
    radiation_parameters = collect_background_radiation(atmosphere, λ, cut_off, target_packets)
    radiationBackground = RadiationContinuum(radiation_parameters...)

    # =============================================================================
    # ATOM
    # =============================================================================
    atom_parameters = collect_atom_data(atmosphere)
    atom = Atom(atom_parameters...)

    n_levels = atom.n_levels
    n_lines = atom.n_lines
    n_transitions = n_levels + n_lines

    # =============================================================================
    # INITIAL POPULATIONS
    # =============================================================================
    populations_LTE = LTE_populations(atom, atmosphere.temperature,
                                            atmosphere.electron_density)
    populations_ZR = zero_radiation_populations(atmosphere, atom)

    plot_populations(populations_LTE, populations_ZR, z)

    # =============================================================================
    # INITIAL TRANSITION RATES
    # =============================================================================
    Bλ = []
    for t=1:n_transitions
        append!(Bλ, [blackbody_lambda(atom.λ[t], atmosphere.temperature)])
    end
    rate_parameters = calculate_transition_rates(atmosphere, atom, Bλ)
    rates = TransitionRates(rate_parameters...)

    # =============================================================================
    # RADIATION
    # =============================================================================

    for level=1:n_levels
        println(level)
        λ = atom.λ[level]
        radiation_parameters = collect_bf_radiation(atmosphere, atom, level, rates, populations_ZR, cut_off, target_packets_bf[level])
        radiation = RadiationContinuum(radiation_parameters...)

        if check == true
            check_radiationContinuum(radiation, atmosphere_size, length(λ))
        end

        if plot == true
            plot_radiation(radiation, λ, level, z)
        end
    end

    line_number = 0
    for l=1:n_levels-1
        for u=l+1:n_levels
            line_number += 1
            λ = atom.λ[n_levels + line_number]
            line_parameters = collect_line_data(atmosphere, atom, u, l)
            line = Line(line_parameters...)
            radiation_parameters = collect_bb_radiation(atmosphere, λ, line, rates, populations_LTE, cut_off, target_packets_bb[line_number])
            radiation = RadiationLine(radiation_parameters...)

            if check == true
                check_line(line, atmosphere_size)
                check_radiationLine(radiation, λ, line, atmosphere_size)
            end

            if plot == true
                plot_radiation(radiation, λ, line, z)
            end
        end
    end

    if check == true
        print("--Run all tests...........")
        check_atmosphere(atmosphere)
        check_radiationContinuum(radiationBackground, atmosphere_size, 1)
        check_populations(populations_LTE, atmosphere_size)
        check_populations(populations_ZR, atmosphere_size)
        check_atom(atom, atmosphere_size)
        check_rates(rates, atmosphere_size, n_levels)
        println("All tests passed.")
    end

    if plot == true
        print("--Plot all parameters.....")
        plot_atmosphere(atmosphere)
        plot_radiationBackground(radiationBackground, z)
        plot_populations(populations_LTE, populations_ZR, z)
        plot_rates(rates, z)
        println("All parameters plotted.\n")
    end

end

run_tests(true, true)
