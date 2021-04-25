include("check_parameters.jl")
include("plot_parameters.jl")

function run_tests(check=true, plot=false)
    target_packets = get_target_packets()
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
    radiationBackground = Radiation(radiation_parameters...)

    # =============================================================================
    # ATOM
    # =============================================================================
    atom_parameters = collect_atom_data(atmosphere)
    atom = Atom(atom_parameters...)

    n_levels = atom.n_levels
    n_lines = atom.n_lines
    n_transitions = n_levels + n_lines

    # =============================================================================
    # LOAD LINES
    # =============================================================================
    lines = []
    for l=1:n_levels-1
        for u=2:n_levels
            line =  Line(collect_line_data(atmosphere, atom, u, l)...)
            append!(lines, [line])
        end
    end

    # =============================================================================
    # INITIAL POPULATIONS
    # =============================================================================
    populations_LTE = LTE_populations(atom, atmosphere.temperature,
                                            atmosphere.electron_density)
    populations_ZR = zero_radiation_populations(atmosphere, atom, lines)

    # =============================================================================
    # INITIAL TRANSITION RATES
    # =============================================================================
    Bλ = blackbody_lambda(atom.λ, atmosphere.temperature)
    rate_parameters = calculate_transition_rates(atmosphere, atom, lines, Bλ)
    rates = TransitionRates(rate_parameters...)

    # =============================================================================
    # RADIATION
    # =============================================================================
    lineRadiations = []

    for l=1:n_levels-1
        for u=2:n_levels
            line_number = sum((n_levels-l+1):(n_levels-1)) + (u - l)
            lineRadiation = LineRadiation(collect_line_radiation_data(lines[line_number], rates, populations_LTE)...)
            append!(lineRadiations, [lineRadiation])
        end
    end

    radiation_parameters = collect_radiation(atmosphere, atom, rates, lines, lineRadiations,
                                             populations_LTE, cut_off, target_packets)

    radiation = Radiation(radiation_parameters...)


    if check == true
        print("--Run all tests...........")
        check_atmosphere(atmosphere)
        check_radiation(radiationBackground, atmosphere_size, 1)
        check_populations(populations_LTE, atmosphere_size)
        check_populations(populations_ZR, atmosphere_size)
        check_atom(atom, atmosphere_size)
        check_rates(rates, atmosphere_size, n_levels)
        for l=1:n_lines
            check_line(lines[l], atmosphere_size)
        end
        check_radiation(radiation, atmosphere_size, atom.nλ)
        println("All tests passed.")
    end

    if plot == true
        print("--Plot all parameters.....")
        plot_atmosphere(atmosphere)
        plot_radiationBackground(radiationBackground, z)
        plot_populations(populations_LTE, populations_ZR, z)
        plot_rates(rates, z)
        plot_radiation(radiation, atom, lines, lineRadiations, z)
        println("All parameters plotted.\n")
    end

end

run_tests(false, true)
