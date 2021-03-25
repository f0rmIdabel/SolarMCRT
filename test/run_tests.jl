include("check_parameters.jl")
include("plot_parameters.jl")

function run_tests(check=true, plot=false)

    target_packets = get_target_packets()
    cut_off = get_cut_off()

    # =============================================================================
    # ATMOSPHERE
    # =============================================================================
    atmosphere_parameters = collect_atmosphere_data()
    atmosphere = Atmosphere(atmosphere_parameters...)

    # =============================================================================
    # BACKGROUND RADIATION
    # =============================================================================
    λ = get_background_λ()
    radiation_parameters = collect_radiation_data(atmosphere, λ, cut_off, target_packets)
    radiationBackground = RadiationBackground(radiation_parameters...)

    # =============================================================================
    # ATOM
    # =============================================================================
    atom_parameters = collect_atom_data(atmosphere)
    atom = Atom(atom_parameters...)

    # =============================================================================
    # INITIAL POPULATIONS
    # =============================================================================
    populations_LTE = LTE_populations(atom, atmosphere.temperature,
                                            atmosphere.electron_density)
    populations_zero_radiation = zero_radiation_populations(atom, atmosphere.temperature,
                                                                  atmosphere.electron_density)

    # =============================================================================
    # INITIAL TRANSITION RATES
    # =============================================================================
    Bλ = blackbody_lambda(atom.λ, atmosphere.temperature)
    rate_parameters = calculate_transition_rates(atom, atmosphere.temperature,
                                                       atmosphere.electron_density, Bλ)
    rates = TransitionRates(rate_parameters...)

    # =============================================================================
    # RADIATION
    # =============================================================================
    radiation_parameters = collect_radiation_data(atmosphere, atom, rates, populations_LTE, cut_off, target_packets)
    radiation = Radiation(radiation_parameters...)

    if check == true
        print("--Run all tests...........")
        atmosphere_size = size(atmosphere.temperature)
        check_atmosphere(atmosphere)
        check_radiationBackground(radiationBackground, atmosphere_size)
        check_atom(atom, atmosphere_size)
        check_populations(populations_LTE, atmosphere_size)
        check_populations(populations_zero_radiation, atmosphere_size)
        check_rates(rates, atmosphere_size)
        check_radiation(radiation, atom, atmosphere_size)
        println("All tests passed.")
    end

    if plot == true
        print("--Plot all parameters.....")
        z = atmosphere.z[1:end-1]
        plot_atmosphere(atmosphere)
        plot_radiationBackground(radiationBackground, z)
        plot_populations(populations_LTE, populations_zero_radiation, z)
        plot_rates(rates, z)
        plot_radiation(radiation, atom, z)
        println("All parameters plotted.\n")
    end
end

run_tests(true, true)
