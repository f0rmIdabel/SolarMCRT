include("check_parameters.jl")
include("plot_parameters.jl")

function run_tests()

    #populations = collect_initial_populations()
    #I_changed_the_LTE_calculation(populations)

    # =============================================================================
    # ATMOSPHERE
    # =============================================================================
    atmosphere_parameters = collect_atmosphere_data()
    atmosphere = Atmosphere(atmosphere_parameters...)

    check_atmosphere(atmosphere)
    plot_atmosphere(atmosphere)

    atmosphere_size = size(atmosphere.temperature)
    # =============================================================================
    # BACKGROUND RADIATION
    # =============================================================================
    λ = get_background_λ()
    radiation_parameters = collect_radiation_data(atmosphere, λ)
    radiationBackground = RadiationBackground(radiation_parameters...)

    check_radiationBackground(radiationBackground, atmosphere_size)
    #plot_radiationBackground(radiationBackground, atmosphere.z)

    # =============================================================================
    # ATOM
    # =============================================================================
    atom_parameters = collect_atom_data(atmosphere)
    atom = Atom(atom_parameters...)

    check_atom(atom, atmosphere_size)

    # =============================================================================
    # INITIAL POPULATIONS
    # =============================================================================
    populations = collect_initial_populations()

    check_populations(populations, atmosphere_size)
    plot_populations(populations)

    # =============================================================================
    # INITIAL TRANSITION RATES
    # =============================================================================
    Bλ = blackbody_lambda(atom.λ, atmosphere.temperature)
    rate_parameters = calculate_transition_rates(atom, atmosphere, populations, Bλ)
    rates = TransitionRates(rate_parameters...)

    check_rates(rates, atmosphere_size)
    #plot_rates(rates, atmosphere.z)

    # =============================================================================
    # RADIATION
    # =============================================================================
    radiation_parameters = collect_radiation_data(atmosphere, atom, rates, populations)
    radiation = Radiation(radiation_parameters...)
    check_radiation(radiation, atom, atmosphere_size)
    #plot_radiation(radiation, atmosphere.z, atom.λ)
end

"""
I am not yet sure if my LTE population function works,
so whenever I change it, I need to update the
atmosphere file and the initial populations file.
"""
function I_changed_the_LTE_calculation(populations)


    h5open(get_initial_populations_path(), "r+") do file
           delete_object(file, "populations")
           write(file, "populations", ustrip.(populations))
    end

    h5open(get_atmosphere_path(), "r+") do file
        delete_object(file, "hydrogen_populations")
        write(file, "hydrogen_populations", ustrip.(populations))
    end
end

run_tests()
