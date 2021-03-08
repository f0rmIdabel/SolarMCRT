include("check_parameters.jl")
include("plot_parameters.jl")

function run_tests(check=true, plot=false)

    # =============================================================================
    # ATMOSPHERE
    # =============================================================================
    atmosphere_parameters = collect_atmosphere_data()
    atmosphere = Atmosphere(atmosphere_parameters...)

    # =============================================================================
    # BACKGROUND RADIATION
    # =============================================================================
    λ = get_background_λ()
    radiation_parameters = collect_radiation_data(atmosphere, λ)
    radiationBackground = RadiationBackground(radiation_parameters...)

    # =============================================================================
    # ATOM
    # =============================================================================
    atom_parameters = collect_atom_data(atmosphere)
    atom = Atom(atom_parameters...)

    # =============================================================================
    # INITIAL POPULATIONS
    # =============================================================================
    populations = collect_initial_populations()

    # =============================================================================
    # INITIAL TRANSITION RATES
    # =============================================================================
    Bλ = blackbody_lambda(atom.λ, atmosphere.temperature)
    rate_parameters = calculate_transition_rates(atom, atmosphere, populations, Bλ)
    rates = TransitionRates(rate_parameters...)

    # =============================================================================
    # RADIATION
    # =============================================================================
    radiation_parameters = collect_radiation_data(atmosphere, atom, rates, populations)
    radiation = Radiation(radiation_parameters...)

    if check == true
        print("--Run all tests...........")
        atmosphere_size = size(atmosphere.temperature)
        check_atmosphere(atmosphere)
        check_radiationBackground(radiationBackground, atmosphere_size)
        check_atom(atom, atmosphere_size)
        check_populations(populations, atmosphere_size)
        check_rates(rates, atmosphere_size)
        check_radiation(radiation, atom, atmosphere_size)
        println("All tests passed.")
    end

    if plot == true
        print("--Plot all parameters.....")
        z = atmosphere.z[1:end-1]
        plot_atmosphere(atmosphere)
        plot_radiationBackground(radiationBackground, z)
        plot_populations(populations, z)
        plot_rates(rates, z)
        plot_radiation(radiation, atom, z)
        println("All parameters plotted.\n")
    end
end

run_tests(true, true)


"""
I am not yet sure if my LTE population function works,
so whenever I change it, I need to update the
atmosphere file and the initial populations file.
"""
function update_population_input()

    # =============================================================================
    # ATMOSPHERE
    # =============================================================================
    atmosphere_parameters = collect_atmosphere_data()
    atmosphere = Atmosphere(atmosphere_parameters...)

    # =============================================================================
    # ATOM
    # =============================================================================
    atom_parameters = collect_atom_data(atmosphere)
    atom = Atom(atom_parameters...)

    # =============================================================================
    # INITIAL POPULATIONS
    # =============================================================================
    old_populations = collect_initial_populations()
    populations = LTE_populations(atom, old_populations, atmosphere.temperature, atmosphere.electron_density)

    """h5open(get_initial_populations_path(), "r+") do file
        delete_object(file, "populations")
        write(file, "populations", ustrip.(populations))
    end

    h5open(get_atmosphere_path(), "r+") do file
        delete_object(file, "hydrogen_populations")
        write(file, "hydrogen_populations", ustrip.(populations))
    end"""

    println(minimum(old_populations))
    println(minimum(populations))
end
