include("lambda_iteration.jl")
include("shift_tools.jl")
#include("feautrier2.jl")

function run()
    println("\n", "="^83, "\n", " "^20,
              "INTEGRAL RT IN A SOLAR ATMOSPHERE MODEL",
              "\n", "="^83)

      # ==================================================================
      # LOAD WAVELENGTHS
      # ==================================================================
      print("\n--Loading wavelengths......................")
      λ = get_λ()[1]
      println(@sprintf(" %d wavelength(s) loaded.", length(λ)))

      # ==================================================================
      # LOAD ATMOSPHERE DATA AND CALCULATE BOUNDARY
      # ==================================================================
      print("--Loading atmosphere data..................")
      atmosphere_parameters = collect_atmosphere_data(λ, false, false)
      atmosphere = Atmosphere(atmosphere_parameters...)
      println(" Full atmosphere loaded.")

      # ==================================================================
      # FEAUTRIER CALCULATION
      # ==================================================================
      print("--Starting λ-iteration.....................")
      lambda_iteration(atmosphere, λ)
      println(" λ-iteration finished.")
end

run()



"""
function run2()
        λ = get_λ()[1]
        atmosphere_parameters = collect_atmosphere_data(λ, false, false)
        atmosphere = Atmosphere(atmosphere_parameters...)

        T = atmosphere.temperature
        z = atmosphere.z
        x = atmosphere.x
        pix = x[2] - x[1]
        z = z[1:end-1]
        T = permutedims(T, [2,3,1])

        T2 = copy(T)


        shift_variable!(T, z, pix, -1.0)
        shift_variable!(T, z, pix, -1.0)

        """
        translate!(T, z, pix, 0.5, π/2)
        T /= 0.5
        translate!(T, z, pix, 0.5, π/2)
        T /= 0.5
        translate!(T, z, pix, 0.5, π/2)
        T /= 0.5
        translate!(T, z, pix, 0.5, π/4)
        T /= 0.5
        translate!(T, z, pix, 0.5, π/8)
        T /= 0.5
        translate!(T, z, pix, 0.5, π/8)
        T /= 0.5
        """
        println(sum(abs.(T .- T2) ./ T2))
end"""
#run2()
