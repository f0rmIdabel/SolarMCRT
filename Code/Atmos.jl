using Unitful

struct Atmosphere
   # Dimensions
   x::Array{<:Unitful.Length, 1}
   y::Array{<:Unitful.Length, 1}
   z::Array{<:Unitful.Length, 1}

   # Local properties
   velocity_x::Array{<:Unitful.Velocity, 3}
   velocity_y::Array{<:Unitful.Velocity, 3}
   velocity_z::Array{<:Unitful.Velocity, 3}
   temperature::Array{<:Unitful.Temperature, 3}
   χ_continuum::Array{<:Unitful.Quantity, 3}
   ε_continuum::Array{Real, 3}
end
