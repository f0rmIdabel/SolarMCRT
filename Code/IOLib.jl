module IOLib


"""
    read_parameter(textfile::AbstractString)

Reads 1D, 2D or 3D spatial data into Arrays. 
Assumes 
"""
function read_parameter(textfile::AbstractString)

    f = open(textfile)
    unit = readline(f)
    shape = split(readline(f)[2:end-1], ",")
    shape = filter!(e->e!="",shape)
    shape = parse.(Int, shape)
    dim = length(shape)
    skip_blankline = readline(f)

    if dim == 1
        array = parse.(Float64, split(readline(f)))

    elseif dim == 2
        array = Array{Float64, 2}(undef, tuple(shape...))

        for j=1:shape[2]
            x_vals = parse.(Float64, split(readline(f)))
            for i=1:shape[1]
                array[i,j] = x_vals[i]
            end
        end

    elseif dim == 3
        array = Array{Float64, 3}(undef, tuple(shape...))
        for k=1:shape[3]
            for j=1:shape[2]
                x_vals = parse.(Float64, split(readline(f)))
                for i=1:shape[1]
                    array[i,j,k] = x_vals[i]
                end
            end
            skip_blankline = readline(f)
        end

    end

    close(f)

    return array
end

end