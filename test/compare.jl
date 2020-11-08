





function name()

    data = h5open("../out/output.hdf5", "r")

    J_mc = h5read(data, )
    J_feautrier =

    diff = abs.(J_mc .- J_feautrier) ./ J_feautrier

    hist(diff)
end
