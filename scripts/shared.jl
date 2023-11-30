
"""
    read_phase_map(filename)

Reads a 2D phase map from a file.
"""
function read_phase_map(filename)

    f = open(filename, "r")
    Nbins = read(f, Int64)
    x_lim = read!(f, Vector{Float64}(undef, 2))
    y_lim = read!(f, Vector{Float64}(undef, 2))
    phase_map = read!(f, Matrix{Float64}(undef, Nbins, Nbins))
    close(f)

    return x_lim, y_lim, phase_map
end

"""
    read_1D_data(filename)

Reads 1D data from a file.
"""
function read_1D_data(filename)

    f = open(filename, "r")
    Nbins = read(f, Int64)
    x_lim = read!(f, Vector{Float64}(undef, 2))
    B = read!(f, Vector{Float64}(undef, Nbins))
    close(f)

    dn = log10(x_lim[2] / x_lim[1]) / Nbins
    ne_bins = [x_lim[1] * 10.0 .^ ((i - 1) * dn) for i = 1:Nbins]

    return ne_bins, B
end