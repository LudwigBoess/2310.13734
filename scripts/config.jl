# path to the package envirnoment
const global dir_path = joinpath(@__DIR__, "..")

# activate the environment
# using Pkg;
# Pkg.activate(dir_path);


# paths to data, plots and maps
const global data_path = dir_path * "/data/"
const global plot_path = dir_path * "/Plots/"
const global map_path = dir_path * "/maps/"