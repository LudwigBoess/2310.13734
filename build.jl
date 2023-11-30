# load standard packages
using Pkg

# change into directory of this file
cd(@__DIR__)

# install non-registered package
Pkg.add(url="https://github.com/LudwigBoess/SpectralCRsUtility.jl")
Pkg.add(url="https://github.com/LudwigBoess/PyPlotUtility.jl")

# activate environment and install packages
Pkg.activate(".")
Pkg.instantiate()

# create folders for maps
mkdir("maps")
mkdir("data")
mkdir("Plots")
