println("allocating cores")
using Distributed, ClusterManagers

# automatically decide if it needs to be run in slurm envirnoment
try
    println("allocating $(ENV["SLURM_NTASKS"]) slurm tasks")
    addprocs_slurm(parse(Int64, ENV["SLURM_NTASKS"]))
catch err
    if isa(err, KeyError)
        println("allocating 4 normal tasks")
        addprocs(4)
    end
end

println("done")
flush(stdout);
flush(stderr);

using GadgetIO, GadgetUnits
using SPHKernels, SPHtoGrid
using Printf
using ProgressMeter
using LinearAlgebra
using Base.Threads
using SpectralCRsUtility


#include("../../../config.jl")

const global center_comov = [247.980, 245.480, 255.290] .* 1.e3
const global radius_limits = [0.0, Inf]

const global data_path = "/gpfs/work/pn36ze/di93son/LocalUniverse/Coma/L5/cr6p20e/maps/"
const global sim_path = "/gpfs/work/pn36ze/di93son/LocalUniverse/Coma/L5/cr6p20e/"

const global snap = 74
const global snap_base = sim_path * "snapdir_$(@sprintf("%03i", snap))/snap_$(@sprintf("%03i", snap))"
const global GU = GadgetPhysical(read_header(snap_base))


include(joinpath(@__DIR__, "..", "allsky", "Bfld.jl"))

function rescale_bfield(data, Bfield_funtion)
    
    B = Matrix{Float32}(undef, 3, size(data["BFLD"], 2))
    
    @threads for i ∈ 1:size(data["BFLD"], 2)
        Bnorm = norm(data["BFLD"][:, i])
        B_ref = Bfield_funtion(data, i)
        
        if iszero(Bnorm)
            B_rescale = 0.0
        else
            B_rescale = B_ref / Bnorm
        end

        for j ∈ 1:3
            B[j, i] = data["BFLD"][j, i] * Float32(B_rescale)
        end
    end

    return B
end


function make_RM_maps(snap, cluster, gpos, side_length, scale)

    println("reading data")
    blocks = ["POS", "MASS", "HSML", "RHO", "U", "BFLD", "VRMS"]


    # default
    image_path = data_path * "$(cluster)_$(scale)_$(@sprintf("%03i", snap))."
    data = read_particles_in_volume(snap_base, blocks, gpos, side_length, use_keys=false)

    println("done")
    flush(stdout)
    flush(stderr)

    # select kernel
    kernel = WendlandC4(Float64, 2)

    h = read_header(snap_base)

    # convert to physical code units for mapping
    pos = data["POS"] .* GU.x_physical
    hsml = data["HSML"] .* GU.x_physical
    rho = data["RHO"] .* GU.rho_physical
    mass = data["MASS"] .* GU.m_physical

    ncm3 = data["RHO"] .* GU.rho_ncm3
    m_cgs = data["MASS"] .* GU.m_cgs

    xy_size = 2side_length
    z_size = 2side_length

    # define mapping parameters
    param = mappingParameters(center=gpos .* GU.x_physical,
        x_size=xy_size * GU.x_physical,
        y_size=xy_size * GU.x_physical,
        z_size=z_size * GU.x_physical,
        Npixels=1024)


    units = "rad/cm^2"
    synch = "RM"
    weights = part_weight_physical(length(ncm3), param)
    flux = false
    reduce_image = false


    # Synch
    println("B sim")
    flush(stdout)
    flush(stderr)

    image_prefix = image_path * "$(synch)_Bsim"
    B = rescale_bfield(data, Bfield_sim)
    RM = @. rotation_measure(ncm3, B[2,:])
    
    map_it(pos, hsml, mass, rho, RM, weights,
        parallel=true,
        projection="xz";
        reduce_image, units,
        kernel, snap, param, image_prefix)

    B = RM = nothing
    GC.gc()

    println("B from FF")
    flush(stdout)
    flush(stderr)

    image_prefix = image_path * "$(synch)_BFF"
    B = rescale_bfield(data, Bfield_FF)
    RM = @. rotation_measure(ncm3, B[2,:])
    
    map_it(pos, hsml, mass, rho, RM, weights,
        parallel=true,
        projection="xz";
        reduce_image, units,
        kernel, snap, param, image_prefix)

    B = RM = nothing
    GC.gc()

    println("B from beta = 50")
    flush(stdout)
    flush(stderr)

    image_prefix = image_path * "$(synch)_beta50"
    B = rescale_bfield(data, Bfield_Beta)
    RM = @. rotation_measure(ncm3, B[2,:])
    
    map_it(pos, hsml, mass, rho, RM, weights,
        parallel=true,
        projection="xz";
        reduce_image, units,
        kernel, snap, param, image_prefix)
    B = RM = nothing
    GC.gc()

    println("B dynamo h")
    flush(stdout)
    flush(stderr)

    image_prefix = image_path * "$(synch)_dyn_h"
    B = rescale_bfield(data, Bfield_dyn_h)
    RM = @. rotation_measure(ncm3, B[2,:])
    
    map_it(pos, hsml, mass, rho, RM, weights,
        parallel=true,
        projection="xz";
        reduce_image, units,
        kernel, snap, param, image_prefix)

    B = RM = nothing
    GC.gc()

    println("B dynamo l")
    flush(stdout)
    flush(stderr)

    image_prefix = image_path * "$(synch)_dyn_l"
    B = rescale_bfield(data, Bfield_dyn_l)
    RM = @. rotation_measure(ncm3, B[2,:])
    
    map_it(pos, hsml, mass, rho, RM, weights,
        parallel=true,
        projection="xz";
        reduce_image, units,
        kernel, snap, param, image_prefix)

    B = RM = nothing
    GC.gc()

    println("B vturb")
    flush(stdout)
    flush(stderr)

    image_prefix = image_path * "$(synch)_01Pturb"
    B = rescale_bfield(data, Bfield_vturb)
    RM = @. rotation_measure(ncm3, B[2,:])
    
    map_it(pos, hsml, mass, rho, RM, weights,
        parallel=true,
        projection="xz";
        reduce_image, units,
        kernel, snap, param, image_prefix)

    B = RM = nothing
    pos = hsml = mass = rho = weights = nothing
    data = nothing
    GC.gc()

end

# coma
gpos = [244977.58, 327853.62, 245989.75]
cluster = "coma"

side_length = 2_500.0

# SB 
make_RM_maps(snap, cluster, gpos, side_length, "5Mpc")