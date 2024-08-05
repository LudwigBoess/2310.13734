println("allocating cores")
using Distributed, ClusterManagers

# automatically decide if it needs to be run in slurm envirnoment
try
    println("allocating $(ENV["SLURM_NTASKS"]) slurm tasks")
    addprocs_slurm(parse(Int64, ENV["SLURM_NTASKS"]))
catch err
    if isa(err, KeyError)
        println("allocating 2 normal tasks")
        addprocs(2)
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

const global data_path = "/gpfs/work/pn68va/di67meg/PaperRepos/SynchWeb/maps/coma/debug/"
#const global sim_path = "/mnt/home/lboess/ceph/LocalUniverse/Coma/L4/cr4p16e_pinj/"
const global sim_path = "/gpfs/work/pn68va/di67meg/LocalUniverse/"

const global snap = 36
const global snap_base = sim_path * "snapdir_$(@sprintf("%03i", snap))/snap_$(@sprintf("%03i", snap))"
const global GU = GadgetPhysical(read_header(snap_base))

function get_gamma(data, h)

    Npart = length(data["CRpC"])
    Fγ = Vector{Float64}(undef, Npart)

    # cr setup 
    Nbins = size(data["CRpN"], 1)
    par = CRMomentumDistributionConfig(0.1, 1.e5, Nbins)
    # construct boundaries 
    bounds = momentum_bin_boundaries(par)

    V = 4π
    d = 1.0

    p = Progress(Npart)
    @threads for i = 1:Npart

        norm = GU.CR_norm .* 10.0 .^ data["CRpN"][:, i]

        slope = Float64.(data["CRpS"][:, i])
        cut = Float64(data["CRpC"][i])

        nH = data["RHO"][i] * GU.rho_ncm3

        Fγ[i] = gamma_flux_pions(norm, slope, cut, bounds, nH, V, d, N_integration_steps=20)

        next!(p)
    end

    Fγ
end


function make_quantity_maps(snap, cluster, gpos, side_length, scale)

    # default
    image_path = data_path * "$(cluster)_$(scale)_$(@sprintf("%03i", snap))."

    @info "reading data"
    h = read_header(snap_base)
    blocks = ["POS", "HSML", "RHO", "U", "MASS",
        "CRpN", "CRpS", "CRpC",
        "CRpP"
    ]

    data = read_particles_in_volume(snap_base, blocks, gpos, side_length, use_keys=true)

    @info "done"

    # select kernel
    kernel = WendlandC4(Float64, 2)

    h = read_header(snap_base)

    # convert to physical code units for mapping
    pos = data["POS"] .* GU.x_physical
    hsml = data["HSML"] .* GU.x_physical
    rho = data["RHO"] .* GU.rho_physical
    mass = data["MASS"] .* GU.m_physical

    xy_size = 2side_length
    z_size = 2side_length

    # define mapping parameters
    param = mappingParameters(center=gpos .* GU.x_physical,
        x_size=xy_size * GU.x_physical,
        y_size=xy_size * GU.x_physical,
        z_size=z_size * GU.x_physical,
        Npixels=2048)

    @info "gamma"
    flush(stdout)
    flush(stderr)

    image_prefix = image_path * "gamma"
    j_gamma = get_gamma(data, GU)

    weights = part_weight_physical(length(j_gamma), param)
    map_it(pos, hsml, mass, rho, j_gamma, weights,
        units="erg/cm^2",
        reduce_image=false,
        parallel=true,
        projection="xz";
        kernel, snap, param, image_prefix)

    j_gamma = nothing
    GC.gc()
    
    @info "CRpP"
    flush(stdout); flush(stderr);

    image_prefix = image_path * "CRpP"
    CRpP = @. data["CRpP"] * GU.P_CR_cgs
    map_it(pos, hsml, mass, rho, CRpP, rho, 
            units="erg/cm^3", reduce_image=true,
            projection="xz"; 
            kernel, snap, param, image_prefix)


    @info "Xcr"
    flush(stdout)
    flush(stderr)

    image_prefix = image_path * "Xcr"
    Pth = @. 2/3 * data["RHO"] * data["U"] * GU.P_th_cgs
    Xcr = CRpP ./ Pth

    map_it(pos, hsml, mass, rho, Xcr, rho,
        units="", calc_mean=false,
        projection="xz";
        kernel, snap, param, image_prefix)

    Xcr = Pth = nothing
    GC.gc()

    data = nothing
    GC.gc()
end


# coma
gpos = [244977.58, 327853.62, 245989.75]
cluster = "coma"

side_length = 0.1 * 100.0e3

make_quantity_maps(snap, cluster, gpos, side_length, "20Mpc")
