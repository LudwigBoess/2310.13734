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

const global data_path = "/gpfs/work/pn36ze/di93son/PaperRepos/2310.13734/maps/zoom_HB/"
const global sim_path = "/gpfs/work/pn36ze/di93son/LocalUniverse/Coma/L5/cr6p20e/"

const global snap = 74
const global snap_base = sim_path * "snapdir_$(@sprintf("%03i", snap))/snap_$(@sprintf("%03i", snap))"
const global GU = GadgetPhysical(read_header(snap_base))


include(joinpath(@__DIR__, "..", "allsky", "Bfld.jl"))
include(joinpath(@__DIR__, "..", "allsky", "synchrotron.jl"))


function make_synch_maps(snap, cluster, gpos, side_length, scale, map_type)

    println("reading data")
    blocks = ["POS", "MASS", "HSML", "RHO", "U", "BFLD", "VRMS", "CReN", "CReS", "CReC"]


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

    rho_cgs = data["RHO"] .* GU.rho_cgs
    m_cgs = data["MASS"] .* GU.m_cgs

    xy_size = 2side_length
    z_size = 2side_length

    # define mapping parameters
    param = mappingParameters(center=gpos .* GU.x_physical,
        x_size=xy_size * GU.x_physical,
        y_size=xy_size * GU.x_physical,
        z_size=z_size * GU.x_physical,
        Npixels=2048)


    units = "erg/s/Hz/cm^2"
    factor = ones(length(m_cgs))
    weights = part_weight_physical(length(rho_cgs), param)
    flux = false
    reduce_image = false


    z_coma = 0.0231
    #nu = 1.4e9 * ( 1 + z_coma )
    #synch = "synch_Inu_1.4GHz"
    nu = 144.0e6 * (1 + z_coma)
    synch = "synch_Inu_144MHz"

    # Synch
    println("B sim")
    flush(stdout)
    flush(stderr)

    image_prefix = image_path * "$(synch)_Bsim"
    j_ν = get_synchrotron(data, nu, Bfield_sim, flux)

    
    map_it(pos, hsml, mass, rho, j_ν, weights,
        parallel=true,
        projection="xz";
        reduce_image, units,
        kernel, snap, param, image_prefix)

    j_ν = nothing
    GC.gc()

    println("B from FF")
    flush(stdout)
    flush(stderr)

    image_prefix = image_path * "$(synch)_BFF"
    j_ν = get_synchrotron(data, nu, Bfield_FF, flux)
    

    map_it(pos, hsml, mass, rho, j_ν, weights,
        parallel=true,
        projection="xz";
        reduce_image, units,
        kernel, snap, param, image_prefix)

    j_ν = nothing
    GC.gc()

    println("B from beta = 50")
    flush(stdout)
    flush(stderr)

    image_prefix = image_path * "$(synch)_beta50"
    j_ν = get_synchrotron(data, nu, Bfield_Beta, flux)
    

    map_it(pos, hsml, mass, rho, j_ν, weights,
        parallel=true,
        projection="xz";
        reduce_image, units,
        kernel, snap, param, image_prefix)
    j_ν = nothing
    GC.gc()

    println("B dynamo h")
    flush(stdout)
    flush(stderr)

    image_prefix = image_path * "$(synch)_dyn_h"
    j_ν = get_synchrotron(data, nu, Bfield_dyn_h, flux)
    

    map_it(pos, hsml, mass, rho, j_ν, weights,
        parallel=true,
        projection="xz";
        reduce_image, units,
        kernel, snap, param, image_prefix)

    j_ν = nothing
    GC.gc()

    println("B dynamo l")
    flush(stdout)
    flush(stderr)

    image_prefix = image_path * "$(synch)_dyn_l"
    j_ν = get_synchrotron(data, nu, Bfield_dyn_l, flux)
    

    map_it(pos, hsml, mass, rho, j_ν, weights,
        parallel=true,
        projection="xz";
        reduce_image, units,
        kernel, snap, param, image_prefix)

    j_ν = nothing
    GC.gc()

    println("B vturb")
    flush(stdout)
    flush(stderr)

    image_prefix = image_path * "$(synch)_01Pturb"
    j_ν = get_synchrotron(data, nu, Bfield_vturb, flux)
    

    map_it(pos, hsml, mass, rho, j_ν, weights,
        parallel=true,
        projection="xz";
        reduce_image, units,
        kernel, snap, param, image_prefix)

    j_ν = nothing
    pos = hsml = mass = rho = weights = nothing
    data = nothing
    GC.gc()

end

function get_CReE(data, GU)

    Npart = length(data["CReC"])
    CReE = Vector{Float64}(undef, Npart)

    # cr setup 
    Nbins = size(data["CReN"], 1)
    par = CRMomentumDistributionConfig(10.0, 1.e5, Nbins)
    bounds = momentum_bin_boundaries(par)

    p = Progress(Npart)

    @threads for i = 1:Npart

        norm = 10.0 .^ data["CReN"][:, i]
        slope = Float64.(data["CReS"][:, i])
        cut = Float64(data["CReC"][i])

        CReE[i] = cr_energy_in_range(norm, slope, cut, 1.0, bounds)

        next!(p)
    end

    # convert to erg
    return CReE .* (GU.E_cgs / GU.x_cgs^3)

end


function make_quantity_maps(snap, cluster, gpos, side_length, scale)

    # default
    image_path = data_path * "$(cluster)_$(scale)_$(@sprintf("%03i", snap))."

    @info "reading data"
    h = read_header(snap_base)
    blocks = ["POS", "HSML", "RHO", "U", "MASS",
        "BFLD",
        "CReP",
        "CReN", "CReS", "CReC",
        "DPP"
    ]

    data = read_particles_in_volume(snap_base, blocks, gpos, side_length, use_keys=false)

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

    @info "CReE"
    flush(stdout)
    flush(stderr)

    image_prefix = image_path * "CReE_gt1GeV"
    CReE = get_CReE(data, GU)

    weights = part_weight_physical(length(CReE), param)
    map_it(pos, hsml, mass, rho, CReE, weights,
        units="erg/cm^2", calc_mean=false,
        reduce_image=false,
        parallel=true,
        projection="xz";
        kernel, snap, param, image_prefix)

    CReE = nothing
    GC.gc()

    @info "B"
    flush(stdout); flush(stderr);

    image_prefix = image_path * "B"
    B = @. sqrt( data["BFLD"][1,:]^2  + data["BFLD"][2,:]^2 + data["BFLD"][3,:]^2 ) * 1.e6
    map_it(pos, hsml, mass, rho, B, rho, 
            units="muG", reduce_image=true,
            projection="xz"; 
            kernel, snap, param, image_prefix)


    B = weights = nothing 
    GC.gc()

    @info "rho"
    flush(stdout); flush(stderr);

    image_prefix = image_path * "rho"
    rho_cgs = data["RHO"] .* GU.rho_cgs
    weights = part_weight_physical(length(rho_cgs), param)

    map_it(pos, hsml, mass, rho, rho_cgs, weights, 
            units="g/cm^2", reduce_image=false,
            projection="xz"; 
            kernel, snap, param, image_prefix)

    rho_cgs = weights = nothing 
    GC.gc()


    @info "D0"
    flush(stdout)
    flush(stderr)

    image_prefix = image_path * "D0_max"
    D0 = Float64.(data["DPP"])
    # map_it(pos, hsml, mass, rho, D0, rho, 
    #         units="1/s", calc_mean=false; 
    #         kernel, snap, param, image_prefix)
    map_it(pos, hsml, mass, rho, D0, rho,
        units="1/s", calc_mean=false,
        projection="xz";
        kernel, snap, param, image_prefix)
    # map_it(pos, hsml, mass, rho, D0, rho, 
    #         units="1/s", calc_mean=false,
    #         projection="yz";  
    #         kernel, snap, param, image_prefix)

    D0 = nothing
    GC.gc()

    data = nothing
    GC.gc()
end


function make_HB_synch_maps(snap, cluster, gpos, side_length, scale, map_type)

    println("reading data")
    blocks = ["POS", "MASS", "HSML", "RHO", "U", "BFLD", "VRMS", "MACH"]


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

    rho_cgs = data["RHO"] .* GU.rho_cgs
    m_cgs = data["MASS"] .* GU.m_cgs

    xy_size = 2side_length
    z_size = 2side_length

    # define mapping parameters
    param = mappingParameters(center=gpos .* GU.x_physical,
        x_size=xy_size * GU.x_physical,
        y_size=xy_size * GU.x_physical,
        z_size=z_size * GU.x_physical,
        Npixels=1024)

    units = "erg/s/Hz/cm^2"
    synch = "synch_Inu_HB_144MHz"
    weights = part_weight_physical(length(rho_cgs), param)
    flux = false
    reduce_image = false


    z_coma = 0.0231
    #nu = 144.0e6 * ( 1 + z_coma )
    nu = 144.0e6 * (1 + z_coma)

    γ_m1 = 5 / 3 - 1
    rho_cgs = @. data["RHO"] * GU.rho_cgs
    T_keV = @. data["U"] * GU.T_eV * 1.e-3
    m_cgs = @. data["MASS"] * GU.m_cgs
    hsml_cgs = @. data["HSML"] * GU.x_cgs
    mach = Float64.(data["MACH"])

    Npart = length(rho)

    B = Vector{Float64}(undef, Npart)

    # Synch
    println("B sim")
    flush(stdout)
    flush(stderr)

    image_prefix = image_path * "$(synch)_Bsim"

    @threads for i = 1:Npart
        B[i] = Bfield_sim(data, i)
    end

    println("Bmax = ", maximum(B) * 1.e6, " muG")
    flush(stdout)
    flush(stderr)

    j_ν = analytic_synchrotron_HB07(rho_cgs, m_cgs, hsml_cgs, B, T_keV, mach, ν0=nu)

    

    map_it(pos, hsml, mass, rho, j_ν, weights,
        parallel=true,
        projection="xz";
        reduce_image, units,
        kernel, snap, param, image_prefix)

    j_ν = nothing
    GC.gc()

    println("B from FF")
    flush(stdout)
    flush(stderr)

    image_prefix = image_path * "$(synch)_BFF"
    @threads for i = 1:Npart
        B[i] = Bfield_FF(data, i)
    end

    println("Bmax = ", maximum(B) * 1.e6, " muG")
    flush(stdout)
    flush(stderr)

    j_ν = analytic_synchrotron_HB07(rho_cgs, m_cgs, hsml_cgs, B, T_keV, mach, ν0=nu)
    

    map_it(pos, hsml, mass, rho, j_ν, weights,
        parallel=true,
        projection="xz";
        reduce_image, units,
        kernel, snap, param, image_prefix)

    j_ν = nothing
    GC.gc()

    println("B from beta = 50")
    flush(stdout)
    flush(stderr)

    image_prefix = image_path * "$(synch)_beta50"
    @threads for i = 1:Npart
        B[i] = Bfield_Beta(data, i)
    end

    println("Bmax = ", maximum(B) * 1.e6, " muG")
    flush(stdout)
    flush(stderr)

    j_ν = analytic_synchrotron_HB07(rho_cgs, m_cgs, hsml_cgs, B, T_keV, mach, ν0=nu)
    

    map_it(pos, hsml, mass, rho, j_ν, weights,
        parallel=true,
        projection="xz";
        reduce_image, units,
        kernel, snap, param, image_prefix)
    j_ν = nothing
    GC.gc()

    println("B dynamo h")
    flush(stdout)
    flush(stderr)

    image_prefix = image_path * "$(synch)_dyn_h"
    @threads for i = 1:Npart
        B[i] = Bfield_dyn_h(data, i)
    end

    println("Bmax = ", maximum(B) * 1.e6, " muG")
    flush(stdout)
    flush(stderr)

    j_ν = analytic_synchrotron_HB07(rho_cgs, m_cgs, hsml_cgs, B, T_keV, mach, ν0=nu)
    

    map_it(pos, hsml, mass, rho, j_ν, weights,
        parallel=true,
        projection="xz";
        reduce_image, units,
        kernel, snap, param, image_prefix)

    j_ν = nothing
    GC.gc()

    println("B dynamo l")
    flush(stdout)
    flush(stderr)

    image_prefix = image_path * "$(synch)_dyn_l"
    @threads for i = 1:Npart
        B[i] = Bfield_dyn_l(data, i)
    end

    println("Bmax = ", maximum(B) * 1.e6, " muG")
    flush(stdout)
    flush(stderr)

    j_ν = analytic_synchrotron_HB07(rho_cgs, m_cgs, hsml_cgs, B, T_keV, mach, ν0=nu)
    

    map_it(pos, hsml, mass, rho, j_ν, weights,
        parallel=true,
        projection="xz";
        reduce_image, units,
        kernel, snap, param, image_prefix)

    j_ν = nothing
    GC.gc()

    println("B vturb")
    flush(stdout)
    flush(stderr)

    image_prefix = image_path * "$(synch)_01Pturb"
    @threads for i = 1:Npart
        B[i] = Bfield_vturb(data, i)
    end

    println("Bmax = ", maximum(B) * 1.e6, " muG")
    flush(stdout)
    flush(stderr)

    j_ν = analytic_synchrotron_HB07(rho_cgs, m_cgs, hsml_cgs, B, T_keV, mach, ν0=nu)
    

    map_it(pos, hsml, mass, rho, j_ν, weights,
        parallel=true,
        projection="xz";
        reduce_image, units,
        kernel, snap, param, image_prefix)

    j_ν = nothing
    pos = hsml = mass = rho = weights = nothing
    data = nothing
    GC.gc()

end

# coma
gpos = [244977.58, 327853.62, 245989.75]
cluster = "coma"

side_length = 0.1 * 100.0e3

# SB 
#make_synch_maps(snap, cluster, gpos, side_length, "20Mpc", "SB")
make_HB_synch_maps(snap, cluster, gpos, side_length, "20Mpc", "SB")

# quantities
#make_quantity_maps(snap, cluster, gpos, side_length, "20Mpc")