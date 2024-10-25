println("allocating cores")
using Distributed, ClusterManagers
# automatically decide if it needs to be run in slurm envirnoment
try
    println("allocating $(ENV["SLURM_NTASKS"]) slurm tasks")
    addprocs_slurm(parse(Int64, ENV["SLURM_NTASKS"]))
catch err
    if isa(err, KeyError)
        N_tasts_ = 2
        println("allocating $N_tasts_ normal tasks")
        addprocs(N_tasts_)
    end
end

println("done")
flush(stdout);
flush(stderr);

println("loading packages")
@everywhere using GadgetIO, GadgetUnits
@everywhere using Printf
@everywhere using SpectralCRsUtility
@everywhere using Base.Threads
@everywhere using Statistics
@everywhere using ProgressMeter
println("done")
flush(stdout);
flush(stderr);

if ARGS[1] == "sim"
    @everywhere const global Bfield_flag = 1
elseif ARGS[1] == "beta"
    @everywhere const global Bfield_flag = 2
elseif ARGS[1] == "vturb"
    @everywhere const global Bfield_flag = 3
elseif ARGS[1] == "ff"
    @everywhere const global Bfield_flag = 4
elseif ARGS[1] == "dyn_l"
    @everywhere const global Bfield_flag = 5
elseif ARGS[1] == "dyn_h"
    @everywhere const global Bfield_flag = 6
end


# @everywhere const snap = 36
# #@everywhere const global sim_path = "/e/ocean2/users/lboess/LocalUniverseZooms/L5/mhd_cr8p24e/"
# @everywhere const global sim_path = "/path/to/simulation/"
@everywhere const snap_base = "/e/ocean3/Local/3072/nonrad_mhd_crs_new/snapdir_000_z=0/snap_000"
@everywhere const global GU = GadgetPhysical(GadgetIO.read_header(snap_base))

@everywhere const map_path = joinpath(@__DIR__, "..", "..", "data", "phase_maps", "box")

# SLOW 1 paper
@everywhere const global center_comov = [247.980, 245.480, 255.290] .* 1.e3
@everywhere const global center = center_comov .* GU.x_physical
@everywhere global const radius_limits = [10_000.0, 240_000.0 * GU.x_physical]

@everywhere include(joinpath(@__DIR__, "bin_1D.jl"))
@everywhere include(joinpath(@__DIR__, "..", "allsky", "Bfld.jl"))

const jobs = RemoteChannel(() -> Channel{Int}(32))
const results = RemoteChannel(() -> Channel{Tuple}(1600))

# phase map settings
@everywhere const x_lim = [1.e-9, 1.0]
@everywhere const Nbins = 100


@everywhere function B_1D_bin_of_subfile(subfile, Bfield_function, Btype)

    println("B 1D (B$Btype): subfile $subfile running on $(nthreads()) threads")
    flush(stdout)
    flush(stderr)

    data = Dict(block => read_block(snap_base * ".$subfile", block, parttype=0)
                for block ∈ ["POS", "MASS", "RHO", "U", "BFLD", "VRMS"])

    ne = data["RHO"] .* GU.rho_ncm3


    B = Vector{Float64}(undef, length(ne))

    # calculate Bfield 
    #P = Progress(length(B))
    @threads for i ∈ 1:length(B)
        B[i] = Bfield_function(data, i)
        #next!(P)
    end

    ne_count, B_sum = bin_1D_log(ne, x_lim, B, show_progress=false, calc_sigma=false)
    j_nu = nothing
    data = nothing
    ne = nothing
    GC.gc()

    return ne_count, B_sum
end



@everywhere function write_B_binning(filename, B, B_sigma)
    Nbins = length(B)
    f = open(filename, "w")
    write(f, Nbins)
    write(f, x_lim)
    write(f, B)
    write(f, B_sigma)
    close(f)
end


@everywhere function B_1D_bin_of_subfile(subfile)
    if Bfield_flag == 1
        return B_1D_bin_of_subfile(subfile, Bfield_sim, "sim")
    elseif Bfield_flag == 2
        return B_1D_bin_of_subfile(subfile, Bfield_Beta, "beta")
    elseif Bfield_flag == 3
        return B_1D_bin_of_subfile(subfile, Bfield_vturb, "vturb")
    elseif Bfield_flag == 4
        return B_1D_bin_of_subfile(subfile, Bfield_FF, "ff")
    elseif Bfield_flag == 5
        return B_1D_bin_of_subfile(subfile, Bfield_dyn_l, "dyn_l")
    elseif Bfield_flag == 6
        return B_1D_bin_of_subfile(subfile, Bfield_dyn_h, "dyn_h")
    end
end


@everywhere function get_B_phase_map_filename()

    if Bfield_flag == 1
        return map_path * "/bin_1D_B_sim.dat"
    elseif Bfield_flag == 2
        return map_path * "/bin_1D_B_beta50.dat"
    elseif Bfield_flag == 3
        return map_path * "/bin_1D_B_Pturb.dat"
    elseif Bfield_flag == 4
        return map_path * "/bin_1D_B_FF.dat"
    elseif Bfield_flag == 5
        return map_path * "/bin_1D_B_dyn_l.dat"
    elseif Bfield_flag == 6
        return map_path * "/bin_1D_B_dyn_h.dat"
    end

end


@everywhere function do_work(jobs, results) # define work function everywhere
    while true
        job_id = take!(jobs)
        put!(results, B_1D_bin_of_subfile(job_id - 1))
    end
end

function make_jobs(n)
    for i in 1:n
        put!(jobs, i)
    end
end


function run_B_1D()

    println("starting workers")

    for p in workers() # start tasks on the workers to process requests in parallel
        remote_do(do_work, p, jobs, results)
    end

    # n = 2048
    #n = 4

    n = GadgetIO.read_header(snap_base).num_files


    sum_ne_count = zeros(Int64, Nbins)
    sum_B = zeros(Float64, Nbins)

    errormonitor(@async make_jobs(n)) # feed the jobs channel with "n" jobs

    println("running")
    flush(stdout)
    flush(stderr)

    @time while n > 0 # print out results

        ne_count, B = take!(results)
        sum_ne_count .+= ne_count
        sum_B .+= B
        n -= 1

        ne_count = B = nothing
        GC.gc()
    end
    flush(stdout)
    flush(stderr)

    B_sigma = σ_1D_quantity(sum_B, sum_ne_count)


    for i = 1:length(sum_ne_count)
        if !iszero(sum_ne_count[i])
            sum_B[i] /= sum_ne_count[i]
        end
    end

    println("B 1D")
    println("    Min:    $(minimum(sum_B))")
    println("    Max:    $(maximum(sum_B))")
    println("    Mean:   $(mean(sum_B))")
    println("    Median: $(median(sum_B))")

    flush(stdout)
    flush(stderr)
    filename = get_B_phase_map_filename()
    write_B_binning(filename, sum_B, B_sigma)
end

run_B_1D()