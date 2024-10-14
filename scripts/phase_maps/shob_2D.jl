println("allocating cores")
using Distributed, ClusterManagers
# automatically decide if it needs to be run in slurm envirnoment
try
    println("allocating $(ENV["SLURM_NTASKS"]) slurm tasks")
    addprocs_slurm(parse(Int64, ENV["SLURM_NTASKS"]))
catch err
    if isa(err, KeyError)
        N_tasts_ = 4
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
@everywhere using PyPlotUtility
println("done")
flush(stdout);
flush(stderr);


# mapping settings
@everywhere const snap = 74
@everywhere const global sim_path = "/gpfs/work/pn36ze/di93son/LocalUniverse/Coma/L5/cr6p20e/"
@everywhere const snap_base = sim_path * "snapdir_$(@sprintf("%03i", snap))/snap_$(@sprintf("%03i", snap))"
@everywhere const map_path = joinpath(@__DIR__, "..", "..", "data", "phase_maps", "zoom_inj")
@everywhere const global GU = GadgetPhysical(GadgetIO.read_header(snap_base))


const jobs = RemoteChannel(() -> Channel{Int}(32))
const results = RemoteChannel(() -> Channel{Matrix}(80000))

# phase map settings
@everywhere const x_lim = log10.([1.e-9, 1.0])
@everywhere const y_lim = [0.0, 90.0]

@everywhere const Nbins = 100

@everywhere function M_phase_map_of_subfile(subfile)

    println("SHOB: subfile $subfile running on $(nthreads()) threads")
    flush(stdout)
    flush(stderr)


    println("\treading data")
    flush(stdout)
    flush(stderr)
    ne = log10.(read_block(snap_base * ".$subfile", "RHO", parttype=0) .* GU.rho_ncm3)

    M = read_block(snap_base * ".$subfile", "MASS", parttype=0) .* GU.m_msun

    # calculate synchrotron emissivity
    θ_B = rad2deg.(read_block(snap_base * ".$subfile", "SHOB", parttype=0))

    θ_B[isnan.(θ_B)] .= -1.0
    θ_B[isinf.(θ_B)] .= -1.0

    println("\tphase map")
    phase_M = bin_2D(ne, θ_B, x_lim, y_lim, M,
                Nbins=Nbins, show_progress=false, calc_mean=false)
    θ_B = nothing
    M = ne = nothing
    GC.gc()

    return phase_M
end


@everywhere function write_phase_map(filename, phase_map)
    Nbins = size(phase_map, 1)
    f = open(filename, "w")
    write(f, Nbins)
    write(f, x_lim)
    write(f, y_lim)
    write(f, phase_map)
    close(f)
end


@everywhere function do_work(jobs, results)
    while true
        job_id = take!(jobs)
        put!(results, M_phase_map_of_subfile(job_id - 1))
    end
end


function make_jobs(n)
    for i in 1:n
        put!(jobs, i)
    end
end


function run_phase_maps(filename)

    println("starting workers")

    for p in workers() # start tasks on the workers to process requests in parallel
        remote_do(do_work, p, jobs, results)
    end

    #n = 2048
    n = 16

    sum_phase_M = zeros(Float64, Nbins, Nbins)

    errormonitor(@async make_jobs(n)) # feed the jobs channel with "n" jobs

    println("running")
    flush(stdout)
    flush(stderr)

    @time while n > 0 # print out results

        phase_M = take!(results)

        # sum up contribution
        sum_phase_M += phase_M
        n -= 1

        phase_M = nothing
    end
    flush(stdout)
    flush(stderr)

    println("Pnu phase map")
    println("    Min:    $(minimum(sum_phase_M))")
    println("    Max:    $(maximum(sum_phase_M))")
    println("    Mean:   $(mean(sum_phase_M))")
    println("    Median: $(median(sum_phase_M))")

    flush(stdout)
    flush(stderr)
    write_phase_map(filename, sum_phase_M)
end

filename = map_path * "/phase_map_shob.dat"
run_phase_maps(filename)
