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
@everywhere using PyPlotUtility
println("done")
flush(stdout);
flush(stderr);

# mapping settings
@everywhere const snap = 36
@everywhere const global sim_path = "path/to/simulation/"
@everywhere const snap_base = sim_path * "snapdir_$(@sprintf("%03i", snap))/snap_$(@sprintf("%03i", snap))"
@everywhere const map_path = joinpath(@__DIR__, "..", "..", "data", "phase_data")
@everywhere const global GU = GadgetPhysical(GadgetIO.read_header(snap_base))


const jobs = RemoteChannel(() -> Channel{Int}(32))
const results = RemoteChannel(() -> Channel{Matrix}(80000))

# phase map settings
@everywhere const x_lim = [1.e-9, 1.0]
@everywhere const y_lim = [10.0, 1.e10]

@everywhere const Nbins = 100

@everywhere function get_phase_map(rho, T, quantity)
    bin_2D_log(rho, T, x_lim, y_lim, quantity,
        Nbins=Nbins, show_progress=false, calc_mean=false)
end

@everywhere function get_CReE(data)

    Npart = length(data["CReC"])
    CReE = Vector{Float64}(undef, Npart)

    # cr setup 
    Nbins = size(data["CReN"], 1)
    par = CRMomentumDistributionConfig(0.1, 1.e5, Nbins)
    bounds = momentum_bin_boundaries(par)

    @threads for i = 1:Npart

        norm = 10.0 .^ data["CReN"][:, i]
        slope = Float64.(data["CReS"][:, i])
        cut = Float64(data["CReC"][i])
        rho    = data["RHO"][i] * GU.rho_physical
        m = data["MASS"][i] * GU.m_physical

        CReE[i] = cr_energy_in_range(norm, slope, cut, rho, bounds) * m

    end

    # convert to erg
    return CReE .* GU.E_cgs
end

@everywhere function CReE_phase_map_of_subfile(subfile)

    println("CReE: subfile $subfile running on $(nthreads()) threads")
    flush(stdout)
    flush(stderr)

    data = Dict(block => read_block(snap_base * ".$subfile", block, parttype=0)
                for block ∈ ["MASS", "RHO", "U",
        "CReN", "CReS", "CReC"])

    println("\ttemperature and density")
    flush(stdout)
    flush(stderr)
    T = data["U"] .* GU.T_K
    ne = data["RHO"] .* GU.rho_ncm3

    # calculate synchrotron emissivity
    CReE = get_CReE(data)

    phase_CReE = get_phase_map(ne, T, CReE)
    CReE = nothing
    data = nothing
    T = ne = nothing
    GC.gc()

    return phase_CReE
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
        put!(results, CReE_phase_map_of_subfile(job_id - 1))
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

    n = 2048
    #n = 4

    sum_phase_CReE = zeros(Float64, Nbins, Nbins)

    errormonitor(@async make_jobs(n)) # feed the jobs channel with "n" jobs

    println("running")
    flush(stdout)
    flush(stderr)

    @time while n > 0 # print out results

        phase_CReE = take!(results)

        # sum up contribution
        sum_phase_CReE += phase_CReE
        n -= 1

        phase_CReE = nothing
    end
    flush(stdout)
    flush(stderr)

    println("Pnu phase map")
    println("    Min:    $(minimum(sum_phase_CReE))")
    println("    Max:    $(maximum(sum_phase_CReE))")
    println("    Mean:   $(mean(sum_phase_CReE))")
    println("    Median: $(median(sum_phase_CReE))")

    flush(stdout)
    flush(stderr)
    write_phase_map(filename, sum_phase_CReE)
end

filename = map_path * "phase_map_CReE_high.dat"
run_phase_maps(filename)
