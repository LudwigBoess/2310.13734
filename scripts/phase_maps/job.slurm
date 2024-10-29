#!/bin/bash
#SBATCH -J phase_jnu_1D                        # name of the job
#SBATCH -o ./log/%x.%j.out                       # output log file with name <job_name>.<job_id>.out
#SBATCH -e ./log/%x.%j.err                       # error log file with name <job_name>.<job_id>.out
#SBATCH -D ./                                # output directory
#SBATCH --nodes=1                            # number of nodes (needs to be 1 currently)
#SBATCH --ntasks-per-node=4                 # number of MPI ranks per node
#SBATCH --cpus-per-task=2                    # number of OpenMP threads per MPI rank
#SBATCH --mem=100G
#SBATCH --time=168:00:00                      # time limit of the run

# set number of OpenMP 
export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK}

#julia shob_2D.jl

#julia synch_jnu_2D.jl sim
#julia synch_jnu_2D.jl beta
#julia synch_jnu_2D.jl vturb
#julia synch_jnu_2D.jl ff
#julia synch_jnu_2D.jl dyn_l
#julia synch_jnu_2D.jl dyn_h

# julia B_ne_1D.jl sim
# julia B_ne_1D.jl beta
#julia B_ne_1D.jl vturb
# julia B_ne_1D.jl ff
# julia B_ne_1D.jl dyn_l
# julia B_ne_1D.jl dyn_h

#julia CReE_2D.jl

# set number of OpenMP 
#export JULIA_NUM_THREADS=1
#julia synch_jnu_2D.jl sim
#julia synch_jnu_2D.jl beta
julia synch_jnu_1D.jl vturb
#julia synch_jnu_2D.jl ff
#julia synch_jnu_2D.jl dyn_l
#julia synch_jnu_2D.jl dyn_h