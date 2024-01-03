#!/bin/bash
#SBATCH -J Pnu_2D      # name of the job
#SBATCH -o ./log/%x.%j.out      # output log file with name <job_name>.<job_id>.out
#SBATCH -e ./log/%x.%j.err      # error log file with name <job_name>.<job_id>.out
#SBATCH -D ./                   # output directory
#SBATCH --nodes=1               # number of nodes (needs to be 1 currently)
#SBATCH --ntasks-per-node=2    # number of MPI ranks per node
#SBATCH --cpus-per-task=8   # number of OpenMP threads per MPI rank
#SBATCH --mem=50G
#SBATCH --time=168:00:00         # time limit of the run

source ~/.bashrc 

# set number of OpenMP 
export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK}

#julia synch_Pnu_2D.jl sim
#julia synch_Pnu_2D.jl beta
#julia synch_Pnu_2D.jl vturb
#julia synch_Pnu_2D.jl ff
#julia synch_Pnu_2D.jl dyn_l
julia synch_Pnu_2D.jl dyn_h