#!/bin/bash
#SBATCH -J allsky
#SBATCH -o ./log/%x.%j.out
#SBATCH -e ./log/%x.%j.err
#SBATCH -D ./
#SBATCH --cluster=c2pap
#SBATCH --partition=c2pap_parallel
#SBATCH --qos=c2pap_parallel
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=4   # number of OpenMP threads per MPI rank
#SBATCH --time=48:00:00
#SBATCH --export=NONE
#SBATCH --account=pn36ze
##SBATCH --mail-type=end
##SBATCH --mail-user=lboess@usm.lmu.de

# set number of OpenMP 
export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK}

julia synch_jnu_1D.jl sim
julia synch_jnu_1D.jl beta
julia synch_jnu_1D.jl vturb
julia synch_jnu_1D.jl ff
julia synch_jnu_1D.jl dyn_l
julia synch_jnu_1D.jl dyn_h
