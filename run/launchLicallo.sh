#!/bin/bash
#SBATCH --job-name=Wavkins_Strat_Asymp_{$1}_M{$2}_kdinf{$3}_kdsup{$4}_kfh{$5}_kfz{$6}
##SBATCH --dependency=afterany:33465627
#SBATCH --dependency=singleton
##SBATCH --nodes=1
##SBATCH --ntasks=2
##SBATCH -c 20
##SBATCH --ntasks-per-node=1
#SBATCH -p x40 # fdr: nodes with 20 cores, x40: nodes with 40 cores, fdr-1t0 node with 32 cores, amd: nodes with 128 cores
#SBATCH --time=20:00:00

module purge


# Run bench
#echo "Bench amd"
#sh run_bench_threads.sh


# Run Stratified_Asymp
echo "type: $1"
echo "M: $2"
echo "kmax: $3"
echo "kdsup: $4"
echo "kfh (if localized) or kf (if isotropic): $5"
echo "kfz (if localized): $6"

let "nthreads = 4 * $2"
echo "number of threads: $nthreads"

julia --threads $nthreads RunCluster_Stratified_Asymp.jl $1 $2 $3 $4 $5 $6
