#!/bin/bash
#SBATCH --job-name=nbody
#SBATCH --output=./benchmark/data/epito/seq_nbody_report_run5.txt
#SBATCH --time=01:00:00
#SBATCH --partition=epito

# Array di potenze di 2 per il numero di corpi
bodies=(2 4 8 16 32 64 128 256 512 1024 2048 4096 8192 16384)

# Array per il numero di thread
threads=(1 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 34 36 38 40 42 44 46 48 50 52 54 56 58 60 62 64 66 68 70 72 74 76 78 80)

seed=5

schedtypes=("static" "dynamic")

# Eseguibile da testare
seq="/beegfs/home/andrea.zito/SCPD-Project/nbody_sequential_ap"
seq_avx="/beegfs/home/andrea.zito/SCPD-Project/nbody_sequential_ap_avx"
seq_omp="/beegfs/home/andrea.zito/SCPD-Project/nbody_omp_ap"

# make clean
make nbody_sequential_ap
make nbody_sequential_ap_avx
make nbody_omp_ap

# Loop sui numeri di corpi
for num_bodies in "${bodies[@]}"; do
  $seq $num_bodies $seed
  # for num_threads in "${threads[@]}"; do
    # export OMP_NUM_THREADS=$num_threads
    # echo "Running with $num_bodies bodies and $num_threads threads"
    # srun --ntasks=1 --cpus-per-task=$num_threads seq_omp $num_bodies
  # done
done

printf "\n"

for num_bodies in "${bodies[@]}"; do
  $seq_avx $num_bodies 
done

printf "\n"

for num_bodies in "${bodies[@]}"; do
  for sched_type in "${schedtypes[@]}"; do
    for num_threads in "${threads[@]}"; do
      $seq_omp $num_bodies $sched_type $((num_bodies / num_threads)) $num_threads $seed
    done
  done
done
