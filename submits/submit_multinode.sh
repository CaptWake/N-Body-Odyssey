#!/bin/bash

body=200000
seed=1
nodes=50
exe="/beegfs/home/andrea.zito/SCPD-Project/nbody_mpi_ap"

for s in ${seed[@]}; do
	for n in ${nodes[@]}; do
        	sbatch <<EOT
#!/bin/bash

#SBATCH --job-name=nbody<3
#SBATCH --partition=broadwell
#SBATCH --time=01:00:00
#SBATCH --nodes="$n"
#SBATCH --ntasks-per-node=1
#SBATCH -o test_mpi_ap_v3_definitve/test"$s"_"$n".txt

mpirun "$exe" "$body" "$s"
EOT
    done
done
