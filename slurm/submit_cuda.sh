#!/bin/bash
seed=(1 2 3 4 5)
cuda_ap="/beegfs/home/andrea.zito/SCPD-Project/nbody_cuda_ap_double"
bodies=(65536)

for s in ${seed[@]}; do
  sbatch <<EOT
#!/bin/bash

#SBATCH --job-name=nbody<3
#SBATCH --partition=gracehopper
#SBATCH --output=test_cuda_double_precision/test_cuda"$s".txt
#SBATCH --nodes=1
#SBATCH --gres=gpu:1
#SBATCH --time=04:00:00

bodies=(${bodies[@]})

for body in \${bodies[@]}; do
    ${cuda_ap} \${body} ${s} gh200
done

EOT
done
