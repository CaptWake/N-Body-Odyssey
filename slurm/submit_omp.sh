#!/bin/bash

#!/bin/bash
#SBATCH --job-name=nbody
#SBATCH --output=test_omp_1.log
#SBATCH --error=test_omp_1.err
#SBATCH --time=06:00:00
#SBATCH --partition=epito
#SBATCH --nodelist=epito01

seed=0

make nbody_sequential_ap
make nbody_omp_ap

threads=(1 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 40 42 44 46 48 50 52 54 56 58 60 62 64 66 68 70 72 74 76 78 80)

ap_seq="/beegfs/home/davide.camino/SCPD-Project/nbody_sequential_ap"
ap_omp="/beegfs/home/davide.camino/SCPD-Project/nbody_omp_ap"

for i in 1 2 3; do
body=65536

echo "sequential"
echo "info run $body 1"
echo "$ap_seq $body $((seed + i))"
$ap_seq $body $((seed + i))

echo "------------------------------"
echo "----         STOP         ----"
echo "------------------------------"

body=65536
echo "strong scalability static $body"
for num_threads in "${threads[@]}"; do
  echo "info run $body $num_threads"
  echo "command: $ap_omp $body static $(($body / num_threads)) $num_threads $((seed + i))"
  $ap_omp $body static $(($body / num_threads)) $num_threads $((seed + i))
done

echo "------------------------------"
echo "----         STOP         ----"
echo "------------------------------"

echo "strong scalability dynamic $body"
for num_threads in "${threads[@]}"; do
  echo "info run $body $num_threads"
  echo "command: $ap_omp $body dynamic $((65536 / num_threads)) $num_threads $((seed + i))"
  $ap_omp $body dynamic $((65536 / num_threads)) $num_threads $((seed + i))
done

echo "------------------------------"
echo "----         STOP         ----"
echo "------------------------------"

body=4096
echo "strong scalability static $body"
for num_threads in "${threads[@]}"; do
  echo "info run $body $num_threads"
  echo "command: $ap_omp $body static $(($body / num_threads)) $num_threads $((seed + i))"
  $ap_omp $body static $(($body / num_threads)) $num_threads $((seed + i))
done

echo "------------------------------"
echo "----         STOP         ----"
echo "------------------------------"

echo "strong scalability dynamic $body"
for num_threads in "${threads[@]}"; do
  echo "info run $body $num_threads"
  echo "command: $ap_omp $body dynamic $(($body / num_threads)) $num_threads $((seed + i))"
  $ap_omp $body dynamic $(($body / num_threads)) $num_threads $((seed + i))
done

echo "------------------------------"
echo "----         STOP         ----"
echo "------------------------------"

body=1024
echo "strong scalability static $body"
for num_threads in "${threads[@]}"; do
  echo "info run $body $num_threads"
  echo "command: $ap_omp $body static $(($body / num_threads)) $num_threads $((seed + i))"
  $ap_omp $body static $(($body / num_threads)) $num_threads $((seed + i))
done

echo "------------------------------"
echo "----         STOP         ----"
echo "------------------------------"

echo "strong scalability dynamic $body"
for num_threads in "${threads[@]}"; do
  echo "info run $body $num_threads"
  echo "command: $ap_omp $body dynamic $(($body / num_threads)) $num_threads $((seed + i))"
  $ap_omp $body dynamic $(($body / num_threads)) $num_threads $((seed + i))
done

echo "------------------------------"
echo "----         STOP         ----"
echo "------------------------------"

body=16384
echo "weak scalability sqrt step static start at $body"
for num_threads in "${threads[@]}"; do
  n_body=$(echo | awk -v t=$num_threads -v b=$body 'BEGIN {print (sqrt(t)*b)}')
  n_body=$(printf "%.0f" "$n_body")
  echo "info run $n_body $num_threads"
  echo "command: $ap_omp $n_body static $((n_body / num_threads)) $num_threads $((seed + i))"
  $ap_omp $n_body static $((n_body / num_threads)) $num_threads $((seed + i))
done

echo "------------------------------"
echo "----         STOP         ----"
echo "------------------------------"

echo "weak scalability sqrt step dynamic start at $body"
for num_threads in "${threads[@]}"; do
  n_body=$(echo | awk -v t=$num_threads -v b=$body 'BEGIN {print (sqrt(t)*b)}')
  n_body=$(printf "%.0f" "$n_body")
  echo "info run $n_body $num_threads"
  echo "command: $ap_omp $n_body dynamic $((n_body / num_threads)) $num_threads $((seed + i))"
  $ap_omp $n_body dynamic $((n_body / num_threads)) $num_threads $((seed + i))
done

echo "------------------------------"
echo "----         STOP         ----"
echo "------------------------------"

body=4096
echo "weak scalability static $body"
for num_threads in "${threads[@]}"; do
  n_body=$(echo | awk -v t=$num_threads -v b=$body 'BEGIN {print (t*b)}')
  echo "info run $n_body $num_threads"
  echo "command: $ap_omp $n_body static $body $num_threads $((seed + i))"
  $ap_omp $n_body static $body $num_threads $((seed + i))
done

echo "------------------------------"
echo "----         STOP         ----"
echo "------------------------------"

echo "weak scalability dynamic $body"
for num_threads in "${threads[@]}"; do
  n_body=$(echo | awk -v t=$num_threads -v b=$body 'BEGIN {print (t*b)}')
  echo "info run $n_body $num_threads"
  echo "command: $ap_omp $n_body dynamic $body $num_threads $((seed + i))"
  $ap_omp $n_body dynamic $body $num_threads $((seed + i))
done

echo "------------------------------"
echo "----         STOP         ----"
echo "------------------------------"

done
