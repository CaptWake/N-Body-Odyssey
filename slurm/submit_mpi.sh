#!/bin/bash

# Carica l'ambiente Spack e OpenMPI
spack load openmpi
make nbody_mpi_ap
spack unload openmpi

#sbatch submit_multinode.sh

