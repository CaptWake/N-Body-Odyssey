#ifndef MPI_AP_H_
#define MPI_AP_H_

#include <cstdint>

#include <iostream> // FOR DEBUG PURPOSES, TODO: REMOVE
#include <random>
#include <string>
#include <vector>

#include "fileIO.h"
#include "mpi.h"
#include "nbody.h"


  // generate random samples
  void mpi_ap(int n_bodies, const float grav_const);

#endif
