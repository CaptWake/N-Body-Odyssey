#ifndef NBODY_H
#define NBODY_H

#include "body.h"
#include <vector>

class NBody {
 public:
  virtual void Update(double dt) = 0;  // Pure virtual method

};

#endif
