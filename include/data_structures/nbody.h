#ifndef NBODY_H
#define NBODY_H

class NBody {
 public:
  virtual void Update(float dt) = 0;  // Pure virtual method
};

#endif
