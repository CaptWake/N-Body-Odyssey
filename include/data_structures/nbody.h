#ifndef NBODY_H
#define NBODY_H

#include "body.h"
#include <vector>

class NBody {
 public:
  virtual void update(double dt) = 0;  // Pure virtual method
};

class SequentialAPNBody : public NBody {
 public:
  SequentialAPNBody() {}
  SequentialAPNBody(std::vector<Body>& bodies, double gravitationalConstant) {
    this->bodies = bodies;
    this->gravitationalConstant = gravitationalConstant;
  }

  // Getters//
  const std::vector<Body>& GetBodies() const { return bodies; }

  // Setters//
  void SetBodies(const std::vector<Body>& bodies) {
    this->bodies = bodies;
  }

  // Update//
  void update(const double dt) override;

  friend std::ostream& operator<<(std::ostream& os, const NBody& nbody);

  void addBody(Body body);

 protected:
  std::vector<Body> bodies;
  double gravitationalConstant;

 private:
  vec3 CalculateGravitationalForce(Body& particle1,
                                   Body& particle2) const;
};

#endif
