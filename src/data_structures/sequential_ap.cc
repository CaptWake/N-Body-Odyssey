//
// Created by maste on 3/21/2024.
//
#include "sequential_ap.h"

#include "data_structures/nbody.h"

// PROBABLY WE CAN SET IT TO PRIVATE
// AND EXPOSE ONLY THE START METHOD ??
// INSTEAD OF USING THE SIMULATION CLASS?
void SequentialAP::Update(const double dt) {
  UpdateAcc();
  UpdatePos(dt);
  UpdateVel(dt);
  std::cout << *this << std::endl;
}

std::ostream& operator<<(std::ostream& os, const SequentialAP& nbody) {
  os << "NBody bodies:\n";
  for (const Body& particle : nbody.GetBodies()) {
    os << particle << "\n";
  }
  return os;
}

void SequentialAP::AddBody(Body body) {
  this->bodies.push_back(body);
}
void SequentialAP::CalcAcc(Body& b1, Body& b2) {
  vec3 distance12 = b2.GetPosition() - b1.GetPosition();
  double factor1 = +G * b2.GetMass() / distance12.length() / distance12.length();
  double factor2 = -G * b1.GetMass() / distance12.length() / distance12.length();
  std::cout<<"p1=" << b1.GetPosition()<<std::endl;
  std::cout<<"p2=" << b2.GetPosition()<<std::endl;
  std::cout<<"dist=" << distance12<<std::endl;
  std::cout << "factor = " << G << "*" << b2.GetMass() << "/" << distance12.length() << "/" << distance12.length() << std::endl;
  std::cout << "acc pre " << b1.GetAcceleration() << std::endl;
  std::cout << "acc add " << unit_vector(distance12) * factor1 << std::endl;
  //std::cout << "acc pre " << b2.getAcceleration() << std::endl;
  b1.SetAcceleration(b1.GetAcceleration() + unit_vector(distance12) * factor1);
  b2.SetAcceleration(b2.GetAcceleration() + unit_vector(distance12) * factor2);
  std::cout << "acc post " << b1.GetAcceleration() << std::endl;

}
void SequentialAP::UpdateAcc() {
  for(auto& body : bodies) {
    body.SetAcceleration(vec3());
  }
  for(int i = 0; i < bodies.size(); i++) {
    for(int j = i+1; j < bodies.size(); j++) {
      Body& body1 = bodies[i];
      Body& body2 = bodies[j];
      CalcAcc(body1, body2);
      std::cout << "acc post body " << i << " = " << body1.GetAcceleration() << std::endl;
      std::cout << "acc post body " << j << " = " << body2.GetAcceleration() << std::endl;
    }
  }
}
void SequentialAP::UpdateVel(const double dt) {
  for(auto& b : bodies) {
    b.SetVelocity(b.GetVelocity() + b.GetAcceleration() * dt);
  }
}
void SequentialAP::UpdatePos(const double dt) {
  for(int i = 0; i < bodies.size(); i++) {
    Body& b = bodies[i];
    std::cout << "pos pre "<< i << " " << b.GetPosition() << std::endl;

    std::cout <<b.GetPosition() <<"+" <<b.GetVelocity()<< "*"<< dt<< "+"  <<b.GetAcceleration()<< "*"<< 0.5<< "*"<< dt << "*"<< dt <<std::endl;
    b.SetPosition(b.GetPosition() + b.GetVelocity() * dt +  b.GetAcceleration() * 0.5 * dt * dt);
    std::cout << "pos post "<< i << " " << b.GetPosition() << std::endl;
    std::cout << "pos post length"<< i << " " << b.GetPosition().length() << std::endl;
  }
}
