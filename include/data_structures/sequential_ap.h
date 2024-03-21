//
// Created by maste on 3/21/2024.
//

#ifndef SEQUENTIALAP_H_
#define SEQUENTIALAP_H_

#include "nbody.h"
class SequentialAP : public NBody {
  public:
   SequentialAP() { G = 1;}
   SequentialAP(std::vector<Body>& bodies, double gravitationalConstant) {
     this->bodies = bodies;
     this->G = gravitationalConstant;
   }

   // Getters//
   const std::vector<Body>& GetBodies() const { return bodies; }

   // Setters//
   void SetBodies(const std::vector<Body>& in_bodies) {
     bodies = in_bodies;
   }

   // Update//
   void Update(double dt) override;

   friend std::ostream& operator<<(std::ostream& os, const SequentialAP& nbody);

   void AddBody(Body body);


  protected:
    std::vector<Body> bodies;
    double G;

  private:
    void CalcAcc(Body& b1, Body& b2);
    void UpdateAcc();
    void UpdateVel(double dt);
    void UpdatePos(double dt);

};

#endif  // SCPD_PROJECT_INCLUDE_DATA_STRUCTURES_SEQUENTIALAP_H_
