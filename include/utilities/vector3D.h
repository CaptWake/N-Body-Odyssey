
// Created by dadi on 3/18/24.
//

#ifndef N_BODY_VECTOR3D_H
#define N_BODY_VECTOR3D_H


#include <ostream>
class vector3D {
 private:

 public:
  double x;
  double y;
  double z;

  vector3D(double x, double y, double z);

  vector3D();

  double getX() const;

  void setX(double x);

  double getY() const;

  void setY(double y);

  double getZ() const;

  void setZ(double z);

  friend std::ostream &operator<<(std::ostream &os, const vector3D &d);

  bool operator==(const vector3D &rhs) const;

  bool operator!=(const vector3D &rhs) const;

  double magnitude() const;

  vector3D operator+(const vector3D &rhs) const;

  vector3D operator-(const vector3D &rhs) const;

  double operator*(const vector3D &rhs) const;

  vector3D operator*(double factor) const;

  vector3D operator^(const vector3D &rhs) const;

  vector3D versor() const;

  std::string toCSV() const;
};


#endif //N_BODY_VECTOR3D_H
