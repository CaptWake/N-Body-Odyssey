//
// Created by dadi on 3/18/24.
//

#include <cmath>
#include "vector3D.h"

vector3D::vector3D(double x, double y, double z) : x(x), y(y), z(z) {}

vector3D::vector3D() {
  x=0;
  y=0;
  z=0;
}

double vector3D::getX() const {
  return x;
}

void vector3D::setX(double x) {
  vector3D::x = x;
}

double vector3D::getY() const {
  return y;
}

void vector3D::setY(double y) {
  vector3D::y = y;
}

double vector3D::getZ() const {
  return z;
}

void vector3D::setZ(double z) {
  vector3D::z = z;
}

std::ostream &operator<<(std::ostream &os, const vector3D &d) {
  os << "x: " << d.x << " y: " << d.y << " z: " << d.z;
  return os;
}

bool vector3D::operator==(const vector3D &rhs) const {
  return x == rhs.x &&
         y == rhs.y &&
         z == rhs.z;
}

bool vector3D::operator!=(const vector3D &rhs) const {
  return !(rhs == *this);
}

double vector3D::magnitude() const {
  return sqrt(x*x + y*y + z*z);
}

vector3D vector3D::operator+(const vector3D &rhs) const {
  return {x + rhs.x, y + rhs.y, z + rhs.z};
}

vector3D vector3D::operator-(const vector3D &rhs) const {
  return {x - rhs.x, y - rhs.y, z - rhs.z};
}

double vector3D::operator*(const vector3D &rhs) const {
  return x * rhs.x + y * rhs.y + z * rhs.z;
}

vector3D vector3D::operator^(const vector3D &rhs) const {
  return {y * rhs.z - z * rhs.y, x * rhs.z - z * rhs.x, x * rhs.y - y * rhs.x};
}

vector3D vector3D::operator*(double factor) const {
  return {x * factor, y * factor, z * factor};
}


vector3D vector3D::versor() const {
  return operator*(1/magnitude());
}

std::string vector3D::toCSV() const {
  std::string str;
  str = std::to_string(x) + "," + std::to_string(y) + "," + std::to_string(z);
  return str;
}

