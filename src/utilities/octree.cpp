#include "octree.h"

octree::octree(std::vector<vec3> points, std::vector<float>& masses) {
  this->masses = masses;
  this->root = build_impl(bbox(points.begin(), points.end()), points.begin(), points.end());
}

vec3 octree::force_at(const vec3 &p, node_id id, float theta) {
  if(id != 0) {
    auto const &n = this->nodes[id];

    if (n.center == p)
      return vec3{0.f, 0.f, 0.f};

    vec3 d = n.center - p;

    auto l = d.x() * d.x() + d.y() * d.y() + d.z() * d.z() + 1e-9f;
    auto d_inv = 1.0f / sqrtf(l);
    auto d_inv3 = d_inv * d_inv * d_inv;

    if (d.length() > n.size * theta)
      return n.mass * d * d_inv3;

  }
  vec3 total_force{0.f, 0.f, 0.f};
  for (auto child_id: this->nodes[id].children) {
    if (child_id != null)
      total_force += force_at(p, child_id, theta);
  }

  return total_force;
}

template<typename Iterator>
inline node_id octree::build_impl(const box &bbox, Iterator begin, Iterator end) {
  if (begin == end) return null;

  node_id result = this->nodes.size();
  this->nodes.emplace_back();

  // Compute bbox length on building time to avoid recompute multiple times
  // when computing the force contributions
  this->nodes[result].size = bbox.max.x() - bbox.min.x();

  if (begin + 1 == end) {
    auto id = begin - this->points.begin();
    this->nodes[result].mass = 1.0f;//this->masses[id];
    this->nodes[result].center = this->points[id];
    return result;
  }

  vec3 center = middle(bbox.min, bbox.max);

  auto bottom = [center](vec3 const& p) { return p.y() < center.y(); };
  auto left = [center](vec3 const& p) { return p.x() < center.x(); };
  auto front = [center](vec3 const& p) { return p.z() < center.z(); };

  // Split the points along Z
  Iterator split_z = std::partition(begin, end, front);

  // Split the points along Y
  Iterator split_y_front = std::partition(begin, split_z, bottom);
  Iterator split_y_back = std::partition(split_z, end, bottom);

  // Split the points along X
  Iterator split_x_front_lower = std::partition(begin, split_y_front, left);
  Iterator split_x_front_upper = std::partition(split_y_front, split_z, left);
  Iterator split_x_back_lower = std::partition(split_z, split_y_back, left);
  Iterator split_x_back_upper = std::partition(split_y_back, end, left);

  /*
     +--------+
    /       / |
   +---+---+  |
   | 3 | 4 |  |
   +---+---+ /
   | 1 | 2 |/
   +---+---+

   */

  // front slice of the cube
  // first quadrant
  this->nodes[result].children[0] = build_impl({ {bbox.min.x(), bbox.min.y(), bbox.min.z()}, {center.x(), center.y(), center.z()} }, begin, split_x_front_lower);
  // second quadrant
  this->nodes[result].children[1] = build_impl({ {center.x(), bbox.min.y(), bbox.min.z()}, {bbox.max.x(), center.y(), center.z()} }, split_x_front_lower, split_y_front);
  // third quadrant
  this->nodes[result].children[2] = build_impl({ {bbox.min.x(), center.y(), bbox.min.z()}, {center.x(), bbox.max.y(), center.z()} }, split_y_front, split_x_front_upper);
  // fourth quadrant
  this->nodes[result].children[3] = build_impl({ {center.x(), center.y(), bbox.min.z()}, {bbox.max.x(), bbox.max.y(), center.z()} }, split_x_front_upper, split_z);

  // back slice of the cube
  this->nodes[result].children[4] = build_impl({ {bbox.min.x(), bbox.min.y(), center.z()}, {center.x(), center.y(), bbox.max.z()} }, split_z, split_x_back_lower);
  this->nodes[result].children[5] = build_impl({ {center.x(), bbox.min.y(), center.z()}, {bbox.max.x(), bbox.max.y(), bbox.max.z()} }, split_x_back_lower, split_y_back);
  this->nodes[result].children[6] = build_impl({ {bbox.min.x(), center.y(), center.z()}, {center.x(), bbox.max.y(), bbox.max.z()} }, split_y_back, split_x_back_upper);
  this->nodes[result].children[7] = build_impl({ {center.x(), center.y(), center.z()}, {bbox.max.x(), bbox.max.y(), bbox.max.z()} },split_x_back_upper, end);

  vec3 sum{0.f, 0.f, 0.f};
  uint64_t i = 0;
  for (auto child_id : this->nodes[result].children) {

    if (child_id != null) {
      this->nodes[result].mass += this->nodes[child_id].mass;
      sum += this->nodes[child_id].center * this->nodes[child_id].mass;
      i++;
    }
  }
  this->nodes[result].center = sum / this->nodes[result].mass ;
  return result;
}