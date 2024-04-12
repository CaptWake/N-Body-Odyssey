#include "octree.h"

#include <immintrin.h>  // AVX intrinsics

octree::octree(std::vector<vec3> points, std::vector<float>& masses) {
  this->masses = masses;
  this->root = build_impl(
      bbox(points.begin(), points.end()), points.begin(), points.end());
}

vec3 octree::force_at(const vec3& p, node_id id, float theta) {
  if (id != 0) {
    auto const& n = this->nodes[id];

    vec3 d = n.center - p;

    auto l = d.x() * d.x() + d.y() * d.y() + d.z() * d.z() + 1e-9f;
    auto d_inv = 1.0f / sqrtf(l);
    auto d_inv3 = d_inv * d_inv * d_inv;

    // TODO: discuss the following fact: we have to set a bool indicating if we
    // are a leaf or not? compute the force if we are a leaf or the theta is
    // respected
    if (n.children[0] == null && n.children[1] == null &&
            n.children[2] == null && n.children[3] == null &&
            n.children[4] == null && n.children[5] == null &&
            n.children[6] == null && n.children[7] == null ||
        d.length() < n.size * theta)
      return n.mass * d * d_inv3;
  }
  vec3 total_force{0.f, 0.f, 0.f};
  for (auto child_id : this->nodes[id].children) {
    if (child_id != null) total_force += force_at(p, child_id, theta);
  }

  return total_force;
}

template <typename Iterator>
inline node_id octree::build_impl(const box& bbox, Iterator begin,
                                  Iterator end) {
  if (begin == end) return null;

  node_id result = this->nodes.size();
  this->nodes.emplace_back();

  // Compute bbox length on building time to avoid recompute multiple times
  // when computing the force contributions
  this->nodes[result].size = bbox.max.x() - bbox.min.x();

  if (begin + 1 == end) {
    auto id = begin - this->points.begin();
    this->nodes[result].mass = 1.0f;  // this->masses[id];
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
  this->nodes[result].children[0] =
      build_impl({{bbox.min.x(), bbox.min.y(), bbox.min.z()},
                  {center.x(), center.y(), center.z()}},
                 begin,
                 split_x_front_lower);
  // second quadrant
  this->nodes[result].children[1] =
      build_impl({{center.x(), bbox.min.y(), bbox.min.z()},
                  {bbox.max.x(), center.y(), center.z()}},
                 split_x_front_lower,
                 split_y_front);
  // third quadrant
  this->nodes[result].children[2] =
      build_impl({{bbox.min.x(), center.y(), bbox.min.z()},
                  {center.x(), bbox.max.y(), center.z()}},
                 split_y_front,
                 split_x_front_upper);
  // fourth quadrant
  this->nodes[result].children[3] =
      build_impl({{center.x(), center.y(), bbox.min.z()},
                  {bbox.max.x(), bbox.max.y(), center.z()}},
                 split_x_front_upper,
                 split_z);

  // back slice of the cube
  this->nodes[result].children[4] =
      build_impl({{bbox.min.x(), bbox.min.y(), center.z()},
                  {center.x(), center.y(), bbox.max.z()}},
                 split_z,
                 split_x_back_lower);
  this->nodes[result].children[5] =
      build_impl({{center.x(), bbox.min.y(), center.z()},
                  {bbox.max.x(), bbox.max.y(), bbox.max.z()}},
                 split_x_back_lower,
                 split_y_back);
  this->nodes[result].children[6] =
      build_impl({{bbox.min.x(), center.y(), center.z()},
                  {center.x(), bbox.max.y(), bbox.max.z()}},
                 split_y_back,
                 split_x_back_upper);
  this->nodes[result].children[7] =
      build_impl({{center.x(), center.y(), center.z()},
                  {bbox.max.x(), bbox.max.y(), bbox.max.z()}},
                 split_x_back_upper,
                 end);

  vec3 sum{0.f, 0.f, 0.f};
  uint64_t i = 0;
  for (auto child_id : this->nodes[result].children) {
    if (child_id != null) {
      this->nodes[result].mass += this->nodes[child_id].mass;
      sum += this->nodes[child_id].center * this->nodes[child_id].mass;
      i++;
    }
  }
  this->nodes[result].center = sum / this->nodes[result].mass;
  return result;
}

void swap(float* p1, float* p2, uint64_t i, uint64_t j) {
  float aux = p1[i];
  p1[i] = p2[j];
  p2[j] = aux;
}

template <class UnaryPred>
uint64_t partition(float* p1, float* p2, float* p3, uint64_t begin,
                   uint64_t end, UnaryPred p) {
  while (begin != end && p(p1[begin])) {
    begin++;
  }

  if (begin == end) return end;

  for (auto i = begin; i != end; ++i) {
    if (p(p1[i])) {
      std::swap(p1[i], p1[begin]);
      std::swap(p2[i], p2[begin]);
      std::swap(p3[i], p3[begin]);
      begin++;
    }
  }
  return begin;
}

inline node_id octreeSOA::build_impl(const boxSOA& bbox, uint64_t begin,
                                     uint64_t end) {
  if (begin >= end) return null;

  node_id result = this->nodes.size();
  this->nodes.emplace_back();

  // Compute bbox length on building time to avoid recompute multiple times
  // when computing the force contributions
  this->nodes[result].size = bbox.maxx - bbox.minx;

  if (begin + 1 == end) {
    this->nodes[result].mass = 1.0f;  // this->masses[id];
    this->nodes[result].cx = this->px[begin];
    this->nodes[result].cy = this->py[begin];
    this->nodes[result].cz = this->pz[begin];
    return result;
  }

  float cx = (bbox.maxx + bbox.minx) / 2.0f;
  float cy = (bbox.maxy + bbox.miny) / 2.0f;
  float cz = (bbox.maxz + bbox.minz) / 2.0f;

  auto bottom = [cy](float y) { return y < cy; };
  auto left = [cx](float x) { return x < cx; };
  auto front = [cz](float z) { return z < cz; };

  uint64_t split_z =
      partition(pz.data(), px.data(), py.data(), begin, end, front);

  // Split the points along Y
  uint64_t split_y_front =
      partition(py.data(), px.data(), pz.data(), begin, split_z, bottom);
  uint64_t split_y_back =
      partition(py.data(), px.data(), pz.data(), split_z, end, bottom);

  // Split the points along X
  uint64_t split_x_front_lower =
      partition(px.data(), py.data(), pz.data(), begin, split_y_front, left);
  uint64_t split_x_front_upper =
      partition(px.data(), py.data(), pz.data(), split_y_front, split_z, left);
  uint64_t split_x_back_lower =
      partition(px.data(), py.data(), pz.data(), split_z, split_y_back, left);
  uint64_t split_x_back_upper =
      partition(px.data(), py.data(), pz.data(), split_y_back, end, left);

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
  this->nodes[result].children[0] =
      build_impl({bbox.minx, bbox.miny, bbox.minz, cx, cy, cz},
                 begin,
                 split_x_front_lower);
  // second quadrant
  this->nodes[result].children[1] =
      build_impl({cx, bbox.miny, bbox.minz, bbox.maxx, cy, cz},
                 split_x_front_lower,
                 split_y_front);
  // third quadrant
  this->nodes[result].children[2] =
      build_impl({bbox.minx, cy, bbox.minz, cx, bbox.maxy, cz},
                 split_y_front,
                 split_x_front_upper);
  // fourth quadrant
  this->nodes[result].children[3] =
      build_impl({cx, cy, bbox.minz, bbox.maxx, bbox.maxy, cz},
                 split_x_front_upper,
                 split_z);

  // back slice of the cube
  this->nodes[result].children[4] =
      build_impl({bbox.minx, bbox.miny, cz, cx, cy, bbox.maxz},
                 split_z,
                 split_x_back_lower);
  this->nodes[result].children[5] =
      build_impl({cx, bbox.miny, cz, bbox.maxx, bbox.maxy, bbox.maxz},
                 split_x_back_lower,
                 split_y_back);
  this->nodes[result].children[6] =
      build_impl({bbox.minx, cy, cz, cx, bbox.maxy, bbox.maxz},
                 split_y_back,
                 split_x_back_upper);
  this->nodes[result].children[7] = build_impl(
      {cx, cy, cz, bbox.maxx, bbox.maxy, bbox.maxz}, split_x_back_upper, end);

  float sumx, sumy, sumz;
  sumx = sumy = sumz = 0.0f;
  for (auto child_id : this->nodes[result].children) {
    if (child_id != null) {
      this->nodes[result].mass += this->nodes[child_id].mass;
      sumx += this->nodes[child_id].cx * this->nodes[child_id].mass;
      sumy += this->nodes[child_id].cy * this->nodes[child_id].mass;
      sumz += this->nodes[child_id].cz * this->nodes[child_id].mass;
    }
  }
  this->nodes[result].cx = sumx / this->nodes[result].mass;
  this->nodes[result].cy = sumy / this->nodes[result].mass;
  this->nodes[result].cz = sumz / this->nodes[result].mass;
  return result;
}
octreeSOA::octreeSOA(std::vector<float> px, std::vector<float> py,
                     std::vector<float> pz, std::vector<float>& masses) {
  this->masses = masses;
  this->px = px;
  this->py = py;
  this->pz = pz;
  this->root = build_impl(
      bbox(px.begin(), px.end(), py.begin(), pz.begin()), 0, px.size());
}
