#ifndef QUADTREE_H_
#define QUADTREE_H_

#include <limits>
#include <algorithm>
#include <memory>
#include <cstdint>
#include <fstream>
#include "vec3.h"

using node_id = std::uint64_t;
using point = vec3;
static constexpr node_id null = node_id(-1);
static constexpr float inf = std::numeric_limits<float>::infinity();
static const float theta = 1.f;

// double BHTreeNode::s_theta = 0.9; -> theta parameter = l/D

struct box {
  point min{inf, inf, inf};
  point max{-inf, -inf, -inf};

  // extending a bounding box with a point
  //
  box & operator |= (point const & p) {
    min[0] = std::min(min.x(), p.x());
    min[1] = std::min(min.y(), p.y());
    min[2] = std::min(min.z(), p.z());
    max[0] = std::max(max.x(), p.x());
    max[1] = std::max(max.y(), p.y());
    max[2] = std::max(max.z(), p.z());
    return *this;
  }
};

// template function that computes the bounding box of a sequence of points
template <typename  Iterator>
box bbox(Iterator begin, Iterator end) {
  box result;
  for (auto it = begin; it != end; ++it)
    result |= *it;
  return result;
}

point middle(point const& p1, point const& p2) {
  return { (p1.x() + p2.x()) / 2.f, (p1.y() + p2.y()) / 2.f, (p1.z() + p2.z()) / 2.f };
}

struct node {
  float mass{0.f};
  float size{0.f};
  vec3 center{0.f, 0.f, 0.f};
  // why not storing children in a unique array like that: children[8] ?
  //node_id children[2][2][2] {{{null, null},{null, null}},
  //    {{null, null},{null, null}}};
  node_id children[8]{null, null, null, null,
      null, null, null, null};
};

struct quadtree {
  box bbox;
  node_id root;
  std::vector<node> nodes;
  std::vector<float> masses;
  std::vector<point> points;
  std::vector<point> velocities;
  //std::vector<std::uint64_t> node_points_begin;
};

template <typename  Iterator>
// [begin, end) is the sequence of points to build the quadtree on.
node_id build_impl(quadtree& tree, box const& bbox, Iterator begin, Iterator end) {
  if (begin == end) return null;

  node_id result = tree.nodes.size();
  tree.nodes.emplace_back();

  tree.nodes[result].size = bbox.max.x() - bbox.min.x();
  if (begin + 1 == end) {
    auto id = begin - tree.points.begin();
    tree.nodes[result].mass = 1.0f;//tree.masses[id];
    tree.nodes[result].center = tree.points[id];
    return result;
  }

  point center = middle(bbox.min, bbox.max);

  auto bottom = [center](point const& p) { return p.y() < center.y(); };
  auto left = [center](point const& p) { return p.x() < center.x(); };
  auto front = [center](point const& p) { return p.z() < center.z(); };

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
  tree.nodes[result].children[0] = build_impl(tree, { {bbox.min.x(), bbox.min.y(), bbox.min.z()}, {center.x(), center.y(), center.z()} }, begin, split_x_front_lower);
  // second quadrant
  tree.nodes[result].children[1] = build_impl(tree, { {center.x(), bbox.min.y(), bbox.min.z()}, {bbox.max.x(), center.y(), center.z()} }, split_x_front_lower, split_y_front);
  // third quadrant
  tree.nodes[result].children[2] = build_impl(tree, { {bbox.min.x(), center.y(), bbox.min.z()}, {center.x(), bbox.max.y(), center.z()} }, split_y_front, split_x_front_upper);
  // fourth quadrant
  tree.nodes[result].children[3] = build_impl(tree, { {center.x(), center.y(), bbox.min.z()}, {bbox.max.x(), bbox.max.y(), center.z()} }, split_x_front_upper, split_z);

  // back slice of the cube
  tree.nodes[result].children[4] = build_impl(tree, { {bbox.min.x(), bbox.min.y(), center.z()}, {center.x(), center.y(), bbox.max.z()} }, split_z, split_x_back_lower);
  tree.nodes[result].children[5] = build_impl(tree, { {center.x(), bbox.min.y(), center.z()}, {bbox.max.x(), bbox.max.y(), bbox.max.z()} }, split_x_back_lower, split_y_back);
  tree.nodes[result].children[6] = build_impl(tree, { {bbox.min.x(), center.y(), center.z()}, {center.x(), bbox.max.y(), bbox.max.z()} }, split_y_back, split_x_back_upper);
  tree.nodes[result].children[7] = build_impl(tree, { {center.x(), center.y(), center.z()}, {bbox.max.x(), bbox.max.y(), bbox.max.z()} },split_x_back_upper, end);

  vec3 sum{0.f, 0.f, 0.f};
  uint64_t i = 0;
  for (auto child_id : tree.nodes[result].children) {

    if (child_id != null) {
      tree.nodes[result].mass += tree.nodes[child_id].mass;
      sum += tree.nodes[child_id].center * tree.nodes[child_id].mass;
      i++;
    }
  }
  tree.nodes[result].center = sum / tree.nodes[result].mass ;
  return result;
}

quadtree build(std::vector<point> points, std::vector<float> masses) {
  quadtree result;
  result.root = build_impl(result, bbox(points.begin(), points.end()), points.begin(), points.end());
  return result;
}

vec3 force_at(quadtree& tree, point const& p, uint64_t id) {
  if(id != 0) {
    auto const &n = tree.nodes[id];

    if (n.center == p)
      return vec3{0.f, 0.f, 0.f};

    vec3 d = n.center - p;

    auto l = d.x() * d.x() + d.y() * d.y() + d.z() * d.z() + 1e-9f;
    auto d_inv = 1.0f / sqrtf(l);
    auto d_inv3 = d_inv * d_inv * d_inv;

    if (d.length() > n.size * 1)
      return n.mass * d * d_inv3;

  }
  vec3 sum{0.f, 0.f, 0.f};
  for (auto child_id: tree.nodes[id].children) {
      if (child_id != null)
        sum += force_at(tree, p, child_id);
  }

  return sum;
};

void update(std::vector<point>& points, std::vector<float>& masses, std::vector<vec3>& velocities, float dt) {
  std::vector<vec3> delta{points.size(), {0,0,0}};
  quadtree tree = build(points, masses);
  for (auto i = 0; i < points.size(); ++i){
    delta[i] = force_at(tree, points[i], 0);
    velocities[i] += delta[i] * dt;
  }

  for (std::size_t i = 0; i < points.size(); ++i) {
    points[i] += velocities[i] * dt;
  }
}

void LoadFromCSVConfiguration(const std::string &filename, std::vector<float>& m, std::vector<vec3>& p, std::vector<vec3>& v) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Cannot find configuration file at: " << filename << std::endl;
    exit(1);
  }

  std::string line;
  uint64_t i = 0;

  // Assuming the header is not present

  std::getline(file, line);
  std::istringstream iss(line);
  std::string token;
  std::getline(iss, token);
  // Get the gravitational constant (first line)
  float G = std::stof(token);

  while (std::getline(file, line)) {
    iss = std::istringstream(line);
    // Read particle data from CSV fields
    float mass, x, y, z, vx, vy, vz;
    // Read and parse comma-separated values
    std::getline(iss, token, ',');
    mass = std::stof(token);

    std::getline(iss, token, ',');
    x = std::stof(token);

    std::getline(iss, token, ',');
    y = std::stof(token);

    std::getline(iss, token, ',');
    z = std::stof(token);

    std::getline(iss, token, ',');
    vx = std::stof(token);

    std::getline(iss, token, ',');
    vy = std::stof(token);

    std::getline(iss, token, ',');
    vz = std::stof(token);

    m.push_back(mass);
    p.emplace_back(x,y,z);
    v.emplace_back(vx, vy, vz);
    i++;
  }
  file.close();
}

#endif
