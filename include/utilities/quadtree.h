#ifndef QUADTREE_H_
#define QUADTREE_H_

#include <limits>
#include <algorithm>
#include <memory>
#include <cstdint>
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
  //std::vector<std::uint64_t> node_points_begin;
};

template <typename  Iterator>
// [begin, end) is the sequence of points to build the quadtree on.
node_id build_impl(quadtree& tree, box const& bbox, Iterator begin, Iterator end) {
  if (begin == end) return null;

  node_id result = tree.nodes.size();
  tree.nodes.emplace_back();

  if (begin + 1 == end) {
    auto id = begin - tree.points.begin();
    tree.nodes[result].mass += tree.masses[id];
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
      sum += tree.nodes[child_id].center;
      i++;
    }
  }
  tree.nodes[result].center = sum / i;
  return result;
}

quadtree build(std::vector<point> points, std::vector<float> masses) {
  quadtree result;
  result.points = std::move(points);
  result.masses = std::move(masses);
  result.root = build_impl(result, bbox(result.points.begin(), result.points.end()), result.points.begin(), result.points.end());
  return result;
}

vec3 force_at(quadtree& tree, point const& p, uint64_t id) {
  auto const & n = tree.nodes[id];

  // TODO fix this error
  if (n.center == p)
    return vec3{0.f, 0.f, 0.f};

  vec3 d = p - n.center;

  auto l = d.length();

  if (n.mass == 1.f || l > n.size * theta)
    return n.mass * d / std::pow(std::max(l, 0.001f), 3.f);

  vec3 sum{0.f, 0.f, 0.f};

  for (auto child_id: tree.nodes[id].children) {
      if (child_id != null)
        sum += force_at(tree, p, child_id);
  }

  return sum;
};


// for (int iteration = 0; iteration < iterations; ++iteration) {

//  for (std::size_t i = 0; i < points.size(); ++i)
//    delta[i] += 2.f * potential * (origin - points[i]);

  // Barnes-Hut algorithm
void update(quadtree& tree, float dt) {
  std::vector<vec3> delta;
  std::vector<float> velocities;
  for (auto i = 0; i < tree.points.size(); ++i)
    delta[i] += force_at(tree, tree.points[i], 0);

  for (std::size_t i = 0; i < tree.points.size(); ++i) {
    velocities[i] += delta[i] * dt;
    points[i] += velocities[i] * dt;
  }
}

/*
Function compute mass distribution() is
24: if new particles == 1 then
25: center of mass = particle.position
26: mass = particle.mass
27: else
28: forall child quadrants with particles do
29: quadrant.compute mass distribution
30: mass += quadrant.mass
31: center of mass = quadrant.mass * quadrant.center of mass
32: center of mass /= mass
*/

// BHTreeNode::GetQuadrant(double x, double y) const
// void BHTreeNode::ComputeMassDistribution()

#endif
