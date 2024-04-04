#ifndef QUADTREE_H_
#define QUADTREE_H_

#include <limits>
#include <algorithm>
#include <memory>
#include <cstdint>

using node_id = std::uint64_t;
static constexpr node_id null = node_id(-1);
static constexpr float inf = std::numeric_limits<float>::infinity();
// double BHTreeNode::s_theta = 0.9; -> theta parameter = l/D


struct point {
  float x, y, z;
  bool operator == (point const & p) const  {
    return x == p.x && y == p.y && z == p.z;
  }
};

struct box {
  point min{inf, inf, inf};
  point max{-inf, -inf, -inf};

  // extending a bounding box with a point
  //
  box & operator |= (point const & p) {
    min.x = std::min(min.x, p.x);
    min.y = std::min(min.y, p.y);
    min.z = std::min(min.z, p.z);
    max.x = std::max(max.x, p.x);
    max.y = std::max(max.y, p.y);
    min.z = std::max(max.z, p.z);
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
  return { (p1.x + p2.x) / 2.f, (p1.y + p2.y) / 2.f, (p1.z + p2.z) / 2.f };
}

struct node {
  node_id children[2][2][2] {{{null, null},{null, null}},
      {{null, null},{null, null}}};
};

struct quadtree {
  box bbox;
  node_id root = 0;
  std::vector<node> nodes;
  std::vector<float> masses;
  std::vector<point> points;
  std::vector<std::uint64_t> node_points_begin;
};

template <typename  Iterator>
node_id build_impl(quadtree& tree, box const& bbox, Iterator begin, Iterator end) {
  if (begin == end) return null;

  node_id result = tree.nodes.size();
  tree.nodes.emplace_back();

  if (std::equal(begin + 1, end, begin)) return result;

  if (begin + 1 == end) return result;

  tree.node_points_begin[result] = (begin - tree.points.begin());
  point center = middle(bbox.min, bbox.max);

  auto bottom = [center](point const& p) { return p.y < center.y; };
  auto left = [center](point const& p) { return p.x < center.x; };
  auto front = [center](point const& p) { return p.z < center.z; };

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

  tree.nodes[result].children[0][0][0] = build_impl(tree, { {bbox.min.x, center.y, bbox.min.z}, {bbox.max.x, center.y, bbox.min.z} }, begin, split_x_front_lower);
  tree.nodes[result].children[0][0][1] = build_impl(tree, { {center.x, bbox.min.y, bbox.min.z}, {bbox.max.x, center.y, center.z} }, split_x_front_lower, split_y_front);
  tree.nodes[result].children[0][1][0] = build_impl(tree, { {center.x, bbox.min.y, center.z}, {bbox.max.x, center.y, bbox.max.z} }, split_y_front, split_x_front_upper);
  tree.nodes[result].children[0][1][1] = build_impl(tree, { {center.x, bbox.min.y, bbox.max.z}, {bbox.max.x, center.y, center.z} }, split_x_front_upper, split_z);
  tree.nodes[result].children[1][0][0] = build_impl(tree, { {bbox.min.x, center.y, bbox.min.z}, {center.x, bbox.max.y, center.z} }, split_z, split_x_back_lower);
  tree.nodes[result].children[1][0][1] = build_impl(tree, { {center.x, center.y, bbox.min.z}, {bbox.max.x, bbox.max.y, center.z} }, split_x_back_lower, split_y_back);
  tree.nodes[result].children[1][1][0] = build_impl(tree, { {center.x, center.y, center.z}, {bbox.max.x, bbox.max.y, bbox.max.z} }, split_y_back, split_x_back_upper);
  tree.nodes[result].children[1][1][1] = build_impl(tree, { {bbox.min.x, center.y, center.z}, {center.x, bbox.max.y, bbox.max.z} }, split_x_back_upper, end);

  return result;
}

quadtree build(std::vector<point> points) {
  quadtree result;
  result.points = std::move(points);
  result.root = build_impl(result, bbox(result.points.begin(), result.points.end()), result.points.begin(), result.points.end());
  result.node_points_begin.push_back(result.points.size());
  return result;
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
