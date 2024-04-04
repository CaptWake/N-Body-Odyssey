#ifndef QUADTREE_H_
#define QUADTREE_H_

#include <limits>
#include <algorithm>
#include <memory>
#include <cstdint>

using node_id = std::uint64_t;
static constexpr node_id null = node_id(-1);
static constexpr float inf = std::numeric_limits<float>::infinity();

struct point {
  float x, y;
};

struct box {
  point min{inf, inf};
  point max{-inf, -inf};

  // extending a bounding box with a point
  //
  box & operator |= (point const & p) {
    min.x = std::min(min.x, p.x);
    min.y = std::min(min.y, p.y);
    max.x = std::max(max.x, p.x);
    max.y = std::max(max.y, p.y);
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
  return { (p1.x + p2.x) / 2.f, (p1.y + p2.y) / 2.f };
}

struct node
{
  node_id children[2][2]{
      {null, null},
      {null, null}
  };
};

struct quadtree {
  box bbox;
  node_id root;
  std::vector<node> nodes;
  std::vector<point> points;
  std::vector<std::uint64_t> node_points_begin;
};

template <typename  Iterator>
node_id build_impl(quadtree& tree, box const& bbox, Iterator begin, Iterator end) {
  if (begin == end) return null;

  node_id result = tree.nodes.size();
  tree.nodes.emplace_back();

  if (begin + 1 == end) return result;

  point center = middle(bbox.min, bbox.max);

  auto bottom = [center](point const& p) { return p.y < center.y; };
  auto left = [center](point const& p) { return p.x < center.x; };

  // Split the points along Y
  Iterator split_y = std::partition(begin, end, bottom);
  // Now, [begin, split_y) is the lower half,
  // and [split_y, end) is the upper half.
  // Split the lower half along X
  Iterator split_x_lower = std::partition(begin, split_y, left);
  Iterator split_x_upper = std::partition(split_y, end, left);

  tree.nodes[result].children[0][0] = build_impl(tree, {bbox.min, center}, begin, split_x_lower);
  tree.nodes[result].children[0][1] = build_impl(tree, { {center.x, bbox.min.y}, {bbox.max.x, center.y} }, split_x_lower, split_y);
  tree.nodes[result].children[1][0] = build_impl(tree, { {bbox.min.x, center.y}, {center.x, bbox.max.y} }, split_y, split_x_upper);
  tree.nodes[result].children[1][1] = build_impl(tree, {center, bbox.max}, split_x_upper, end);

  return result;
}

quadtree build(std::vector<point> points) {
  quadtree result;
  result.points = std::move(points);
  result.root = build_impl(result, bbox(result.points.begin(), result.points.end()), result.points.begin(), result.points.end());
  result.node_points_begin.push_back(result.points.size());
  return result;
}

#endif
