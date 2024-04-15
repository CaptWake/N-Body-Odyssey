#ifndef OCTREE_H_
#define OCTREE_H_

#include <immintrin.h>

#include <algorithm>
#include <cstdint>
#include <limits>
#include <memory>
#include <vector>

#include "vec3.h"

using node_id = std::uint64_t;
static constexpr node_id null = node_id(-1);
static constexpr float inf = std::numeric_limits<float>::infinity();

struct box {
  vec3 min{inf, inf, inf};
  vec3 max{-inf, -inf, -inf};

  // extending a bounding box with a point
  //
  box& operator|=(vec3 const& p) {
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
template <typename Iterator>
box bbox(Iterator begin, Iterator end) {
  box result;
  for (auto it = begin; it != end; ++it) result |= *it;
  return result;
}

struct node {
  float mass{0.f};
  float size{0.f};
  vec3 center{0.f, 0.f, 0.f};
  bool is_leaf = false;
  node_id children[8]{null, null, null, null, null, null, null, null};
};

inline vec3 middle(vec3 const& p1, vec3 const& p2) {
  return {(p1.x() + p2.x()) / 2.f,
          (p1.y() + p2.y()) / 2.f,
          (p1.z() + p2.z()) / 2.f};
}

class octree {
 public:
  octree(std::vector<vec3> points, std::vector<float>& masses);
  vec3 force_at(vec3 const& p, node_id id, float theta);

 private:
  box _bbox;
  node_id root;
  std::vector<node> nodes;
  std::vector<float> masses;
  std::vector<vec3> points;
  template <typename Iterator>
  // [begin, end) is the sequence of points to build the quadtree on.
  node_id build_impl(box const& bbox, Iterator begin, Iterator end);
};

struct boxSOA {
  float minx = inf;
  float miny = inf;
  float minz = inf;
  float maxx = -inf;
  float maxy = -inf;
  float maxz = -inf;

  // extending a bounding box with a point
  //
  void Update(float x, float y, float z) {
    minx = std::min(minx, x);
    miny = std::min(miny, y);
    minz = std::min(minz, z);
    maxx = std::max(maxx, x);
    maxy = std::max(maxy, y);
    maxz = std::max(maxz, z);
  }
};

template <typename Iterator>
boxSOA bbox(Iterator beginx, Iterator endx, Iterator beginy, Iterator beginz) {
  boxSOA result;
  for (auto itx = beginx; itx != endx; ++itx) {
    result.Update(*itx, *beginy, *beginz);
    ++beginy;
    ++beginz;
  }
  return result;
}

struct nodeSOA {
  float mass{0.f};
  float size{0.f};
  float cx;
  float cy;
  float cz;
  bool is_leaf = false;
  node_id children[8]{null, null, null, null, null, null, null, null};
};

class octreeSOA {
 public:
  octreeSOA(float* px, float* py, float* pz, float* masses, uint64_t n_bodies);
  void force_at(float px, float py, float pz, const node_id* nodes, float theta,
                __m256& fx, __m256& fy, __m256& fz);

 private:
  boxSOA _bbox;
  node_id root;
  std::vector<nodeSOA> nodes;
  float *masses, *px, *py, *pz;
  // [begin, end) is the sequence of points to build the quadtree on.
  node_id build_impl(boxSOA const& bbox, uint64_t begin, uint64_t end);
};

#endif
