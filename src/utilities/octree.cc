#include "octree.h"
#include "avx.h"

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

    if (n.is_leaf || sqrtf(l - 1e-9f) < n.size * theta)
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
    this->nodes[result].is_leaf = true;
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

// Function to apply a predicate operation on AVX register containing distances
__m256 predicate_operation(__m256 distances, __m256 threshold) {
  // Perform comparison operation
  __m256 comparison = _mm256_cmp_ps(distances, threshold, _CMP_LT_OQ);

  // Create a floating-point mask where true elements are 1.0 and false elements are 0.0
  __m256 mask = _mm256_blendv_ps(_mm256_setzero_ps(), _mm256_set1_ps(1.0f), comparison);

  return mask;
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
    this->nodes[result].is_leaf = true;
    return result;
  }

  float cx = (bbox.maxx + bbox.minx) / 2.0f;
  float cy = (bbox.maxy + bbox.miny) / 2.0f;
  float cz = (bbox.maxz + bbox.minz) / 2.0f;

  auto bottom = [cy](float y) { return y < cy; };
  auto left = [cx](float x) { return x < cx; };
  auto front = [cz](float z) { return z < cz; };

  uint64_t split_z =
      partition(pz, px, py, begin, end, front);

  // Split the points along Y
  uint64_t split_y_front =
      partition(py, px, pz, begin, split_z, bottom);
  uint64_t split_y_back =
      partition(py, px, pz, split_z, end, bottom);

  // Split the points along X
  uint64_t split_x_front_lower =
      partition(px, py, pz, begin, split_y_front, left);
  uint64_t split_x_front_upper =
      partition(px, py, pz, split_y_front, split_z, left);
  uint64_t split_x_back_lower =
      partition(px, py, pz, split_z, split_y_back, left);
  uint64_t split_x_back_upper =
      partition(px, py, pz, split_y_back, end, left);

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

  auto children = this->nodes[result].children;

  auto tx = [&](auto cid) { return (cid != null) ? this->nodes[cid].cx : 0; };
  auto ty = [&](auto cid) { return (cid != null) ? this->nodes[cid].cy : 0; };
  auto tz = [&](auto cid) { return (cid != null) ? this->nodes[cid].cz : 0; };
  auto tm = [&](auto cid) { return (cid != null) ? this->nodes[cid].mass : 0; };

  auto m256_x = _mm256_set_ps(tx(children[0]),
                 tx(children[1]),
                 tx(children[2]),
                 tx(children[3]),
                 tx(children[4]),
                 tx(children[5]),
                 tx(children[6]),
                 tx(children[7]));

  auto m256_y =  _mm256_set_ps(ty(children[0]),
                   ty(children[1]),
                   ty(children[2]),
                   ty(children[3]),
                   ty(children[4]),
                   ty(children[5]),
                   ty(children[6]),
                   ty(children[7]));

  auto m256_z = _mm256_set_ps(tz(children[0]),
                   tz(children[1]),
                   tz(children[2]),
                   tz(children[3]),
                   tz(children[4]),
                   tz(children[5]),
                   tz(children[6]),
                   tz(children[7]));


  auto m256_m = _mm256_set_ps(tm(children[0]),
                              tm(children[1]),
                              tm(children[2]),
                              tm(children[3]),
                              tm(children[4]),
                              tm(children[5]),
                              tm(children[6]),
                              tm(children[7]));

  this->nodes[result].mass = hsum_avx(m256_m);
  this->nodes[result].cx = hsum_avx(_mm256_mul_ps(m256_x, m256_m)) / this->nodes[result].mass;
  this->nodes[result].cy = hsum_avx(_mm256_mul_ps(m256_y, m256_m)) / this->nodes[result].mass;
  this->nodes[result].cz = hsum_avx(_mm256_mul_ps(m256_z, m256_m)) / this->nodes[result].mass;

  return result;
}
octreeSOA::octreeSOA(float *px, float *py, float *pz,
                     float *masses, uint64_t n_bodies) {
  this->masses = masses;
  this->px = (float *)malloc(n_bodies * sizeof(float));
  this->py = (float *)malloc(n_bodies * sizeof(float));
  this->pz = (float *)malloc(n_bodies * sizeof(float));

  for (uint64_t i = 0; i < n_bodies; ++i) {
    this->px[i] = px[i];
    this->py[i] = py[i];
    this->pz[i] = pz[i];
  }

  this->root = build_impl(
      bbox(this->px, this->px+n_bodies, this->py, this->pz), 0, n_bodies);
}

void octreeSOA::force_at(const float px,
                         const float py,
                         const float pz,
                         const node_id nodes[8],
                         float theta,
                         __m256& fx,
                         __m256& fy,
                         __m256& fz) {

  auto tx = [&](auto cid) { return (cid != null) ? this->nodes[cid].cx : px; };
  auto ty = [&](auto cid) { return (cid != null) ? this->nodes[cid].cy : py; };
  auto tz = [&](auto cid) { return (cid != null) ? this->nodes[cid].cz : pz; };
  auto ts = [&](auto cid) { return (cid != null) ? this->nodes[cid].size : 0; };
  auto tm = [&](auto cid) { return (cid != null) ? this->nodes[cid].mass : 0; };
  auto tl = [&](auto cid) { return (cid != null) ? this->nodes[cid].is_leaf : 0; };

  const float softening = 1e-9f;
  auto m256_s = _mm256_broadcast_ss(&softening);
  auto m256_sx = _mm256_broadcast_ss(&px);
  auto m256_sy = _mm256_broadcast_ss(&py);
  auto m256_sz = _mm256_broadcast_ss(&pz);

  auto m256_dx = _mm256_set_ps(tx(nodes[7]),
                               tx(nodes[6]),
                               tx(nodes[5]),
                               tx(nodes[4]),
                               tx(nodes[3]),
                               tx(nodes[2]),
                               tx(nodes[1]),
                               tx(nodes[0]));

  auto m256_dy = _mm256_set_ps(ty(nodes[7]),
                               ty(nodes[6]),
                               ty(nodes[5]),
                               ty(nodes[4]),
                               ty(nodes[3]),
                               ty(nodes[2]),
                               ty(nodes[1]),
                               ty(nodes[0]));

  auto m256_dz = _mm256_set_ps(tz(nodes[7]),
                               tz(nodes[6]),
                               tz(nodes[5]),
                               tz(nodes[4]),
                               tz(nodes[3]),
                               tz(nodes[2]),
                               tz(nodes[1]),
                               tz(nodes[0]));

  auto dx = _mm256_sub_ps(m256_dx, m256_sx);
  auto dy = _mm256_sub_ps(m256_dy, m256_sy);
  auto dz = _mm256_sub_ps(m256_dz, m256_sz);

  const __m256 D = _mm256_fmadd_ps(dx, dx, _mm256_fmadd_ps(dy, dy, _mm256_fmadd_ps(dz, dz, m256_s)));
  const __m256 D_inv = _mm256_rsqrt_ps(D);
  const __m256 D_inv3 = _mm256_mul_ps(D_inv,_mm256_mul_ps(D_inv, D_inv));


  const __m256 D_norm = _mm256_sqrt_ps(D);

  auto m256_m = _mm256_set_ps(tm(nodes[7]),
                              tm(nodes[6]),
                              tm(nodes[5]),
                              tm(nodes[4]),
                              tm(nodes[3]),
                              tm(nodes[2]),
                              tm(nodes[1]),
                              tm(nodes[0]));

  auto m256_size = _mm256_set_ps(ts(nodes[7]),
                                 ts(nodes[6]),
                                 ts(nodes[5]),
                                 ts(nodes[4]),
                                 ts(nodes[3]),
                                 ts(nodes[2]),
                                 ts(nodes[1]),
                                 ts(nodes[0]));

  auto m256_theta = _mm256_set1_ps(theta);
  auto m256_threshold = _mm256_mul_ps(m256_size, m256_theta);
  auto mask = predicate_operation(D_norm, m256_threshold);
  auto mask2 = _mm256_set_ps(tl(nodes[7]),
                             tl(nodes[6]),
                             tl(nodes[5]),
                             tl(nodes[4]),
                             tl(nodes[3]),
                             tl(nodes[2]),
                             tl(nodes[1]),
                             tl(nodes[0]));

  auto final_mask = _mm256_or_ps(mask, mask2);
  fx = _mm256_fmadd_ps(D_inv3, _mm256_mul_ps(m256_m, _mm256_mul_ps(final_mask, dx)), fx);
  fy = _mm256_fmadd_ps(D_inv3, _mm256_mul_ps(m256_m, _mm256_mul_ps(final_mask, dy)), fy);
  fz = _mm256_fmadd_ps(D_inv3, _mm256_mul_ps(m256_m, _mm256_mul_ps(final_mask, dz)), fz);

  for (uint_fast8_t i = 0; i < 8; ++i) {
    auto c = nodes[i];
    if (c != null && !this->nodes[c].is_leaf && D_norm[i] >= this->nodes[c].size * theta) {
      force_at(px, py, pz, this->nodes[c].children, theta, fx, fy, fz);
    }
  }
}

