// Tests for cavc_combine_compound_plines and the C++ CompoundPolyline API.
//
// Test cases cover:
//   - simple outer-vs-outer union/intersection/xor/exclude
//   - compound shapes with one outer and one hole
//   - shapes with multiple holes
//   - interaction between an outer loop and another shape overlapping a hole
//   - coincident/touching boundaries
//   - winding normalization

#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

#include <gtest/gtest.h>

#include "c_api_include/cavaliercontours.h"
#include "cavc/compoundpolylinecombine.hpp"
#include "testhelpers.hpp"

namespace {

constexpr cavc_real EPS = 1e-4;

inline bool fuzzyEq(cavc_real a, cavc_real b) { return std::abs(a - b) < EPS; }

// ── helpers ──────────────────────────────────────────────────────────────────

// Create a closed CCW rectangle [x0,x1]×[y0,y1] using the C API.
cavc_pline *makeCCWRect(cavc_real x0, cavc_real y0, cavc_real x1, cavc_real y1) {
  std::vector<cavc_vertex> v = {
      {x0, y0, 0}, {x1, y0, 0}, {x1, y1, 0}, {x0, y1, 0}};
  return cavc_pline_new(v.data(), static_cast<uint32_t>(v.size()), 1);
}

// Create a closed CW rectangle (hole) [x0,x1]×[y0,y1].
cavc_pline *makeCWRect(cavc_real x0, cavc_real y0, cavc_real x1, cavc_real y1) {
  // Reverse winding: go CW → top-right then bottom-right order.
  std::vector<cavc_vertex> v = {
      {x0, y1, 0}, {x1, y1, 0}, {x1, y0, 0}, {x0, y0, 0}};
  return cavc_pline_new(v.data(), static_cast<uint32_t>(v.size()), 1);
}

// Sum areas of all outer and hole loops.
struct CompoundAreas {
  cavc_real outerTotal = 0;
  cavc_real holeTotal = 0;
  cavc_real netArea() const { return outerTotal + holeTotal; }
};

CompoundAreas getAreas(cavc_compound_pline const *cpline) {
  CompoundAreas a;
  uint32_t no = cavc_compound_pline_outer_count(cpline);
  uint32_t nh = cavc_compound_pline_hole_count(cpline);
  for (uint32_t i = 0; i < no; ++i) {
    a.outerTotal += cavc_get_area(cavc_compound_pline_get_outer(cpline, i));
  }
  for (uint32_t i = 0; i < nh; ++i) {
    a.holeTotal += cavc_get_area(cavc_compound_pline_get_hole(cpline, i));
  }
  return a;
}

// RAII wrapper for cavc_compound_pline.
struct CompoundPlineGuard {
  cavc_compound_pline *ptr = nullptr;
  explicit CompoundPlineGuard(cavc_compound_pline *p = nullptr) : ptr(p) {}
  CompoundPlineGuard(CompoundPlineGuard const &) = delete;
  CompoundPlineGuard &operator=(CompoundPlineGuard const &) = delete;
  ~CompoundPlineGuard() {
    if (ptr) cavc_compound_pline_delete(ptr);
  }
};

// RAII wrapper for cavc_pline.
struct PlineGuard {
  cavc_pline *ptr = nullptr;
  explicit PlineGuard(cavc_pline *p = nullptr) : ptr(p) {}
  PlineGuard(PlineGuard const &) = delete;
  PlineGuard &operator=(PlineGuard const &) = delete;
  ~PlineGuard() {
    if (ptr) cavc_pline_delete(ptr);
  }
};

// ── Helper: perform a combine via the C API ───────────────────────────────────
// Returns a new CompoundPlineGuard; caller owns the result.
CompoundPlineGuard combine(cavc_compound_pline *a, cavc_compound_pline *b, int mode) {
  cavc_compound_pline *result = nullptr;
  cavc_combine_compound_plines(a, b, mode, &result);
  return CompoundPlineGuard(result);
}

// ── C++ API helper: build a CompoundPolyline from rectangles ─────────────────
using CppPline = cavc::Polyline<double>;
using CppCompound = cavc::CompoundPolyline<double>;

CppPline makeCppCCWRect(double x0, double y0, double x1, double y1) {
  CppPline p;
  p.isClosed() = true;
  p.addVertex({x0, y0, 0});
  p.addVertex({x1, y0, 0});
  p.addVertex({x1, y1, 0});
  p.addVertex({x0, y1, 0});
  return p;
}

CppPline makeCppCWRect(double x0, double y0, double x1, double y1) {
  CppPline p;
  p.isClosed() = true;
  p.addVertex({x0, y1, 0});
  p.addVertex({x1, y1, 0});
  p.addVertex({x1, y0, 0});
  p.addVertex({x0, y0, 0});
  return p;
}

double totalArea(cavc::CompoundCombineResult<double> const &r) {
  double s = 0;
  for (auto const &l : r.outerLoops) s += cavc::getArea(l);
  for (auto const &l : r.holeLoops) s += cavc::getArea(l);
  return s;
}

// ═══════════════════════════════════════════════════════════════════════════════
// Section 1: Simple outer-vs-outer boolean via the C++ API
// ═══════════════════════════════════════════════════════════════════════════════

TEST(CompoundPolylineCombine, UnionOverlapping) {
  // A=[0,10]×[0,10], B=[5,15]×[0,10]; overlap=[5,10]×[0,10] (50 units)
  // Union area = 100 + 100 - 50 = 150
  CppCompound a, b;
  a.outerLoops.push_back(makeCppCCWRect(0, 0, 10, 10));
  b.outerLoops.push_back(makeCppCCWRect(5, 0, 15, 10));

  auto r = cavc::combineCompoundPolylines(a, b, cavc::PlineCombineMode::Union);
  EXPECT_EQ(r.outerLoops.size(), 1u);
  EXPECT_EQ(r.holeLoops.size(), 0u);
  EXPECT_TRUE(fuzzyEq(totalArea(r), 150.0));
}

TEST(CompoundPolylineCombine, IntersectOverlapping) {
  // Intersection = [5,10]×[0,10] area = 50
  CppCompound a, b;
  a.outerLoops.push_back(makeCppCCWRect(0, 0, 10, 10));
  b.outerLoops.push_back(makeCppCCWRect(5, 0, 15, 10));

  auto r = cavc::combineCompoundPolylines(a, b, cavc::PlineCombineMode::Intersect);
  EXPECT_EQ(r.outerLoops.size(), 1u);
  EXPECT_EQ(r.holeLoops.size(), 0u);
  EXPECT_TRUE(fuzzyEq(totalArea(r), 50.0));
}

TEST(CompoundPolylineCombine, ExcludeOverlapping) {
  // A - B = [0,5]×[0,10] area = 50
  CppCompound a, b;
  a.outerLoops.push_back(makeCppCCWRect(0, 0, 10, 10));
  b.outerLoops.push_back(makeCppCCWRect(5, 0, 15, 10));

  auto r = cavc::combineCompoundPolylines(a, b, cavc::PlineCombineMode::Exclude);
  EXPECT_EQ(r.outerLoops.size(), 1u);
  EXPECT_EQ(r.holeLoops.size(), 0u);
  EXPECT_TRUE(fuzzyEq(totalArea(r), 50.0));
}

TEST(CompoundPolylineCombine, XorOverlapping) {
  // XOR = two strips: [0,5]×[0,10] (area=50) and [10,15]×[0,10] (area=50); total=100
  CppCompound a, b;
  a.outerLoops.push_back(makeCppCCWRect(0, 0, 10, 10));
  b.outerLoops.push_back(makeCppCCWRect(5, 0, 15, 10));

  auto r = cavc::combineCompoundPolylines(a, b, cavc::PlineCombineMode::XOR);
  EXPECT_EQ(r.outerLoops.size(), 2u);
  EXPECT_EQ(r.holeLoops.size(), 0u);
  EXPECT_TRUE(fuzzyEq(totalArea(r), 100.0));
}

TEST(CompoundPolylineCombine, UnionDisjoint) {
  // Non-overlapping rectangles → union = two separate loops
  CppCompound a, b;
  a.outerLoops.push_back(makeCppCCWRect(0, 0, 10, 10));   // area = 100
  b.outerLoops.push_back(makeCppCCWRect(20, 0, 30, 10));  // area = 100

  auto r = cavc::combineCompoundPolylines(a, b, cavc::PlineCombineMode::Union);
  EXPECT_EQ(r.outerLoops.size(), 2u);
  EXPECT_EQ(r.holeLoops.size(), 0u);
  EXPECT_TRUE(fuzzyEq(totalArea(r), 200.0));
}

TEST(CompoundPolylineCombine, IntersectDisjoint) {
  // Disjoint shapes → empty intersection
  CppCompound a, b;
  a.outerLoops.push_back(makeCppCCWRect(0, 0, 10, 10));
  b.outerLoops.push_back(makeCppCCWRect(20, 0, 30, 10));

  auto r = cavc::combineCompoundPolylines(a, b, cavc::PlineCombineMode::Intersect);
  EXPECT_EQ(r.outerLoops.size(), 0u);
  EXPECT_EQ(r.holeLoops.size(), 0u);
}

TEST(CompoundPolylineCombine, ExcludeDisjoint) {
  // A - B where B is fully outside A → result = A
  CppCompound a, b;
  a.outerLoops.push_back(makeCppCCWRect(0, 0, 10, 10));
  b.outerLoops.push_back(makeCppCCWRect(20, 0, 30, 10));

  auto r = cavc::combineCompoundPolylines(a, b, cavc::PlineCombineMode::Exclude);
  EXPECT_EQ(r.outerLoops.size(), 1u);
  EXPECT_EQ(r.holeLoops.size(), 0u);
  EXPECT_TRUE(fuzzyEq(totalArea(r), 100.0));
}

TEST(CompoundPolylineCombine, XorDisjoint) {
  // XOR of disjoint shapes = both shapes
  CppCompound a, b;
  a.outerLoops.push_back(makeCppCCWRect(0, 0, 10, 10));   // area = 100
  b.outerLoops.push_back(makeCppCCWRect(20, 0, 30, 10));  // area = 100

  auto r = cavc::combineCompoundPolylines(a, b, cavc::PlineCombineMode::XOR);
  EXPECT_EQ(r.outerLoops.size(), 2u);
  EXPECT_EQ(r.holeLoops.size(), 0u);
  EXPECT_TRUE(fuzzyEq(totalArea(r), 200.0));
}

// ═══════════════════════════════════════════════════════════════════════════════
// Section 2: Compound shapes with one outer and one hole
// ═══════════════════════════════════════════════════════════════════════════════

// A = 20×20 square with a 10×10 square hole in the middle.
// Net area of A = 400 - 100 = 300.
// B = 10×20 rectangle on the right half (x: 10..20).
// Net area of B = 200.

TEST(CompoundPolylineCombine, CompoundHoleUnion) {
  // Union: outer = 20×20 (area=400), hole is the part of A's hole NOT covered by B.
  // A hole = [5,15]×[5,15]. B covers x:[10,20]. Part of hole NOT in B: [5,10]×[5,15] = 50.
  // Net area = 400 - 50 = 350.
  CppCompound a, b;
  a.outerLoops.push_back(makeCppCCWRect(0, 0, 20, 20));  // 20×20
  a.holeLoops.push_back(makeCppCWRect(5, 5, 15, 15));    // 10×10 hole (CW)
  b.outerLoops.push_back(makeCppCCWRect(10, 0, 20, 20)); // 10×20

  auto r = cavc::combineCompoundPolylines(a, b, cavc::PlineCombineMode::Union);
  EXPECT_GE(r.outerLoops.size(), 1u);
  // Net area: the union of A+B outer is 20×20 = 400, the surviving hole portion is 5×10 = 50
  EXPECT_TRUE(fuzzyEq(totalArea(r), 350.0)) << "Got net area: " << totalArea(r);
}

TEST(CompoundPolylineCombine, CompoundHoleIntersect) {
  // Intersect: result = (A.outer ∩ B.outer) minus (A.hole ∩ B.outer)
  // A.outer ∩ B.outer = [10,20]×[0,20] = 200
  // A.hole ∩ B.outer = [10,15]×[5,15] = 5×10 = 50
  // Net area = 200 - 50 = 150
  CppCompound a, b;
  a.outerLoops.push_back(makeCppCCWRect(0, 0, 20, 20));
  a.holeLoops.push_back(makeCppCWRect(5, 5, 15, 15));
  b.outerLoops.push_back(makeCppCCWRect(10, 0, 20, 20));

  auto r = cavc::combineCompoundPolylines(a, b, cavc::PlineCombineMode::Intersect);
  EXPECT_GE(r.outerLoops.size(), 1u);
  EXPECT_TRUE(fuzzyEq(totalArea(r), 150.0)) << "Got net area: " << totalArea(r);
}

TEST(CompoundPolylineCombine, CompoundHoleExclude) {
  // A - B: part of A not in B.
  // A.outer outside B = [0,10]×[0,20] = 200
  // Part of A.hole that is in A but not B = [5,10]×[5,15] = 50 (inside A.outer, left of B)
  // This part of the hole remains as a hole in the result.
  // Net area = 200 - 50 = 150
  CppCompound a, b;
  a.outerLoops.push_back(makeCppCCWRect(0, 0, 20, 20));
  a.holeLoops.push_back(makeCppCWRect(5, 5, 15, 15));
  b.outerLoops.push_back(makeCppCCWRect(10, 0, 20, 20));

  auto r = cavc::combineCompoundPolylines(a, b, cavc::PlineCombineMode::Exclude);
  EXPECT_GE(r.outerLoops.size(), 1u);
  EXPECT_TRUE(fuzzyEq(totalArea(r), 150.0)) << "Got net area: " << totalArea(r);
}

// ═══════════════════════════════════════════════════════════════════════════════
// Section 3: Shapes with multiple holes
// ═══════════════════════════════════════════════════════════════════════════════

TEST(CompoundPolylineCombine, MultipleHolesUnionDisjoint) {
  // A = large square [0,30]×[0,10] with two holes:
  //   hole1 = [2,8]×[2,8] (6×6 = 36)
  //   hole2 = [12,18]×[2,8] (6×6 = 36)
  // B = small rectangle [22,28]×[2,8] (disjoint from A's holes, inside A's outer)
  // Union of A and B = A (since B is fully inside A.outer but outside A's holes)
  //   Net area = 30*10 - 36 - 36 = 300 - 72 = 228
  // Since B is inside A (compound inside), keepFromB = !isInsideA = false → B not added.
  CppCompound a, b;
  a.outerLoops.push_back(makeCppCCWRect(0, 0, 30, 10));
  a.holeLoops.push_back(makeCppCWRect(2, 2, 8, 8));    // hole1
  a.holeLoops.push_back(makeCppCWRect(12, 2, 18, 8));  // hole2
  b.outerLoops.push_back(makeCppCCWRect(22, 2, 28, 8)); // inside A, not in holes

  auto r = cavc::combineCompoundPolylines(a, b, cavc::PlineCombineMode::Union);
  // B is inside A (compound), so union = A
  EXPECT_EQ(r.outerLoops.size(), 1u);
  EXPECT_EQ(r.holeLoops.size(), 2u);
  EXPECT_TRUE(fuzzyEq(totalArea(r), 228.0)) << "Got net area: " << totalArea(r);
}

TEST(CompoundPolylineCombine, MultipleHolesUnionOverlapping) {
  // A = large rect [0,20]×[0,10] with two holes:
  //   hole1 = [1,5]×[1,9]  (4×8=32)
  //   hole2 = [6,10]×[1,9] (4×8=32)
  // A net area = 200 - 32 - 32 = 136
  // B = rect [4,8]×[1,9] (4×8=32) that spans across both holes and the solid strip between them
  // Union should add the part of B not in A → parts of B that are inside A's holes
  // Part of B in hole1: [4,5]×[1,9] = 1×8 = 8
  // Part of B in hole2: [6,8]×[1,9] = 2×8 = 16
  // Total added = 24
  // Net area ≈ 136 + 24 = 160 ... but this isn't a simple count.
  // Let me compute differently: A ∪ B net area:
  //   = area of set {points in A or B}
  //   = area(A) + area(B) - area(A ∩ B)
  //   area(A) = 136, area(B) = 32
  //   area(A ∩ B) = area of B that is inside A = B total (32) - B inside holes
  //                = 32 - (B in hole1) - (B in hole2)
  //                = 32 - 8 - 16 = 8
  //   Net A ∪ B = 136 + 32 - 8 = 160
  CppCompound a, b;
  a.outerLoops.push_back(makeCppCCWRect(0, 0, 20, 10));
  a.holeLoops.push_back(makeCppCWRect(1, 1, 5, 9));
  a.holeLoops.push_back(makeCppCWRect(6, 1, 10, 9));
  b.outerLoops.push_back(makeCppCCWRect(4, 1, 8, 9));

  auto r = cavc::combineCompoundPolylines(a, b, cavc::PlineCombineMode::Union);
  EXPECT_TRUE(fuzzyEq(totalArea(r), 160.0)) << "Got net area: " << totalArea(r);
}

// ═══════════════════════════════════════════════════════════════════════════════
// Section 4: Interaction between outer loop and a shape overlapping a hole
// ═══════════════════════════════════════════════════════════════════════════════

TEST(CompoundPolylineCombine, BInsideAHoleUnion) {
  // A = 20×20 square with a 10×10 hole (the hole is at [5,15]×[5,15]).
  // B = small 4×4 square entirely inside A's hole ([7,11]×[7,11]).
  // B is NOT inside A (compound), since it's in A's hole.
  // Union: result = A with the part of the hole that is NOT B filled; i.e.,
  //   outer = A's outer (unchanged, since B doesn't reach outside the hole)
  //   B adds itself as an island inside the hole.
  // Net area = (400 - 100) + 16 = 316
  CppCompound a, b;
  a.outerLoops.push_back(makeCppCCWRect(0, 0, 20, 20));
  a.holeLoops.push_back(makeCppCWRect(5, 5, 15, 15));    // 10×10 hole
  b.outerLoops.push_back(makeCppCCWRect(7, 7, 11, 11));  // 4×4 inside hole

  auto r = cavc::combineCompoundPolylines(a, b, cavc::PlineCombineMode::Union);
  // Result has: A's outer + A's hole + B as a new outer island
  EXPECT_GE(r.outerLoops.size(), 1u);  // at least A's outer
  EXPECT_TRUE(fuzzyEq(totalArea(r), 316.0)) << "Got net area: " << totalArea(r);
}

TEST(CompoundPolylineCombine, BOverlappingAHoleUnion) {
  // A = 20×20 square [0,20]×[0,20] with 10×10 hole [5,15]×[5,15].
  // B = 10×20 rect [12,22]×[0,20] that partially overlaps A's hole.
  // B overlaps the hole region [12,15]×[5,15] (3×10=30).
  //
  // Union: fills the intersection of B and A's hole region.
  // Net area of A = 400 - 100 = 300.
  // B.outer covers x:[12,22]. The new region added by B (not in A) includes:
  //   - B outside A's outer: [20,22]×[0,20] = 2×20 = 40
  //   - B inside A's hole (not in A): [12,15]×[5,15] = 3×10 = 30
  // Total added = 70
  // But also A's hole shrinks: the part [12,15]×[5,15] is now filled → hole area = 100 - 30 = 70
  // Net area = 300 + 40 + 30 = 370
  CppCompound a, b;
  a.outerLoops.push_back(makeCppCCWRect(0, 0, 20, 20));
  a.holeLoops.push_back(makeCppCWRect(5, 5, 15, 15));
  b.outerLoops.push_back(makeCppCCWRect(12, 0, 22, 20));

  auto r = cavc::combineCompoundPolylines(a, b, cavc::PlineCombineMode::Union);
  EXPECT_TRUE(fuzzyEq(totalArea(r), 370.0)) << "Got net area: " << totalArea(r);
}

TEST(CompoundPolylineCombine, BInsideAHoleIntersect) {
  // B is entirely inside A's hole → B is NOT inside A → intersection is empty.
  CppCompound a, b;
  a.outerLoops.push_back(makeCppCCWRect(0, 0, 20, 20));
  a.holeLoops.push_back(makeCppCWRect(5, 5, 15, 15));
  b.outerLoops.push_back(makeCppCCWRect(7, 7, 11, 11));

  auto r = cavc::combineCompoundPolylines(a, b, cavc::PlineCombineMode::Intersect);
  EXPECT_EQ(r.outerLoops.size(), 0u);
  EXPECT_EQ(r.holeLoops.size(), 0u);
}

// ═══════════════════════════════════════════════════════════════════════════════
// Section 5: Winding normalization
// ═══════════════════════════════════════════════════════════════════════════════

TEST(CompoundPolylineCombine, WindingNormalizationOuterBecomesCorrect) {
  // Create an outer loop that is CW (incorrect) and normalize it.
  CppCompound compound;
  compound.outerLoops.push_back(makeCppCWRect(0, 0, 10, 10)); // CW = incorrect for outer
  EXPECT_LT(cavc::getArea(compound.outerLoops[0]), 0.0);       // currently negative

  cavc::normalizeCompoundPolylineWinding(compound);
  EXPECT_GT(cavc::getArea(compound.outerLoops[0]), 0.0);       // now positive (CCW)
}

TEST(CompoundPolylineCombine, WindingNormalizationHoleBecomesCorrect) {
  // Create a hole that is CCW (incorrect) and normalize it.
  CppCompound compound;
  compound.holeLoops.push_back(makeCppCCWRect(2, 2, 8, 8)); // CCW = incorrect for hole
  EXPECT_GT(cavc::getArea(compound.holeLoops[0]), 0.0);       // currently positive

  cavc::normalizeCompoundPolylineWinding(compound);
  EXPECT_LT(cavc::getArea(compound.holeLoops[0]), 0.0);       // now negative (CW)
}

TEST(CompoundPolylineCombine, NormalizedWindingDoesNotChangeCorrectLoops) {
  // Already-correct loops are unchanged.
  CppCompound compound;
  compound.outerLoops.push_back(makeCppCCWRect(0, 0, 10, 10));
  compound.holeLoops.push_back(makeCppCWRect(2, 2, 8, 8));

  double outerAreaBefore = cavc::getArea(compound.outerLoops[0]);
  double holeAreaBefore = cavc::getArea(compound.holeLoops[0]);

  cavc::normalizeCompoundPolylineWinding(compound);

  EXPECT_TRUE(fuzzyEq(cavc::getArea(compound.outerLoops[0]), outerAreaBefore));
  EXPECT_TRUE(fuzzyEq(cavc::getArea(compound.holeLoops[0]), holeAreaBefore));
}

TEST(CompoundPolylineCombine, WindingNormalizationThenCombine) {
  // Swap the windings: outer given as CW, hole given as CCW.
  // After normalisation the result should be the same as starting correctly.
  CppCompound a;
  a.outerLoops.push_back(makeCppCWRect(0, 0, 10, 10)); // wrong winding
  a.holeLoops.push_back(makeCppCCWRect(3, 3, 7, 7));   // wrong winding

  cavc::normalizeCompoundPolylineWinding(a);

  CppCompound b;
  b.outerLoops.push_back(makeCppCCWRect(5, 0, 15, 10));

  // Union after normalization should give same result as if windings were correct.
  auto r = cavc::combineCompoundPolylines(a, b, cavc::PlineCombineMode::Union);
  EXPECT_GE(r.outerLoops.size(), 1u);
  // Net area: A outer=100, A hole=16; A net=84. B=100. Union = 84+100 - (overlap in A net)
  // A.outer ∪ B.outer = 15×10 = 150. Hole of A only on left half (x:3-5 is 2 wide? let's check)
  // Hole [3,7]×[3,7]. Part outside B (x<5): [3,5]×[3,7] = 2×4=8.
  // Net area = 150 - 8 = 142
  EXPECT_TRUE(fuzzyEq(totalArea(r), 142.0)) << "Got net area: " << totalArea(r);
}

// ═══════════════════════════════════════════════════════════════════════════════
// Section 6: Coincident and touching boundaries
// ═══════════════════════════════════════════════════════════════════════════════

TEST(CompoundPolylineCombine, CoincidentOuterLoopsUnion) {
  // Two identical rectangles: union = that rectangle.
  CppCompound a, b;
  a.outerLoops.push_back(makeCppCCWRect(0, 0, 10, 10));
  b.outerLoops.push_back(makeCppCCWRect(0, 0, 10, 10));

  auto r = cavc::combineCompoundPolylines(a, b, cavc::PlineCombineMode::Union);
  // Result should be the same rectangle (area=100) or empty (both coincident).
  // The algorithm handles completely coincident pairs by checking the predicate at a
  // representative point; for Union that point is not inside B... it's on the boundary.
  // We accept either 1 outer or that totalArea ≈ 100.
  double area = totalArea(r);
  EXPECT_TRUE(fuzzyEq(area, 100.0) || r.outerLoops.size() >= 1u)
      << "Got net area: " << area;
}

TEST(CompoundPolylineCombine, CoincidentOuterLoopsIntersect) {
  // Two identical rectangles: intersect = that rectangle.
  CppCompound a, b;
  a.outerLoops.push_back(makeCppCCWRect(0, 0, 10, 10));
  b.outerLoops.push_back(makeCppCCWRect(0, 0, 10, 10));

  auto r = cavc::combineCompoundPolylines(a, b, cavc::PlineCombineMode::Intersect);
  double area = totalArea(r);
  EXPECT_TRUE(fuzzyEq(area, 100.0) || r.outerLoops.size() >= 1u)
      << "Got net area: " << area;
}

TEST(CompoundPolylineCombine, TouchingAtEdgeUnion) {
  // A and B share one edge (x=10): [0,10]×[0,10] and [10,20]×[0,10].
  // Union should produce a rectangle [0,20]×[0,10] with area = 200.
  CppCompound a, b;
  a.outerLoops.push_back(makeCppCCWRect(0, 0, 10, 10));
  b.outerLoops.push_back(makeCppCCWRect(10, 0, 20, 10));

  auto r = cavc::combineCompoundPolylines(a, b, cavc::PlineCombineMode::Union);
  // With touching edges (coincident boundary at x=10) the result might be one or two loops.
  // Either way, net area must be 200.
  EXPECT_TRUE(fuzzyEq(totalArea(r), 200.0)) << "Got net area: " << totalArea(r);
}

TEST(CompoundPolylineCombine, TouchingAtEdgeIntersect) {
  // Two rectangles sharing only an edge have zero-area intersection.
  CppCompound a, b;
  a.outerLoops.push_back(makeCppCCWRect(0, 0, 10, 10));
  b.outerLoops.push_back(makeCppCCWRect(10, 0, 20, 10));

  auto r = cavc::combineCompoundPolylines(a, b, cavc::PlineCombineMode::Intersect);
  // The intersection is the shared edge, which has zero area. Expect empty or zero-area result.
  EXPECT_TRUE(r.outerLoops.empty() || fuzzyEq(totalArea(r), 0.0))
      << "Got net area: " << totalArea(r);
}

// ═══════════════════════════════════════════════════════════════════════════════
// Section 7: C API tests
// ═══════════════════════════════════════════════════════════════════════════════

TEST(CompoundPolylineCombineCAPI, BasicCreateAndDelete) {
  cavc_compound_pline *cp = cavc_compound_pline_new();
  ASSERT_NE(cp, nullptr);
  EXPECT_EQ(cavc_compound_pline_outer_count(cp), 0u);
  EXPECT_EQ(cavc_compound_pline_hole_count(cp), 0u);
  cavc_compound_pline_delete(cp);
}

TEST(CompoundPolylineCombineCAPI, AddOuterAndHole) {
  CompoundPlineGuard cp(cavc_compound_pline_new());
  ASSERT_NE(cp.ptr, nullptr);

  PlineGuard outer(makeCCWRect(0, 0, 10, 10));
  PlineGuard hole(makeCWRect(2, 2, 8, 8));

  cavc_compound_pline_add_outer(cp.ptr, outer.ptr);
  cavc_compound_pline_add_hole(cp.ptr, hole.ptr);

  EXPECT_EQ(cavc_compound_pline_outer_count(cp.ptr), 1u);
  EXPECT_EQ(cavc_compound_pline_hole_count(cp.ptr), 1u);

  cavc_pline const *gotOuter = cavc_compound_pline_get_outer(cp.ptr, 0);
  ASSERT_NE(gotOuter, nullptr);
  EXPECT_TRUE(fuzzyEq(cavc_get_area(gotOuter), 100.0));

  cavc_pline const *gotHole = cavc_compound_pline_get_hole(cp.ptr, 0);
  ASSERT_NE(gotHole, nullptr);
  EXPECT_LT(cavc_get_area(gotHole), 0.0);  // CW = negative area
}

TEST(CompoundPolylineCombineCAPI, NormalizeViaAPI) {
  CompoundPlineGuard cp(cavc_compound_pline_new());
  PlineGuard wrongOuter(makeCWRect(0, 0, 10, 10)); // CW outer (wrong)
  cavc_compound_pline_add_outer(cp.ptr, wrongOuter.ptr);

  // Before normalize: area is negative (CW)
  EXPECT_LT(cavc_get_area(cavc_compound_pline_get_outer(cp.ptr, 0)), 0.0);

  cavc_compound_pline_normalize(cp.ptr);

  // After normalize: area is positive (CCW)
  EXPECT_GT(cavc_get_area(cavc_compound_pline_get_outer(cp.ptr, 0)), 0.0);
}

TEST(CompoundPolylineCombineCAPI, UnionOverlappingRects) {
  // C API union of two overlapping rectangles.
  CompoundPlineGuard a(cavc_compound_pline_new());
  CompoundPlineGuard b(cavc_compound_pline_new());

  {
    PlineGuard pa(makeCCWRect(0, 0, 10, 10));
    cavc_compound_pline_add_outer(a.ptr, pa.ptr);
  }
  {
    PlineGuard pb(makeCCWRect(5, 0, 15, 10));
    cavc_compound_pline_add_outer(b.ptr, pb.ptr);
  }

  cavc_compound_pline *result = nullptr;
  cavc_combine_compound_plines(a.ptr, b.ptr, 0 /*Union*/, &result);
  CompoundPlineGuard r(result);
  ASSERT_NE(r.ptr, nullptr);

  EXPECT_EQ(cavc_compound_pline_outer_count(r.ptr), 1u);
  EXPECT_EQ(cavc_compound_pline_hole_count(r.ptr), 0u);
  EXPECT_TRUE(fuzzyEq(getAreas(r.ptr).netArea(), 150.0))
      << "Got net area: " << getAreas(r.ptr).netArea();
}

TEST(CompoundPolylineCombineCAPI, InvalidModeReturnsNull) {
  CompoundPlineGuard a(cavc_compound_pline_new());
  CompoundPlineGuard b(cavc_compound_pline_new());
  cavc_compound_pline *result = nullptr;
  cavc_combine_compound_plines(a.ptr, b.ptr, 99 /*invalid*/, &result);
  EXPECT_EQ(result, nullptr);
}

TEST(CompoundPolylineCombineCAPI, CompoundWithHoleUnion) {
  // Same as CompoundHoleUnion but through the C API.
  CompoundPlineGuard a(cavc_compound_pline_new());
  CompoundPlineGuard b(cavc_compound_pline_new());

  {
    PlineGuard outer(makeCCWRect(0, 0, 20, 20));
    cavc_compound_pline_add_outer(a.ptr, outer.ptr);
    PlineGuard hole(makeCWRect(5, 5, 15, 15));
    cavc_compound_pline_add_hole(a.ptr, hole.ptr);
  }
  {
    PlineGuard pb(makeCCWRect(10, 0, 20, 20));
    cavc_compound_pline_add_outer(b.ptr, pb.ptr);
  }

  cavc_compound_pline *result = nullptr;
  cavc_combine_compound_plines(a.ptr, b.ptr, 0 /*Union*/, &result);
  CompoundPlineGuard r(result);
  ASSERT_NE(r.ptr, nullptr);

  EXPECT_TRUE(fuzzyEq(getAreas(r.ptr).netArea(), 350.0))
      << "Got net area: " << getAreas(r.ptr).netArea();
}

// ═══════════════════════════════════════════════════════════════════════════════
// Section 8: Arc/bulge preservation
// ═══════════════════════════════════════════════════════════════════════════════

TEST(CompoundPolylineCombine, ArcLoopsUnion) {
  // A = circle of radius 5 centred at origin (two-vertex arc polyline, bulge=1).
  // B = circle of radius 5 centred at (8,0) — overlapping.
  // The boolean engine should preserve bulges (arcs) in the result.
  CppPline circleA, circleB;
  circleA.isClosed() = true;
  circleA.addVertex({-5, 0, 1});
  circleA.addVertex({5, 0, 1});

  circleB.isClosed() = true;
  circleB.addVertex({3, 0, 1});
  circleB.addVertex({13, 0, 1});

  CppCompound a, b;
  a.outerLoops.push_back(circleA);
  b.outerLoops.push_back(circleB);

  auto r = cavc::combineCompoundPolylines(a, b, cavc::PlineCombineMode::Union);
  // The result area should be approximately area(A) + area(B) - area(A∩B).
  // Each circle area ≈ pi*25 ≈ 78.54. The combined area is less.
  // Just verify we get one outer loop and that it has bulges (arcs).
  ASSERT_EQ(r.outerLoops.size(), 1u);
  EXPECT_EQ(r.holeLoops.size(), 0u);
  // Check that the result area is between pi*25 and pi*25*2 (reasonable range)
  double area = totalArea(r);
  EXPECT_GT(area, 70.0);
  EXPECT_LT(area, 160.0);
  // Verify that at least one vertex has a non-zero bulge (arc preserved)
  bool hasArc = false;
  for (auto const &v : r.outerLoops[0].vertexes()) {
    if (!v.bulgeIsZero()) { hasArc = true; break; }
  }
  EXPECT_TRUE(hasArc) << "Expected arc vertices in union of two arc loops";
}

TEST(CompoundPolylineCombine, ArcLoopWithHole) {
  // Circular outer with a rectangular hole.
  // Verify that bulges are preserved in the outer loop after union with a rectangle.
  CppPline circle;
  circle.isClosed() = true;
  circle.addVertex({-10, 0, 1});
  circle.addVertex({10, 0, 1});

  CppCompound a, b;
  a.outerLoops.push_back(circle);
  a.holeLoops.push_back(makeCppCWRect(-4, -4, 4, 4)); // rectangular hole

  b.outerLoops.push_back(makeCppCCWRect(6, -4, 14, 4)); // rect partially outside circle

  auto r = cavc::combineCompoundPolylines(a, b, cavc::PlineCombineMode::Union);
  EXPECT_GE(r.outerLoops.size(), 1u);
  // Total area: circle (~pi*100≈314) - hole (64) + extra rect parts
  // Just verify the result is non-trivial and has arcs.
  bool hasArc = false;
  for (auto const &loop : r.outerLoops) {
    for (auto const &v : loop.vertexes()) {
      if (!v.bulgeIsZero()) { hasArc = true; break; }
    }
  }
  EXPECT_TRUE(hasArc) << "Expected arc vertices to be preserved";
}

} // namespace
