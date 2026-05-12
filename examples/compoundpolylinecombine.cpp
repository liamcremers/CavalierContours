// Example: compound-path boolean operations.
//
// This example demonstrates combineCompoundPolylines() for shapes that have
// both outer loops and hole loops.  Run the program and inspect the printed
// areas to verify the results.

#include "cavc/compoundpolylinecombine.hpp"
#include <iostream>

using namespace cavc;

// Build a CCW rectangle loop (positive area / outer).
static Polyline<double> ccwRect(double x0, double y0, double x1, double y1) {
  Polyline<double> p;
  p.isClosed() = true;
  p.addVertex(x0, y0, 0);
  p.addVertex(x1, y0, 0);
  p.addVertex(x1, y1, 0);
  p.addVertex(x0, y1, 0);
  return p;
}

// Build a CW rectangle loop (negative area / hole).
static Polyline<double> cwRect(double x0, double y0, double x1, double y1) {
  Polyline<double> p;
  p.isClosed() = true;
  p.addVertex(x0, y1, 0);
  p.addVertex(x1, y1, 0);
  p.addVertex(x1, y0, 0);
  p.addVertex(x0, y0, 0);
  return p;
}

static void printResult(char const *label, CompoundCombineResult<double> const &r) {
  double outerArea = 0, holeArea = 0;
  for (auto const &l : r.outerLoops) outerArea += getArea(l);
  for (auto const &l : r.holeLoops) holeArea += getArea(l);
  std::cout << label << "\n"
            << "  outer loops : " << r.outerLoops.size() << "  (area sum = " << outerArea << ")\n"
            << "  hole  loops : " << r.holeLoops.size() << "  (area sum = " << holeArea << ")\n"
            << "  net area    : " << (outerArea + holeArea) << "\n\n";
}

int main(int argc, char *argv[]) {
  (void)argc;
  (void)argv;

  // ── Example 1: simple outer-vs-outer ──────────────────────────────────────
  // A = [0,10]×[0,10] (area = 100), B = [5,15]×[0,10] (area = 100)
  {
    CompoundPolyline<double> a, b;
    a.outerLoops.push_back(ccwRect(0, 0, 10, 10));
    b.outerLoops.push_back(ccwRect(5, 0, 15, 10));

    printResult("Simple outer-vs-outer union (expected net area=150):",
                combineCompoundPolylines(a, b, PlineCombineMode::Union));
    printResult("Simple outer-vs-outer intersect (expected net area=50):",
                combineCompoundPolylines(a, b, PlineCombineMode::Intersect));
    printResult("Simple outer-vs-outer exclude A-B (expected net area=50):",
                combineCompoundPolylines(a, b, PlineCombineMode::Exclude));
    printResult("Simple outer-vs-outer XOR (expected net area=100, 2 outer loops):",
                combineCompoundPolylines(a, b, PlineCombineMode::XOR));
  }

  // ── Example 2: compound path with a hole ──────────────────────────────────
  // A = 20×20 square with a 10×10 square hole in the middle.
  // Net area of A = 400 - 100 = 300.
  // B = right half of A's bounding box: [10,20]×[0,20] (area = 200).
  {
    CompoundPolyline<double> a, b;
    a.outerLoops.push_back(ccwRect(0, 0, 20, 20));
    a.holeLoops.push_back(cwRect(5, 5, 15, 15)); // 10×10 hole

    b.outerLoops.push_back(ccwRect(10, 0, 20, 20)); // right half

    // After normalize (already correct here, but shown for documentation):
    normalizeCompoundPolylineWinding(a);
    normalizeCompoundPolylineWinding(b);

    printResult("Compound-hole union (expected net area=350):",
                combineCompoundPolylines(a, b, PlineCombineMode::Union));
    printResult("Compound-hole intersect (expected net area=150):",
                combineCompoundPolylines(a, b, PlineCombineMode::Intersect));
    printResult("Compound-hole exclude A-B (expected net area=150):",
                combineCompoundPolylines(a, b, PlineCombineMode::Exclude));
  }

  // ── Example 3: B entirely inside A's hole ─────────────────────────────────
  // B is a small square fully within A's hole → B is NOT "inside A" (compound).
  // Union should add B as an island inside the hole.
  {
    CompoundPolyline<double> a, b;
    a.outerLoops.push_back(ccwRect(0, 0, 20, 20));
    a.holeLoops.push_back(cwRect(5, 5, 15, 15));  // 10×10 hole
    b.outerLoops.push_back(ccwRect(7, 7, 11, 11)); // 4×4 inside hole

    printResult("B inside A's hole - union (expected net area=316):",
                combineCompoundPolylines(a, b, PlineCombineMode::Union));
    printResult("B inside A's hole - intersect (expected empty):",
                combineCompoundPolylines(a, b, PlineCombineMode::Intersect));
  }

  // ── Example 4: winding normalization ──────────────────────────────────────
  // Same as example 2 but with loops given in the wrong winding order.
  {
    CompoundPolyline<double> a;
    a.outerLoops.push_back(cwRect(0, 0, 20, 20));  // wrong: CW outer
    a.holeLoops.push_back(ccwRect(5, 5, 15, 15));  // wrong: CCW hole

    std::cout << "Before normalization:\n"
              << "  outer area = " << getArea(a.outerLoops[0]) << " (expected negative)\n"
              << "  hole  area = " << getArea(a.holeLoops[0]) << " (expected positive)\n\n";

    normalizeCompoundPolylineWinding(a);

    std::cout << "After normalization:\n"
              << "  outer area = " << getArea(a.outerLoops[0]) << " (expected positive)\n"
              << "  hole  area = " << getArea(a.holeLoops[0]) << " (expected negative)\n\n";
  }

  return 0;
}
