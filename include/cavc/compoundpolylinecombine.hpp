#ifndef CAVC_COMPOUNDPOLYLINECOMBINE_HPP
#define CAVC_COMPOUNDPOLYLINECOMBINE_HPP

// Compound-path boolean operations built on top of the closed-polyline boolean engine.
//
// A CompoundPolyline represents a filled 2-D region composed of:
//   - outer loops: CCW winding (positive area)
//   - hole  loops: CW  winding (negative area)
//
// A point p is "inside" a CompoundPolyline when the sum of winding numbers from all its loops is
// non-zero (the standard nonzero fill rule). With CCW outers (+1) and CW holes (-1) this
// naturally handles holes that subtract from the filled area.
//
// Supported boolean modes: Union, Exclude (A−B), Intersect, XOR.
// All use the same PlineCombineMode enum as combinePolylines.
//
// Usage:
//   CompoundPolyline<double> a, b;
//   a.outerLoops.push_back(outerLoop);
//   a.holeLoops.push_back(holeLoop);
//   normalizeCompoundPolylineWinding(a);   // ensure correct winding
//   auto result = combineCompoundPolylines(a, b, PlineCombineMode::Union);
//   // result.outerLoops / result.holeLoops

#include "plinesegment.hpp"
#include "polyline.hpp"
#include "polylinecombine.hpp"
#include <algorithm>
#include <unordered_map>
#include <vector>

namespace cavc {

/// A compound polyline: a filled 2-D region with outer loops (CCW) and hole loops (CW).
template <typename Real> struct CompoundPolyline {
  std::vector<Polyline<Real>> outerLoops;
  std::vector<Polyline<Real>> holeLoops;
};

/// Result of combining two CompoundPolylines.
template <typename Real> struct CompoundCombineResult {
  std::vector<Polyline<Real>> outerLoops;
  std::vector<Polyline<Real>> holeLoops;
};

/// Normalize winding in-place: outer loops → CCW (positive area), hole loops → CW (negative area).
/// Call this before combining if the loop windings are not already known to be correct.
template <typename Real> void normalizeCompoundPolylineWinding(CompoundPolyline<Real> &compound) {
  for (auto &loop : compound.outerLoops) {
    if (loop.isClosed() && getArea(loop) < Real(0)) {
      invertDirection(loop);
    }
  }
  for (auto &loop : compound.holeLoops) {
    if (loop.isClosed() && getArea(loop) > Real(0)) {
      invertDirection(loop);
    }
  }
}

// ─── Implementation details ────────────────────────────────────────────────────────────────────

namespace internal {

/// Returns true when the nonzero-winding fill rule places `pt` inside the compound path.
/// Correctly handles outer loops (+1) and hole loops (-1).
template <typename Real>
bool isInsideCompound(std::vector<Polyline<Real>> const &outerLoops,
                      std::vector<Polyline<Real>> const &holeLoops, Vector2<Real> const &pt) {
  int winding = 0;
  for (auto const &loop : outerLoops) {
    winding += getWindingNumber(loop, pt);
  }
  for (auto const &loop : holeLoops) {
    winding += getWindingNumber(loop, pt);
  }
  return winding != 0;
}

/// Slice `loop` at all intersection points accumulated in `rawIntrMap`
/// (mapping segment-start-index → list of intersection positions), keeping slices whose last-
/// segment midpoint satisfies `keepPred`. Materialised open polylines are appended to `result`.
/// If `invertSlices` is true each kept slice is direction-inverted before insertion.
template <typename Real, typename Pred>
void sliceLoopMultiPair(Polyline<Real> const &loop,
                        std::unordered_map<std::size_t, std::vector<Vector2<Real>>> const &rawIntrMap,
                        Pred &&keepPred, bool invertSlices, Real joinThreshold,
                        std::vector<Polyline<Real>> &result) {
  if (rawIntrMap.empty() || !loop.isClosed() || loop.size() < 2) {
    return;
  }

  // Build a sorted copy of rawIntrMap (sorted by distance from segment start vertex).
  std::unordered_map<std::size_t, std::vector<Vector2<Real>>> intrMap;
  intrMap.reserve(rawIntrMap.size());
  for (auto const &kvp : rawIntrMap) {
    Vector2<Real> const startPos = loop[kvp.first].pos();
    auto pts = kvp.second;
    std::sort(pts.begin(), pts.end(), [&](Vector2<Real> const &a, Vector2<Real> const &b) {
      return distSquared(a, startPos) < distSquared(b, startPos);
    });
    // Deduplicate (multiple pairs might contribute the same intersection point).
    pts.erase(std::unique(pts.begin(), pts.end(),
                          [](Vector2<Real> const &a, Vector2<Real> const &b) {
                            return fuzzyEqual(a, b, utils::realPrecision<Real>());
                          }),
              pts.end());
    if (!pts.empty()) {
      intrMap.emplace(kvp.first, std::move(pts));
    }
  }
  if (intrMap.empty()) {
    return;
  }

  // Helper: materialise a slice view and push to result if the predicate passes.
  auto maybePushSlice = [&](std::optional<PlineSliceViewData<Real>> const &viewOpt) {
    if (!viewOpt) {
      return;
    }
    // Predicate point: midpoint of the last segment of the slice (same convention as
    // the existing sliceAtIntersects / slicePredicatePoint).
    PlineVertex<Real> lastSegStart = viewOpt->firstVertex(loop);
    PlineVertex<Real> lastSegEnd = viewOpt->lastVertex(loop);
    viewOpt->visitSegments(loop, [&](PlineVertex<Real> const &v1, PlineVertex<Real> const &v2) {
      lastSegStart = v1;
      lastSegEnd = v2;
      return true;
    });
    if (!keepPred(segMidpoint(lastSegStart, lastSegEnd))) {
      return;
    }
    Polyline<Real> slicePline = viewOpt->toPolyline(loop, joinThreshold);
    if (invertSlices) {
      invertDirection(slicePline);
    }
    result.push_back(std::move(slicePline));
  };

  for (auto const &kvp : intrMap) {
    std::size_t const sIndex = kvp.first;
    std::vector<Vector2<Real>> const &intrsList = kvp.second;

    PlineVertex<Real> const &firstSegStartVertex = loop[sIndex];
    std::size_t const nextIndex = utils::nextWrappingIndex(sIndex, loop);
    PlineVertex<Real> const &firstSegEndVertex = loop[nextIndex];

    // Build sub-slices between consecutive intersection points within the same segment.
    if (intrsList.size() > 1) {
      SplitResult<Real> firstSplit =
          splitAtPoint(firstSegStartVertex, firstSegEndVertex, intrsList[0]);
      auto prevVertex = firstSplit.splitVertex;
      for (std::size_t i = 1; i < intrsList.size(); ++i) {
        SplitResult<Real> split = splitAtPoint(prevVertex, firstSegEndVertex, intrsList[i]);
        prevVertex = split.splitVertex;

        if (fuzzyEqual(split.updatedStart.pos(), split.splitVertex.pos(),
                       utils::realPrecision<Real>())) {
          continue;
        }
        maybePushSlice(PlineSliceViewData<Real>::createOnSingleSegment(loop, sIndex,
                                                                        split.updatedStart,
                                                                        intrsList[i]));
      }
    }

    // Build the wrap-around slice: from the last intersection point on this segment, traverse the
    // loop forward until the first intersection point of the next intersected segment is found.
    std::size_t index = nextIndex;
    std::size_t loopCount = 0;
    std::size_t const maxLoopCount = loop.size();
    while (true) {
      if (loopCount++ > maxLoopCount) {
        CAVC_ASSERT(false, "Bug detected in sliceLoopMultiPair: infinite loop guard triggered");
        break;
      }
      auto nextIntr = intrMap.find(index);
      if (nextIntr != intrMap.end()) {
        maybePushSlice(PlineSliceViewData<Real>::createFromSlicePoints(
            loop, intrsList.back(), sIndex, nextIntr->second[0], index));
        break;
      }
      index = utils::nextWrappingIndex(index, loop);
    }
  }
}

/// One pass of the compound boolean: collect kept slices from all A-loops and all B-loops.
/// Loops with no pairwise intersections are checked with a single representative point.
/// `invertBSlices` reverses B-loop slice directions (used for Exclude / XOR).
template <typename Real, typename PredA, typename PredB>
void collectCompoundSlices(CompoundPolyline<Real> const &a, CompoundPolyline<Real> const &b,
                           PredA &&keepFromA, PredB &&keepFromB, bool invertBSlices,
                           Real joinThreshold, std::vector<Polyline<Real>> &allSlices) {
  // Flatten A-loops and B-loops into pointer arrays.
  std::vector<Polyline<Real> const *> aLoops;
  aLoops.reserve(a.outerLoops.size() + a.holeLoops.size());
  for (auto const &l : a.outerLoops) {
    aLoops.push_back(&l);
  }
  for (auto const &l : a.holeLoops) {
    aLoops.push_back(&l);
  }

  std::vector<Polyline<Real> const *> bLoops;
  bLoops.reserve(b.outerLoops.size() + b.holeLoops.size());
  for (auto const &l : b.outerLoops) {
    bLoops.push_back(&l);
  }
  for (auto const &l : b.holeLoops) {
    bLoops.push_back(&l);
  }

  if (aLoops.empty() || bLoops.empty()) {
    // One compound path has no loops; handle each loop individually.
    for (auto *la : aLoops) {
      if (la->size() < 2 || !la->isClosed()) {
        continue;
      }
      auto midpt = segMidpoint((*la)[0], (*la)[1]);
      if (keepFromA(midpt)) {
        Polyline<Real> copy = *la;
        copy.isClosed() = false;
        copy.addVertex(copy[0]);
        allSlices.push_back(std::move(copy));
      }
    }
    for (auto *lb : bLoops) {
      if (lb->size() < 2 || !lb->isClosed()) {
        continue;
      }
      auto midpt = segMidpoint((*lb)[0], (*lb)[1]);
      if (keepFromB(midpt)) {
        Polyline<Real> copy = *lb;
        copy.isClosed() = false;
        copy.addVertex(copy[0]);
        if (invertBSlices) {
          invertDirection(copy);
        }
        allSlices.push_back(std::move(copy));
      }
    }
    return;
  }

  // Per-loop accumulated intersection maps and intersection flags.
  using IntrMap = std::unordered_map<std::size_t, std::vector<Vector2<Real>>>;
  std::vector<IntrMap> aIntrMaps(aLoops.size());
  std::vector<IntrMap> bIntrMaps(bLoops.size());
  std::vector<bool> aHasIntr(aLoops.size(), false);
  std::vector<bool> bHasIntr(bLoops.size(), false);

  // Build spatial indices for all A-loops (reused across pairs).
  std::vector<StaticSpatialIndex<Real>> aSpatialIdx;
  aSpatialIdx.reserve(aLoops.size());
  for (auto const *la : aLoops) {
    if (la->size() > 1) {
      aSpatialIdx.push_back(createApproxSpatialIndex(*la));
    } else {
      // Placeholder; pairs involving this loop will be skipped below.
      aSpatialIdx.push_back(StaticSpatialIndex<Real>(1));
    }
  }

  // Find all pairwise intersections and accumulate per-loop.
  for (std::size_t i = 0; i < aLoops.size(); ++i) {
    auto const *la = aLoops[i];
    if (la->size() < 2) {
      continue;
    }
    for (std::size_t j = 0; j < bLoops.size(); ++j) {
      auto const *lb = bLoops[j];
      if (lb->size() < 2) {
        continue;
      }

      auto pairInfo = processForCombine(*la, *lb, aSpatialIdx[i]);
      if (!pairInfo.anyIntersects()) {
        continue;
      }

      if (pairInfo.completelyCoincident()) {
        // For completely coincident pairs, test a representative point with both predicates and
        // add the loop once if either accepts it.
        if (la->size() >= 2) {
          auto midpt = segMidpoint((*la)[0], (*la)[1]);
          if (keepFromA(midpt)) {
            Polyline<Real> copy = *la;
            copy.isClosed() = false;
            copy.addVertex(copy[0]);
            allSlices.push_back(std::move(copy));
          }
        }
        continue;
      }

      aHasIntr[i] = true;
      bHasIntr[j] = true;

      // Accumulate non-coincident intersection points.
      for (auto const &intr : pairInfo.nonCoincidentIntersects) {
        aIntrMaps[i][intr.sIndex1].push_back(intr.pos);
        bIntrMaps[j][intr.sIndex2].push_back(intr.pos);
      }

      // Treat coincident-slice endpoints as regular intersection points so the boundary is
      // cut at the right places (the midpoint predicate handles keep/discard correctly).
      for (auto const &slice : pairInfo.coincidentSlices) {
        aIntrMaps[i][slice.startPointOnA.sIndex].push_back(slice.startPointOnA.pos);
        aIntrMaps[i][slice.endPointOnA.sIndex].push_back(slice.endPointOnA.pos);
        bIntrMaps[j][slice.startPointOnB.sIndex].push_back(slice.startPointOnB.pos);
        bIntrMaps[j][slice.endPointOnB.sIndex].push_back(slice.endPointOnB.pos);
      }
    }
  }

  // Slice each A-loop that has intersections.
  for (std::size_t i = 0; i < aLoops.size(); ++i) {
    if (!aHasIntr[i]) {
      continue;
    }
    std::size_t beforeSize = allSlices.size();
    sliceLoopMultiPair(*aLoops[i], aIntrMaps[i], keepFromA, false, joinThreshold, allSlices);
    // If no slices were produced despite having intersections, fall back to whole-loop test
    // (handles tangent / touching-only cases).
    if (allSlices.size() == beforeSize && aLoops[i]->size() >= 2) {
      auto midpt = segMidpoint((*aLoops[i])[0], (*aLoops[i])[1]);
      if (keepFromA(midpt)) {
        Polyline<Real> copy = *aLoops[i];
        copy.isClosed() = false;
        copy.addVertex(copy[0]);
        allSlices.push_back(std::move(copy));
      }
    }
  }

  // Slice each B-loop that has intersections.
  for (std::size_t j = 0; j < bLoops.size(); ++j) {
    if (!bHasIntr[j]) {
      continue;
    }
    std::size_t beforeSize = allSlices.size();
    sliceLoopMultiPair(*bLoops[j], bIntrMaps[j], keepFromB, invertBSlices, joinThreshold,
                       allSlices);
    if (allSlices.size() == beforeSize && bLoops[j]->size() >= 2) {
      auto midpt = segMidpoint((*bLoops[j])[0], (*bLoops[j])[1]);
      if (keepFromB(midpt)) {
        Polyline<Real> copy = *bLoops[j];
        copy.isClosed() = false;
        copy.addVertex(copy[0]);
        if (invertBSlices) {
          invertDirection(copy);
        }
        allSlices.push_back(std::move(copy));
      }
    }
  }

  // Handle A-loops with no intersections (whole-loop predicate check).
  for (std::size_t i = 0; i < aLoops.size(); ++i) {
    if (aHasIntr[i] || aLoops[i]->size() < 2 || !aLoops[i]->isClosed()) {
      continue;
    }
    auto midpt = segMidpoint((*aLoops[i])[0], (*aLoops[i])[1]);
    if (keepFromA(midpt)) {
      Polyline<Real> copy = *aLoops[i];
      copy.isClosed() = false;
      copy.addVertex(copy[0]);
      allSlices.push_back(std::move(copy));
    }
  }

  // Handle B-loops with no intersections.
  for (std::size_t j = 0; j < bLoops.size(); ++j) {
    if (bHasIntr[j] || bLoops[j]->size() < 2 || !bLoops[j]->isClosed()) {
      continue;
    }
    auto midpt = segMidpoint((*bLoops[j])[0], (*bLoops[j])[1]);
    if (keepFromB(midpt)) {
      Polyline<Real> copy = *bLoops[j];
      copy.isClosed() = false;
      copy.addVertex(copy[0]);
      if (invertBSlices) {
        invertDirection(copy);
      }
      allSlices.push_back(std::move(copy));
    }
  }
}

} // namespace internal

// ─── Public API ───────────────────────────────────────────────────────────────────────────────

/// Combine two compound polylines using the given boolean mode.
///
/// Input loops should be normalised (outer = CCW, holes = CW). Call
/// normalizeCompoundPolylineWinding() on each input if this is not already guaranteed.
///
/// The result is classified by area sign: positive area → outerLoops, negative area → holeLoops.
///
/// Supports arcs/bulges: arc vertices are preserved through the boolean; no conversion to line
/// segments is performed.
template <typename Real>
CompoundCombineResult<Real> combineCompoundPolylines(CompoundPolyline<Real> const &a,
                                                     CompoundPolyline<Real> const &b,
                                                     PlineCombineMode mode) {
  using namespace internal;
  Real const joinThreshold = utils::sliceJoinThreshold<Real>();

  // Compound inside predicates (nonzero winding rule, accounting for all loops).
  auto isInsideA = [&](Vector2<Real> const &pt) {
    return isInsideCompound(a.outerLoops, a.holeLoops, pt);
  };
  auto isInsideB = [&](Vector2<Real> const &pt) {
    return isInsideCompound(b.outerLoops, b.holeLoops, pt);
  };

  std::vector<Polyline<Real>> allSlices;

  switch (mode) {
  case PlineCombineMode::Union:
    // Keep A boundary segments outside B; keep B boundary segments outside A.
    collectCompoundSlices(a, b, [&](auto const &pt) { return !isInsideB(pt); },
                          [&](auto const &pt) { return !isInsideA(pt); }, false, joinThreshold,
                          allSlices);
    break;

  case PlineCombineMode::Intersect:
    // Keep A boundary segments inside B; keep B boundary segments inside A.
    collectCompoundSlices(a, b, [&](auto const &pt) { return isInsideB(pt); },
                          [&](auto const &pt) { return isInsideA(pt); }, false, joinThreshold,
                          allSlices);
    break;

  case PlineCombineMode::Exclude:
    // A − B: keep A outside B; keep B inside A (inverted so it forms the inner boundary).
    collectCompoundSlices(a, b, [&](auto const &pt) { return !isInsideB(pt); },
                          [&](auto const &pt) { return isInsideA(pt); }, true, joinThreshold,
                          allSlices);
    break;

  case PlineCombineMode::XOR:
    // XOR = (A − B) ∪ (B − A), processed in two passes.
    // Pass 1: A − B
    collectCompoundSlices(a, b, [&](auto const &pt) { return !isInsideB(pt); },
                          [&](auto const &pt) { return isInsideA(pt); }, true, joinThreshold,
                          allSlices);
    // Pass 2: B − A
    collectCompoundSlices(b, a, [&](auto const &pt) { return !isInsideA(pt); },
                          [&](auto const &pt) { return isInsideB(pt); }, true, joinThreshold,
                          allSlices);
    break;
  }

  // Stitch open slices into closed polylines.
  auto closedLoops =
      stitchOrderedSlicesIntoClosedPolylines(allSlices, StitchFirstAvailable(), joinThreshold);

  // Classify by area sign: positive (CCW) → outer, negative (CW) → hole.
  CompoundCombineResult<Real> result;
  for (auto &loop : closedLoops) {
    if (getArea(loop) >= Real(0)) {
      result.outerLoops.push_back(std::move(loop));
    } else {
      result.holeLoops.push_back(std::move(loop));
    }
  }
  return result;
}

} // namespace cavc

#endif // CAVC_COMPOUNDPOLYLINECOMBINE_HPP
