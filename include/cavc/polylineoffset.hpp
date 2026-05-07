#ifndef CAVC_POLYLINEOFFSET_HPP
#define CAVC_POLYLINEOFFSET_HPP
#include <algorithm>
#include <cstdint>
#include <cmath>
#include <functional>

#include "internal/plinesliceview.hpp"
#include "polyline.hpp"
#include "polylineintersects.hpp"
#include <limits>
#include <map>
#include <unordered_map>
#include <vector>

// This header has functions for offsetting polylines

namespace cavc {
enum class OffsetJoinType { Round = 0, Miter = 1, Bevel = 2 };

enum class OffsetEndCapType { Round = 0, Square = 1, Butt = 2 };

template <typename Real> struct ParallelOffsetOptions {
  bool hasSelfIntersects = false;
  OffsetJoinType joinType = OffsetJoinType::Round;
  OffsetEndCapType endCapType = OffsetEndCapType::Round;
  Real miterLimit = Real(4);
};

namespace internal {
/// Represents a raw polyline offset segment.
template <typename Real> struct PlineOffsetSegment {
  PlineVertex<Real> v1;
  PlineVertex<Real> v2;
  Vector2<Real> origV2Pos;
  bool collapsedArc;
};

/// Creates all the raw polyline offset segments.
template <typename Real>
std::vector<PlineOffsetSegment<Real>> createUntrimmedOffsetSegments(Polyline<Real> const &pline,
                                                                    Real offset) {
  std::size_t segmentCount = pline.isClosed() ? pline.size() : pline.size() - 1;

  std::vector<PlineOffsetSegment<Real>> result;
  result.reserve(segmentCount);

  auto lineVisitor = [&](PlineVertex<Real> const &v1, PlineVertex<Real> const &v2) {
    result.emplace_back();
    PlineOffsetSegment<Real> &seg = result.back();
    seg.collapsedArc = false;
    seg.origV2Pos = v2.pos();
    Vector2<Real> edge = v2.pos() - v1.pos();
    Vector2<Real> offsetV = offset * safeUnitPerp(edge);
    seg.v1.pos() = v1.pos() + offsetV;
    seg.v1.bulge() = v1.bulge();
    seg.v2.pos() = v2.pos() + offsetV;
    seg.v2.bulge() = v2.bulge();
  };

  auto arcVisitor = [&](PlineVertex<Real> const &v1, PlineVertex<Real> const &v2) {
    auto arc = arcRadiusAndCenter(v1, v2);
    Real offs = v1.bulgeIsNeg() ? offset : -offset;
    Real radiusAfterOffset = arc.radius + offs;
    Vector2<Real> v1ToCenter = v1.pos() - arc.center;
    safeNormalize(v1ToCenter);
    Vector2<Real> v2ToCenter = v2.pos() - arc.center;
    safeNormalize(v2ToCenter);

    result.emplace_back();
    PlineOffsetSegment<Real> &seg = result.back();
    seg.origV2Pos = v2.pos();
    seg.v1.pos() = offs * v1ToCenter + v1.pos();
    seg.v2.pos() = offs * v2ToCenter + v2.pos();
    seg.v2.bulge() = v2.bulge();

    if (radiusAfterOffset < utils::realThreshold<Real>()) {
      // collapsed arc, offset arc start and end points towards arc center and turn into line
      // handles case where offset vertexes are equal and simplifies path for clipping algorithm
      seg.collapsedArc = true;
      seg.v1.bulge() = Real(0);
    } else {
      seg.collapsedArc = false;
      seg.v1.bulge() = v1.bulge();
    }
  };

  auto offsetVisitor = [&](PlineVertex<Real> const &v1, PlineVertex<Real> const &v2) {
    if (v1.bulgeIsZero()) {
      lineVisitor(v1, v2);
    } else {
      arcVisitor(v1, v2);
    }
  };

  for (std::size_t i = 1; i < pline.size(); ++i) {
    offsetVisitor(pline[i - 1], pline[i]);
  }

  if (pline.isClosed()) {
    offsetVisitor(pline.lastVertex(), pline[0]);
  }

  return result;
}

template <typename Real> bool falseIntersect(Real t) { return t < 0.0 || t > 1.0; }

// Gets the bulge to describe the arc going from start point to end point with the given arc center
// and curve orientation, if orientation is negative then bulge is negative otherwise it is positive
template <typename Real>
Real bulgeForConnection(Vector2<Real> const &arcCenter, Vector2<Real> const &sp,
                        Vector2<Real> const &ep, bool isCCW) {
  Real a1 = angle(arcCenter, sp);
  Real a2 = angle(arcCenter, ep);
  Real absSweepAngle = std::abs(utils::deltaAngle(a1, a2));
  Real absBulge = std::tan(absSweepAngle / Real(4));
  if (isCCW) {
    return absBulge;
  }

  return -absBulge;
}

template <typename Real>
void lineToLineJoin(PlineOffsetSegment<Real> const &s1, PlineOffsetSegment<Real> const &s2,
                    bool connectionArcsAreCCW, OffsetJoinType joinType, Real miterLimit,
                    Polyline<Real> &result) {
  const auto &v1 = s1.v1;
  const auto &v2 = s1.v2;
  const auto &u1 = s2.v1;
  const auto &u2 = s2.v2;
  CAVC_ASSERT(v1.bulgeIsZero() && u1.bulgeIsZero(), "both segs should be lines");
  CAVC_ASSERT(joinType == OffsetJoinType::Round || !s1.collapsedArc,
              "non-round join with collapsed arc is not supported");
  CAVC_ASSERT(joinType == OffsetJoinType::Round || !s2.collapsedArc,
              "non-round join with collapsed arc is not supported");

  auto connectUsingArc = [&] {
    auto const &arcCenter = s1.origV2Pos;
    auto const &sp = v2.pos();
    auto const &ep = u1.pos();
    Real bulge = bulgeForConnection(arcCenter, sp, ep, connectionArcsAreCCW);
    addOrReplaceIfSamePos(result, PlineVertex<Real>(sp, bulge));
    addOrReplaceIfSamePos(result, PlineVertex<Real>(ep, Real(0)));
  };

  auto miterRatio = [&](Vector2<Real> const &miter_point) {
    Real offset_dist = length(v2.pos() - s1.origV2Pos);
    if (offset_dist <= utils::realThreshold<Real>()) {
      return std::numeric_limits<Real>::infinity();
    }

    Real miter_dist = length(miter_point - s1.origV2Pos);
    return miter_dist / offset_dist;
  };

  if (s1.collapsedArc || s2.collapsedArc) {
    // connecting to/from collapsed arc, always connect using arc
    connectUsingArc();
  } else {
    auto intrResult = intrLineSeg2LineSeg2(v1.pos(), v2.pos(), u1.pos(), u2.pos());

    switch (intrResult.intrType) {
    case LineSeg2LineSeg2IntrType::None:
      if (joinType == OffsetJoinType::Round) {
        // Parallel offset segments of a collapsed loop should connect with a half circle instead of
        // a straight chord.
        connectUsingArc();
      } else {
        addOrReplaceIfSamePos(result, PlineVertex<Real>(v2.pos(), Real(0)));
        addOrReplaceIfSamePos(result, u1);
      }
      break;
    case LineSeg2LineSeg2IntrType::True:
      addOrReplaceIfSamePos(result, PlineVertex<Real>(intrResult.point, Real(0)));
      break;
    case LineSeg2LineSeg2IntrType::Coincident:
      addOrReplaceIfSamePos(result, PlineVertex<Real>(v2.pos(), Real(0)));
      break;
    case LineSeg2LineSeg2IntrType::False:
      if (joinType == OffsetJoinType::Round) {
        if (intrResult.t0 > Real(1) && falseIntersect(intrResult.t1)) {
          // extend and join the lines together using an arc
          connectUsingArc();
        } else {
          addOrReplaceIfSamePos(result, PlineVertex<Real>(v2.pos(), Real(0)));
          addOrReplaceIfSamePos(result, u1);
        }
      } else if (joinType == OffsetJoinType::Miter && intrResult.t0 > Real(1) &&
                 intrResult.t1 < Real(0) && miterRatio(intrResult.point) <= miterLimit) {
        addOrReplaceIfSamePos(result, PlineVertex<Real>(intrResult.point, Real(0)));
      } else {
        addOrReplaceIfSamePos(result, PlineVertex<Real>(v2.pos(), Real(0)));
        addOrReplaceIfSamePos(result, u1);
      }
      break;
    }
  }
}

template <typename Real>
Vector2<Real> joinSegmentTangent(PlineOffsetSegment<Real> const &segment, bool atEndPoint) {
  if (segment.v1.bulgeIsZero()) {
    return segment.v2.pos() - segment.v1.pos();
  }

  Vector2<Real> pointOnSeg = atEndPoint ? segment.v2.pos() : segment.v1.pos();
  return segTangentVector(segment.v1, segment.v2, pointOnSeg);
}

template <typename Real>
bool tryComputeMiterPointForArcJoin(PlineOffsetSegment<Real> const &s1,
                                    PlineOffsetSegment<Real> const &s2, Vector2<Real> &miterPoint) {
  Vector2<Real> t1 = joinSegmentTangent(s1, true);
  Vector2<Real> t2 = joinSegmentTangent(s2, false);
  if (fuzzyZero(t1) || fuzzyZero(t2)) {
    return false;
  }

  normalize(t1);
  normalize(t2);
  Vector2<Real> p = s1.v2.pos();
  Vector2<Real> q = s2.v1.pos();
  Vector2<Real> r = t1;
  Vector2<Real> s = -t2;
  Real denom = perpDot(r, s);
  if (std::abs(denom) <= utils::realThreshold<Real>()) {
    return false;
  }

  Vector2<Real> qp = q - p;
  Real t = perpDot(qp, s) / denom;
  Real u = perpDot(qp, r) / denom;
  Real eps = utils::realPrecision<Real>();
  if (t < -eps || u < -eps) {
    return false;
  }

  miterPoint = p + t * r;
  return true;
}

template <typename Real>
void trimPreviousJoinSegmentAtPoint(PlineOffsetSegment<Real> const &s1, Vector2<Real> const &point,
                                    Polyline<Real> &result) {
  CAVC_ASSERT(result.size() > 0, "join result must contain previous segment start");

  if (!s1.v1.bulgeIsZero()) {
    PlineVertex<Real> &prevVertex = result.lastVertex();
    if (!prevVertex.bulgeIsZero() && !fuzzyEqual(prevVertex.pos(), s1.v2.pos())) {
      auto prevArc = arcRadiusAndCenter(prevVertex, s1.v2);
      Real prevArcStartAngle = angle(prevArc.center, prevVertex.pos());
      Real updatedPrevTheta = utils::deltaAngle(prevArcStartAngle, angle(prevArc.center, point));
      if ((updatedPrevTheta > Real(0)) == prevVertex.bulgeIsPos()) {
        prevVertex.bulge() = std::tan(updatedPrevTheta / Real(4));
      }
    }
  }

  addOrReplaceIfSamePos(result, PlineVertex<Real>(point, Real(0)));
}

template <typename Real>
PlineVertex<Real> createTrimmedNextJoinStart(PlineOffsetSegment<Real> const &s2,
                                             Vector2<Real> const &point) {
  if (s2.v1.bulgeIsZero()) {
    return PlineVertex<Real>(point, Real(0));
  }

  auto arc = arcRadiusAndCenter(s2.v1, s2.v2);
  Real a = angle(arc.center, point);
  Real endAngle = angle(arc.center, s2.v2.pos());
  Real theta = utils::deltaAngle(a, endAngle);
  if ((theta > Real(0)) == s2.v1.bulgeIsPos()) {
    return PlineVertex<Real>(point, std::tan(theta / Real(4)));
  }

  return PlineVertex<Real>(point, s2.v1.bulge());
}

template <typename Real>
bool miterLimitAllowsPoint(PlineOffsetSegment<Real> const &s1, Vector2<Real> const &point,
                           Real miterLimit) {
  Real offsetDist = length(s1.v2.pos() - s1.origV2Pos);
  if (offsetDist <= utils::realThreshold<Real>()) {
    return false;
  }

  Real miterDist = length(point - s1.origV2Pos);
  return miterDist / offsetDist <= miterLimit + utils::realPrecision<Real>();
}

template <typename Real>
bool tryAddExactNonRoundArcJoin(PlineOffsetSegment<Real> const &s1,
                                PlineOffsetSegment<Real> const &s2, Real miterLimit,
                                Polyline<Real> &result) {
  const bool s1IsLine = s1.v1.bulgeIsZero();
  const bool s2IsLine = s2.v1.bulgeIsZero();
  Real const eps = utils::realPrecision<Real>();

  auto betterCandidate = [&](bool &found, Vector2<Real> &best, Vector2<Real> const &candidate) {
    if (!miterLimitAllowsPoint(s1, candidate, miterLimit)) {
      return;
    }

    if (!found || distSquared(candidate, s1.origV2Pos) < distSquared(best, s1.origV2Pos)) {
      best = candidate;
      found = true;
    }
  };

  bool found = false;
  Vector2<Real> best = Vector2<Real>::zero();

  if (s1IsLine && !s2IsLine) {
    auto arc = arcRadiusAndCenter(s2.v1, s2.v2);
    auto addIfValid = [&](Real t) {
      if (t < -eps) {
        return;
      }

      Vector2<Real> candidate = pointFromParametric(s1.v1.pos(), s1.v2.pos(), t);
      if (pointWithinArcSweepAngle(arc.center, s2.v1.pos(), s2.v2.pos(), s2.v1.bulge(),
                                   candidate, eps)) {
        betterCandidate(found, best, candidate);
      }
    };

    auto intrResult = intrLineSeg2Circle2(s1.v1.pos(), s1.v2.pos(), arc.radius, arc.center);
    if (intrResult.numIntersects > 0) {
      addIfValid(intrResult.t0);
    }
    if (intrResult.numIntersects > 1) {
      addIfValid(intrResult.t1);
    }

    if (found) {
      addOrReplaceIfSamePos(result, createTrimmedNextJoinStart(s2, best));
      return true;
    }

    return false;
  }

  if (!s1IsLine && s2IsLine) {
    auto arc = arcRadiusAndCenter(s1.v1, s1.v2);
    auto addIfValid = [&](Real t) {
      if (t > Real(1) + eps) {
        return;
      }

      Vector2<Real> candidate = pointFromParametric(s2.v1.pos(), s2.v2.pos(), t);
      if (pointWithinArcSweepAngle(arc.center, s1.v1.pos(), s1.v2.pos(), s1.v1.bulge(),
                                   candidate, eps)) {
        betterCandidate(found, best, candidate);
      }
    };

    auto intrResult = intrLineSeg2Circle2(s2.v1.pos(), s2.v2.pos(), arc.radius, arc.center);
    if (intrResult.numIntersects > 0) {
      addIfValid(intrResult.t0);
    }
    if (intrResult.numIntersects > 1) {
      addIfValid(intrResult.t1);
    }

    if (found) {
      trimPreviousJoinSegmentAtPoint(s1, best, result);
      return true;
    }

    return false;
  }

  if (!s1IsLine && !s2IsLine) {
    auto arc1 = arcRadiusAndCenter(s1.v1, s1.v2);
    auto arc2 = arcRadiusAndCenter(s2.v1, s2.v2);
    auto addIfValid = [&](Vector2<Real> const &candidate) {
      if (pointWithinArcSweepAngle(arc1.center, s1.v1.pos(), s1.v2.pos(), s1.v1.bulge(),
                                   candidate, eps) &&
          pointWithinArcSweepAngle(arc2.center, s2.v1.pos(), s2.v2.pos(), s2.v1.bulge(),
                                   candidate, eps)) {
        betterCandidate(found, best, candidate);
      }
    };

    auto intrResult = intrCircle2Circle2(arc1.radius, arc1.center, arc2.radius, arc2.center);
    switch (intrResult.intrType) {
    case Circle2Circle2IntrType::NoIntersect:
    case Circle2Circle2IntrType::Coincident:
      break;
    case Circle2Circle2IntrType::OneIntersect:
      addIfValid(intrResult.point1);
      break;
    case Circle2Circle2IntrType::TwoIntersects:
      addIfValid(intrResult.point1);
      addIfValid(intrResult.point2);
      break;
    }

    if (found) {
      trimPreviousJoinSegmentAtPoint(s1, best, result);
      addOrReplaceIfSamePos(result, createTrimmedNextJoinStart(s2, best));
      return true;
    }
  }

  return false;
}

template <typename Real>
void nonRoundArcJoin(PlineOffsetSegment<Real> const &s1, PlineOffsetSegment<Real> const &s2,
                     OffsetJoinType joinType, Real miterLimit, Polyline<Real> &result) {
  CAVC_ASSERT(joinType != OffsetJoinType::Round, "use round join functions for round joins");

  auto connectUsingBevel = [&] {
    addOrReplaceIfSamePos(result, PlineVertex<Real>(s1.v2.pos(), Real(0)));
    addOrReplaceIfSamePos(result, s2.v1);
  };

  if (joinType == OffsetJoinType::Bevel || s1.collapsedArc || s2.collapsedArc) {
    connectUsingBevel();
    return;
  }

  CAVC_ASSERT(joinType == OffsetJoinType::Miter, "unsupported non-round join type");
  if (tryAddExactNonRoundArcJoin(s1, s2, miterLimit, result)) {
    return;
  }

  Vector2<Real> miterPoint;
  if (!tryComputeMiterPointForArcJoin(s1, s2, miterPoint)) {
    connectUsingBevel();
    return;
  }

  if (!miterLimitAllowsPoint(s1, miterPoint, miterLimit)) {
    connectUsingBevel();
    return;
  }

  // A tangent-line miter point for an arc-involved join is generally not on the offset circle.
  // Using it as an arc endpoint corrupts the bulge representation and can generate unstable
  // spikes/self-intersections after slicing. Preserve arc endpoints and connect them with straight
  // miter bridge legs; line segments are still extended/trimmed directly to the miter point.
  if (!s1.v1.bulgeIsZero()) {
    addOrReplaceIfSamePos(result, PlineVertex<Real>(s1.v2.pos(), Real(0)));
  }

  addOrReplaceIfSamePos(result, PlineVertex<Real>(miterPoint, Real(0)));

  if (!s2.v1.bulgeIsZero()) {
    addOrReplaceIfSamePos(result, s2.v1);
  }
}

template <typename Real>
Vector2<Real> openPolylineEndpointTangent(Polyline<Real> const &pline, bool atStart) {
  CAVC_ASSERT(!pline.isClosed(), "endpoint tangent is only for open polyline");
  CAVC_ASSERT(pline.size() >= 2, "open polyline must have at least 2 vertexes");

  auto tangentAlongSegment = [](PlineVertex<Real> const &segStart, PlineVertex<Real> const &segEnd,
                                bool atSegmentStartPoint) {
    if (segStart.bulgeIsZero()) {
      Vector2<Real> d = segEnd.pos() - segStart.pos();
      safeNormalize(d);
      return d;
    }

    auto arc = arcRadiusAndCenter(segStart, segEnd);
    Vector2<Real> radial =
        atSegmentStartPoint ? (segStart.pos() - arc.center) : (segEnd.pos() - arc.center);
    safeNormalize(radial);
    Vector2<Real> tangent = unitPerp(radial);
    if (segStart.bulgeIsNeg()) {
      tangent = -tangent;
    }
    return tangent;
  };

  if (atStart) {
    return tangentAlongSegment(pline[0], pline[1], true);
  }

  std::size_t const prev_index = pline.size() - 2;
  return tangentAlongSegment(pline[prev_index], pline.lastVertex(), false);
}

template <typename Real>
Vector2<Real> openPolylineEndCapCircleCenter(Polyline<Real> const &pline, Real offset, bool atStart,
                                             OffsetEndCapType endCapType) {
  CAVC_ASSERT(!pline.isClosed(), "end cap circle center is only for open polyline");
  Vector2<Real> center = atStart ? pline[0].pos() : pline.lastVertex().pos();
  if (endCapType == OffsetEndCapType::Square) {
    Real extension = std::abs(offset);
    Vector2<Real> tangent = openPolylineEndpointTangent(pline, atStart);
    if (atStart) {
      center = center - extension * tangent;
    } else {
      center = center + extension * tangent;
    }
  }
  return center;
}

template <typename Real>
void applySquareEndCapToRawOpenOffset(Polyline<Real> const &originalPline, Real offset,
                                      Polyline<Real> &rawOffsetPline) {
  CAVC_ASSERT(!originalPline.isClosed(), "square end cap helper is only for open polyline");
  if (rawOffsetPline.size() < 2) {
    return;
  }

  Real extension = std::abs(offset);
  Vector2<Real> start_tangent = openPolylineEndpointTangent(originalPline, true);
  Vector2<Real> end_tangent = openPolylineEndpointTangent(originalPline, false);

  // Preserve arc geometry at endpoints by adding explicit cap line segments.
  if (rawOffsetPline[0].bulgeIsZero()) {
    rawOffsetPline[0].pos() = rawOffsetPline[0].pos() - extension * start_tangent;
  } else {
    PlineVertex<Real> cap_start = rawOffsetPline[0];
    cap_start.pos() = cap_start.pos() - extension * start_tangent;
    cap_start.bulge() = Real(0);
    rawOffsetPline.vertexes().insert(rawOffsetPline.vertexes().begin(), cap_start);
  }

  if (rawOffsetPline[rawOffsetPline.size() - 2].bulgeIsZero()) {
    rawOffsetPline.lastVertex().pos() = rawOffsetPline.lastVertex().pos() + extension * end_tangent;
  } else {
    rawOffsetPline.lastVertex().bulge() = Real(0);
    PlineVertex<Real> cap_end = rawOffsetPline.lastVertex();
    cap_end.pos() = cap_end.pos() + extension * end_tangent;
    cap_end.bulge() = Real(0);
    rawOffsetPline.addVertex(cap_end);
  }
}

template <typename Real>
void lineToArcJoin(PlineOffsetSegment<Real> const &s1, PlineOffsetSegment<Real> const &s2,
                   bool connectionArcsAreCCW, Polyline<Real> &result) {

  const auto &v1 = s1.v1;
  const auto &v2 = s1.v2;
  const auto &u1 = s2.v1;
  const auto &u2 = s2.v2;
  CAVC_ASSERT(v1.bulgeIsZero() && !u1.bulgeIsZero(),
              "first seg should be arc, second seg should be line");

  auto connectUsingArc = [&] {
    auto const &arcCenter = s1.origV2Pos;
    auto const &sp = v2.pos();
    auto const &ep = u1.pos();
    Real bulge = bulgeForConnection(arcCenter, sp, ep, connectionArcsAreCCW);
    // bulge > 1 means arc sweep > 180°; clamp to avoid degenerate geometry
    bulge = std::max(Real(-1), std::min(Real(1), bulge));
    addOrReplaceIfSamePos(result, PlineVertex<Real>(sp, bulge));
    addOrReplaceIfSamePos(result, u1);
  };

  const auto arc = arcRadiusAndCenter(u1, u2);

  auto processIntersect = [&](Real t, Vector2<Real> const &intersect) {
    const bool trueSegIntersect = !falseIntersect(t);
    const bool trueArcIntersect =
        pointWithinArcSweepAngle(arc.center, u1.pos(), u2.pos(), u1.bulge(), intersect,
                                 utils::realPrecision<Real>());
    if (trueSegIntersect && trueArcIntersect) {
      // trim at intersect
      Real a = angle(arc.center, intersect);
      Real arcEndAngle = angle(arc.center, u2.pos());
      Real theta = utils::deltaAngle(a, arcEndAngle);
      // ensure the sign matches (may get flipped if intersect is at the very end of the arc, in
      // which case we do not want to update the bulge)
      if ((theta > Real(0)) == u1.bulgeIsPos()) {
        addOrReplaceIfSamePos(result, PlineVertex<Real>(intersect, std::tan(theta / Real(4))));
      } else {
        addOrReplaceIfSamePos(result, PlineVertex<Real>(intersect, u1.bulge()));
      }
    } else if (t > Real(1) && !trueArcIntersect) {
      connectUsingArc();
    } else if (s1.collapsedArc) {
      // collapsed arc connecting to arc, connect using arc
      connectUsingArc();
    } else {
      // connect using line
      addOrReplaceIfSamePos(result, PlineVertex<Real>(v2.pos(), Real(0)));
      addOrReplaceIfSamePos(result, u1);
    }
  };

  auto intrResult = intrLineSeg2Circle2(v1.pos(), v2.pos(), arc.radius, arc.center);
  if (intrResult.numIntersects == 0) {
    connectUsingArc();
  } else if (intrResult.numIntersects == 1) {
    processIntersect(intrResult.t0, pointFromParametric(v1.pos(), v2.pos(), intrResult.t0));
  } else {
    CAVC_ASSERT(intrResult.numIntersects == 2, "should have 2 intersects here");
    // always use intersect closest to original point
    Vector2<Real> i1 = pointFromParametric(v1.pos(), v2.pos(), intrResult.t0);
    Real dist1 = distSquared(i1, s1.origV2Pos);
    Vector2<Real> i2 = pointFromParametric(v1.pos(), v2.pos(), intrResult.t1);
    Real dist2 = distSquared(i2, s1.origV2Pos);

    if (dist1 < dist2) {
      processIntersect(intrResult.t0, i1);
    } else {
      processIntersect(intrResult.t1, i2);
    }
  }
}

template <typename Real>
void arcToLineJoin(PlineOffsetSegment<Real> const &s1, PlineOffsetSegment<Real> const &s2,
                   bool connectionArcsAreCCW, Polyline<Real> &result) {

  const auto &v1 = s1.v1;
  const auto &v2 = s1.v2;
  const auto &u1 = s2.v1;
  const auto &u2 = s2.v2;
  CAVC_ASSERT(!v1.bulgeIsZero() && u1.bulgeIsZero(),
              "first seg should be line, second seg should be arc");

  auto connectUsingArc = [&] {
    auto const &arcCenter = s1.origV2Pos;
    auto const &sp = v2.pos();
    auto const &ep = u1.pos();
    Real bulge = bulgeForConnection(arcCenter, sp, ep, connectionArcsAreCCW);
    // bulge > 1 means arc sweep > 180°; clamp to avoid degenerate geometry
    bulge = std::max(Real(-1), std::min(Real(1), bulge));
    addOrReplaceIfSamePos(result, PlineVertex<Real>(sp, bulge));
    addOrReplaceIfSamePos(result, u1);
  };

  const auto arc = arcRadiusAndCenter(v1, v2);

  auto processIntersect = [&](Real t, Vector2<Real> const &intersect) {
    const bool trueSegIntersect = !falseIntersect(t);
    const bool trueArcIntersect =
        pointWithinArcSweepAngle(arc.center, v1.pos(), v2.pos(), v1.bulge(), intersect,
                                 utils::realPrecision<Real>());
    if (trueSegIntersect && trueArcIntersect) {
      PlineVertex<Real> &prevVertex = result.lastVertex();

      if (!prevVertex.bulgeIsZero() && !fuzzyEqual(prevVertex.pos(), v2.pos())) {
        // modify previous bulge and trim at intersect
        Real a = angle(arc.center, intersect);
        auto prevArc = arcRadiusAndCenter(prevVertex, v2);
        Real prevArcStartAngle = angle(prevArc.center, prevVertex.pos());
        Real updatedPrevTheta = utils::deltaAngle(prevArcStartAngle, a);

        // ensure the sign matches (may get flipped if intersect is at the very end of the arc, in
        // which case we do not want to update the bulge)
        if ((updatedPrevTheta > Real(0)) == prevVertex.bulgeIsPos()) {
          result.lastVertex().bulge() = std::tan(updatedPrevTheta / Real(4));
        }
      }

      addOrReplaceIfSamePos(result, PlineVertex<Real>(intersect, Real(0)));

    } else {
      connectUsingArc();
    }
  };

  auto intrResult = intrLineSeg2Circle2(u1.pos(), u2.pos(), arc.radius, arc.center);
  if (intrResult.numIntersects == 0) {
    connectUsingArc();
  } else if (intrResult.numIntersects == 1) {
    processIntersect(intrResult.t0, pointFromParametric(u1.pos(), u2.pos(), intrResult.t0));
  } else {
    CAVC_ASSERT(intrResult.numIntersects == 2, "should have 2 intersects here");
    const auto &origPoint = s2.collapsedArc ? u1.pos() : s1.origV2Pos;
    Vector2<Real> i1 = pointFromParametric(u1.pos(), u2.pos(), intrResult.t0);
    Real dist1 = distSquared(i1, origPoint);
    Vector2<Real> i2 = pointFromParametric(u1.pos(), u2.pos(), intrResult.t1);
    Real dist2 = distSquared(i2, origPoint);

    if (dist1 < dist2) {
      processIntersect(intrResult.t0, i1);
    } else {
      processIntersect(intrResult.t1, i2);
    }
  }
}

template <typename Real>
void arcToArcJoin(PlineOffsetSegment<Real> const &s1, PlineOffsetSegment<Real> const &s2,
                  bool connectionArcsAreCCW, Polyline<Real> &result) {

  const auto &v1 = s1.v1;
  const auto &v2 = s1.v2;
  const auto &u1 = s2.v1;
  const auto &u2 = s2.v2;
  CAVC_ASSERT(!v1.bulgeIsZero() && !u1.bulgeIsZero(), "both segs should be arcs");

  const auto arc1 = arcRadiusAndCenter(v1, v2);
  const auto arc2 = arcRadiusAndCenter(u1, u2);

  auto connectUsingArc = [&] {
    auto const &arcCenter = s1.origV2Pos;
    auto const &sp = v2.pos();
    auto const &ep = u1.pos();
    Real bulge = bulgeForConnection(arcCenter, sp, ep, connectionArcsAreCCW);
    // bulge > 1 means arc sweep > 180°; clamp to avoid degenerate geometry
    bulge = std::max(Real(-1), std::min(Real(1), bulge));
    addOrReplaceIfSamePos(result, PlineVertex<Real>(sp, bulge));
    addOrReplaceIfSamePos(result, u1);
  };

  auto processIntersect = [&](Vector2<Real> const &intersect) {
    const bool trueArcIntersect1 =
        pointWithinArcSweepAngle(arc1.center, v1.pos(), v2.pos(), v1.bulge(), intersect,
                                 utils::realPrecision<Real>());
    const bool trueArcIntersect2 =
        pointWithinArcSweepAngle(arc2.center, u1.pos(), u2.pos(), u1.bulge(), intersect,
                                 utils::realPrecision<Real>());

    if (trueArcIntersect1 && trueArcIntersect2) {
      PlineVertex<Real> &prevVertex = result.lastVertex();
      if (!prevVertex.bulgeIsZero() && !fuzzyEqual(prevVertex.pos(), v2.pos())) {
        // modify previous bulge and trim at intersect
        Real a1 = angle(arc1.center, intersect);
        auto prevArc = arcRadiusAndCenter(prevVertex, v2);
        Real prevArcStartAngle = angle(prevArc.center, prevVertex.pos());
        Real updatedPrevTheta = utils::deltaAngle(prevArcStartAngle, a1);

        // ensure the sign matches (may get flipped if intersect is at the very end of the arc, in
        // which case we do not want to update the bulge)
        if ((updatedPrevTheta > Real(0)) == prevVertex.bulgeIsPos()) {
          result.lastVertex().bulge() = std::tan(updatedPrevTheta / Real(4));
        }
      }

      // add the vertex at our current trim/join point
      Real a2 = angle(arc2.center, intersect);
      Real endAngle = angle(arc2.center, u2.pos());
      Real theta = utils::deltaAngle(a2, endAngle);

      // ensure the sign matches (may get flipped if intersect is at the very end of the arc, in
      // which case we do not want to update the bulge)
      if ((theta > Real(0)) == u1.bulgeIsPos()) {
        addOrReplaceIfSamePos(result, PlineVertex<Real>(intersect, std::tan(theta / Real(4))));
      } else {
        addOrReplaceIfSamePos(result, PlineVertex<Real>(intersect, u1.bulge()));
      }

    } else {
      connectUsingArc();
    }
  };

  const auto intrResult = intrCircle2Circle2(arc1.radius, arc1.center, arc2.radius, arc2.center);
  switch (intrResult.intrType) {
  case Circle2Circle2IntrType::NoIntersect:
    connectUsingArc();
    break;
  case Circle2Circle2IntrType::OneIntersect:
    processIntersect(intrResult.point1);
    break;
  case Circle2Circle2IntrType::TwoIntersects: {
    Real dist1 = distSquared(intrResult.point1, s1.origV2Pos);
    Real dist2 = distSquared(intrResult.point2, s1.origV2Pos);
    if (dist1 < dist2) {
      processIntersect(intrResult.point1);
    } else {
      processIntersect(intrResult.point2);
    }
  } break;
  case Circle2Circle2IntrType::Coincident:
    // same constant arc radius and center, just add the vertex (nothing to trim/extend)
    addOrReplaceIfSamePos(result, u1);
    break;
  }
}

template <typename Real>
void offsetCircleIntersectsWithPline(Polyline<Real> const &pline, Real offset,
                                     Vector2<Real> const &circleCenter,
                                     StaticSpatialIndex<Real> const &spatialIndex,
                                     std::vector<std::pair<std::size_t, Vector2<Real>>> &output,
                                     std::vector<std::size_t> &queryResults,
                                     std::vector<std::size_t> &queryStack) {

  const Real circleRadius = std::abs(offset);

  queryResults.clear();
  spatialIndex.query(circleCenter.x() - circleRadius, circleCenter.y() - circleRadius,
                     circleCenter.x() + circleRadius, circleCenter.y() + circleRadius, queryResults,
                     queryStack);

  auto validLineSegIntersect = [](Real t) {
    return !falseIntersect(t) && std::abs(t) > utils::realPrecision<Real>();
  };

  auto validArcSegIntersect = [](Vector2<Real> const &arcCenter, Vector2<Real> const &arcStart,
                                 Vector2<Real> const &arcEnd, Real bulge,
                                 Vector2<Real> const &intrPoint) {
    return !fuzzyEqual(arcStart, intrPoint, utils::realPrecision<Real>()) &&
           pointWithinArcSweepAngle(arcCenter, arcStart, arcEnd, bulge, intrPoint,
                                    utils::realPrecision<Real>());
  };

  for (std::size_t sIndex : queryResults) {
    PlineVertex<Real> const &v1 = pline[sIndex];
    PlineVertex<Real> const &v2 = pline[sIndex + 1];
    if (v1.bulgeIsZero()) {
      IntrLineSeg2Circle2Result<Real> intrResult =
          intrLineSeg2Circle2(v1.pos(), v2.pos(), circleRadius, circleCenter);
      if (intrResult.numIntersects == 0) {
        continue;
      } else if (intrResult.numIntersects == 1) {
        if (validLineSegIntersect(intrResult.t0)) {
          output.emplace_back(sIndex, pointFromParametric(v1.pos(), v2.pos(), intrResult.t0));
        }
      } else {
        CAVC_ASSERT(intrResult.numIntersects == 2, "should be two intersects here");
        if (validLineSegIntersect(intrResult.t0)) {
          output.emplace_back(sIndex, pointFromParametric(v1.pos(), v2.pos(), intrResult.t0));
        }
        if (validLineSegIntersect(intrResult.t1)) {
          output.emplace_back(sIndex, pointFromParametric(v1.pos(), v2.pos(), intrResult.t1));
        }
      }
    } else {
      auto arc = arcRadiusAndCenter(v1, v2);
      IntrCircle2Circle2Result<Real> intrResult =
          intrCircle2Circle2(arc.radius, arc.center, circleRadius, circleCenter);
      switch (intrResult.intrType) {
      case Circle2Circle2IntrType::NoIntersect:
        break;
      case Circle2Circle2IntrType::OneIntersect:
        if (validArcSegIntersect(arc.center, v1.pos(), v2.pos(), v1.bulge(), intrResult.point1)) {
          output.emplace_back(sIndex, intrResult.point1);
        }
        break;
      case Circle2Circle2IntrType::TwoIntersects:
        if (validArcSegIntersect(arc.center, v1.pos(), v2.pos(), v1.bulge(), intrResult.point1)) {
          output.emplace_back(sIndex, intrResult.point1);
        }
        if (validArcSegIntersect(arc.center, v1.pos(), v2.pos(), v1.bulge(), intrResult.point2)) {
          output.emplace_back(sIndex, intrResult.point2);
        }
        break;
      case Circle2Circle2IntrType::Coincident:
        break;
      }
    }
  }
}

template <typename Real>
void offsetLineIntersectsWithPline(Polyline<Real> const &pline, Vector2<Real> const &linePoint,
                                   Vector2<Real> const &lineNormal,
                                   std::vector<std::pair<std::size_t, Vector2<Real>>> &output) {
  if (pline.size() < 2) {
    return;
  }

  Vector2<Real> normalizedNormal = lineNormal;
  CAVC_ASSERT(!fuzzyZero(normalizedNormal), "line normal must be non-zero");
  normalize(normalizedNormal);
  Vector2<Real> lineDirection = unitPerp(normalizedNormal);

  auto validLineSegIntersect = [](Real t) {
    return !falseIntersect(t) && std::abs(t) > utils::realPrecision<Real>();
  };

  auto validArcSegIntersect = [](Vector2<Real> const &arcCenter, Vector2<Real> const &arcStart,
                                 Vector2<Real> const &arcEnd, Real bulge,
                                 Vector2<Real> const &intrPoint) {
    return !fuzzyEqual(arcStart, intrPoint, utils::realPrecision<Real>()) &&
           pointWithinArcSweepAngle(arcCenter, arcStart, arcEnd, bulge, intrPoint,
                                    utils::realPrecision<Real>());
  };

  std::size_t const segCount = pline.isClosed() ? pline.size() : pline.size() - 1;
  for (std::size_t sIndex = 0; sIndex < segCount; ++sIndex) {
    std::size_t const nextIndex =
        pline.isClosed() ? utils::nextWrappingIndex(sIndex, pline) : sIndex + 1;
    PlineVertex<Real> const &v1 = pline[sIndex];
    PlineVertex<Real> const &v2 = pline[nextIndex];

    if (v1.bulgeIsZero()) {
      Vector2<Real> segDelta = v2.pos() - v1.pos();
      Real const denom = dot(segDelta, normalizedNormal);
      if (std::abs(denom) <= utils::realThreshold<Real>()) {
        continue;
      }

      Real const t = dot(linePoint - v1.pos(), normalizedNormal) / denom;
      if (validLineSegIntersect(t)) {
        output.emplace_back(sIndex, pointFromParametric(v1.pos(), v2.pos(), t));
      }
      continue;
    }

    auto arc = arcRadiusAndCenter(v1, v2);
    Real signedDistance = dot(arc.center - linePoint, normalizedNormal);
    Real absDistance = std::abs(signedDistance);
    if (absDistance > arc.radius + utils::realPrecision<Real>()) {
      continue;
    }

    Vector2<Real> projectedPoint = arc.center - signedDistance * normalizedNormal;

    auto addArcIntersectIfValid = [&](Vector2<Real> const &intersectPoint) {
      if (validArcSegIntersect(arc.center, v1.pos(), v2.pos(), v1.bulge(), intersectPoint)) {
        output.emplace_back(sIndex, intersectPoint);
      }
    };

    if (utils::fuzzyEqual(absDistance, arc.radius, utils::realPrecision<Real>())) {
      addArcIntersectIfValid(projectedPoint);
      continue;
    }

    Real offsetSquared = arc.radius * arc.radius - absDistance * absDistance;
    if (offsetSquared < Real(0)) {
      offsetSquared = Real(0);
    }
    Real offsetAlongLine = std::sqrt(offsetSquared);
    addArcIntersectIfValid(projectedPoint + offsetAlongLine * lineDirection);
    addArcIntersectIfValid(projectedPoint - offsetAlongLine * lineDirection);
  }
}

/// Function to test if a point is a valid distance from the original polyline.
template <typename Real, std::size_t N>
bool pointValidForOffset(Polyline<Real> const &pline, Real offset,
                         StaticSpatialIndex<Real, N> const &spatialIndex,
                         Vector2<Real> const &point, std::vector<std::size_t> &queryStack,
                         Real offsetTol = utils::offsetThreshold<Real>()) {
  const Real absOffset = std::abs(offset) - offsetTol;
  const Real minDist = absOffset * absOffset;

  bool pointValid = true;

  auto visitor = [&](std::size_t i) {
    std::size_t j = utils::nextWrappingIndex(i, pline.vertexes());
    auto closestPoint = closestPointOnSeg(pline[i], pline[j], point, utils::realPrecision<Real>());
    Real dist = distSquared(closestPoint, point);
    pointValid = dist > minDist;
    return pointValid;
  };

  spatialIndex.visitQuery(point.x() - absOffset, point.y() - absOffset, point.x() + absOffset,
                          point.y() + absOffset, visitor, queryStack);
  return pointValid;
}

template <typename Real>
Real slicePathLength(Polyline<Real> const &source, PlineSliceViewData<Real> const &slice) {
  Real result = Real(0);
  slice.visitSegments(source, [&](PlineVertex<Real> const &v1, PlineVertex<Real> const &v2) {
    result += segLength(v1, v2);
    return true;
  });
  return result;
}

template <typename Real> struct PointPairHash {
  std::size_t operator()(std::pair<Real, Real> const &value) const noexcept {
    std::size_t seed = std::hash<Real>{}(value.first);
    std::size_t hash2 = std::hash<Real>{}(value.second);
    return seed ^ (hash2 + 0x9e3779b9u + (seed << 6) + (seed >> 2));
  }
};

template <typename Real, typename PointValidF, typename PointValidAtRawIndexF,
          typename SegmentIntersectsF, typename RawSegmentIntersectsF>
bool offsetSliceIsValid(Polyline<Real> const &rawOffsetPline, PlineSliceViewData<Real> const &slice,
                        PointValidF &&pointValid, PointValidAtRawIndexF &&pointValidAtRawIndex,
                        SegmentIntersectsF &&segmentIntersectsOrig,
                        RawSegmentIntersectsF &&rawSegmentIntersectsOrig) {
  if (slice.endIndexOffset == 0) {
    PlineVertex<Real> const &sliceStart = slice.updatedStart;
    PlineVertex<Real> sliceEnd(slice.endPoint, Real(0));
    if (!pointValid(sliceStart.pos()) || !pointValid(sliceEnd.pos())) {
      return false;
    }

    Vector2<Real> midpoint = segMidpoint(sliceStart, sliceEnd);
    return pointValid(midpoint) && !segmentIntersectsOrig(sliceStart, sliceEnd);
  }

  std::size_t const nextIndex = utils::nextWrappingIndex(slice.startIndex, rawOffsetPline);
  Vector2<Real> startSegMidpoint = segMidpoint(slice.updatedStart, rawOffsetPline[nextIndex]);
  if (!pointValid(startSegMidpoint)) {
    return false;
  }

  std::size_t endIndex = slice.startIndex + slice.endIndexOffset;
  if (endIndex >= rawOffsetPline.size()) {
    endIndex -= rawOffsetPline.size();
  }
  PlineVertex<Real> endSegStart = rawOffsetPline[endIndex];
  endSegStart.bulge() = slice.updatedEndBulge;
  Vector2<Real> endSegMidpoint = segMidpoint(endSegStart, PlineVertex<Real>(slice.endPoint, Real(0)));
  if (!pointValid(endSegMidpoint)) {
    return false;
  }

  Real const posEqualEps = utils::realPrecision<Real>();
  bool const startIsRawVertex =
      fuzzyEqual(slice.updatedStart.pos(), rawOffsetPline[slice.startIndex].pos(), posEqualEps);
  bool const startPointValid =
      startIsRawVertex ? pointValidAtRawIndex(slice.startIndex) : pointValid(slice.updatedStart.pos());
  bool const startSegIntersectsOrig =
      startIsRawVertex && utils::fuzzyEqual(slice.updatedStart.bulge(),
                                            rawOffsetPline[slice.startIndex].bulge(), posEqualEps)
          ? rawSegmentIntersectsOrig(slice.startIndex)
          : segmentIntersectsOrig(slice.updatedStart, rawOffsetPline[nextIndex]);

  if (!startPointValid || startSegIntersectsOrig) {
    return false;
  }

  std::size_t segIndex = nextIndex;
  for (std::size_t i = 1; i < slice.endIndexOffset; ++i) {
    if (!pointValidAtRawIndex(segIndex) || rawSegmentIntersectsOrig(segIndex)) {
      return false;
    }
    segIndex = utils::nextWrappingIndex(segIndex, rawOffsetPline);
  }

  if (!pointValidAtRawIndex(endIndex) || segmentIntersectsOrig(endSegStart, {slice.endPoint, Real(0)})) {
    return false;
  }

  return pointValid(slice.endPoint);
}

/// Creates the raw offset polyline.
template <typename Real>
Polyline<Real> createRawOffsetPline(Polyline<Real> const &pline, Real offset,
                                    ParallelOffsetOptions<Real> const &options) {

  Polyline<Real> result;
  if (pline.size() < 2) {
    return result;
  }

  std::vector<PlineOffsetSegment<Real>> rawOffsets = createUntrimmedOffsetSegments(pline, offset);
  if (rawOffsets.size() == 0) {
    return result;
  }

  // detect single collapsed arc segment (this may be removed in the future if invalid segments are
  // tracked in join functions to be pruned at slice creation)
  if (rawOffsets.size() == 1 && rawOffsets[0].collapsedArc) {
    return result;
  }

  result.vertexes().reserve(pline.size());
  result.isClosed() = pline.isClosed();

  const bool connectionArcsAreCCW = offset < Real(0);

  auto joinResultVisitor = [&options, connectionArcsAreCCW](PlineOffsetSegment<Real> const &s1,
                                                            PlineOffsetSegment<Real> const &s2,
                                                            Polyline<Real> &p_result) {
    const bool s1IsLine = s1.v1.bulgeIsZero();
    const bool s2IsLine = s2.v1.bulgeIsZero();
    if (s1IsLine && s2IsLine) {
      if (options.joinType == OffsetJoinType::Round || (!s1.collapsedArc && !s2.collapsedArc)) {
        internal::lineToLineJoin(s1, s2, connectionArcsAreCCW, options.joinType, options.miterLimit,
                                 p_result);
      } else {
        internal::nonRoundArcJoin(s1, s2, options.joinType, options.miterLimit, p_result);
      }
    } else if (s1IsLine) {
      if (options.joinType == OffsetJoinType::Round) {
        internal::lineToArcJoin(s1, s2, connectionArcsAreCCW, p_result);
      } else {
        internal::nonRoundArcJoin(s1, s2, options.joinType, options.miterLimit, p_result);
      }
    } else if (s2IsLine) {
      if (options.joinType == OffsetJoinType::Round) {
        internal::arcToLineJoin(s1, s2, connectionArcsAreCCW, p_result);
      } else {
        internal::nonRoundArcJoin(s1, s2, options.joinType, options.miterLimit, p_result);
      }
    } else {
      if (options.joinType == OffsetJoinType::Round) {
        internal::arcToArcJoin(s1, s2, connectionArcsAreCCW, p_result);
      } else {
        internal::nonRoundArcJoin(s1, s2, options.joinType, options.miterLimit, p_result);
      }
    }
  };

  result.addVertex(rawOffsets[0].v1);

  // join first two segments and determine if first vertex was replaced (to know how to handle last
  // two segment joins for closed polyline)
  if (rawOffsets.size() > 1) {

    auto const &seg01 = rawOffsets[0];
    auto const &seg12 = rawOffsets[1];
    joinResultVisitor(seg01, seg12, result);
  }
  const bool firstVertexReplaced = result.size() == 1;

  for (std::size_t i = 2; i < rawOffsets.size(); ++i) {
    const auto &seg1 = rawOffsets[i - 1];
    const auto &seg2 = rawOffsets[i];
    joinResultVisitor(seg1, seg2, result);
  }

  if (pline.isClosed() && result.size() > 1) {
    // joining segments at vertex indexes (n, 0) and (0, 1)
    const auto &s1 = rawOffsets.back();
    const auto &s2 = rawOffsets[0];

    // temp polyline to capture results of joining (to avoid mutating result)
    Polyline<Real> closingPartResult;
    closingPartResult.addVertex(result.lastVertex());
    joinResultVisitor(s1, s2, closingPartResult);

    // update last vertexes
    result.lastVertex() = closingPartResult[0];
    for (std::size_t i = 1; i < closingPartResult.size(); ++i) {
      result.addVertex(closingPartResult[i]);
    }
    result.vertexes().pop_back();

    // update first vertex (only if it has not already been updated/replaced)
    if (!firstVertexReplaced) {
      const Vector2<Real> &updatedFirstPos = closingPartResult.lastVertex().pos();
      if (result[0].bulgeIsZero()) {
        // just update position
        result[0].pos() = updatedFirstPos;
      } else if (result.size() > 1) {
        // update position and bulge
        const auto arc = arcRadiusAndCenter(result[0], result[1]);
        const Real a1 = angle(arc.center, updatedFirstPos);
        const Real a2 = angle(arc.center, result[1].pos());
        const Real updatedTheta = utils::deltaAngle(a1, a2);
        if ((updatedTheta < Real(0) && result[0].bulgeIsPos()) ||
            (updatedTheta > Real(0) && result[0].bulgeIsNeg())) {
          // first vertex not valid, just update its position to be removed later
          result[0].pos() = updatedFirstPos;
        } else {
          // update position and bulge
          result[0].pos() = updatedFirstPos;
          result[0].bulge() = std::tan(updatedTheta / Real(4));
        }
      }
    }

    // must do final singularity prune between first and second vertex after joining curves (n, 0)
    // and (0, 1)
    if (result.size() > 1) {
      if (fuzzyEqual(result[0].pos(), result[1].pos(), utils::realPrecision<Real>())) {
        result.vertexes().erase(result.vertexes().begin());
      }
    }
  } else {
    internal::addOrReplaceIfSamePos(result, rawOffsets.back().v2);
  }

  if (!pline.isClosed() && options.endCapType == OffsetEndCapType::Square && result.size() > 1) {
    internal::applySquareEndCapToRawOpenOffset(pline, offset, result);
  }

  // if due to joining of segments we are left with only 1 vertex then return no raw offset (empty
  // polyline)
  if (result.size() == 1) {
    result.vertexes().clear();
  }

  return result;
}

template <typename Real> using OpenPolylineSlice = OffsetSliceRef<Real>;

/// Slices a raw offset polyline at all of its self intersects.
template <typename Real>
std::vector<OpenPolylineSlice<Real>>
slicesFromRawOffset(Polyline<Real> const &originalPline, Polyline<Real> const &rawOffsetPline,
                    Real offset, bool enforceMinDistance = true,
                    bool skipOrigIntersectionCheck = false) {
  CAVC_ASSERT(originalPline.isClosed(), "use dual slice at intersects for open polylines");

  std::vector<OpenPolylineSlice<Real>> result;
  if (rawOffsetPline.size() < 2) {
    return result;
  }

  StaticSpatialIndex<Real> origPlineSpatialIndex = createApproxSpatialIndex(originalPline);
  StaticSpatialIndex<Real> rawOffsetPlineSpatialIndex = createApproxSpatialIndex(rawOffsetPline);

  std::vector<PlineIntersect<Real>> selfIntersects;
  allSelfIntersects(rawOffsetPline, selfIntersects, rawOffsetPlineSpatialIndex);

  std::vector<std::size_t> queryStack;
  queryStack.reserve(8);
  std::unordered_map<std::pair<Real, Real>, std::uint8_t, PointPairHash<Real>>
      intersectPointValidCache;
  if (selfIntersects.size() != 0) {
    intersectPointValidCache.reserve(2 * selfIntersects.size());
  }
  std::vector<std::int8_t> rawVertexPointValidCache(rawOffsetPline.size(), -1);
  std::vector<std::int8_t> rawSegmentIntersectsOrigCache(rawOffsetPline.size(), -1);
  auto pointValid = [&](Vector2<Real> const &p) {
    if (!enforceMinDistance) {
      return true;
    }

    return pointValidForOffset(originalPline, offset, origPlineSpatialIndex, p, queryStack);
  };
  auto intersectPointValid = [&](Vector2<Real> const &p) {
    auto const key = std::make_pair(p.x(), p.y());
    auto const iter = intersectPointValidCache.find(key);
    if (iter != intersectPointValidCache.end()) {
      return iter->second != 0;
    }

    bool const valid = pointValid(p);
    intersectPointValidCache.emplace(key, valid ? 1 : 0);
    return valid;
  };
  auto pointValidAtRawIndex = [&](std::size_t index) {
    std::int8_t &cached = rawVertexPointValidCache[index];
    if (cached != -1) {
      return cached != 0;
    }

    bool const valid = pointValid(rawOffsetPline[index].pos());
    cached = valid ? 1 : 0;
    return valid;
  };

  if (selfIntersects.size() == 0) {
    if (!pointValidAtRawIndex(0)) {
      return result;
    }
    result.push_back(
        {std::numeric_limits<std::size_t>::max(), PlineSliceViewData<Real>::fromEntirePolyline(rawOffsetPline)});
    return result;
  }

  // using unordered_map rather than map for performance (as is used in
  // dualSliceAtIntersectsForOffset) since all slices will stitch together to form closed
  // loops/polylines so later when slices are stitched together the order that slices are visited
  // does not matter
  std::unordered_map<std::size_t, std::vector<Vector2<Real>>> intersectsLookup;
  intersectsLookup.reserve(2 * selfIntersects.size());

  for (PlineIntersect<Real> const &si : selfIntersects) {
    intersectsLookup[si.sIndex1].push_back(si.pos);
    intersectsLookup[si.sIndex2].push_back(si.pos);
  }

  // sort intersects by distance from start vertex
  for (auto &kvp : intersectsLookup) {
    Vector2<Real> startPos = rawOffsetPline[kvp.first].pos();
    auto cmp = [&](Vector2<Real> const &si1, Vector2<Real> const &si2) {
      return distSquared(si1, startPos) < distSquared(si2, startPos);
    };
    std::sort(kvp.second.begin(), kvp.second.end(), cmp);
  }

  auto intersectsOrigPline = [&](PlineVertex<Real> const &v1, PlineVertex<Real> const &v2) {
    if (skipOrigIntersectionCheck) {
      return false;
    }

    AABB<Real> approxBB = createFastApproxBoundingBox(v1, v2);
    bool hasIntersect = false;
    auto visitor = [&](std::size_t i) {
      using namespace internal;
      std::size_t j = utils::nextWrappingIndex(i, originalPline);
      IntrPlineSegsResult<Real> intrResult =
          intrPlineSegs(v1, v2, originalPline[i], originalPline[j]);
      hasIntersect = intrResult.intrType != PlineSegIntrType::NoIntersect;
      return !hasIntersect;
    };

    origPlineSpatialIndex.visitQuery(approxBB.xMin, approxBB.yMin, approxBB.xMax, approxBB.yMax,
                                     visitor, queryStack);

    return hasIntersect;
  };
  auto rawSegmentIntersectsOrig = [&](std::size_t index) {
    std::int8_t &cached = rawSegmentIntersectsOrigCache[index];
    if (cached != -1) {
      return cached != 0;
    }

    std::size_t const nextIndex = utils::nextWrappingIndex(index, rawOffsetPline);
    bool const intersects = intersectsOrigPline(rawOffsetPline[index], rawOffsetPline[nextIndex]);
    cached = intersects ? 1 : 0;
    return intersects;
  };
  auto rawVertexSegmentIntersectsOrig = [&](PlineVertex<Real> const &v1, std::size_t endIndex) {
    std::size_t prevIndex = utils::prevWrappingIndex(endIndex, rawOffsetPline);
    if (fuzzyEqual(v1.pos(), rawOffsetPline[prevIndex].pos(), utils::realPrecision<Real>()) &&
        utils::fuzzyEqual(v1.bulge(), rawOffsetPline[prevIndex].bulge(),
                          utils::realPrecision<Real>())) {
      return rawSegmentIntersectsOrig(prevIndex);
    }

    return intersectsOrigPline(v1, rawOffsetPline[endIndex]);
  };

  for (auto const &kvp : intersectsLookup) {
    // start index for the slice we're about to build
    std::size_t sIndex = kvp.first;
    // self intersect list for this start index
    std::vector<Vector2<Real>> const &siList = kvp.second;

    const auto &startVertex = rawOffsetPline[sIndex];
    std::size_t nextIndex = utils::nextWrappingIndex(sIndex, rawOffsetPline);
    const auto &endVertex = rawOffsetPline[nextIndex];

    if (siList.size() != 1) {
      // build all the segments between the N intersects in siList (N > 1), skipping the first
      // segment (to be processed at the end)
      SplitResult<Real> firstSplit = splitAtPoint(startVertex, endVertex, siList[0]);
      auto prevVertex = firstSplit.splitVertex;
      for (std::size_t i = 1; i < siList.size(); ++i) {
        SplitResult<Real> split = splitAtPoint(prevVertex, endVertex, siList[i]);
        // update prevVertex for next loop iteration
        prevVertex = split.splitVertex;
        // skip if they're ontop of each other
        if (fuzzyEqual(split.updatedStart.pos(), split.splitVertex.pos(),
                       utils::realPrecision<Real>())) {
          continue;
        }

        // test start point
        if (!pointValid(split.updatedStart.pos())) {
          continue;
        }

        // test end point
        if (!pointValid(split.splitVertex.pos())) {
          continue;
        }

        // test mid point
        auto midpoint = segMidpoint(split.updatedStart, split.splitVertex);
        if (!pointValid(midpoint)) {
          continue;
        }

        // test intersection with original polyline
        if (intersectsOrigPline(split.updatedStart, split.splitVertex)) {
          continue;
        }

        auto slice = PlineSliceViewData<Real>::createOnSingleSegment(
            rawOffsetPline, sIndex, split.updatedStart, split.splitVertex.pos());
        if (slice) {
          result.push_back({sIndex, *slice});
        }
      }
    }

    std::size_t index = nextIndex;
    bool isValidSlice = intersectPointValid(siList.back());
    SplitResult<Real> sliceStartSplit = splitAtPoint(startVertex, endVertex, siList.back());
    PlineVertex<Real> currLastVertex = sliceStartSplit.splitVertex;
    std::size_t vertexCount = 1;
    std::size_t loopCount = 0;
    const std::size_t maxLoopCount = rawOffsetPline.size();
    while (true) {
      if (loopCount++ > maxLoopCount) {
        CAVC_ASSERT(false, "Bug detected, should never loop this many times!");
        break;
      }
      if (!pointValidAtRawIndex(index) || rawVertexSegmentIntersectsOrig(currLastVertex, index)) {
        isValidSlice = false;
        break;
      }

      if (fuzzyEqual(currLastVertex.pos(), rawOffsetPline[index].pos(),
                     utils::realPrecision<Real>())) {
        currLastVertex.bulge() = rawOffsetPline[index].bulge();
      } else {
        currLastVertex = rawOffsetPline[index];
        vertexCount += 1;
      }

      auto nextIntr = intersectsLookup.find(index);
      if (nextIntr != intersectsLookup.end()) {
        Vector2<Real> const &intersectPos = nextIntr->second[0];
        if (!intersectPointValid(intersectPos)) {
          isValidSlice = false;
          break;
        }

        std::size_t l_nextIndex = utils::nextWrappingIndex(index, rawOffsetPline);
        SplitResult<Real> l_split =
            splitAtPoint(currLastVertex, rawOffsetPline[l_nextIndex], intersectPos);
        PlineVertex<Real> sliceEndVertex(intersectPos, Real(0));
        if (!pointValid(segMidpoint(l_split.updatedStart, sliceEndVertex))) {
          isValidSlice = false;
          break;
        }

        auto slice = PlineSliceViewData<Real>::createFromSlicePoints(
            rawOffsetPline, siList.back(), sIndex, intersectPos, index);
        if (isValidSlice && slice) {
          if (vertexCount > 1 &&
              fuzzyEqual(slice->firstPoint(rawOffsetPline), slice->lastPoint(rawOffsetPline),
                         utils::realPrecision<Real>()) &&
              slicePathLength(rawOffsetPline, *slice) <= Real(1e-2)) {
            isValidSlice = false;
          }

          if (isValidSlice) {
            result.push_back({sIndex, *slice});
          }
        }
        break;
      }
      index = utils::nextWrappingIndex(index, rawOffsetPline);
    }
  }

  return result;
}

/// Slices a raw offset polyline at all of its self intersects and intersects with its dual.
template <typename Real>
std::vector<OpenPolylineSlice<Real>>
dualSliceAtIntersectsForOffset(Polyline<Real> const &originalPline,
                               Polyline<Real> const &rawOffsetPline,
                               Polyline<Real> const &dualRawOffsetPline, Real offset,
                               ParallelOffsetOptions<Real> const &options,
                               bool enforceMinDistance = true,
                               bool skipOrigIntersectionCheck = false) {
  std::vector<OpenPolylineSlice<Real>> result;
  if (rawOffsetPline.size() < 2) {
    return result;
  }

  StaticSpatialIndex<Real> origPlineSpatialIndex = createApproxSpatialIndex(originalPline);
  StaticSpatialIndex<Real> rawOffsetPlineSpatialIndex = createApproxSpatialIndex(rawOffsetPline);

  std::vector<PlineIntersect<Real>> selfIntersects;
  allSelfIntersects(rawOffsetPline, selfIntersects, rawOffsetPlineSpatialIndex);

  PlineIntersectsResult<Real> dualIntersects;
  findIntersects(rawOffsetPline, dualRawOffsetPline, rawOffsetPlineSpatialIndex, dualIntersects);

  // using map rather than unordered map since we want to construct the slices in vertex index order
  // and we do so by looping through all intersects (required later when slices are stitched
  // together, because slices may not all form closed loops/polylines so must go in order of
  // indexes to ensure longest sitched results are formed)
  std::map<std::size_t, std::vector<Vector2<Real>>> intersectsLookup;
  auto addIntersect = [&](std::size_t sIndex, Vector2<Real> const &pos) {
    intersectsLookup[sIndex].push_back(pos);
  };
  std::vector<std::size_t> queryStack;
  queryStack.reserve(8);
  std::unordered_map<std::pair<Real, Real>, std::uint8_t, PointPairHash<Real>>
      intersectPointValidCache;
  std::vector<std::int8_t> rawVertexPointValidCache(rawOffsetPline.size(), -1);
  std::vector<std::int8_t> rawSegmentIntersectsOrigCache(rawOffsetPline.size(), -1);
  auto pointValid = [&](Vector2<Real> const &p) {
    if (!enforceMinDistance) {
      return true;
    }
    if (!originalPline.isClosed()) {
      auto const key = std::make_pair(p.x(), p.y());
      auto const iter = intersectPointValidCache.find(key);
      if (iter != intersectPointValidCache.end()) {
        return iter->second != 0;
      }

      bool const valid =
          pointValidForOffset(originalPline, offset, origPlineSpatialIndex, p, queryStack);
      intersectPointValidCache.emplace(key, valid ? 1 : 0);
      return valid;
    }

    return pointValidForOffset(originalPline, offset, origPlineSpatialIndex, p, queryStack);
  };
  auto pointValidAtRawIndex = [&](std::size_t index) {
    std::int8_t &cached = rawVertexPointValidCache[index];
    if (cached != -1) {
      return cached != 0;
    }

    bool const valid = pointValid(rawOffsetPline[index].pos());
    cached = valid ? 1 : 0;
    return valid;
  };
  if (!originalPline.isClosed()) {
    intersectPointValidCache.reserve(2 * selfIntersects.size() + dualIntersects.intersects.size() +
                                     2 * dualIntersects.coincidentIntersects.size() +
                                     rawOffsetPline.size());
  }

  if (!originalPline.isClosed()) {
    // find intersects between raw offset polyline and end-cap clip geometry
    std::vector<std::pair<std::size_t, Vector2<Real>>> intersects;
    intersects.reserve(rawOffsetPline.size());
    if (options.endCapType == OffsetEndCapType::Butt) {
      Vector2<Real> start_tangent = internal::openPolylineEndpointTangent(originalPline, true);
      Vector2<Real> end_tangent = internal::openPolylineEndpointTangent(originalPline, false);
      internal::offsetLineIntersectsWithPline(rawOffsetPline, originalPline[0].pos(), start_tangent,
                                              intersects);
      internal::offsetLineIntersectsWithPline(rawOffsetPline, originalPline.lastVertex().pos(),
                                              end_tangent, intersects);
    } else {
      Vector2<Real> start_circle_center =
          internal::openPolylineEndCapCircleCenter(originalPline, offset, true, options.endCapType);
      Vector2<Real> end_circle_center = internal::openPolylineEndCapCircleCenter(
          originalPline, offset, false, options.endCapType);
      std::vector<std::size_t> circleQueryResults;
      circleQueryResults.reserve(rawOffsetPline.size());
      internal::offsetCircleIntersectsWithPline(rawOffsetPline, offset, start_circle_center,
                                                rawOffsetPlineSpatialIndex, intersects,
                                                circleQueryResults, queryStack);
      internal::offsetCircleIntersectsWithPline(rawOffsetPline, offset, end_circle_center,
                                                rawOffsetPlineSpatialIndex, intersects,
                                                circleQueryResults, queryStack);
    }
    for (auto const &pair : intersects) {
      addIntersect(pair.first, pair.second);
    }
  }

  for (PlineIntersect<Real> const &si : selfIntersects) {
    addIntersect(si.sIndex1, si.pos);
    addIntersect(si.sIndex2, si.pos);
  }

  for (PlineIntersect<Real> const &intr : dualIntersects.intersects) {
    addIntersect(intr.sIndex1, intr.pos);
  }

  for (PlineCoincidentIntersect<Real> const &intr : dualIntersects.coincidentIntersects) {
    addIntersect(intr.sIndex1, intr.point1);
    addIntersect(intr.sIndex1, intr.point2);
  }

  if (intersectsLookup.size() == 0) {
    if (!pointValidAtRawIndex(0)) {
      return result;
    }
    result.push_back(
        {std::numeric_limits<std::size_t>::max(), PlineSliceViewData<Real>::fromEntirePolyline(rawOffsetPline)});
    return result;
  }

  // sort intersects by distance from start vertex
  for (auto &kvp : intersectsLookup) {
    Vector2<Real> startPos = rawOffsetPline[kvp.first].pos();
    auto cmp = [&](Vector2<Real> const &si1, Vector2<Real> const &si2) {
      return distSquared(si1, startPos) < distSquared(si2, startPos);
    };
    std::sort(kvp.second.begin(), kvp.second.end(), cmp);
  }

  auto intersectsOrigPline = [&](PlineVertex<Real> const &v1, PlineVertex<Real> const &v2) {
    if (skipOrigIntersectionCheck) {
      return false;
    }

    AABB<Real> approxBB = createFastApproxBoundingBox(v1, v2);
    bool intersects = false;
    auto visitor = [&](std::size_t i) {
      using namespace internal;
      std::size_t j = utils::nextWrappingIndex(i, originalPline);
      IntrPlineSegsResult<Real> intrResult =
          intrPlineSegs(v1, v2, originalPline[i], originalPline[j]);
      intersects = intrResult.intrType != PlineSegIntrType::NoIntersect;
      return !intersects;
    };

    origPlineSpatialIndex.visitQuery(approxBB.xMin, approxBB.yMin, approxBB.xMax, approxBB.yMax,
                                     visitor, queryStack);

    return intersects;
  };
  auto rawSegmentIntersectsOrig = [&](std::size_t index) {
    std::int8_t &cached = rawSegmentIntersectsOrigCache[index];
    if (cached != -1) {
      return cached != 0;
    }

    std::size_t const nextIndex = utils::nextWrappingIndex(index, rawOffsetPline);
    bool const intersects = intersectsOrigPline(rawOffsetPline[index], rawOffsetPline[nextIndex]);
    cached = intersects ? 1 : 0;
    return intersects;
  };

  auto sliceIsValid = [&](PlineSliceViewData<Real> const &slice) {
    return offsetSliceIsValid(rawOffsetPline, slice, pointValid, pointValidAtRawIndex,
                              intersectsOrigPline, rawSegmentIntersectsOrig);
  };
  auto maybeAppendSlice = [&](std::size_t sIndex,
                              std::optional<PlineSliceViewData<Real>> const &slice) {
    if (!slice || !sliceIsValid(*slice)) {
      return;
    }

    result.push_back({sIndex, *slice});
  };

  if (!originalPline.isClosed()) {
    // build first open polyline that ends at the first intersect since we will not wrap back to
    // capture it as in the case of a closed polyline
    auto iter = intersectsLookup.begin();
    if (iter != intersectsLookup.end()) {
      auto slice = PlineSliceViewData<Real>::createFromSlicePoints(
          rawOffsetPline, rawOffsetPline[0].pos(), 0, iter->second[0], iter->first);
      maybeAppendSlice(0, slice);
    }
  }

  for (auto const &kvp : intersectsLookup) {
    // start index for the slice we're about to build
    std::size_t sIndex = kvp.first;
    // self intersect list for this start index
    std::vector<Vector2<Real>> const &siList = kvp.second;

    const auto &startVertex = rawOffsetPline[sIndex];
    std::size_t nextIndex = utils::nextWrappingIndex(sIndex, rawOffsetPline);
    const auto &endVertex = rawOffsetPline[nextIndex];

    if (siList.size() != 1) {
      // build all the segments between the N intersects in siList (N > 1), skipping the first
      // segment (to be processed at the end)
      SplitResult<Real> firstSplit = splitAtPoint(startVertex, endVertex, siList[0]);
      auto prevVertex = firstSplit.splitVertex;
      for (std::size_t i = 1; i < siList.size(); ++i) {
        SplitResult<Real> split = splitAtPoint(prevVertex, endVertex, siList[i]);
        // update prevVertex for next loop iteration
        prevVertex = split.splitVertex;
        // skip if they're ontop of each other
        if (fuzzyEqual(split.updatedStart.pos(), split.splitVertex.pos(),
                       utils::realPrecision<Real>())) {
          continue;
        }

        // test start point
        if (!pointValid(split.updatedStart.pos())) {
          continue;
        }

        // test end point
        if (!pointValid(split.splitVertex.pos())) {
          continue;
        }

        // test mid point
        auto midpoint = segMidpoint(split.updatedStart, split.splitVertex);
        if (!pointValid(midpoint)) {
          continue;
        }

        // test intersection with original polyline
        if (intersectsOrigPline(split.updatedStart, split.splitVertex)) {
          continue;
        }

        auto slice = PlineSliceViewData<Real>::createOnSingleSegment(
            rawOffsetPline, sIndex, split.updatedStart, split.splitVertex.pos());
        maybeAppendSlice(sIndex, slice);
      }
    }

    auto nextIntr = intersectsLookup.lower_bound(nextIndex);
    if (nextIntr != intersectsLookup.end()) {
      auto slice = PlineSliceViewData<Real>::createFromSlicePoints(
          rawOffsetPline, siList.back(), sIndex, nextIntr->second[0], nextIntr->first);
      maybeAppendSlice(sIndex, slice);
      continue;
    }

    if (originalPline.isClosed()) {
      auto wrapIntr = intersectsLookup.begin();
      if (wrapIntr != intersectsLookup.end()) {
        auto slice = PlineSliceViewData<Real>::createFromSlicePoints(
            rawOffsetPline, siList.back(), sIndex, wrapIntr->second[0], wrapIntr->first);
        maybeAppendSlice(sIndex, slice);
      }
      continue;
    }

    auto tailSlice = PlineSliceViewData<Real>::createFromSlicePoints(
        rawOffsetPline, siList.back(), sIndex, rawOffsetPline.lastVertex().pos(),
        rawOffsetPline.size() - 1);
    maybeAppendSlice(sIndex, tailSlice);
  }

  return result;
}

/// Stitches raw offset polyline slices together, discarding any that are not valid.
template <typename Real>
std::vector<Polyline<Real>>
stitchOffsetSlicesTogether(Polyline<Real> const &sourcePline,
                           std::vector<OpenPolylineSlice<Real>> const &slices, bool closedPolyline,
                           std::size_t origMaxIndex,
                           Real joinThreshold = utils::sliceJoinThreshold<Real>()) {
  std::vector<Polyline<Real>> result;
  if (slices.size() == 0) {
    return result;
  }

  if (slices.size() == 1) {
    result.emplace_back(slices[0].viewData.toPolyline(sourcePline, joinThreshold));
    if (closedPolyline &&
        fuzzyEqual(result[0][0].pos(), result[0].lastVertex().pos(), joinThreshold)) {
      result[0].isClosed() = true;
      result[0].vertexes().pop_back();
    }

    return result;
  }

  // load spatial index with all start points
  StaticSpatialIndex<Real> spatialIndex(slices.size());
  std::vector<Vector2<Real>> sliceStartPoints(slices.size());
  std::vector<Vector2<Real>> sliceEndPoints(slices.size());

  for (std::size_t i = 0; i < slices.size(); ++i) {
    auto const &slice = slices[i];
    auto const point = slice.viewData.firstPoint(sourcePline);
    sliceStartPoints[i] = point;
    sliceEndPoints[i] = slice.viewData.lastPoint(sourcePline);
    spatialIndex.add(point.x() - joinThreshold, point.y() - joinThreshold,
                     point.x() + joinThreshold, point.y() + joinThreshold);
  }

  spatialIndex.finish();

  std::vector<std::uint8_t> visitedIndexes(slices.size(), 0);
  std::vector<std::size_t> queryResults;
  queryResults.reserve(8);
  std::vector<std::size_t> queryStack;
  queryStack.reserve(8);
  for (std::size_t i = 0; i < slices.size(); ++i) {
    if (visitedIndexes[i] != 0) {
      continue;
    }

    visitedIndexes[i] = 1;

    Polyline<Real> currPline;
    std::size_t currIndex = i;
    auto const &initialStartPoint = sliceStartPoints[i];
    std::size_t loopCount = 0;
    const std::size_t maxLoopCount = slices.size();
    while (true) {
      if (loopCount++ > maxLoopCount) {
        CAVC_ASSERT(false, "Bug detected, should never loop this many times!");
        // break to avoid infinite loop
        break;
      }
      const std::size_t currLoopStartIndex = slices[currIndex].intrStartIndex;
      auto const &currSlice = slices[currIndex].viewData;
      auto const &currEndPoint = sliceEndPoints[currIndex];
      currSlice.appendTo(currPline, sourcePline, joinThreshold);
      queryResults.clear();
      auto queryVisitor = [&](std::size_t index) {
        if (visitedIndexes[index] == 0) {
          queryResults.push_back(index);
        }
        return true;
      };
      spatialIndex.visitQuery(currEndPoint.x() - joinThreshold, currEndPoint.y() - joinThreshold,
                              currEndPoint.x() + joinThreshold, currEndPoint.y() + joinThreshold,
                              queryVisitor, queryStack);

      auto indexDistAndEqualInitial = [&](std::size_t index) {
        auto const &slice = slices[index];
        std::size_t indexDist;
        if (currLoopStartIndex <= slice.intrStartIndex) {
          indexDist = slice.intrStartIndex - currLoopStartIndex;
        } else {
          // forward wrapping distance (distance to end + distance to index)
          indexDist = origMaxIndex - currLoopStartIndex + slice.intrStartIndex;
        }

        bool equalToInitial =
            fuzzyEqual(sliceEndPoints[index], initialStartPoint, joinThreshold);

        return std::make_pair(indexDist, equalToInitial);
      };

      std::sort(queryResults.begin(), queryResults.end(),
                [&](std::size_t index1, std::size_t index2) {
                  auto distAndEqualInitial1 = indexDistAndEqualInitial(index1);
                  auto distAndEqualInitial2 = indexDistAndEqualInitial(index2);
                  if (distAndEqualInitial1.first == distAndEqualInitial2.first) {
                    // index distances are equal, compare on position being equal to initial start
                    // (testing index1 < index2, we want the longest closed loop possible)
                    return distAndEqualInitial1.second < distAndEqualInitial2.second;
                  }

                  return distAndEqualInitial1.first < distAndEqualInitial2.first;
                });

      if (queryResults.size() == 0) {
        // we're done
        if (currPline.size() > 1) {
          if (closedPolyline &&
              fuzzyEqual(currPline[0].pos(), currPline.lastVertex().pos(), joinThreshold)) {
            currPline.vertexes().pop_back();
            currPline.isClosed() = true;
          }
          result.emplace_back(std::move(currPline));
        }
        break;
      }

      // else continue stitching
      visitedIndexes[queryResults[0]] = 1;
      currIndex = queryResults[0];
    }
  }

  return result;
}

template <typename Real> bool isSimpleClosedPolyline(Polyline<Real> const &candidate) {
  if (!candidate.isClosed() || candidate.size() < 3) {
    return false;
  }

  auto spatialIndex = createApproxSpatialIndex(candidate);
  std::vector<PlineIntersect<Real>> intersects;
  allSelfIntersects(candidate, intersects, spatialIndex);
  return intersects.empty();
}

template <typename Real>
std::vector<Polyline<Real>> filterSimpleClosedLoops(std::vector<Polyline<Real>> const &candidates) {
  std::vector<Polyline<Real>> filtered;
  filtered.reserve(candidates.size());
  for (auto const &candidate : candidates) {
    if (!isSimpleClosedPolyline(candidate)) {
      continue;
    }

    filtered.push_back(candidate);
  }

  return filtered;
}

template <typename Real> struct OffsetResultQualityThresholds {
  Real minPathLength;
  Real minClosedAbsArea;
  Real minRelaxedClosedAbsArea;
};

template <typename Real> Real finiteNonNegative(Real value) {
  if (!std::isfinite(value) || value < Real(0)) {
    return Real(0);
  }

  return value;
}

template <typename Real> Real polylineLengthScale(Polyline<Real> const &pline) {
  auto const extents = getExtents(pline);
  Real bboxDiagonal = Real(0);
  if (std::isfinite(extents.xMin) && std::isfinite(extents.yMin) &&
      std::isfinite(extents.xMax) && std::isfinite(extents.yMax)) {
    bboxDiagonal = std::hypot(extents.xMax - extents.xMin, extents.yMax - extents.yMin);
  }

  return std::max({Real(1), finiteNonNegative(bboxDiagonal),
                   finiteNonNegative(getPathLength(pline))});
}

template <typename Real>
OffsetResultQualityThresholds<Real> offsetResultQualityThresholds(Polyline<Real> const &reference,
                                                                  Real offset) {
  Real const lengthScale = std::max(polylineLengthScale(reference), std::abs(offset));
  // Keep the floor tied to runtime tolerances so small-coordinate workloads can still produce
  // usable open offsets after callers tighten the epsilon config.
  Real const tolerancePathFloor =
      std::max({utils::realPrecision<Real>() * Real(10), utils::sliceJoinThreshold<Real>(),
                utils::realThreshold<Real>() * Real(100)});
  Real const minPathLength =
      std::max(tolerancePathFloor, lengthScale * utils::realThreshold<Real>() * Real(100));
  Real const referenceAbsArea = std::abs(getArea(reference));
  Real const offsetAbsAreaScale = std::abs(offset * offset);
  Real const minClosedAbsArea =
      std::max(minPathLength * minPathLength, referenceAbsArea * Real(1e-8));
  Real const minRelaxedClosedAbsArea = std::max(minClosedAbsArea, offsetAbsAreaScale);
  return {minPathLength, minClosedAbsArea, minRelaxedClosedAbsArea};
}

template <typename Real>
bool hasOnlyEndpointTouchIntersects(Polyline<Real> const &candidate) {
  if (!candidate.isClosed() || candidate.size() < 3) {
    return false;
  }

  auto spatialIndex = createApproxSpatialIndex(candidate);
  std::vector<PlineIntersect<Real>> intersects;
  globalSelfIntersects(candidate, intersects, spatialIndex);
  Real const epsilon = utils::realPrecision<Real>();

  for (auto const &intr : intersects) {
    std::size_t const endIndex1 = utils::nextWrappingIndex(intr.sIndex1, candidate);
    std::size_t const endIndex2 = utils::nextWrappingIndex(intr.sIndex2, candidate);
    bool const onSeg1Endpoint = fuzzyEqual(candidate[intr.sIndex1].pos(), intr.pos, epsilon) ||
                                fuzzyEqual(candidate[endIndex1].pos(), intr.pos, epsilon);
    bool const onSeg2Endpoint = fuzzyEqual(candidate[intr.sIndex2].pos(), intr.pos, epsilon) ||
                                fuzzyEqual(candidate[endIndex2].pos(), intr.pos, epsilon);
    if (!(onSeg1Endpoint && onSeg2Endpoint)) {
      return false;
    }
  }

  return true;
}

template <typename Real>
std::vector<Polyline<Real>> filterClosedLoopsAllowingEndpointTouches(
    std::vector<Polyline<Real>> const &candidates,
    OffsetResultQualityThresholds<Real> const &qualityThresholds) {
  std::vector<Polyline<Real>> filtered;
  filtered.reserve(candidates.size());
  for (auto const &candidate : candidates) {
    if (!candidate.isClosed() || candidate.size() < 3) {
      continue;
    }
    if (getPathLength(candidate) <= qualityThresholds.minPathLength ||
        std::abs(getArea(candidate)) < qualityThresholds.minClosedAbsArea) {
      continue;
    }
    if (!hasOnlyEndpointTouchIntersects(candidate)) {
      continue;
    }

    filtered.push_back(candidate);
  }

  return filtered;
}

template <typename Real> bool offsetResultVertexesAreFinite(Polyline<Real> const &candidate) {
  for (auto const &vertex : candidate.vertexes()) {
    if (!std::isfinite(vertex.x()) || !std::isfinite(vertex.y()) ||
        !std::isfinite(vertex.bulge())) {
      return false;
    }
  }

  return true;
}

template <typename Real>
bool isUsableOpenOffsetPolyline(Polyline<Real> const &candidate,
                                OffsetResultQualityThresholds<Real> const &qualityThresholds) {
  return !candidate.isClosed() && candidate.size() >= 2 &&
         getPathLength(candidate) > qualityThresholds.minPathLength &&
         offsetResultVertexesAreFinite(candidate);
}

template <typename Real>
std::vector<Polyline<Real>> filterUsableOpenOffsetPolylines(
    std::vector<Polyline<Real>> const &candidates,
    OffsetResultQualityThresholds<Real> const &qualityThresholds) {
  std::vector<Polyline<Real>> filtered;
  filtered.reserve(candidates.size());
  for (auto const &candidate : candidates) {
    if (isUsableOpenOffsetPolyline(candidate, qualityThresholds)) {
      filtered.push_back(candidate);
    }
  }

  return filtered;
}

template <typename Real>
std::vector<Polyline<Real>> keepDominantOpenOffsetPolyline(
    std::vector<Polyline<Real>> const &candidates,
    OffsetResultQualityThresholds<Real> const &qualityThresholds) {
  if (candidates.size() < 2) {
    return candidates;
  }

  std::size_t bestIndex = 0;
  Real bestLength = getPathLength(candidates[0]);
  Real secondBestLength = Real(0);
  for (std::size_t i = 1; i < candidates.size(); ++i) {
    Real const candidateLength = getPathLength(candidates[i]);
    if (candidateLength > bestLength) {
      secondBestLength = bestLength;
      bestLength = candidateLength;
      bestIndex = i;
    } else if (candidateLength > secondBestLength) {
      secondBestLength = candidateLength;
    }
  }

  Real const dominanceThreshold = std::max(qualityThresholds.minPathLength * Real(10),
                                           secondBestLength * Real(4));
  if (bestLength <= dominanceThreshold) {
    return candidates;
  }

  std::vector<Polyline<Real>> result;
  result.reserve(1);
  result.push_back(candidates[bestIndex]);
  return result;
}


template <typename Real>
std::vector<Polyline<Real>> filterClosedLoopsWithMinimumAbsArea(
    std::vector<Polyline<Real>> const &candidates,
    OffsetResultQualityThresholds<Real> const &qualityThresholds) {
  std::vector<Polyline<Real>> filtered;
  filtered.reserve(candidates.size());
  for (auto const &candidate : candidates) {
    if (!candidate.isClosed() || candidate.size() < 3) {
      continue;
    }
    if (!offsetResultVertexesAreFinite(candidate)) {
      continue;
    }
    if (getPathLength(candidate) <= qualityThresholds.minPathLength) {
      continue;
    }
    if (std::abs(getArea(candidate)) + utils::realPrecision<Real>() <
        qualityThresholds.minClosedAbsArea) {
      continue;
    }

    filtered.push_back(candidate);
  }

  return filtered;
}

template <typename Real>
std::vector<Polyline<Real>>
stitchSlicesIntoSimpleClosedLoops(Polyline<Real> const &sourcePline,
                                  std::vector<OpenPolylineSlice<Real>> const &slices,
                                  std::size_t origMaxIndex,
                                  Real joinThreshold = utils::sliceJoinThreshold<Real>());

template <typename Real>
std::vector<Polyline<Real>> recoverOpenOffsetPolylinesFromRelaxedSlices(
    Polyline<Real> const &cleaned, Polyline<Real> const &rawOffset,
    Polyline<Real> const &dualRawOffset, Real offset,
    ParallelOffsetOptions<Real> const &options) {
  CAVC_ASSERT(!cleaned.isClosed(), "relaxed open-offset recovery requires an open polyline");

  auto const qualityThresholds = offsetResultQualityThresholds(cleaned, offset);
  auto recoverFromSlices = [&](bool skipOrigIntersectionCheck) {
    auto relaxedSlices = dualSliceAtIntersectsForOffset(cleaned, rawOffset, dualRawOffset, offset,
                                                        options, false, skipOrigIntersectionCheck);
    auto relaxedResult =
        stitchOffsetSlicesTogether(rawOffset, relaxedSlices, cleaned.isClosed(), rawOffset.size() - 1);
    return keepDominantOpenOffsetPolyline(
        filterUsableOpenOffsetPolylines(relaxedResult, qualityThresholds), qualityThresholds);
  };

  auto recovered = recoverFromSlices(false);
  if (!recovered.empty()) {
    return recovered;
  }

  return recoverFromSlices(true);
}

template <typename Real>
std::vector<Polyline<Real>> recoverClosedOffsetLoopsFromRelaxedSlices(
    Polyline<Real> const &cleaned, Polyline<Real> const &rawOffset, Real offset) {
  CAVC_ASSERT(cleaned.isClosed(), "relaxed closed-loop recovery requires a closed polyline");
  if (rawOffset.size() < 2) {
    return {};
  }

  auto qualityThresholds = offsetResultQualityThresholds(cleaned, offset);
  qualityThresholds.minClosedAbsArea = qualityThresholds.minRelaxedClosedAbsArea;

  auto recoverFromSlices = [&](bool skipOrigIntersectionCheck) {
    auto relaxedSlices =
        slicesFromRawOffset(cleaned, rawOffset, offset, false, skipOrigIntersectionCheck);
    auto relaxedStitched =
        stitchOffsetSlicesTogether(rawOffset, relaxedSlices, true, rawOffset.size() - 1);

    auto simpleRecovered = filterClosedLoopsWithMinimumAbsArea(
        filterSimpleClosedLoops(relaxedStitched), qualityThresholds);
    if (!simpleRecovered.empty()) {
      return simpleRecovered;
    }

    auto graphRecovered = stitchSlicesIntoSimpleClosedLoops(rawOffset, relaxedSlices, rawOffset.size() - 1);
    auto graphFiltered = filterClosedLoopsWithMinimumAbsArea(graphRecovered, qualityThresholds);
    if (!graphFiltered.empty()) {
      return graphFiltered;
    }

    return filterClosedLoopsWithMinimumAbsArea(
        filterClosedLoopsAllowingEndpointTouches(relaxedStitched, qualityThresholds),
        qualityThresholds);
  };

  // First keep the original-intersection rejection enabled. If that still fails, relax only that
  // final check, while retaining the non-round raw join geometry and area/path/finite filters.
  auto recovered = recoverFromSlices(false);
  if (!recovered.empty()) {
    return recovered;
  }

  return recoverFromSlices(true);
}

template <typename Real>
std::vector<Polyline<Real>> filterCollapsedLineLoops(std::vector<Polyline<Real>> const &candidates) {
  if (candidates.size() != 1) {
    return {};
  }

  if (candidates[0].size() != 2) {
    return {};
  }

  return candidates;
}

template <typename Real>
std::vector<Polyline<Real>>
stitchSlicesIntoSimpleClosedLoops(Polyline<Real> const &sourcePline,
                                  std::vector<OpenPolylineSlice<Real>> const &slices,
                                  std::size_t origMaxIndex,
                                  Real joinThreshold) {
  std::vector<Polyline<Real>> result;
  if (slices.empty()) {
    return result;
  }

  // This search is only a salvage path for fragmented closed offsets; keep it bounded so pathological
  // slice graphs do not explode exponentially. Callers fall back to returning no loops rather than
  // invalid geometry when the search budget is exceeded.
  constexpr std::size_t maxSearchSlices = 24;
  if (slices.size() > maxSearchSlices) {
    return result;
  }

  StaticSpatialIndex<Real> spatialIndex(slices.size());
  for (auto const &slice : slices) {
    auto const &point = slice.viewData.firstPoint(sourcePline);
    spatialIndex.add(point.x() - joinThreshold, point.y() - joinThreshold,
                     point.x() + joinThreshold, point.y() + joinThreshold);
  }
  spatialIndex.finish();

  auto forwardDistance = [&](std::size_t from, std::size_t to) {
    std::size_t const maxIndex = std::numeric_limits<std::size_t>::max();
    std::size_t const fromIndex = slices[from].intrStartIndex;
    std::size_t const toIndex = slices[to].intrStartIndex;
    if (fromIndex == maxIndex || toIndex == maxIndex) {
      return std::size_t{0};
    }
    if (fromIndex <= toIndex) {
      return toIndex - fromIndex;
    }
    return origMaxIndex - fromIndex + toIndex;
  };

  std::vector<std::vector<std::size_t>> nextIndexes(slices.size());
  std::vector<std::size_t> queryResults;
  queryResults.reserve(8);
  std::vector<std::size_t> queryStack;
  queryStack.reserve(8);

  for (std::size_t i = 0; i < slices.size(); ++i) {
    auto const &currEndPoint = slices[i].viewData.lastPoint(sourcePline);
    queryResults.clear();
    auto queryVisitor = [&](std::size_t index) {
      if (index != i) {
        queryResults.push_back(index);
      }
      return true;
    };
    spatialIndex.visitQuery(currEndPoint.x() - joinThreshold, currEndPoint.y() - joinThreshold,
                            currEndPoint.x() + joinThreshold, currEndPoint.y() + joinThreshold,
                            queryVisitor, queryStack);
    std::sort(queryResults.begin(), queryResults.end(),
              [&](std::size_t left, std::size_t right) {
                std::size_t const leftDist = forwardDistance(i, left);
                std::size_t const rightDist = forwardDistance(i, right);
                if (leftDist != rightDist) {
                  return leftDist > rightDist;
                }
                if (slices[left].viewData.vertexCount() != slices[right].viewData.vertexCount()) {
                  return slices[left].viewData.vertexCount() > slices[right].viewData.vertexCount();
                }
                return left < right;
              });
    nextIndexes[i] = std::move(queryResults);
  }

  std::vector<std::uint8_t> usedIndexes(slices.size(), 0);
  std::vector<std::uint8_t> localUsed(slices.size(), 0);
  std::vector<std::size_t> path;
  std::vector<std::size_t> acceptedPath;
  Polyline<Real> acceptedLoop;
  auto materializePath = [&](std::vector<std::size_t> const &slicePath) {
    Polyline<Real> pline;
    for (std::size_t sliceIndex : slicePath) {
      slices[sliceIndex].viewData.appendTo(pline, sourcePline, joinThreshold);
    }
    return pline;
  };

  auto dfs = [&](auto &&self, std::size_t currIndex) -> bool {
    Polyline<Real> candidate = materializePath(path);
    if (candidate.size() > 2 &&
        fuzzyEqual(candidate[0].pos(), candidate.lastVertex().pos(), joinThreshold)) {
      candidate.isClosed() = true;
      candidate.vertexes().pop_back();
      if (isSimpleClosedPolyline(candidate) && getPathLength(candidate) > Real(1e-2)) {
        acceptedPath = path;
        acceptedLoop = std::move(candidate);
        return true;
      }
    }

    for (std::size_t nextIndex : nextIndexes[currIndex]) {
      if (usedIndexes[nextIndex] != 0 || localUsed[nextIndex] != 0) {
        continue;
      }

      localUsed[nextIndex] = 1;
      path.push_back(nextIndex);
      if (self(self, nextIndex)) {
        return true;
      }
      path.pop_back();
      localUsed[nextIndex] = 0;
    }

    return false;
  };

  for (std::size_t i = 0; i < slices.size(); ++i) {
    if (usedIndexes[i] != 0) {
      continue;
    }

    path.clear();
    path.push_back(i);
    std::fill(localUsed.begin(), localUsed.end(), false);
    localUsed[i] = 1;
    acceptedPath.clear();
    acceptedLoop = Polyline<Real>();

    if (dfs(dfs, i)) {
      for (std::size_t index : acceptedPath) {
        usedIndexes[index] = 1;
      }
      result.emplace_back(std::move(acceptedLoop));
    }
  }

  return result;
}
} // namespace internal

/// Creates the paralell offset polylines to the polyline given.
template <typename Real>
std::vector<Polyline<Real>> parallelOffset(Polyline<Real> const &pline, Real offset,
                                           ParallelOffsetOptions<Real> const &options = {}) {
  using namespace internal;
  if (options.joinType == OffsetJoinType::Miter) {
    CAVC_ASSERT(options.miterLimit >= Real(1), "miterLimit must be >= 1");
  }

  Polyline<Real> cleaned = removeRedundant(pline, utils::realPrecision<Real>());
  if (cleaned.size() < 2) {
    return std::vector<Polyline<Real>>();
  }
  auto rawOffset = createRawOffsetPline(cleaned, offset, options);
  if (rawOffset.size() < 2) {
    return std::vector<Polyline<Real>>();
  }
  auto const qualityThresholds = offsetResultQualityThresholds(cleaned, offset);
  if (cleaned.isClosed() && !options.hasSelfIntersects) {
    auto slices = slicesFromRawOffset(cleaned, rawOffset, offset);
    auto result =
        stitchOffsetSlicesTogether(rawOffset, slices, cleaned.isClosed(), rawOffset.size() - 1);
    auto filteredResult = filterSimpleClosedLoops(result);
    if (!filteredResult.empty()) {
      return filteredResult;
    }

    auto collapsedLineResult = filterCollapsedLineLoops(result);
    if (!collapsedLineResult.empty()) {
      return collapsedLineResult;
    }

    auto rescuedResult = stitchSlicesIntoSimpleClosedLoops(rawOffset, slices, rawOffset.size() - 1);
    if (!rescuedResult.empty()) {
      return rescuedResult;
    }

    auto endpointTouchResult = filterClosedLoopsAllowingEndpointTouches(result, qualityThresholds);
    if (!endpointTouchResult.empty()) {
      return endpointTouchResult;
    }

    auto relaxedRecoveredResult =
        recoverClosedOffsetLoopsFromRelaxedSlices(cleaned, rawOffset, offset);
    if (!relaxedRecoveredResult.empty()) {
      return relaxedRecoveredResult;
    }

    return std::vector<Polyline<Real>>();
  }

  // not closed polyline or has self intersects, must apply dual clipping
  auto dualRawOffset = createRawOffsetPline(cleaned, -offset, options);
  bool const enforceMinDistance =
      !cleaned.isClosed() || options.joinType == OffsetJoinType::Round;
  auto slices = dualSliceAtIntersectsForOffset(cleaned, rawOffset, dualRawOffset, offset, options,
                                               enforceMinDistance);
  auto result =
      stitchOffsetSlicesTogether(rawOffset, slices, cleaned.isClosed(), rawOffset.size() - 1);
  if (!cleaned.isClosed()) {
    auto filteredOpenResult = filterUsableOpenOffsetPolylines(result, qualityThresholds);
    if (options.joinType != OffsetJoinType::Round) {
      filteredOpenResult =
          keepDominantOpenOffsetPolyline(filteredOpenResult, qualityThresholds);
    }
    if (!filteredOpenResult.empty()) {
      return filteredOpenResult;
    }

    if (options.joinType != OffsetJoinType::Round) {
      auto relaxedOpenResult = recoverOpenOffsetPolylinesFromRelaxedSlices(
          cleaned, rawOffset, dualRawOffset, offset, options);
      if (!relaxedOpenResult.empty()) {
        return relaxedOpenResult;
      }
    }

    return {};
  }

  auto filteredResult = filterSimpleClosedLoops(result);
  if (!filteredResult.empty()) {
    return filteredResult;
  }

  auto collapsedLineResult = filterCollapsedLineLoops(result);
  if (!collapsedLineResult.empty()) {
    return collapsedLineResult;
  }

  auto rescuedResult = stitchSlicesIntoSimpleClosedLoops(rawOffset, slices, rawOffset.size() - 1);
  if (!rescuedResult.empty()) {
    return rescuedResult;
  }

  auto endpointTouchResult = filterClosedLoopsAllowingEndpointTouches(result, qualityThresholds);
  if (!endpointTouchResult.empty()) {
    return endpointTouchResult;
  }

  if (cleaned.isClosed() && options.joinType != OffsetJoinType::Round) {
    auto relaxedRecoveredResult =
        recoverClosedOffsetLoopsFromRelaxedSlices(cleaned, rawOffset, offset);
    if (!relaxedRecoveredResult.empty()) {
      return relaxedRecoveredResult;
    }
  }

  return std::vector<Polyline<Real>>();
}

} // namespace cavc
#endif // CAVC_POLYLINEOFFSET_HPP
