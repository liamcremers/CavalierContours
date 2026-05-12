#include "cavaliercontours.h"
#include "cavc/compoundpolylinecombine.hpp"
#include "cavc/polylinecombine.hpp"
#include "cavc/polylineoffset.hpp"
#include "cavc/polylineoffsetislands.hpp"
#include <cmath>
#include <exception>
#include <iostream>
#include <limits>
#include <memory>
#include <span>

#define CAVC_BEGIN_TRY_CATCH try {

#define CAVC_END_TRY_CATCH                                                                         \
  }                                                                                                \
  catch (std::exception const &ex) {                                                               \
    std::cerr << "Unexpected exception thrown in " << __FUNCTION__ << ": " << ex.what();           \
    std::terminate();                                                                              \
  }

struct cavc_pline {
  cavc::Polyline<cavc_real> data;
  cavc_pline() = default;
  cavc_pline(cavc::Polyline<cavc_real> &&p_data) noexcept : data(std::move(p_data)) {}
};

struct cavc_pline_list {
  std::vector<std::unique_ptr<cavc_pline>> data;
};

struct cavc_spatial_index {
  cavc::StaticSpatialIndex<cavc_real> data;
  uint32_t item_count;
  cavc_spatial_index(cavc::StaticSpatialIndex<cavc_real> &&p_data, uint32_t p_item_count) noexcept
      : data(std::move(p_data)), item_count(p_item_count) {}
};

struct cavc_offset_loop_topology {
  std::vector<cavc_offset_loop_topology_node> nodes;
};

struct cavc_compound_pline {
  std::vector<std::unique_ptr<cavc_pline>> outerLoops;
  std::vector<std::unique_ptr<cavc_pline>> holeLoops;

  cavc::CompoundPolyline<cavc_real> toCpp() const {
    cavc::CompoundPolyline<cavc_real> result;
    result.outerLoops.reserve(outerLoops.size());
    result.holeLoops.reserve(holeLoops.size());
    for (auto const &l : outerLoops) {
      result.outerLoops.push_back(l->data);
    }
    for (auto const &l : holeLoops) {
      result.holeLoops.push_back(l->data);
    }
    return result;
  }
};

static cavc_tolerances to_api_tolerances(cavc::utils::EpsilonConfig<cavc_real> const &config) {
  return {config.realThreshold, config.realPrecision, config.sliceJoinThreshold,
          config.offsetThreshold};
}

static cavc::utils::EpsilonConfig<cavc_real> to_cpp_tolerances(cavc_tolerances const &config) {
  return {config.real_threshold, config.real_precision, config.slice_join_threshold,
          config.offset_threshold};
}

static cavc::ClosedPolylineWinding to_cpp_closed_winding(cavc_closed_winding winding) {
  switch (winding) {
  case CAVC_WINDING_KEEP:
    return cavc::ClosedPolylineWinding::Keep;
  case CAVC_WINDING_COUNTER_CLOCKWISE:
    return cavc::ClosedPolylineWinding::CounterClockwise;
  case CAVC_WINDING_CLOCKWISE:
    return cavc::ClosedPolylineWinding::Clockwise;
  }
  CAVC_ASSERT(false, "Unhandled closed winding enum");
  std::terminate();
}

static bool is_valid_offset_join_type(cavc_offset_join_type join_type) {
  switch (join_type) {
  case CAVC_OFFSET_JOIN_ROUND:
  case CAVC_OFFSET_JOIN_MITER:
  case CAVC_OFFSET_JOIN_BEVEL:
    return true;
  }
  return false;
}

static bool is_valid_offset_end_cap_type(cavc_offset_end_cap_type end_cap_type) {
  switch (end_cap_type) {
  case CAVC_OFFSET_END_CAP_ROUND:
  case CAVC_OFFSET_END_CAP_SQUARE:
  case CAVC_OFFSET_END_CAP_BUTT:
    return true;
  }
  return false;
}

static bool is_valid_parallel_offset_options(cavc_parallel_offset_options const &options) {
  if (!is_valid_offset_join_type(options.join_type) ||
      !is_valid_offset_end_cap_type(options.end_cap_type)) {
    return false;
  }

  if (options.join_type == CAVC_OFFSET_JOIN_MITER &&
      (!std::isfinite(options.miter_limit) || options.miter_limit < cavc_real(1))) {
    return false;
  }

  return true;
}

static cavc::OffsetJoinType to_cpp_offset_join_type(cavc_offset_join_type join_type) {
  switch (join_type) {
  case CAVC_OFFSET_JOIN_ROUND:
    return cavc::OffsetJoinType::Round;
  case CAVC_OFFSET_JOIN_MITER:
    return cavc::OffsetJoinType::Miter;
  case CAVC_OFFSET_JOIN_BEVEL:
    return cavc::OffsetJoinType::Bevel;
  }
  CAVC_ASSERT(false, "Unhandled offset join enum");
  std::terminate();
}

static cavc::OffsetEndCapType to_cpp_offset_end_cap_type(cavc_offset_end_cap_type end_cap_type) {
  switch (end_cap_type) {
  case CAVC_OFFSET_END_CAP_ROUND:
    return cavc::OffsetEndCapType::Round;
  case CAVC_OFFSET_END_CAP_SQUARE:
    return cavc::OffsetEndCapType::Square;
  case CAVC_OFFSET_END_CAP_BUTT:
    return cavc::OffsetEndCapType::Butt;
  }
  CAVC_ASSERT(false, "Unhandled offset end cap enum");
  std::terminate();
}

static cavc::ParallelOffsetOptions<cavc_real>
to_cpp_parallel_offset_options(cavc_parallel_offset_options const &options) {
  CAVC_ASSERT(is_valid_parallel_offset_options(options), "invalid parallel offset options");
  cavc::ParallelOffsetOptions<cavc_real> result;
  result.hasSelfIntersects = options.may_have_self_intersects != 0;
  result.joinType = to_cpp_offset_join_type(options.join_type);
  result.endCapType = to_cpp_offset_end_cap_type(options.end_cap_type);
  result.miterLimit = options.miter_limit;
  return result;
}

// helper to move vector of plines to cavc_pline_list
static void move_to_list(std::vector<cavc::Polyline<cavc_real>> &&plines, cavc_pline_list *list) {
  list->data.reserve(plines.size());

  for (std::size_t i = 0; i < plines.size(); ++i) {
    list->data.push_back(std::make_unique<cavc_pline>(std::move(plines[i])));
  }
}

// helper to copy vertex data to cavc_pline
static void copy_to_pline(cavc_pline *api_pline, cavc_vertex const *vertex_data,
                          uint32_t vertex_count) {
  api_pline->data.vertexes().clear();
  api_pline->data.vertexes().reserve(vertex_count);

  // Use std::span for safe buffer access
  if (vertex_data != nullptr && vertex_count > 0) {
    std::span<const cavc_vertex> vertices_span(vertex_data, vertex_count);
    for (const auto &vertex : vertices_span) {
      api_pline->data.addVertex(vertex.x, vertex.y, vertex.bulge);
    }
  }
}

// helper to copy to vertex_data from cavc_pline
static void copy_to_vertex_data(cavc_pline const *api_pline, cavc_vertex *vertex_data) {
  auto const &pline = api_pline->data;
  uint32_t vertex_count = static_cast<uint32_t>(pline.size());

  if (vertex_data != nullptr && vertex_count > 0) {
    std::span<cavc_vertex> output_span(vertex_data, vertex_count);
    for (uint32_t i = 0; i < vertex_count; ++i) {
      auto const &v = pline[i];
      output_span[i] = cavc_vertex{v.x(), v.y(), v.bulge()};
    }
  }
}

static std::vector<std::size_t> query_spatial_index(cavc_spatial_index const *spatial_index,
                                                    cavc_real min_x, cavc_real min_y,
                                                    cavc_real max_x, cavc_real max_y) {
  std::vector<std::size_t> results;
  spatial_index->data.query(min_x, min_y, max_x, max_y, results);
  return results;
}

static cavc_offset_loop_role to_api_offset_loop_role(cavc::OffsetLoopRole role) {
  switch (role) {
  case cavc::OffsetLoopRole::Outer:
    return CAVC_OFFSET_LOOP_ROLE_OUTER;
  case cavc::OffsetLoopRole::Hole:
    return CAVC_OFFSET_LOOP_ROLE_HOLE;
  }
  CAVC_ASSERT(false, "Unhandled offset loop role");
  return CAVC_OFFSET_LOOP_ROLE_OUTER;
}

static cavc_offset_loop_topology_node
to_api_topology_node(cavc::OffsetLoopTopologyNode<cavc_real> const &node) {
  CAVC_ASSERT(node.sourceIndex <= std::numeric_limits<uint32_t>::max(),
              "source index exceeds uint32_t");

  uint32_t parent_index = CAVC_OFFSET_LOOP_NO_PARENT;
  if (node.parentIndex != cavc::kNoParentOffsetLoop) {
    CAVC_ASSERT(node.parentIndex <= std::numeric_limits<uint32_t>::max(),
                "parent index exceeds uint32_t");
    parent_index = static_cast<uint32_t>(node.parentIndex);
  }

  return {to_api_offset_loop_role(node.role), static_cast<uint32_t>(node.sourceIndex),
          parent_index};
}

static std::vector<cavc::OffsetLoop<cavc_real>> to_cpp_offset_loops(cavc_pline const *const *loops,
                                                                    uint32_t loop_count) {
  CAVC_ASSERT(loop_count == 0 || loops != nullptr, "non-zero loop count requires loop pointer");

  std::vector<cavc::OffsetLoop<cavc_real>> result;
  result.reserve(loop_count);

  for (uint32_t i = 0; i < loop_count; ++i) {
    CAVC_ASSERT(loops[i] != nullptr, "null loop pointer not allowed");
    auto const &loop = loops[i]->data;
    CAVC_ASSERT(loop.isClosed(), "offset topology requires closed loops");
    CAVC_ASSERT(loop.size() > 1, "offset topology requires loops with at least 2 vertices");

    result.push_back({static_cast<std::size_t>(i), loop, cavc::createApproxSpatialIndex(loop)});
  }

  return result;
}

cavc_pline *cavc_pline_new(const cavc_vertex *vertex_data, uint32_t vertex_count, int is_closed) {
  CAVC_BEGIN_TRY_CATCH
  cavc_pline *result = new cavc_pline();
  if (vertex_data) {
    copy_to_pline(result, vertex_data, vertex_count);
  } else {
    result->data.vertexes().reserve(vertex_count);
  }
  result->data.isClosed() = is_closed;
  return result;
  CAVC_END_TRY_CATCH
}

void cavc_pline_delete(cavc_pline *polyline) {
  CAVC_BEGIN_TRY_CATCH
  delete polyline;
  CAVC_END_TRY_CATCH
}

uint32_t cavc_pline_capacity(cavc_pline const *pline) {
  CAVC_ASSERT(pline, "null pline not allowed");
  CAVC_BEGIN_TRY_CATCH
  return static_cast<uint32_t>(pline->data.vertexes().capacity());
  CAVC_END_TRY_CATCH
}

void cavc_pline_set_capacity(cavc_pline *pline, uint32_t size) {
  CAVC_ASSERT(pline, "null pline not allowed");
  CAVC_BEGIN_TRY_CATCH
  pline->data.vertexes().reserve(size);
  CAVC_END_TRY_CATCH
}

uint32_t cavc_pline_vertex_count(cavc_pline const *pline) {
  CAVC_ASSERT(pline, "null pline not allowed");
  CAVC_BEGIN_TRY_CATCH
  return static_cast<uint32_t>(pline->data.size());
  CAVC_END_TRY_CATCH
}

void cavc_pline_vertex_data(cavc_pline const *pline, cavc_vertex *vertex_data) {
  CAVC_ASSERT(pline, "null pline not allowed");
  CAVC_BEGIN_TRY_CATCH
  copy_to_vertex_data(pline, vertex_data);
  CAVC_END_TRY_CATCH
}

int cavc_pline_is_closed(cavc_pline const *pline) {
  CAVC_ASSERT(pline, "null pline not allowed");
  CAVC_BEGIN_TRY_CATCH
  return pline->data.isClosed();
  CAVC_END_TRY_CATCH
}

void cavc_pline_set_vertex_data(cavc_pline *pline, const cavc_vertex *vertex_data,
                                uint32_t vertex_count) {
  CAVC_ASSERT(pline, "null pline not allowed");
  CAVC_BEGIN_TRY_CATCH
  copy_to_pline(pline, vertex_data, vertex_count);
  CAVC_END_TRY_CATCH
}

void cavc_pline_add_vertex(cavc_pline *pline, cavc_vertex vertex) {
  CAVC_ASSERT(pline, "null pline not allowed");
  CAVC_BEGIN_TRY_CATCH
  pline->data.addVertex(vertex.x, vertex.y, vertex.bulge);
  CAVC_END_TRY_CATCH
}

void cavc_pline_remove_range(cavc_pline *pline, uint32_t start_index, uint32_t count) {
  CAVC_ASSERT(pline, "null pline not allowed");
  CAVC_ASSERT(start_index < pline->data.size(), "start_index is out of vertexes range");
  CAVC_ASSERT(start_index + count <= pline->data.size(), "count is out of vertexes range");
  CAVC_BEGIN_TRY_CATCH
  auto &vertexes = pline->data.vertexes();
  auto start_it = vertexes.begin() + static_cast<std::ptrdiff_t>(start_index);
  vertexes.erase(start_it, start_it + static_cast<std::ptrdiff_t>(count));
  CAVC_END_TRY_CATCH
}

void cavc_pline_clear(cavc_pline *pline) {
  CAVC_ASSERT(pline, "null pline not allowed");
  CAVC_BEGIN_TRY_CATCH
  pline->data.vertexes().clear();
  CAVC_END_TRY_CATCH
}

void cavc_pline_set_is_closed(cavc_pline *pline, int is_closed) {
  CAVC_ASSERT(pline, "null pline not allowed");
  CAVC_BEGIN_TRY_CATCH
  pline->data.isClosed() = is_closed;
  CAVC_END_TRY_CATCH
}

void cavc_pline_invert_direction(cavc_pline *pline) {
  CAVC_ASSERT(pline, "null pline not allowed");
  CAVC_BEGIN_TRY_CATCH
  cavc::invertDirection(pline->data);
  CAVC_END_TRY_CATCH
}

void cavc_pline_prune_singularities(cavc_pline const *pline, cavc_real epsilon,
                                    cavc_pline **output) {
  CAVC_ASSERT(pline, "null pline not allowed");
  CAVC_ASSERT(output, "null output not allowed");
  CAVC_BEGIN_TRY_CATCH
  *output = new cavc_pline(cavc::pruneSingularities(pline->data, epsilon));
  CAVC_END_TRY_CATCH
}

void cavc_pline_remove_redundant(cavc_pline *pline, cavc_real epsilon) {
  CAVC_ASSERT(pline, "null pline not allowed");
  CAVC_BEGIN_TRY_CATCH
  pline->data = cavc::removeRedundant(pline->data, epsilon);
  CAVC_END_TRY_CATCH
}

void cavc_pline_convert_arcs_to_lines(cavc_pline const *pline, cavc_real error,
                                      cavc_pline **output) {
  CAVC_ASSERT(pline, "null pline not allowed");
  CAVC_ASSERT(output, "null output not allowed");
  CAVC_BEGIN_TRY_CATCH
  *output = new cavc_pline(cavc::convertArcsToLines(pline->data, error));
  CAVC_END_TRY_CATCH
}

void cavc_pline_normalize(cavc_pline const *pline, cavc_real epsilon,
                          cavc_closed_winding closed_winding, cavc_pline **output) {
  CAVC_ASSERT(pline, "null pline not allowed");
  CAVC_ASSERT(output, "null output not allowed");
  CAVC_BEGIN_TRY_CATCH
  auto winding = to_cpp_closed_winding(closed_winding);
  *output = new cavc_pline(cavc::normalizePolyline(pline->data, epsilon, winding));
  CAVC_END_TRY_CATCH
}

void cavc_pline_list_delete(cavc_pline_list *pline_list) {
  CAVC_BEGIN_TRY_CATCH
  delete pline_list;
  CAVC_END_TRY_CATCH
}

uint32_t cavc_pline_list_count(cavc_pline_list const *pline_list) {
  CAVC_ASSERT(pline_list, "null pline_list not allowed");
  CAVC_BEGIN_TRY_CATCH
  return static_cast<uint32_t>(pline_list->data.size());
  CAVC_END_TRY_CATCH
}

cavc_pline *cavc_pline_list_get(cavc_pline_list const *pline_list, uint32_t index) {
  CAVC_ASSERT(pline_list, "null pline_list not allowed");
  CAVC_ASSERT(index < pline_list->data.size(), "index is out of vertexes range");
  CAVC_BEGIN_TRY_CATCH
  return pline_list->data[index].get();
  CAVC_END_TRY_CATCH
}

cavc_pline *cavc_pline_list_release(cavc_pline_list *pline_list, uint32_t index) {
  CAVC_ASSERT(pline_list, "null pline_list not allowed");
  CAVC_ASSERT(index < pline_list->data.size(), "index is out of vertexes range");
  CAVC_BEGIN_TRY_CATCH
  cavc_pline *target = pline_list->data[index].release();
  pline_list->data.erase(pline_list->data.begin() + static_cast<std::ptrdiff_t>(index));
  return target;
  CAVC_END_TRY_CATCH
}

// cavc_spatial_index APIs
// -------------------------
cavc_spatial_index *cavc_spatial_index_create(cavc_pline const *pline) {
  CAVC_ASSERT(pline, "null pline not allowed");
  CAVC_ASSERT(pline->data.size() > 1,
              "need at least 2 vertexes to form segments for spatial index");
  CAVC_BEGIN_TRY_CATCH
  uint32_t segment_count =
      static_cast<uint32_t>(pline->data.isClosed() ? pline->data.size() : pline->data.size() - 1);
  return new cavc_spatial_index(cavc::createApproxSpatialIndex(pline->data), segment_count);
  CAVC_END_TRY_CATCH
}

void cavc_spatial_index_delete(cavc_spatial_index *spatial_index) {
  CAVC_BEGIN_TRY_CATCH
  delete spatial_index;
  CAVC_END_TRY_CATCH
}

uint32_t cavc_spatial_index_item_count(cavc_spatial_index const *spatial_index) {
  CAVC_ASSERT(spatial_index, "null spatial_index not allowed");
  CAVC_BEGIN_TRY_CATCH
  return spatial_index->item_count;
  CAVC_END_TRY_CATCH
}

uint32_t cavc_spatial_index_query_count(cavc_spatial_index const *spatial_index, cavc_real min_x,
                                        cavc_real min_y, cavc_real max_x, cavc_real max_y) {
  CAVC_ASSERT(spatial_index, "null spatial_index not allowed");
  CAVC_BEGIN_TRY_CATCH
  auto results = query_spatial_index(spatial_index, min_x, min_y, max_x, max_y);
  return static_cast<uint32_t>(results.size());
  CAVC_END_TRY_CATCH
}

void cavc_spatial_index_query(cavc_spatial_index const *spatial_index, cavc_real min_x,
                              cavc_real min_y, cavc_real max_x, cavc_real max_y,
                              uint32_t *results_out) {
  CAVC_ASSERT(spatial_index, "null spatial_index not allowed");
  CAVC_ASSERT(results_out, "null results_out not allowed");
  CAVC_BEGIN_TRY_CATCH
  auto results = query_spatial_index(spatial_index, min_x, min_y, max_x, max_y);
  for (std::size_t i = 0; i < results.size(); ++i) {
    results_out[i] = static_cast<uint32_t>(results[i]);
  }
  CAVC_END_TRY_CATCH
}

cavc_offset_loop_topology *cavc_offset_loop_topology_build(cavc_pline const *const *ccw_loops,
                                                           uint32_t ccw_loop_count,
                                                           cavc_pline const *const *cw_loops,
                                                           uint32_t cw_loop_count,
                                                           cavc_real boundary_epsilon) {
  CAVC_ASSERT(boundary_epsilon >= cavc_real(0), "boundary_epsilon must be >= 0");
  CAVC_BEGIN_TRY_CATCH
  cavc::OffsetLoopSet<cavc_real> loop_set;
  loop_set.ccwLoops = to_cpp_offset_loops(ccw_loops, ccw_loop_count);
  loop_set.cwLoops = to_cpp_offset_loops(cw_loops, cw_loop_count);

  auto cpp_topology = cavc::buildOffsetLoopTopology(loop_set, boundary_epsilon);

  auto *result = new cavc_offset_loop_topology();
  result->nodes.reserve(cpp_topology.size());
  for (auto const &node : cpp_topology) {
    result->nodes.push_back(to_api_topology_node(node));
  }

  return result;
  CAVC_END_TRY_CATCH
}

void cavc_offset_loop_topology_delete(cavc_offset_loop_topology *topology) {
  CAVC_BEGIN_TRY_CATCH
  delete topology;
  CAVC_END_TRY_CATCH
}

uint32_t cavc_offset_loop_topology_count(cavc_offset_loop_topology const *topology) {
  CAVC_ASSERT(topology, "null topology not allowed");
  CAVC_BEGIN_TRY_CATCH
  CAVC_ASSERT(topology->nodes.size() <= std::numeric_limits<uint32_t>::max(),
              "node count exceeds uint32_t");
  return static_cast<uint32_t>(topology->nodes.size());
  CAVC_END_TRY_CATCH
}

cavc_offset_loop_topology_node
cavc_offset_loop_topology_get(cavc_offset_loop_topology const *topology, uint32_t index) {
  CAVC_ASSERT(topology, "null topology not allowed");
  CAVC_ASSERT(index < topology->nodes.size(), "index is out of topology node range");
  CAVC_BEGIN_TRY_CATCH
  return topology->nodes[index];
  CAVC_END_TRY_CATCH
}

cavc_parallel_offset_options cavc_parallel_offset_default_options(void) {
  CAVC_BEGIN_TRY_CATCH
  return cavc_parallel_offset_options{0, CAVC_OFFSET_JOIN_ROUND, CAVC_OFFSET_END_CAP_ROUND, 4.0};
  CAVC_END_TRY_CATCH
}

void cavc_parallel_offset(cavc_pline const *pline, cavc_real delta, cavc_pline_list **output,
                          cavc_parallel_offset_options options) {
  CAVC_ASSERT(pline, "null pline not allowed");
  CAVC_ASSERT(output, "null output not allowed");
  CAVC_BEGIN_TRY_CATCH
  *output = new cavc_pline_list();
  if (!std::isfinite(delta) || !is_valid_parallel_offset_options(options)) {
    return;
  }

  auto cpp_options = to_cpp_parallel_offset_options(options);
  auto results = cavc::parallelOffset(pline->data, delta, cpp_options);
  move_to_list(std::move(results), *output);
  CAVC_END_TRY_CATCH
}

void cavc_combine_plines(cavc_pline const *pline_a, cavc_pline const *pline_b, int combine_mode,
                         cavc_pline_list **remaining, cavc_pline_list **subtracted) {
  CAVC_ASSERT(pline_a, "null pline_a not allowed");
  CAVC_ASSERT(pline_b, "null pline_b not allowed");
  CAVC_ASSERT(combine_mode >= 0 && combine_mode <= 3, "combine_mode must be 0, 1, 2, or 3");
  CAVC_BEGIN_TRY_CATCH
  cavc::PlineCombineMode mode;
  switch (combine_mode) {
  case 0:
    mode = cavc::PlineCombineMode::Union;
    break;
  case 1:
    mode = cavc::PlineCombineMode::Exclude;
    break;
  case 2:
    mode = cavc::PlineCombineMode::Intersect;
    break;
  case 3:
    mode = cavc::PlineCombineMode::XOR;
    break;
  default:
    mode = cavc::PlineCombineMode::Union;
    break;
  }
  auto results = cavc::combinePolylines(pline_a->data, pline_b->data, mode);

  *remaining = new cavc_pline_list();
  *subtracted = new cavc_pline_list();
  move_to_list(std::move(results.remaining), *remaining);
  move_to_list(std::move(results.subtracted), *subtracted);
  CAVC_END_TRY_CATCH
}

cavc_real cavc_get_path_length(cavc_pline const *pline) {
  CAVC_ASSERT(pline, "null pline not allowed");
  CAVC_BEGIN_TRY_CATCH
  return cavc::getPathLength(pline->data);
  CAVC_END_TRY_CATCH
}

cavc_real cavc_get_area(cavc_pline const *pline) {
  CAVC_ASSERT(pline, "null pline not allowed");
  CAVC_BEGIN_TRY_CATCH
  return cavc::getArea(pline->data);
  CAVC_END_TRY_CATCH
}

int cavc_get_winding_number(cavc_pline const *pline, cavc_point point) {
  CAVC_ASSERT(pline, "null pline not allowed");
  CAVC_BEGIN_TRY_CATCH
  return cavc::getWindingNumber(pline->data, cavc::Vector2<cavc_real>(point.x, point.y));
  CAVC_END_TRY_CATCH
}

cavc_point_containment cavc_get_point_containment(cavc_pline const *pline, cavc_point point,
                                                  cavc_real boundary_epsilon) {
  CAVC_ASSERT(pline, "null pline not allowed");
  CAVC_BEGIN_TRY_CATCH
  auto containment = cavc::getPointContainment(
      pline->data, cavc::Vector2<cavc_real>(point.x, point.y), boundary_epsilon);
  switch (containment) {
  case cavc::PointContainment::Outside:
    return CAVC_POINT_OUTSIDE;
  case cavc::PointContainment::Inside:
    return CAVC_POINT_INSIDE;
  case cavc::PointContainment::OnBoundary:
    return CAVC_POINT_ON_BOUNDARY;
  }
  CAVC_ASSERT(false, "Unhandled point containment enum");
  return CAVC_POINT_OUTSIDE;
  CAVC_END_TRY_CATCH
}

void cavc_get_extents(cavc_pline const *pline, cavc_real *min_x, cavc_real *min_y, cavc_real *max_x,
                      cavc_real *max_y) {
  CAVC_ASSERT(pline, "null pline not allowed");
  CAVC_BEGIN_TRY_CATCH
  auto result = cavc::getExtents(pline->data);
  *min_x = result.xMin;
  *min_y = result.yMin;
  *max_x = result.xMax;
  *max_y = result.yMax;
  CAVC_END_TRY_CATCH
}

void cavc_get_closest_point(cavc_pline const *pline, cavc_point input_point,
                            uint32_t *closest_start_index, cavc_point *closest_point,
                            cavc_real *distance) {
  CAVC_ASSERT(pline, "null pline not allowed");
  CAVC_ASSERT(pline->data.size() != 0, "empty pline not allowed");
  CAVC_BEGIN_TRY_CATCH
  cavc::ClosestPoint<cavc_real> closestPoint(
      pline->data, cavc::Vector2<cavc_real>(input_point.x, input_point.y));
  *closest_start_index = static_cast<uint32_t>(closestPoint.index());
  *closest_point = cavc_point{closestPoint.point().x(), closestPoint.point().y()};
  *distance = closestPoint.distance();
  CAVC_END_TRY_CATCH
}

void cavc_get_tolerances(cavc_tolerances *tolerances_out) {
  CAVC_ASSERT(tolerances_out, "null tolerances_out not allowed");
  CAVC_BEGIN_TRY_CATCH
  *tolerances_out = to_api_tolerances(cavc::utils::getEpsilonConfig<cavc_real>());
  CAVC_END_TRY_CATCH
}

void cavc_set_tolerances(cavc_tolerances const *tolerances) {
  CAVC_ASSERT(tolerances, "null tolerances not allowed");
  CAVC_BEGIN_TRY_CATCH
  cavc::utils::setEpsilonConfig<cavc_real>(to_cpp_tolerances(*tolerances));
  CAVC_END_TRY_CATCH
}

void cavc_reset_tolerances(void) {
  CAVC_BEGIN_TRY_CATCH
  cavc::utils::resetEpsilonConfig<cavc_real>();
  CAVC_END_TRY_CATCH
}

// cavc_compound_pline APIs
// -------------------------

cavc_compound_pline *cavc_compound_pline_new(void) {
  CAVC_BEGIN_TRY_CATCH
  return new cavc_compound_pline{};
  CAVC_END_TRY_CATCH
}

void cavc_compound_pline_delete(cavc_compound_pline *cpline) {
  CAVC_BEGIN_TRY_CATCH
  delete cpline;
  CAVC_END_TRY_CATCH
}

void cavc_compound_pline_add_outer(cavc_compound_pline *cpline, cavc_pline const *pline) {
  CAVC_ASSERT(cpline, "null cpline not allowed");
  CAVC_ASSERT(pline, "null pline not allowed");
  CAVC_BEGIN_TRY_CATCH
  cpline->outerLoops.push_back(std::make_unique<cavc_pline>(cavc::Polyline<cavc_real>(pline->data)));
  CAVC_END_TRY_CATCH
}

void cavc_compound_pline_add_hole(cavc_compound_pline *cpline, cavc_pline const *pline) {
  CAVC_ASSERT(cpline, "null cpline not allowed");
  CAVC_ASSERT(pline, "null pline not allowed");
  CAVC_BEGIN_TRY_CATCH
  cpline->holeLoops.push_back(std::make_unique<cavc_pline>(cavc::Polyline<cavc_real>(pline->data)));
  CAVC_END_TRY_CATCH
}

uint32_t cavc_compound_pline_outer_count(cavc_compound_pline const *cpline) {
  CAVC_ASSERT(cpline, "null cpline not allowed");
  CAVC_BEGIN_TRY_CATCH
  return static_cast<uint32_t>(cpline->outerLoops.size());
  CAVC_END_TRY_CATCH
}

uint32_t cavc_compound_pline_hole_count(cavc_compound_pline const *cpline) {
  CAVC_ASSERT(cpline, "null cpline not allowed");
  CAVC_BEGIN_TRY_CATCH
  return static_cast<uint32_t>(cpline->holeLoops.size());
  CAVC_END_TRY_CATCH
}

cavc_pline const *cavc_compound_pline_get_outer(cavc_compound_pline const *cpline,
                                                uint32_t index) {
  CAVC_ASSERT(cpline, "null cpline not allowed");
  CAVC_ASSERT(index < cpline->outerLoops.size(), "index out of range");
  CAVC_BEGIN_TRY_CATCH
  return cpline->outerLoops[index].get();
  CAVC_END_TRY_CATCH
}

cavc_pline const *cavc_compound_pline_get_hole(cavc_compound_pline const *cpline,
                                               uint32_t index) {
  CAVC_ASSERT(cpline, "null cpline not allowed");
  CAVC_ASSERT(index < cpline->holeLoops.size(), "index out of range");
  CAVC_BEGIN_TRY_CATCH
  return cpline->holeLoops[index].get();
  CAVC_END_TRY_CATCH
}

void cavc_compound_pline_normalize(cavc_compound_pline *cpline) {
  CAVC_ASSERT(cpline, "null cpline not allowed");
  CAVC_BEGIN_TRY_CATCH
  for (auto &l : cpline->outerLoops) {
    if (l->data.isClosed() && cavc::getArea(l->data) < cavc_real(0)) {
      cavc::invertDirection(l->data);
    }
  }
  for (auto &l : cpline->holeLoops) {
    if (l->data.isClosed() && cavc::getArea(l->data) > cavc_real(0)) {
      cavc::invertDirection(l->data);
    }
  }
  CAVC_END_TRY_CATCH
}

void cavc_combine_compound_plines(cavc_compound_pline const *a, cavc_compound_pline const *b,
                                  int combine_mode, cavc_compound_pline **result_out) {
  CAVC_ASSERT(a, "null a not allowed");
  CAVC_ASSERT(b, "null b not allowed");
  CAVC_ASSERT(result_out, "null result_out not allowed");
  CAVC_BEGIN_TRY_CATCH
  cavc::PlineCombineMode mode;
  switch (combine_mode) {
  case 0:
    mode = cavc::PlineCombineMode::Union;
    break;
  case 1:
    mode = cavc::PlineCombineMode::Exclude;
    break;
  case 2:
    mode = cavc::PlineCombineMode::Intersect;
    break;
  case 3:
    mode = cavc::PlineCombineMode::XOR;
    break;
  default:
    *result_out = nullptr;
    return;
  }

  auto cppA = a->toCpp();
  auto cppB = b->toCpp();
  auto cppResult = cavc::combineCompoundPolylines(cppA, cppB, mode);

  auto *result = new cavc_compound_pline{};
  for (auto &loop : cppResult.outerLoops) {
    result->outerLoops.push_back(std::make_unique<cavc_pline>(std::move(loop)));
  }
  for (auto &loop : cppResult.holeLoops) {
    result->holeLoops.push_back(std::make_unique<cavc_pline>(std::move(loop)));
  }
  *result_out = result;
  CAVC_END_TRY_CATCH
}
