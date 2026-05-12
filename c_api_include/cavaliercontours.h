#ifndef CAVALIERCONTOURS_HPP
#define CAVALIERCONTOURS_HPP

#include <stdint.h>

#ifdef CAVC_STATIC_LIB
#define CAVC_API
#else
#if defined _WIN32 || defined __CYGWIN__
#ifdef CAVC_EXPORTS
#ifdef __GNUC__
#define CAVC_API __attribute__((dllexport))
#else
#define CAVC_API __declspec(dllexport)
#endif
#else
#ifdef __GNUC__
#define CAVC_API __attribute__((dllimport))
#else
#define CAVC_API __declspec(dllimport)
#endif
#endif
#else
#if __GNUC__ >= 4
#define CAVC_API __attribute__((visibility("default")))
#else
#define CAVC_API
#endif
#endif
#endif

#ifndef CAVC_REAL
#define CAVC_REAL double
#endif

typedef CAVC_REAL cavc_real;

#ifdef __cplusplus
extern "C" {
#endif

typedef struct cavc_pline cavc_pline;

typedef struct cavc_pline_list cavc_pline_list;

typedef struct cavc_spatial_index cavc_spatial_index;

typedef struct cavc_offset_loop_topology cavc_offset_loop_topology;

typedef struct cavc_vertex {
  cavc_real x;
  cavc_real y;
  cavc_real bulge;
} cavc_vertex;

typedef struct cavc_point {
  cavc_real x;
  cavc_real y;
} cavc_point;

typedef struct cavc_tolerances {
  cavc_real real_threshold;
  cavc_real real_precision;
  cavc_real slice_join_threshold;
  cavc_real offset_threshold;
} cavc_tolerances;

typedef enum cavc_point_containment {
  CAVC_POINT_OUTSIDE = 0,
  CAVC_POINT_INSIDE = 1,
  CAVC_POINT_ON_BOUNDARY = 2
} cavc_point_containment;

typedef enum cavc_closed_winding {
  CAVC_WINDING_KEEP = 0,
  CAVC_WINDING_COUNTER_CLOCKWISE = 1,
  CAVC_WINDING_CLOCKWISE = 2
} cavc_closed_winding;

typedef enum cavc_offset_join_type {
  CAVC_OFFSET_JOIN_ROUND = 0,
  CAVC_OFFSET_JOIN_MITER = 1,
  CAVC_OFFSET_JOIN_BEVEL = 2
} cavc_offset_join_type;

typedef enum cavc_offset_end_cap_type {
  CAVC_OFFSET_END_CAP_ROUND = 0,
  CAVC_OFFSET_END_CAP_SQUARE = 1,
  CAVC_OFFSET_END_CAP_BUTT = 2
} cavc_offset_end_cap_type;

typedef enum cavc_offset_loop_role {
  CAVC_OFFSET_LOOP_ROLE_OUTER = 0,
  CAVC_OFFSET_LOOP_ROLE_HOLE = 1
} cavc_offset_loop_role;

#define CAVC_OFFSET_LOOP_NO_PARENT UINT32_MAX

typedef struct cavc_offset_loop_topology_node {
  cavc_offset_loop_role role;
  uint32_t source_index;
  uint32_t parent_index;
} cavc_offset_loop_topology_node;

typedef struct cavc_parallel_offset_options {
  int may_have_self_intersects;
  cavc_offset_join_type join_type;
  cavc_offset_end_cap_type end_cap_type;
  cavc_real miter_limit;
} cavc_parallel_offset_options;

// Functions for working with cavc_pline

// Create/alloc a new cavc_pline initialized with data given. vertex_data may be null (vertex_count
// will still be used to set capacity), if is_closed is 0 then polyline will be open, otherwise it
// will be closed.
CAVC_API cavc_pline *cavc_pline_new(cavc_vertex const *vertex_data, uint32_t vertex_count,
                                    int is_closed);

// Delete/free a cavc_pline that was created with cavc_pline_new or released from a cavc_pline_list.
// NOTE: Do not call this function on cavc_pline's that are owned by a cavc_pline_list!
CAVC_API void cavc_pline_delete(cavc_pline *pline);

// Returns the current capacity of pline.
CAVC_API uint32_t cavc_pline_capacity(cavc_pline const *pline);

// Reserves memory for size vertexes to be stored in the pline if size is greater than current
// capacity, otherwise does nothing.
CAVC_API void cavc_pline_set_capacity(cavc_pline *pline, uint32_t size);

// Returns the current vertex count of a cavc_pline.
CAVC_API uint32_t cavc_pline_vertex_count(cavc_pline const *pline);

// Returns the vertex data of a cavc_pline by filling vertex_data as an array. NOTE: vertex_data
// must be large enough to hold the total vertex count! Call cavc_pline_vertex_count first to
// determine size required.
CAVC_API void cavc_pline_vertex_data(cavc_pline const *pline, cavc_vertex *vertex_data);

// Returns whether the cavc_pline is closed or not.
CAVC_API int cavc_pline_is_closed(cavc_pline const *pline);

// Sets the vertex data of a cavc_pline. Must pass in vertex_count to indicate how many vertexes to
// copy from the vertex_data array.
CAVC_API void cavc_pline_set_vertex_data(cavc_pline *pline, cavc_vertex const *vertex_data,
                                         uint32_t vertex_count);

// Adds a vertex to the cavc_pline.
CAVC_API void cavc_pline_add_vertex(cavc_pline *pline, cavc_vertex vertex);

// Removes a range of vertexes from the cavc_pline, starting at start_index and removing count
// number of vertexes. No bounds checking is performed (ensure start_index + count <=
// pline_vertex_count).
CAVC_API void cavc_pline_remove_range(cavc_pline *pline, uint32_t start_index, uint32_t count);

// Clears the cavc_pline (capacity is left unchanged).
CAVC_API void cavc_pline_clear(cavc_pline *pline);

// Sets the polyline to be closed or open, if is_closed is 0 then polyline will be open, otherwise
// it will be closed.
CAVC_API void cavc_pline_set_is_closed(cavc_pline *pline, int is_closed);

// Invert the direction of the polyline (preserving geometric shape).
CAVC_API void cavc_pline_invert_direction(cavc_pline *pline);

// Create a new polyline with repeated positions removed from the input.
CAVC_API void cavc_pline_prune_singularities(cavc_pline const *pline, cavc_real epsilon,
                                             cavc_pline **output);

// Remove redundant vertexes from the polyline in place. Redundant vertexes can come from repeated
// positions, collinear line segments, or concentric arc segments that can be merged without
// changing the represented curve.
CAVC_API void cavc_pline_remove_redundant(cavc_pline *pline, cavc_real epsilon);

// Create a new polyline with all arc segments approximated by line segments.
CAVC_API void cavc_pline_convert_arcs_to_lines(cavc_pline const *pline, cavc_real error,
                                               cavc_pline **output);

// Create a normalized polyline by pruning singularities and optionally normalizing winding.
// Winding normalization only applies to closed polylines.
CAVC_API void cavc_pline_normalize(cavc_pline const *pline, cavc_real epsilon,
                                   cavc_closed_winding closed_winding, cavc_pline **output);

// Functions for working with cavc_pline_list

// Delete/free a cavc_pline_list, this will also delete all elements in the list.
CAVC_API void cavc_pline_list_delete(cavc_pline_list *pline_list);

// Get the element count of the cavc_pline_list.
CAVC_API uint32_t cavc_pline_list_count(cavc_pline_list const *pline_list);

// Get a cavc_pline from the cavc_pline_list. No bounds checking is performed (ensure index <
// cavc_pline_list count). NOTE: The cavc_pline is still owned by the list!
CAVC_API cavc_pline *cavc_pline_list_get(cavc_pline_list const *pline_list, uint32_t index);

// Release a cavc_pline from the cavc_pline_list's ownership, returning the cavc_pline that was
// released and removing it from the list. NOTE: cavc_pline_delete must now be called on the
// released cavc_pline!
CAVC_API cavc_pline *cavc_pline_list_release(cavc_pline_list *pline_list, uint32_t index);

// Functions for working with cavc_spatial_index

// Build an approximate segment spatial index from the input polyline. Polyline must have at least
// 2 vertexes.
CAVC_API cavc_spatial_index *cavc_spatial_index_create(cavc_pline const *pline);

// Delete/free a cavc_spatial_index.
CAVC_API void cavc_spatial_index_delete(cavc_spatial_index *spatial_index);

// Returns the number of indexed segments in cavc_spatial_index.
CAVC_API uint32_t cavc_spatial_index_item_count(cavc_spatial_index const *spatial_index);

// Query the spatial index and return how many segment indexes overlap the given bounding box.
CAVC_API uint32_t cavc_spatial_index_query_count(cavc_spatial_index const *spatial_index,
                                                 cavc_real min_x, cavc_real min_y, cavc_real max_x,
                                                 cavc_real max_y);

// Query the spatial index and write overlapping segment indexes to results_out.
// results_out must have enough storage for cavc_spatial_index_query_count(...) values.
CAVC_API void cavc_spatial_index_query(cavc_spatial_index const *spatial_index, cavc_real min_x,
                                       cavc_real min_y, cavc_real max_x, cavc_real max_y,
                                       uint32_t *results_out);

// Build stable topology ordering and parent-child mapping for offset loops.
// The input loops are split by role (ccw=outer, cw=hole).
// - output order is stable geometry order (area-desc + bbox key)
// - source_index is within the corresponding input role array
// - parent_index is index in the returned topology order, or CAVC_OFFSET_LOOP_NO_PARENT
// All input loops must be closed and contain at least 2 vertices.
// boundary_epsilon must be >= 0.
CAVC_API cavc_offset_loop_topology *
cavc_offset_loop_topology_build(cavc_pline const *const *ccw_loops, uint32_t ccw_loop_count,
                                cavc_pline const *const *cw_loops, uint32_t cw_loop_count,
                                cavc_real boundary_epsilon);

// Delete/free a cavc_offset_loop_topology.
CAVC_API void cavc_offset_loop_topology_delete(cavc_offset_loop_topology *topology);

// Get total node count from cavc_offset_loop_topology.
CAVC_API uint32_t cavc_offset_loop_topology_count(cavc_offset_loop_topology const *topology);

// Get a topology node at index.
CAVC_API cavc_offset_loop_topology_node
cavc_offset_loop_topology_get(cavc_offset_loop_topology const *topology, uint32_t index);

// Algorithm functions

// Returns default parallel offset options:
// - may_have_self_intersects = 0
// - join_type = CAVC_OFFSET_JOIN_ROUND
// - end_cap_type = CAVC_OFFSET_END_CAP_ROUND
// - miter_limit = 4.0
CAVC_API cavc_parallel_offset_options cavc_parallel_offset_default_options(void);

// Generates the parallel offset of a polyline. delta is the offset delta, output is filled with
// the result. Invalid offset options (unknown join/end-cap enum, non-finite delta, or
// miter_limit < 1 when join_type is miter) produce an allocated empty output list instead of
// silently falling back to another join/end-cap type.
// Join support:
// - round: the original rounded raw-offset join behavior for lines and arcs
// - miter: intersects adjacent raw-offset tangents and clips to miter_limit; when the miter limit
//   is exceeded, the local join is beveled as standard miter-limit behavior
// - bevel: connects adjacent raw-offset endpoints directly with a bevel segment
// - arc-involved non-round joins use endpoint tangents; collapsed-arc joins degrade locally to
//   bevel because no stable tangent-preserving miter exists for the collapsed segment
// End cap support for open polylines:
// - round: endpoint-circle clipping
// - square: endpoint-tangent extension by abs(delta), including line/arc start/end segments
// - butt: endpoint-normal line clipping
// These end caps operate on the returned open offset curve. They do not generate a closed outline
// cap face.
CAVC_API void cavc_parallel_offset(cavc_pline const *pline, cavc_real delta,
                                   cavc_pline_list **output, cavc_parallel_offset_options options);

// Combines two non-self intersecting closed polylines, pline_a and pline_b.
// For union combine_mode = 0
// For exclude combine_mode = 1
// For intersect combine_mode = 2
// For XOR combine_mode = 3
// If combine_mode is any other value then no output parameters are filled.
// remaining is filled with the closed polylines that remain after combining, subtracted is filled
// with closed polylines that represent subtracted space (in the case of islands after a union or
// exclude).
// If pline_a or pline_b is an open polyline or has self intersects then the result is undefined.
CAVC_API void cavc_combine_plines(cavc_pline const *pline_a, cavc_pline const *pline_b,
                                  int combine_mode, cavc_pline_list **remaining,
                                  cavc_pline_list **subtracted);

// Returns the path length of the cavc_pline given. If pline vertex count is less than 2 then 0 is
// returned.
CAVC_API cavc_real cavc_get_path_length(cavc_pline const *pline);

// Returns the signed area of the cavc_pline given. If pline is open or vertex count is less than 2
// then 0 is returned. If pline goes clockwise then a negative signed area is returned, otherwise a
// positive signed area is returned.
CAVC_API cavc_real cavc_get_area(cavc_pline const *pline);

// Compute the winding number of the 2d point relative to the pline. If pline is open or vertex
// count is less than 2 then 0 is returned. For more on winding number see:
// https://en.wikipedia.org/wiki/Winding_number
CAVC_API int cavc_get_winding_number(cavc_pline const *pline, cavc_point point);

// Classify a point against a polyline as Outside/Inside/OnBoundary.
// For open polylines, only Outside or OnBoundary are returned.
CAVC_API cavc_point_containment cavc_get_point_containment(cavc_pline const *pline,
                                                           cavc_point point,
                                                           cavc_real boundary_epsilon);

// Compute the axis aligned extents of the pline, results are written to min_x, min_y, max_x, and
// max_y output parameters. If pline is empty then min_x and min_y are filled with positive infinity
// and max_x and max_y are filled with negative infinity.
CAVC_API void cavc_get_extents(cavc_pline const *pline, cavc_real *min_x, cavc_real *min_y,
                               cavc_real *max_x, cavc_real *max_y);

// Finds the closest point on the polyline to the point given. closest_start_index is filled with
// the starting vertex index of the segment that the closest point is on. closest_point is filled
// with the closest point on the polyline. distance is filled with the distance from the point given
// to the closest point. pline must not be empty.
CAVC_API void cavc_get_closest_point(cavc_pline const *pline, cavc_point input_point,
                                     uint32_t *closest_start_index, cavc_point *closest_point,
                                     cavc_real *distance);

// Get the current geometry tolerances used by algorithms.
CAVC_API void cavc_get_tolerances(cavc_tolerances *tolerances_out);

// Override geometry tolerances used by algorithms.
CAVC_API void cavc_set_tolerances(cavc_tolerances const *tolerances);

// Reset geometry tolerances to default values.
CAVC_API void cavc_reset_tolerances(void);

// Functions for working with cavc_compound_pline
// A compound polyline represents a filled 2-D region with outer loops (CCW) and hole loops (CW).
// A point is "inside" the compound path when the sum of winding numbers across all loops is
// non-zero (the standard nonzero fill rule).

typedef struct cavc_compound_pline cavc_compound_pline;

// Create a new, empty compound polyline. Must be freed with cavc_compound_pline_delete.
CAVC_API cavc_compound_pline *cavc_compound_pline_new(void);

// Delete/free a compound polyline.
CAVC_API void cavc_compound_pline_delete(cavc_compound_pline *cpline);

// Add an outer (CCW) loop to the compound polyline. A copy of pline is stored.
CAVC_API void cavc_compound_pline_add_outer(cavc_compound_pline *cpline, cavc_pline const *pline);

// Add a hole (CW) loop to the compound polyline. A copy of pline is stored.
CAVC_API void cavc_compound_pline_add_hole(cavc_compound_pline *cpline, cavc_pline const *pline);

// Return the number of outer loops.
CAVC_API uint32_t cavc_compound_pline_outer_count(cavc_compound_pline const *cpline);

// Return the number of hole loops.
CAVC_API uint32_t cavc_compound_pline_hole_count(cavc_compound_pline const *cpline);

// Get an outer loop by index. The returned pointer is owned by the compound polyline.
// No bounds checking is performed (ensure index < cavc_compound_pline_outer_count).
CAVC_API cavc_pline const *cavc_compound_pline_get_outer(cavc_compound_pline const *cpline,
                                                         uint32_t index);

// Get a hole loop by index. The returned pointer is owned by the compound polyline.
// No bounds checking is performed (ensure index < cavc_compound_pline_hole_count).
CAVC_API cavc_pline const *cavc_compound_pline_get_hole(cavc_compound_pline const *cpline,
                                                        uint32_t index);

// Normalize winding in-place: outer loops are made CCW (positive area), hole loops are made CW
// (negative area). Call this before combining if winding may be incorrect.
CAVC_API void cavc_compound_pline_normalize(cavc_compound_pline *cpline);

// Combine two compound polylines using the given boolean mode.
// combine_mode: 0=Union, 1=Exclude (a-b), 2=Intersect, 3=XOR.
// If combine_mode is any other value, *result_out is set to NULL.
// The caller owns the returned cavc_compound_pline and must free it with
// cavc_compound_pline_delete.
// All input loops must be closed and have at least 2 vertices; outer loops should be CCW and hole
// loops CW (call cavc_compound_pline_normalize if unsure).
CAVC_API void cavc_combine_compound_plines(cavc_compound_pline const *a,
                                           cavc_compound_pline const *b, int combine_mode,
                                           cavc_compound_pline **result_out);

#ifdef __cplusplus
}
#endif

#endif // CAVALIERCONTOURS_HPP
