import vtk
import os
import math
from collections import defaultdict
from time import perf_counter

# Refinement version: exact VTK clipping + permissive retry/skip +
# local cross-section shape, small-opening, and overlap QC.

# Refinement: every local cut can retry at nearby positions along its
# associated centerline branch. In exploratory/permissive mode, a cut that
# still looks geometrically bad is skipped instead of aborting the full run.


# ============================================================
# USER SETTINGS
# ============================================================
File_location = r"C:\Users\Admin\Documents\3D_Slicer_Files\JL_9-gen-test"
# Full, untrimmed centerline tree.
FULL_CENTERLINE_FILE = r"C:\Users\Admin\Documents\3D_Slicer_Files\JL_4-gen-test\Network_model_V2.vtp"

# Centerline tree after trimming to the extent you want to keep.
TRIMMED_CENTERLINE_FILE = File_location + "\JL_9-gen-centerlines.vtp"

# Original pulmonary artery surface to prune.
SURFACE_FILE = r"C:\Users\Admin\Documents\3D_Slicer_Files\JL_4-gen-test\JL_Full_Tree_Meshed.vtp"

# Exact-clipped PA tree with the new outlet rings left OPEN.
# This is the most useful QC output because you can inspect the new
# plane-cut boundary before any cap is added.
OUTPUT_OPEN_FILE = File_location + "\Pulmonary_Artery_Surface_Pruned_Exact_Open.vtp"

# Exact-clipped PA tree with flat triangulated caps added at the new cuts.
OUTPUT_CAPPED_FILE = File_location + "\Pulmonary_Artery_Surface_Pruned_Exact_Capped.vtp"

# Matching SimVascular face-name file written beside the FINAL capped model.
# Format is intentionally the same as the supplied .vtp.facenames example:
#   set gPolyDataFaceNames(1) Wall
#   set gPolyDataFaceNames(2) noname_2
#   ...
OUTPUT_CAPPED_FACENAMES = OUTPUT_CAPPED_FILE + ".facenames"

# Debug copy of the final boundary loops that are capped after ALL clipping
# and mesh-quality cleanup has finished.
FINAL_CAP_LOOPS_DEBUG_FILE = File_location + "\Pulmonary_Artery_Final_Cap_Loops_Debug.vtp"

# Build caps only after the final endpoint-box and quality-cleanup passes.
# This is deliberately different from the older workflow, which capped the
# pre-cleanup cut loops and could therefore create caps that no longer matched
# the actual final opening.
FINAL_CAP_WALL_FACE_ID = 1

# Face 2 is reserved for the inlet, regardless of its Sidedness.
FINAL_CAP_INLET_NAME = "INLET"

# The trimmed centerline should carry PointData Sidedness as a vtkStringArray
# (or another array with string-like values) containing Left, Right, or any
# third category. Any other/missing value is called Other.
SIDEDNESS_ARRAY = "Sidedness"
LEFT_CAP_PREFIX = "Left"
RIGHT_CAP_PREFIX = "Right"
OTHER_CAP_PREFIX = "Other"
# False: still run with older centerline files; unmatched caps become Other#.
REQUIRE_SIDEDNESS_ARRAY = False

# The loop closest to the accepted inlet cut MUST be close enough to really
# be that new inlet opening. Do not silently rename some distant outlet INLET.
INLET_CAP_MAX_DISTANCE_RADII = 3.0

# Radial-ring caps avoid the long polygon diagonals / tiny ear-clipping slivers
# that can appear with vtkContourTriangulator on large irregular loops.
# 0 = simple center fan; 1-3 adds concentric interior rings and usually gives
# much more regular planar triangles.
FINAL_CAP_INTERIOR_RINGS = 2

# Warn if a generated cap triangle is extremely skinny by simple edge ratio.
FINAL_CAP_EDGE_RATIO_WARN = 25.0

# Debug output: the intended circular search/QC disk for every cut.
CUTS_DEBUG_FILE = File_location  + "\Pulmonary_Artery_Cuts_Debug.vtp"

# Debug output: the ACTUAL plane/surface intersection loop selected for
# every cut. These should look like clean oval/circular outlet rings.
CUT_LOOPS_DEBUG_FILE = File_location + "\Pulmonary_Artery_Exact_Cut_Loops_Debug.vtp"

# Debug output: after all exact cuts, but BEFORE connectivity removes the
# detached distal vessel pieces.
SEVERED_DEBUG_FILE = File_location + "\Pulmonary_Artery_Exact_Severed_Debug.vtp"

# Debug output: intended cutting circles for EVERY failed/rejected retry
# position. Load this together with the centerlines/surface to see exactly
# where the script was unhappy without stopping the rest of the run.
FAILED_CUT_ATTEMPTS_DEBUG_FILE = File_location  + "\Pulmonary_Artery_Failed_Cut_Attempts_Debug.vtp"

# Rejected inlet contours and their *confirmed* 2-D segment intersection
# locations. View these files together with the source PA surface in ParaView.
INLET_REJECTED_CONTOURS_DEBUG_FILE = File_location + "\Pulmonary_Artery_Inlet_Rejected_Contours_Debug.vtp"
INLET_INTERSECTION_POINTS_DEBUG_FILE = File_location + "\Pulmonary_Artery_Inlet_Intersection_Points_Debug.vtp"

# Validate the inlet on the original triangulated surface BEFORE Stage 1.
# If none of the nine planes is plausible, fail early and write inlet debug
# instead of wasting time on all pruning/outlet cuts. Stage 3 rechecks the
# inlet on the actual post-pruning geometry before applying it.
ENABLE_EARLY_INLET_PREFLIGHT = True
INLET_PREFLIGHT_FAIL_FAST = True



# ============================================================
# CENTERLINE / CUT SETTINGS
# ============================================================

# In the supplied trimmed centerline file, CellIds map back to the
# original full centerline cells.
CELL_ID_ARRAY = "CellIds"

# Radius array on the full centerline.
RADIUS_ARRAY = "Radius"

# NEW inlet plane: halfway along the CENTERLINE ARC LENGTH from the IsInlet
# point to the FIRST anatomical bifurcation (not halfway by point index).
ENABLE_INLET_CUT = True
INLET_CUT_PATH_FRACTION = 0.40

# If the midpoint lies in a poor local mesh region, retry small shifts along
# the same pre-bifurcation trunk; never move directly onto the bifurcation.
INLET_RETRY_PATH_FRACTIONS = (0.45, 0.55, 0.40, 0.60, 0.35, 0.65, 0.30, 0.70)

# Use a point safely on the retained/branchward side of the cut for every
# downstream connectivity operation, not the original IsInlet coordinate.
INLET_RETAINED_ANCHOR_PATH_FRACTION = 0.80

# Do not write an incorrectly labeled model if the new inlet was not cut.
REQUIRE_SUCCESSFUL_INLET_CUT = True

# Used as a QC radius for deciding whether the plane intersection loop
# selected by VTK is plausibly the vessel surrounding this centerline.
# The exact cut itself is made from the true plane/surface intersection;
# it is NOT approximated by this circle.
CUT_RADIUS_FACTOR = 1.75

# When a side branch disappears from the trimmed centerline, move the
# plane slightly INTO the removed branch so we are not cutting directly
# through the bifurcation crotch.
#
# 0.0 = exactly at the branch junction
# 0.35 = 35% of the way from the junction to the next point in the
#        removed centerline branch
SIDE_BRANCH_OFFSET_FRACTION = 0.25

# If a proposed cut is geometrically awkward, retry at nearby locations
# along the SAME centerline cell before giving up. Offsets are in local
# vessel radii relative to the original plane position.
#
# For terminal and true IsOutlet cuts, positive offsets move UPSTREAM toward
# the parent vessel. Those are intentionally tried first because an ugly
# half-branch/compound cut is usually improved by moving toward the parent.
#
# For side-branch pruning cuts, positive offsets move farther INTO the branch
# being removed, which remains preferable near a bifurcation crotch.
TERMINAL_RETRY_RADIUS_OFFSETS = (
    0.25,
    0.50,
    0.75,
    1.00,
    1.50,
    2.00,
    3.00,
    5.00,
)

SIDE_BRANCH_RETRY_RADIUS_OFFSETS = (
    0.25,
    0.50,
    0.75,
    1.00,
    1.50,
    2.00,
    -0.25,
    -1.50,
)

# True IsOutlet retries are intentionally NOT scaled by endpoint Radius.
# Some centerline Radius values collapse near old terminal tips, making
# radius-based retries microscopic.  Instead, retry at fixed fractions of
# the entire final centerline cell length, measured from the IsOutlet back
# toward its parent branch.  Only positions farther upstream than the base
# cut are retained.
OUTLET_RETRY_PATH_FRACTIONS = (
    0.10,
    0.20,
    0.30,
    0.40,
    0.55,
    0.70,
    0.85,
    0.90,
)

# Never move an alternate plane beyond this fraction of the source polyline
# cell, which helps avoid wandering into the next anatomical junction.
CUT_RETRY_MAX_PATH_FRACTION = 0.90

# Suspicious loops are retryable for ALL cut types. This now includes the
# synthetic terminal cuts because those can also land in the middle of a
# compound/branching region and produce the mangled openings we want to avoid.
RETRY_SUSPICIOUS_SIDE_BRANCH_LOOPS = True
RETRY_SUSPICIOUS_TERMINAL_LOOPS = True
RETRY_SUSPICIOUS_OUTLET_LOOPS = True

# -------------------------------------------------------------------------
# Cross-section SHAPE / AREA quality control
# -------------------------------------------------------------------------
#
# A global average cap area is not ideal for a pulmonary tree because a
# proximal PA and a distal PA are supposed to have very different areas.
# Instead, normalize each proposed loop by the LOCAL centerline Radius:
#
#       expected area = pi * Radius^2
#
# This catches both abnormally tiny and abnormally huge cuts without
# penalizing legitimate vessel-size differences throughout the tree.
ENABLE_LOOP_SHAPE_QC = True

# Broad limits are intentional. We only want to reject obvious outliers.
MIN_LOOP_AREA_TO_EXPECTED_RATIO = 0.30
MAX_LOOP_AREA_TO_EXPECTED_RATIO = 3.50

# Circularity = 4*pi*A / P^2. A circle is 1.0; a normal oval is lower but
# should still be well above a strongly branched/lopsided contour.
MIN_LOOP_CIRCULARITY = 0.40

# Shape radii are measured from the loop's own point-centroid (not from the
# centerline), making this sensitive to a one-sided protrusion while allowing
# an ordinary oval. Values much above ~2 usually indicate a branch-like arm.
MAX_LOOP_RADIAL_ASYMMETRY = 2.25

# Large offset of the selected contour from the centerline is another sign
# that the plane picked up the wrong/compound region. Expressed in local radii.
MAX_LOOP_CENTROID_OFFSET_TO_RADIUS = 1.50

# -------------------------------------------------------------------------
# Small-opening / overlapping-contour quality control
# -------------------------------------------------------------------------
#
# Radius can collapse artificially at a terminal cap, so the ordinary
# A/(pi*r^2) test can sometimes accept a physically tiny opening simply
# because the centerline Radius at that exact point is also tiny.  Each cut
# therefore carries a robust reference radius taken from the WHOLE source
# centerline cell.  The equivalent radius of the proposed opening must be at
# least this fraction of that branch-scale reference radius.
MIN_OPENING_EFFECTIVE_RADIUS_TO_REFERENCE = 0.45

# A clean vtkCutter contour contains each geometric line segment once.
# Coincident/overlapping surface triangles can contribute duplicate cutter
# segments, making the total raw segment length appreciably larger than the
# perimeter of the stripped contour.
MAX_LOOP_SEGMENT_REDUNDANCY_RATIO = 1.10

# A valid vessel cross-section should be a simple non-self-intersecting loop.
# Any crossing of non-adjacent projected contour edges is treated as evidence
# of overlapping/self-intersecting elements and triggers an upstream retry.
MAX_LOOP_SELF_INTERSECTIONS = 0

# Non-adjacent contour vertices that collapse onto essentially the same
# location are another indication of a self-touching/overlapping opening.
# Tolerance is expressed as a fraction of the branch reference radius.
OVERLAP_VERTEX_TOLERANCE_FACTOR = 1.0e-4

# Use a MUCH tighter tolerance for deciding whether contour EDGES really
# intersect than for flagging non-adjacent vertices as suspiciously close.
# A 2-D orientation determinant has units of length^2, so it cannot be
# compared to a length tolerance directly. The geometry below compares
# perpendicular distances, with an independent small spatial tolerance.
LOOP_INTERSECTION_TOLERANCE_RADIUS_FACTOR = 1.0e-8
LOOP_INTERSECTION_TOLERANCE_ABSOLUTE = 1.0e-10

# Exploratory mode requested for refinement: if all candidate locations for
# an individual cut fail, skip ONLY that cut and continue processing the
# remaining tree. The skipped/failed attempt circles are written to a debug
# VTP so they can be inspected afterward.
SKIP_FAILED_PRUNING_CUTS = True
SKIP_FAILED_OUTLET_CUTS = True

# Keep normal runs readable.  Set True only when diagnosing a particular cut.
VERBOSE_CUT_LOG = False
PROGRESS_PRINT_EVERY = 50

# -------------------------------------------------------------------------
# Performance settings
# -------------------------------------------------------------------------
# Exact cuts used to clean + triangulate the ENTIRE ~600k-cell surface after
# every individual outlet. That is safe but extremely expensive for ~1000
# outlets. The local cut itself is still exact; we now defer the global weld/
# triangulation and do it in batches. Set to 1 to recover the old behavior.
EXACT_CUT_GLOBAL_CLEAN_INTERVAL = 25

# The final endpoint safety sweep similarly used to run vtkCleanPolyData and
# vtkTriangleFilter after every one of ~1000 boxes. The box clipping remains
# exact, but global cleanup is deferred to this interval and once at the end.
# Set to 1 for the old behavior.
END_BOX_GLOBAL_CLEAN_INTERVAL = 100

# Old versions built a temporary flat cap for every exact cut even though the
# final workflow discards those caps and rebuilds better radial-ring caps only
# after ALL cleanup. Leave this False for a substantial amount of wasted work
# to be skipped.
BUILD_LEGACY_PER_CUT_CAPS = False

# If the retained centerline terminates at a pruning junction, make one
# transverse cut using the last retained branch direction. Positive values
# move the plane slightly upstream into the retained parent branch.
#
# This applies ONLY to the synthetic pruning-terminal cuts created by
# comparing the full centerline with the trimmed centerline.
TERMINAL_OFFSET_FRACTION = 0.0

# Separate offset for the TRUE IsOutlet points on the trimmed centerline.
# The IsOutlet point is the terminal centerline point, and the cutting-plane
# normal is defined by the vector from the immediately upstream point toward
# that IsOutlet point.
#
# 0.0 = plane exactly through the IsOutlet point
# 0.5 = halfway from the IsOutlet point back toward the previous point
# 1.0 = through the immediately previous centerline point
#
# A modest positive value is usually preferable because it removes the old
# distal tip/cap and leaves a real vessel cross-section for the exact VTK cut.
OUTLET_OFFSET_FRACTION = 0.75

# For IsOutlet cuts only: after vtkCutter finds the closest plane/surface
# contour, require that contour to be reasonably local to the centerline
# endpoint. This prevents an outlet that was already removed by the pruning
# stage from accidentally selecting a distant vessel crossed by the same
# infinite plane. If no local contour exists, that IsOutlet is simply skipped.
OUTLET_LOCAL_SEARCH_FACTOR = 2.5

# Coordinate tolerance used to group coincident centerline endpoints that
# use different VTK point IDs at the same anatomical bifurcation.
COORDINATE_TOLERANCE = 1e-5

# Small geometric tolerance.
PLANE_EPS = 1e-10

# Debug circle resolution only. This does NOT control the actual cut.
DEBUG_CIRCLE_RESOLUTION = 64

# If the selected plane/surface intersection loop extends much farther
# from the centerline than expected, it is probably passing through a
# bifurcation/compound region instead of a clean daughter-vessel section.
# Side-branch cuts will automatically retry farther downstream when
# RETRY_SUSPICIOUS_SIDE_BRANCH_LOOPS is True.
MAX_LOOP_TO_EXPECTED_RADIUS_RATIO = 1.75
FAIL_ON_SUSPICIOUS_LOOP = False

# Tiny tolerance used when welding the newly interpolated VTK clip points
# back to the untouched surface.
MERGE_TOLERANCE = 1e-8

# Add flat caps from the exact cut loops.
WRITE_CAPPED_OUTPUT = True

# After generating new geometry, old GlobalNodeID/GlobalElementID values
# can no longer remain globally valid because clipping creates new points
# and triangles. If True, replace/add them with fresh sequential IDs.
REGENERATE_GLOBAL_IDS = True



# ============================================================
# FINAL MESH-QUALITY / SPIKE CLEANUP SETTINGS
# ============================================================
#
# After the centerline-driven pruning/capping workflow is complete, run the
# VTK equivalent of ParaView's "Mesh Quality" filter on the final surface.
# Connected groups of triangles whose Quality exceeds the threshold are
# treated as suspect spike / sliver regions.
#
# For each suspect group:
#   1) find the nearest location on the TRIMMED centerline,
#   2) use the local centerline tangent to orient a vtkBox,
#   3) size the box from the local centerline Radius,
#   4) clip AWAY the inside of that box,
#   5) keep only the surface still connected to the pulmonary-artery inlet.
#
# Because the box spans the whole local vessel cross-section, the upstream
# box face acts like a clean transverse plane cut and the distal spike/piece
# becomes disconnected and is discarded.
ENABLE_MESH_QUALITY_CLEANUP = True

# SimVascular meshing failure threshold requested for this workflow.
MESH_QUALITY_THRESHOLD = 1.0e5

# ParaView/VTK Mesh Quality supports several triangle measures.  "AspectRatio"
# is used here because very skinny/sliver triangles become very large.
# If your ParaView Mesh Quality panel is using a different measure, change
# this string to one of the supported names in configure_mesh_quality_filter().
MESH_QUALITY_TRIANGLE_MEASURE = "AspectRatio"

# Mesh-quality cleanup is ALWAYS run on the OPEN surface.  The capped model
# is written only as an optional diagnostic/output and can never feed back
# into OUTPUT_QUALITY_CLEANED_FILE.

# Oriented removal box dimensions.  These are HALF-widths expressed as
# multiples of the local centerline Radius.
#
# radial half-width must exceed the vessel radius so the box fully severs the
# vessel rather than merely punching a rectangular side hole.
QUALITY_BOX_RADIAL_HALF_WIDTH_FACTOR = 1.75

# The box is centered at the nearest centerline location.  Therefore this
# axial half-length also controls how far toward the parent vessel the new
# retained opening is moved away from the bad elements.
QUALITY_BOX_AXIAL_HALF_LENGTH_FACTOR = 1.25

# Multiple disconnected bad-cell clusters can belong to the same spike.
# If their proposed box centers are closer than this many local radii, merge
# them into a single box to avoid repeatedly chopping the same location.
QUALITY_BOX_MERGE_DISTANCE_FACTOR = 1.25

# -------------------------------------------------------------------------
# FINAL END-BOX SWEEP
# -------------------------------------------------------------------------
# After all centerline-driven exact cuts are finished, make one final
# oriented Box Clip pass at EVERY intended terminal opening. This includes
# removed side-branch ends, synthetic pruning-terminal ends, true IsOutlet
# ends, and (when present) points marked by an IsTerminalEnd PointData array.
#
# The box is deliberately ONE-SIDED with respect to the desired cut plane:
# its upstream face sits at the accepted/planned cut location and the box
# extends only downstream. This removes overlapping triangles, tiny residual
# caps, spikes, and distal fragments without unnecessarily moving the opening
# farther toward the parent vessel.
ENABLE_FINAL_END_BOX_SWEEP = False

# Requested: endpoint boxes are 2x the current mesh-quality box size.
END_BOX_SIZE_MULTIPLIER = 2.0
END_BOX_RADIAL_HALF_WIDTH_FACTOR = (
    QUALITY_BOX_RADIAL_HALF_WIDTH_FACTOR * END_BOX_SIZE_MULTIPLIER
)
END_BOX_AXIAL_HALF_LENGTH_FACTOR = (
    QUALITY_BOX_AXIAL_HALF_LENGTH_FACTOR * END_BOX_SIZE_MULTIPLIER
)

# Prefer the robust branch-scale Radius when an endpoint Radius collapses near
# an old cap. This prevents microscopic boxes on otherwise normal branches.
END_BOX_USE_REFERENCE_RADIUS = True

# If the nearest surviving surface is farther than this many robust radii
# from the desired cut plane, assume the end was already removed and skip it.
END_BOX_LOCALITY_FACTOR = 3.0

# Optional explicit point marker. If absent, the script simply uses the
# terminal definitions already inferred from the full-vs-trimmed tree.
TERMINAL_END_ARRAY = "IsTerminalEnd"

# Debug surface for the endpoint boxes actually attempted.
END_BOXES_DEBUG_FILE = File_location  + "\Pulmonary_Artery_Final_End_Boxes_Debug.vtp"

# Open surface immediately after the endpoint-box sweep and before the final
# Mesh Quality cleanup. Useful for isolating what the endpoint boxes changed.
END_BOX_OPEN_OUTPUT_FILE = File_location  + "\Pulmonary_Artery_Surface_End_Box_Cleaned_Open.vtp"

# A single cleanup pass is intentionally conservative.  The script reports
# any Quality > threshold cells that remain afterward so they can be inspected
# before deciding whether a second pass is desirable.
MESH_QUALITY_CLEANUP_PASSES = 2

# Final open output after bad-element boxes have been clipped away.  This is
# intentionally OPEN at any quality-cleanup locations.
OUTPUT_QUALITY_CLEANED_FILE = File_location  + "\Pulmonary_Artery_Surface_Quality_Cleaned_Open.vtp"

# Debug: all cells above MESH_QUALITY_THRESHOLD BEFORE cleanup.
BAD_QUALITY_CELLS_DEBUG_FILE = File_location + "\Pulmonary_Artery_Bad_Quality_Cells_Debug.vtp"

# Debug: oriented boxes actually used for quality cleanup.
QUALITY_BOXES_DEBUG_FILE = File_location + "\Pulmonary_Artery_Quality_Clip_Boxes_Debug.vtp"

# Debug: any cells still above threshold AFTER the cleanup pass(es).
REMAINING_BAD_QUALITY_CELLS_DEBUG_FILE = File_location + "\Pulmonary_Artery_Remaining_Bad_Quality_Cells_Debug.vtp"


# ============================================================
# BASIC VTK I/O
# ============================================================

def read_polydata(filename):

    if not os.path.exists(filename):
        raise RuntimeError(
            "File does not exist:\n{}".format(filename)
        )

    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(filename)
    reader.Update()

    poly = vtk.vtkPolyData()
    poly.DeepCopy(reader.GetOutput())

    if poly.GetNumberOfPoints() == 0:
        raise RuntimeError(
            "File contains no points:\n{}".format(filename)
        )

    return poly


def write_polydata(polydata, filename):

    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(filename)
    writer.SetInputData(polydata)

    if writer.Write() != 1:
        raise RuntimeError(
            "Failed to write:\n{}".format(filename)
        )


# ============================================================
# VECTOR HELPERS
# ============================================================

def vec_add(a, b):

    return (
        a[0] + b[0],
        a[1] + b[1],
        a[2] + b[2]
    )


def vec_sub(a, b):

    return (
        a[0] - b[0],
        a[1] - b[1],
        a[2] - b[2]
    )


def vec_scale(a, s):

    return (
        a[0] * s,
        a[1] * s,
        a[2] * s
    )


def dot(a, b):

    return (
        a[0] * b[0]
        + a[1] * b[1]
        + a[2] * b[2]
    )


def cross(a, b):

    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0]
    )


def magnitude(a):

    return math.sqrt(
        a[0] * a[0]
        + a[1] * a[1]
        + a[2] * a[2]
    )


def normalize(a):

    m = magnitude(a)

    if m <= PLANE_EPS:
        raise RuntimeError(
            "Encountered a zero-length centerline direction."
        )

    return (
        a[0] / m,
        a[1] / m,
        a[2] / m
    )


def point_tuple(polydata, point_id):

    p = polydata.GetPoint(point_id)

    return (
        float(p[0]),
        float(p[1]),
        float(p[2])
    )


def distance(a, b):

    return magnitude(
        vec_sub(a, b)
    )


def distance_from_axis(point, center, axis):

    rel = vec_sub(point, center)

    axial = dot(rel, axis)

    radial_vector = vec_sub(
        rel,
        vec_scale(axis, axial)
    )

    return magnitude(radial_vector)


# ============================================================
# CLEAN / TRIANGULATE HELPERS
# ============================================================

def clean_polydata(polydata, tolerance=MERGE_TOLERANCE):

    cleaner = vtk.vtkCleanPolyData()
    cleaner.SetInputData(polydata)
    cleaner.PointMergingOn()
    cleaner.ToleranceIsAbsoluteOn()
    cleaner.SetAbsoluteTolerance(tolerance)
    cleaner.Update()

    out = vtk.vtkPolyData()
    out.DeepCopy(cleaner.GetOutput())

    return out


def triangulate_polydata(polydata):

    tri = vtk.vtkTriangleFilter()
    tri.SetInputData(polydata)
    tri.PassLinesOff()
    tri.PassVertsOff()
    tri.Update()

    out = vtk.vtkPolyData()
    out.DeepCopy(tri.GetOutput())

    return out


# ============================================================
# COORDINATE KEY
#
# The full centerline branches may use separate point IDs at a
# bifurcation even though the endpoint coordinates are identical.
# ============================================================

def coordinate_key(point):

    tol = COORDINATE_TOLERANCE

    return (
        int(round(point[0] / tol)),
        int(round(point[1] / tol)),
        int(round(point[2] / tol))
    )


# ============================================================
# CENTERLINE CELL ENDPOINT HELPERS
# ============================================================

def get_endpoint_entry(centerline, cell_id, endpoint_local_id):

    cell = centerline.GetCell(cell_id)

    if endpoint_local_id == 0:
        point_id = cell.GetPointId(0)
    else:
        point_id = cell.GetPointId(
            cell.GetNumberOfPoints() - 1
        )

    return {
        "cell_id": int(cell_id),
        "endpoint_local_id": int(endpoint_local_id),
        "point_id": int(point_id)
    }


def get_interior_neighbor(centerline, entry):

    cell = centerline.GetCell(
        entry["cell_id"]
    )

    n = cell.GetNumberOfPoints()

    if n < 2:
        raise RuntimeError(
            "Centerline cell {} has fewer than 2 points."
            .format(entry["cell_id"])
        )

    if entry["endpoint_local_id"] == 0:
        local_neighbor = 1
    else:
        local_neighbor = n - 2

    neighbor_id = cell.GetPointId(
        local_neighbor
    )

    return int(neighbor_id)


# ============================================================
# CENTERLINE PATH SAMPLING FOR SIDE-BRANCH RETRIES
#
# A removed centerline cell is a polyline running from one anatomical
# junction to the next. These helpers order that polyline starting at the
# pruning junction and let us place retry planes a physical distance
# downstream along the branch instead of only a fraction of the first edge.
# ============================================================

def ordered_cell_point_ids_from_endpoint(centerline, entry):

    cell = centerline.GetCell(entry["cell_id"])

    ids = [
        int(cell.GetPointId(i))
        for i in range(cell.GetNumberOfPoints())
    ]

    if entry["endpoint_local_id"] == 0:
        return ids

    return list(reversed(ids))


def polyline_total_length(centerline, point_ids):

    total = 0.0

    for i in range(len(point_ids) - 1):
        total += distance(
            point_tuple(centerline, point_ids[i]),
            point_tuple(centerline, point_ids[i + 1])
        )

    return total


def robust_cell_reference_radius(centerline, cell_id, radius_array):

    cell = centerline.GetCell(int(cell_id))

    if cell is None or cell.GetNumberOfPoints() == 0:
        return None

    values = []

    for i in range(cell.GetNumberOfPoints()):
        point_id = int(cell.GetPointId(i))
        r = float(radius_array.GetTuple1(point_id))

        if math.isfinite(r) and r > 0.0:
            values.append(r)

    if len(values) == 0:
        return None

    # Median is deliberately used instead of the endpoint value.  It is
    # stable when Radius collapses near an old cap but does not force the
    # distal opening to be as large as the proximal end of a tapering vessel.
    values.sort()
    mid = len(values) // 2

    if len(values) % 2 == 1:
        return float(values[mid])

    return 0.5 * float(values[mid - 1] + values[mid])


def sample_centerline_cell_at_distance(
    centerline,
    entry,
    distance_from_junction,
    radius_array
):

    point_ids = ordered_cell_point_ids_from_endpoint(
        centerline,
        entry
    )

    if len(point_ids) < 2:
        return None

    target = max(0.0, float(distance_from_junction))
    travelled = 0.0

    for i in range(len(point_ids) - 1):

        id0 = point_ids[i]
        id1 = point_ids[i + 1]

        p0 = point_tuple(centerline, id0)
        p1 = point_tuple(centerline, id1)

        segment_vector = vec_sub(p1, p0)
        segment_length = magnitude(segment_vector)

        if segment_length <= PLANE_EPS:
            continue

        if travelled + segment_length >= target:

            t = (target - travelled) / segment_length

            if t < 0.0:
                t = 0.0
            elif t > 1.0:
                t = 1.0

            center = vec_add(
                p0,
                vec_scale(segment_vector, t)
            )

            normal = normalize(segment_vector)

            r0 = float(radius_array.GetTuple1(id0))
            r1 = float(radius_array.GetTuple1(id1))

            valid0 = math.isfinite(r0) and r0 > 0.0
            valid1 = math.isfinite(r1) and r1 > 0.0

            if valid0 and valid1:
                radius = r0 + t * (r1 - r0)
            elif valid0:
                radius = r0
            elif valid1:
                radius = r1
            else:
                radius = None

            return {
                "center": center,
                "normal": normal,
                "radius": radius,
                "distance_from_junction": target
            }

        travelled += segment_length

    return None


def build_cut_retry_candidates(
    centerline,
    endpoint_entry,
    base_cut,
    radius_array,
    flip_sample_normal=False
):

    point_ids = ordered_cell_point_ids_from_endpoint(
        centerline,
        endpoint_entry
    )

    total_length = polyline_total_length(
        centerline,
        point_ids
    )

    if total_length <= PLANE_EPS:
        return []

    endpoint = point_tuple(
        centerline,
        endpoint_entry["point_id"]
    )

    # Every base cut used here lies on the first centerline edge measured
    # from this endpoint, so Euclidean distance is also its path distance.
    base_distance = distance(
        base_cut["center"],
        endpoint
    )

    max_distance = (
        total_length
        * CUT_RETRY_MAX_PATH_FRACTION
    )

    base_radius = float(base_cut["radius"])

    if not math.isfinite(base_radius) or base_radius <= 0.0:
        return []

    candidates = []
    used_distances = [base_distance]

    if base_cut.get("type") == "terminal":
        retry_offsets = TERMINAL_RETRY_RADIUS_OFFSETS
    else:
        retry_offsets = SIDE_BRANCH_RETRY_RADIUS_OFFSETS

    for radius_offset in retry_offsets:

        target_distance = (
            base_distance
            + float(radius_offset) * base_radius
        )

        # Stay on this branch cell. Zero is allowed because it corresponds
        # exactly to the anatomical endpoint/junction.
        if target_distance < 0.0:
            target_distance = 0.0

        if target_distance > max_distance:
            target_distance = max_distance

        if any(
            abs(target_distance - d) <= PLANE_EPS
            for d in used_distances
        ):
            continue

        sample = sample_centerline_cell_at_distance(
            centerline,
            endpoint_entry,
            target_distance,
            radius_array
        )

        if sample is None:
            continue

        if (
            sample["radius"] is None
            or not math.isfinite(sample["radius"])
            or sample["radius"] <= 0.0
        ):
            continue

        candidate = dict(base_cut)
        candidate["center"] = sample["center"]

        sampled_normal = sample["normal"]

        if flip_sample_normal:
            sampled_normal = vec_scale(
                sampled_normal,
                -1.0
            )

        candidate["normal"] = sampled_normal
        candidate["radius"] = sample["radius"]
        candidate["cut_radius"] = (
            sample["radius"]
            * CUT_RADIUS_FACTOR
        )
        candidate["retry_distance_from_anchor"] = target_distance
        candidate["retry_radius_offset"] = float(radius_offset)

        candidates.append(candidate)
        used_distances.append(target_distance)

    return candidates


def build_outlet_retry_candidates(
    centerline,
    endpoint_entry,
    base_cut,
    radius_array
):

    point_ids = ordered_cell_point_ids_from_endpoint(
        centerline,
        endpoint_entry
    )

    total_length = polyline_total_length(
        centerline,
        point_ids
    )

    if total_length <= PLANE_EPS:
        return []

    endpoint = point_tuple(
        centerline,
        endpoint_entry["point_id"]
    )

    base_distance = distance(
        base_cut["center"],
        endpoint
    )

    max_distance = total_length * CUT_RETRY_MAX_PATH_FRACTION
    candidates = []
    used_distances = [base_distance]

    for path_fraction in OUTLET_RETRY_PATH_FRACTIONS:

        target_distance = min(
            float(path_fraction) * total_length,
            max_distance
        )

        # These retries are specifically intended to move farther toward the
        # parent branch than the originally requested IsOutlet cut.
        if target_distance <= base_distance + PLANE_EPS:
            continue

        if any(
            abs(target_distance - d) <= PLANE_EPS
            for d in used_distances
        ):
            continue

        sample = sample_centerline_cell_at_distance(
            centerline,
            endpoint_entry,
            target_distance,
            radius_array
        )

        if sample is None:
            continue

        if (
            sample["radius"] is None
            or not math.isfinite(sample["radius"])
            or sample["radius"] <= 0.0
        ):
            continue

        candidate = dict(base_cut)
        candidate["center"] = sample["center"]

        # Sampling proceeds from the IsOutlet back toward its parent, so the
        # sampled tangent is upstream.  Flip it to preserve the outward cut
        # normal used by the clipping logic.
        candidate["normal"] = vec_scale(sample["normal"], -1.0)
        candidate["radius"] = sample["radius"]
        candidate["cut_radius"] = sample["radius"] * CUT_RADIUS_FACTOR
        candidate["retry_distance_from_anchor"] = target_distance
        candidate["retry_path_fraction"] = float(path_fraction)

        candidates.append(candidate)
        used_distances.append(target_distance)

    return candidates


def build_side_branch_retry_candidates(
    full_centerline,
    removed_entry,
    base_cut,
    radius_array
):

    # The endpoint is the pruning junction and the ordered cell direction
    # points INTO the removed branch, which is also the desired plane normal.
    return build_cut_retry_candidates(
        full_centerline,
        removed_entry,
        base_cut,
        radius_array,
        flip_sample_normal=False
    )


# ============================================================
# FIND RETAINED FULL-CENTERLINE CELL IDS
# ============================================================

def get_retained_cell_ids(trimmed_centerline, full_centerline):

    arr = trimmed_centerline.GetCellData().GetArray(
        CELL_ID_ARRAY
    )

    if arr is None:
        raise RuntimeError(
            "Trimmed centerline is missing CellData array "
            "'{}'. This script expects the trimmed file to "
            "preserve the original full-tree cell IDs."
            .format(CELL_ID_ARRAY)
        )

    retained = set()

    for i in range(arr.GetNumberOfTuples()):

        cell_id = int(
            round(arr.GetTuple1(i))
        )

        if (
            cell_id < 0
            or cell_id >= full_centerline.GetNumberOfCells()
        ):
            raise RuntimeError(
                "Trimmed CellIds contains {}, but the full "
                "centerline has only {} cells."
                .format(
                    cell_id,
                    full_centerline.GetNumberOfCells()
                )
            )

        retained.add(cell_id)

    return retained


# ============================================================
# BUILD FULL-TREE JUNCTION MAP
# ============================================================

def build_endpoint_map(centerline):

    endpoint_map = defaultdict(list)

    for cell_id in range(
        centerline.GetNumberOfCells()
    ):

        cell = centerline.GetCell(cell_id)

        if cell.GetNumberOfPoints() < 2:
            continue

        for endpoint_local_id in (0, -1):

            entry = get_endpoint_entry(
                centerline,
                cell_id,
                endpoint_local_id
            )

            p = point_tuple(
                centerline,
                entry["point_id"]
            )

            key = coordinate_key(p)

            endpoint_map[key].append(entry)

    return endpoint_map


# ============================================================
# BUILD CUT DEFINITIONS
#
# Same centerline-pruning logic as the previous working script.
#
# CASE 1: exactly one retained incident centerline cell
#         -> retained tree terminates here -> one terminal cut
#
# CASE 2: two or more retained incident cells + removed branch
#         -> retained tree continues -> cut only removed branch
# ============================================================

def build_cut_definitions(full_centerline, retained_cell_ids):

    endpoint_map = build_endpoint_map(
        full_centerline
    )

    radius_array = (
        full_centerline
        .GetPointData()
        .GetArray(RADIUS_ARRAY)
    )

    if radius_array is None:
        raise RuntimeError(
            "Full centerline is missing PointData array "
            "'{}'.".format(RADIUS_ARRAY)
        )

    cuts = []
    boundary_junction_count = 0

    for key, entries in endpoint_map.items():

        kept_entries = [
            e for e in entries
            if e["cell_id"] in retained_cell_ids
        ]

        removed_entries = [
            e for e in entries
            if e["cell_id"] not in retained_cell_ids
        ]

        if (
            len(kept_entries) == 0
            or len(removed_entries) == 0
        ):
            continue

        boundary_junction_count += 1

        junction = point_tuple(
            full_centerline,
            entries[0]["point_id"]
        )

        # ----------------------------------------------------
        # CASE 1: retained tree terminates here.
        # ----------------------------------------------------
        if len(kept_entries) == 1:

            kept = kept_entries[0]

            interior_id = get_interior_neighbor(
                full_centerline,
                kept
            )

            interior = point_tuple(
                full_centerline,
                interior_id
            )

            # Normal points from retained parent toward the
            # downstream terminal junction.
            normal = normalize(
                vec_sub(junction, interior)
            )

            cut_center = vec_add(
                junction,
                vec_scale(
                    vec_sub(interior, junction),
                    TERMINAL_OFFSET_FRACTION
                )
            )

            r0 = float(
                radius_array.GetTuple1(
                    kept["point_id"]
                )
            )

            r1 = float(
                radius_array.GetTuple1(
                    interior_id
                )
            )

            valid_radii = [
                r for r in (r0, r1)
                if math.isfinite(r) and r > 0.0
            ]

            if len(valid_radii) == 0:
                raise RuntimeError(
                    "Invalid radius at terminal pruning "
                    "junction associated with full cell {}."
                    .format(kept["cell_id"])
                )

            radius = max(valid_radii)

            reference_radius = robust_cell_reference_radius(
                full_centerline,
                kept["cell_id"],
                radius_array
            )

            if reference_radius is None:
                reference_radius = radius

            base_cut = {
                "type": "terminal",
                "center": cut_center,
                "junction": junction,
                "normal": normal,
                "radius": radius,
                "reference_radius": reference_radius,
                "cut_radius": radius * CUT_RADIUS_FACTOR,
                "source_cell_id": kept["cell_id"],
                "source_endpoint_local_id": kept["endpoint_local_id"],
                "removed_cell_ids": [
                    e["cell_id"]
                    for e in removed_entries
                ]
            }

            # If this terminal plane happens to overlap a previous cut or
            # otherwise fails geometrically, try nearby positions on the
            # retained parent cell. sample_centerline_cell_at_distance()
            # points from the terminal junction upstream, so flip its normal
            # to preserve the desired outward/downstream cut orientation.
            base_cut["retry_candidates"] = build_cut_retry_candidates(
                full_centerline,
                kept,
                base_cut,
                radius_array,
                flip_sample_normal=True
            )

            cuts.append(base_cut)

        # ----------------------------------------------------
        # CASE 2: retained tree continues; remove side branch.
        # ----------------------------------------------------
        else:

            for removed in removed_entries:

                interior_id = get_interior_neighbor(
                    full_centerline,
                    removed
                )

                interior = point_tuple(
                    full_centerline,
                    interior_id
                )

                branch_vector = vec_sub(
                    interior,
                    junction
                )

                normal = normalize(
                    branch_vector
                )

                cut_center = vec_add(
                    junction,
                    vec_scale(
                        branch_vector,
                        SIDE_BRANCH_OFFSET_FRACTION
                    )
                )

                r0 = float(
                    radius_array.GetTuple1(
                        removed["point_id"]
                    )
                )

                r1 = float(
                    radius_array.GetTuple1(
                        interior_id
                    )
                )

                valid_radii = [
                    r for r in (r0, r1)
                    if math.isfinite(r) and r > 0.0
                ]

                if len(valid_radii) == 0:
                    raise RuntimeError(
                        "Invalid radius for removed full "
                        "centerline cell {}."
                        .format(removed["cell_id"])
                    )

                radius = min(valid_radii)

                reference_radius = robust_cell_reference_radius(
                    full_centerline,
                    removed["cell_id"],
                    radius_array
                )

                if reference_radius is None:
                    reference_radius = radius

                base_cut = {
                    "type": "side_branch",
                    "center": cut_center,
                    "junction": junction,
                    "normal": normal,
                    "radius": radius,
                    "reference_radius": reference_radius,
                    "cut_radius": radius * CUT_RADIUS_FACTOR,
                    "source_cell_id": removed["cell_id"],
                    "source_endpoint_local_id": removed["endpoint_local_id"],
                    "removed_cell_ids": [
                        removed["cell_id"]
                    ]
                }

                # Keep the original plane as attempt 1, but pre-compute
                # progressively deeper locations on the removed branch in
                # case the first plane intersects the bifurcation complex.
                base_cut["retry_candidates"] = (
                    build_side_branch_retry_candidates(
                        full_centerline,
                        removed,
                        base_cut,
                        radius_array
                    )
                )

                cuts.append(base_cut)

    return cuts, boundary_junction_count


# ============================================================
# BUILD CUT DEFINITIONS FOR TRUE IsOutlet POINTS
#
# These are separate from the synthetic pruning-terminal cuts above.
# Every IsOutlet == 1 point in the TRIMMED centerline is treated as a
# desired final surface outlet. The plane normal points from the immediately
# upstream centerline point toward the IsOutlet point, and the plane can be
# shifted upstream independently with OUTLET_OFFSET_FRACTION.
# ============================================================

def build_outlet_cut_definitions(trimmed_centerline):

    pd = trimmed_centerline.GetPointData()

    is_outlet = pd.GetArray("IsOutlet")

    if is_outlet is None:
        raise RuntimeError(
            "Trimmed centerline is missing PointData array 'IsOutlet'."
        )

    radius_array = pd.GetArray(RADIUS_ARRAY)

    if radius_array is None:
        raise RuntimeError(
            "Trimmed centerline is missing PointData array "
            "'{}'.".format(RADIUS_ARRAY)
        )

    cuts = []
    skipped_topology = 0

    for outlet_id in range(trimmed_centerline.GetNumberOfPoints()):

        if int(round(is_outlet.GetTuple1(outlet_id))) != 1:
            continue

        # Find a polyline cell for which this IsOutlet point is a true
        # endpoint. That lets us identify the immediately upstream point
        # without making assumptions about cell orientation.
        cell_ids = vtk.vtkIdList()
        trimmed_centerline.GetPointCells(outlet_id, cell_ids)

        endpoint_candidates = []

        for j in range(cell_ids.GetNumberOfIds()):

            cell_id = int(cell_ids.GetId(j))
            cell = trimmed_centerline.GetCell(cell_id)

            n = cell.GetNumberOfPoints()

            if n < 2:
                continue

            first_id = int(cell.GetPointId(0))
            last_id = int(cell.GetPointId(n - 1))

            if outlet_id == first_id:
                interior_id = int(cell.GetPointId(1))
                endpoint_candidates.append((cell_id, interior_id, 0))

            elif outlet_id == last_id:
                interior_id = int(cell.GetPointId(n - 2))
                endpoint_candidates.append((cell_id, interior_id, -1))

        if len(endpoint_candidates) == 0:
            print(
                "WARNING: IsOutlet point {} is not an endpoint of any "
                "polyline cell. Skipping it.".format(outlet_id)
            )
            skipped_topology += 1
            continue

        if len(endpoint_candidates) > 1:
            print(
                "WARNING: IsOutlet point {} belongs to {} endpoint cells. "
                "Using the first one.".format(
                    outlet_id,
                    len(endpoint_candidates)
                )
            )

        source_cell_id, interior_id, source_endpoint_local_id = endpoint_candidates[0]

        outlet = point_tuple(
            trimmed_centerline,
            outlet_id
        )

        interior = point_tuple(
            trimmed_centerline,
            interior_id
        )

        # Outward/downstream direction of the final centerline segment.
        normal = normalize(
            vec_sub(outlet, interior)
        )

        # Positive OUTLET_OFFSET_FRACTION moves the plane upstream, away
        # from the old terminal tip/cap and into the vessel lumen.
        cut_center = vec_add(
            outlet,
            vec_scale(
                vec_sub(interior, outlet),
                OUTLET_OFFSET_FRACTION
            )
        )

        r0 = float(
            radius_array.GetTuple1(outlet_id)
        )

        r1 = float(
            radius_array.GetTuple1(interior_id)
        )

        valid_radii = [
            r for r in (r0, r1)
            if math.isfinite(r) and r > 0.0
        ]

        if len(valid_radii) == 0:
            print(
                "WARNING: IsOutlet point {} has no valid local Radius. "
                "Skipping it.".format(outlet_id)
            )
            skipped_topology += 1
            continue

        # The radius is only used for QC/locality checks. Using the larger
        # of the endpoint and immediately-upstream values is deliberately
        # conservative, especially where Radius collapses near an old cap.
        radius = max(valid_radii)

        reference_radius = robust_cell_reference_radius(
            trimmed_centerline,
            source_cell_id,
            radius_array
        )

        if reference_radius is None:
            reference_radius = radius

        base_cut = {
            "type": "outlet",
            "center": cut_center,
            "junction": outlet,
            "normal": normal,
            "radius": radius,
            "reference_radius": reference_radius,
            "cut_radius": radius * CUT_RADIUS_FACTOR,
            "source_cell_id": source_cell_id,
            "source_endpoint_local_id": source_endpoint_local_id,
            "source_point_id": outlet_id,
            "removed_cell_ids": []
        }

        endpoint_entry = {
            "cell_id": source_cell_id,
            "endpoint_local_id": source_endpoint_local_id,
            "point_id": outlet_id
        }

        # Alternate positions move both closer to the true endpoint and
        # farther upstream along the same final centerline cell. The sampled
        # cell direction points upstream, so flip it to keep the cut normal
        # outward/downstream.
        base_cut["retry_candidates"] = build_outlet_retry_candidates(
            trimmed_centerline,
            endpoint_entry,
            base_cut,
            radius_array
        )

        cuts.append(base_cut)

    return cuts, skipped_topology


# ============================================================
# INLET CUT: HALF THE CENTERLINE ARC LENGTH TO FIRST BIFURCATION
# ============================================================

def trace_inlet_to_first_bifurcation(centerline):
    """Follow connected centerline cells from IsInlet to the first branch.

    Cell endpoints are grouped by coordinate, because VMTK centerlines can
    store geometrically identical junction coordinates with distinct point
    IDs. A bifurcation has at least THREE distinct incident centerline cells.
    """
    inlet_array = centerline.GetPointData().GetArray("IsInlet")
    if inlet_array is None:
        raise RuntimeError("Trimmed centerline is missing PointData 'IsInlet'.")

    inlet_ids = [
        i for i in range(centerline.GetNumberOfPoints())
        if int(round(inlet_array.GetTuple1(i))) == 1
    ]
    if len(inlet_ids) != 1:
        raise RuntimeError(
            "Expected exactly one IsInlet point, but found {}."
            .format(len(inlet_ids))
        )

    inlet_id = int(inlet_ids[0])
    inlet_xyz = point_tuple(centerline, inlet_id)
    endpoints = build_endpoint_map(centerline)
    entries = [
        e for e in endpoints.get(coordinate_key(inlet_xyz), [])
        if e["point_id"] == inlet_id
    ]
    if len(entries) != 1:
        raise RuntimeError(
            "IsInlet must be at one unique centerline cell endpoint; "
            "found {} matching endpoint cells."
            .format(len(entries))
        )

    entry = entries[0]
    visited = set()
    segments = []
    first_branch_xyz = None
    first_branch_cell_ids = []

    while True:
        cell_id = entry["cell_id"]
        if cell_id in visited:
            raise RuntimeError("Loop in the inlet-to-first-bifurcation centerline graph.")
        visited.add(cell_id)

        point_ids = ordered_cell_point_ids_from_endpoint(centerline, entry)
        if len(point_ids) < 2:
            raise RuntimeError("Inlet trunk contains a degenerate centerline cell.")

        for j in range(len(point_ids) - 1):
            id0 = point_ids[j]
            id1 = point_ids[j + 1]
            p0 = point_tuple(centerline, id0)
            p1 = point_tuple(centerline, id1)
            if distance(p0, p1) > PLANE_EPS:
                segments.append((cell_id, id0, id1))

        distal_xyz = point_tuple(centerline, point_ids[-1])
        next_entries = endpoints.get(coordinate_key(distal_xyz), [])
        incident_cells = sorted(set(e["cell_id"] for e in next_entries))

        if len(incident_cells) >= 3:
            first_branch_xyz = distal_xyz
            first_branch_cell_ids = incident_cells
            break

        remaining = [
            e for e in next_entries
            if e["cell_id"] != cell_id and e["cell_id"] not in visited
        ]
        if len(remaining) != 1:
            raise RuntimeError(
                "Could not follow the inlet trunk to its first bifurcation. "
                "At {} the junction has {} incident cells and {} possible "
                "continuations."
                .format(distal_xyz, len(incident_cells), len(remaining))
            )
        entry = remaining[0]

    if not segments:
        raise RuntimeError("No nonzero-length segments on the inlet trunk.")

    total_length = sum(
        distance(point_tuple(centerline, p0), point_tuple(centerline, p1))
        for _, p0, p1 in segments
    )
    if total_length <= PLANE_EPS:
        raise RuntimeError("Zero inlet-to-first-bifurcation path length.")

    return {
        "inlet_point_id": inlet_id,
        "inlet_xyz": inlet_xyz,
        "branch_xyz": first_branch_xyz,
        "branch_cell_ids": first_branch_cell_ids,
        "segments": segments,
        "total_length": total_length,
        "visited_cell_ids": sorted(visited),
    }


def sample_inlet_trunk(centerline, trace, fraction):
    """Interpolate position/radius/tangent using true cumulative arc length."""
    radius_array = centerline.GetPointData().GetArray(RADIUS_ARRAY)
    if radius_array is None:
        raise RuntimeError("Trimmed centerline is missing PointData '{}'".format(RADIUS_ARRAY))

    fraction = float(fraction)
    if not (0.0 < fraction < 1.0):
        raise ValueError("Inlet trunk sampling fraction must be between 0 and 1.")

    target = trace["total_length"] * fraction
    travelled = 0.0
    segments = trace["segments"]

    for segment_index, (cell_id, id0, id1) in enumerate(segments):
        p0 = point_tuple(centerline, id0)
        p1 = point_tuple(centerline, id1)
        length = distance(p0, p1)
        if length <= PLANE_EPS:
            continue

        if travelled + length >= target or segment_index == len(segments) - 1:
            t = min(1.0, max(0.0, (target - travelled) / length))
            forward = normalize(vec_sub(p1, p0))
            r0 = float(radius_array.GetTuple1(id0))
            r1 = float(radius_array.GetTuple1(id1))
            valid0 = math.isfinite(r0) and r0 > 0.0
            valid1 = math.isfinite(r1) and r1 > 0.0
            if valid0 and valid1:
                radius = (1.0-t)*r0 + t*r1
            elif valid0:
                radius = r0
            elif valid1:
                radius = r1
            else:
                raise RuntimeError(
                    "Inlet trunk has no valid Radius at path fraction {}."
                    .format(fraction)
                )
            return {
                "center": vec_add(p0, vec_scale(vec_sub(p1, p0), t)),
                "forward": forward,
                "radius": radius,
                "source_cell_id": int(cell_id),
                "path_fraction": fraction,
            }
        travelled += length

    raise RuntimeError("Could not interpolate the inlet centerline trunk.")


def build_inlet_cut_definition(trimmed_centerline):
    """Make a normal pointing from the first bifurcation TOWARD old inlet.

    The exact cut implementation retains the NEGATIVE side of its plane;
    therefore this reversed normal keeps the bifurcation/pulmonary tree and
    discards the unwanted original inlet stub.
    """
    trace = trace_inlet_to_first_bifurcation(trimmed_centerline)
    base = sample_inlet_trunk(
        trimmed_centerline, trace, INLET_CUT_PATH_FRACTION
    )
    retained = sample_inlet_trunk(
        trimmed_centerline, trace, INLET_RETAINED_ANCHOR_PATH_FRACTION
    )
    cut = {
        "type": "inlet",
        "center": base["center"],
        "junction": trace["inlet_xyz"],
        "normal": vec_scale(base["forward"], -1.0),
        "radius": base["radius"],
        "reference_radius": base["radius"],
        "cut_radius": base["radius"] * CUT_RADIUS_FACTOR,
        "source_cell_id": base["source_cell_id"],
        "source_endpoint_local_id": 0,
        "source_point_id": trace["inlet_point_id"],
        "removed_cell_ids": [],
        "retained_anchor": retained["center"],
        "first_bifurcation": trace["branch_xyz"],
        "inlet_path_length": trace["total_length"],
        "retry_path_fraction": float(INLET_CUT_PATH_FRACTION),
    }

    retries = []
    for fraction in INLET_RETRY_PATH_FRACTIONS:
        if abs(float(fraction) - INLET_CUT_PATH_FRACTION) < 1.0e-12:
            continue
        sample = sample_inlet_trunk(trimmed_centerline, trace, fraction)
        alternate = dict(cut)
        alternate.update({
            "center": sample["center"],
            "normal": vec_scale(sample["forward"], -1.0),
            "radius": sample["radius"],
            "reference_radius": sample["radius"],
            "cut_radius": sample["radius"] * CUT_RADIUS_FACTOR,
            "source_cell_id": sample["source_cell_id"],
            "retry_path_fraction": float(fraction),
            "retry_distance_from_anchor": fraction * trace["total_length"],
        })
        retries.append(alternate)
    cut["retry_candidates"] = retries

    return cut, trace


# ============================================================
# SIDEDNESS / AUTOMATIC CAP NAMES
# ============================================================

def get_centerline_sidedness_array(centerline):
    sidedness = centerline.GetPointData().GetAbstractArray(SIDEDNESS_ARRAY)
    if sidedness is None:
        if REQUIRE_SIDEDNESS_ARRAY:
            raise RuntimeError(
                "Trimmed centerline missing PointData '{}' array."
                .format(SIDEDNESS_ARRAY)
            )
        print(
            "WARNING: trimmed centerline has no '{}' PointData. "
            "Non-inlet caps will be named Other#."
            .format(SIDEDNESS_ARRAY)
        )
    return sidedness


def sidedness_at_point(sidedness_array, point_id):
    if sidedness_array is None:
        return "Other"

    # vtkStringArray, vtkVariantArray and string-compatible VTK arrays.
    try:
        value = sidedness_array.GetVariantValue(int(point_id)).ToString()
    except (AttributeError, TypeError, ValueError):
        try:
            value = str(sidedness_array.GetValue(int(point_id)))
        except (AttributeError, TypeError, ValueError):
            return "Other"

    text = str(value).strip().casefold()
    if text == "left":
        return "Left"
    if text == "right":
        return "Right"
    return "Other"


def build_sidedness_point_locator(trimmed_centerline):
    locator = vtk.vtkStaticPointLocator()
    locator.SetDataSet(trimmed_centerline)
    locator.BuildLocator()
    return locator


def classify_loop_sidedness(
    centroid, trimmed_centerline, sidedness_array, locator,
    accepted_cut=None
):
    """Prefer the correct source IsOutlet point when the loop matches a cut.

    For branch-pruning cuts, synthetic terminal cuts, and unmatched boundary
    loops, use the closest trimmed-centerline point in 3D instead.
    """
    if sidedness_array is None:
        return "Other"

    if accepted_cut is not None and accepted_cut.get("type") == "outlet":
        point_id = accepted_cut.get("source_point_id", -1)
        if 0 <= point_id < trimmed_centerline.GetNumberOfPoints():
            return sidedness_at_point(sidedness_array, point_id)

    closest_point_id = int(locator.FindClosestPoint(centroid))
    if closest_point_id < 0:
        return "Other"
    return sidedness_at_point(sidedness_array, closest_point_id)


# ============================================================
# DEBUG CIRCLES FOR THE INTENDED CUT PLANES
# ============================================================

def add_debug_circle(center, normal, radius, points, lines):

    if abs(normal[0]) < 0.9:
        reference = (1.0, 0.0, 0.0)
    else:
        reference = (0.0, 1.0, 0.0)

    u = normalize(
        cross(normal, reference)
    )

    v = normalize(
        cross(normal, u)
    )

    ids = []

    for i in range(DEBUG_CIRCLE_RESOLUTION):

        theta = (
            2.0
            * math.pi
            * float(i)
            / float(DEBUG_CIRCLE_RESOLUTION)
        )

        p = vec_add(
            center,
            vec_add(
                vec_scale(
                    u,
                    radius * math.cos(theta)
                ),
                vec_scale(
                    v,
                    radius * math.sin(theta)
                )
            )
        )

        ids.append(
            points.InsertNextPoint(p)
        )

    polyline = vtk.vtkPolyLine()

    polyline.GetPointIds().SetNumberOfIds(
        DEBUG_CIRCLE_RESOLUTION + 1
    )

    for i in range(DEBUG_CIRCLE_RESOLUTION):
        polyline.GetPointIds().SetId(i, ids[i])

    polyline.GetPointIds().SetId(
        DEBUG_CIRCLE_RESOLUTION,
        ids[0]
    )

    lines.InsertNextCell(polyline)


def write_cut_debug_file(cuts, filename):

    points = vtk.vtkPoints()
    lines = vtk.vtkCellArray()

    cut_type = vtk.vtkIntArray()
    cut_type.SetName("CutType")

    source_cell = vtk.vtkIntArray()
    source_cell.SetName("SourceCellId")

    cut_radius = vtk.vtkDoubleArray()
    cut_radius.SetName("ExpectedCutRadius")

    attempt_number = vtk.vtkIntArray()
    attempt_number.SetName("AttemptNumber")

    planned_cut_index = vtk.vtkIntArray()
    planned_cut_index.SetName("PlannedCutIndex")

    for cut in cuts:

        add_debug_circle(
            cut["center"],
            cut["normal"],
            cut["cut_radius"],
            points,
            lines
        )

        if cut["type"] == "side_branch":
            cut_type.InsertNextValue(1)
        elif cut["type"] == "terminal":
            cut_type.InsertNextValue(2)
        elif cut["type"] == "inlet":
            cut_type.InsertNextValue(5)
        else:
            # 3 = true IsOutlet cut from the trimmed centerline
            cut_type.InsertNextValue(3)

        source_cell.InsertNextValue(
            int(cut["source_cell_id"])
        )

        cut_radius.InsertNextValue(
            float(cut["cut_radius"])
        )

        attempt_number.InsertNextValue(
            int(cut.get("attempt_number", 1))
        )

        planned_cut_index.InsertNextValue(
            int(cut.get("planned_cut_index", -1))
        )

    debug = vtk.vtkPolyData()
    debug.SetPoints(points)
    debug.SetLines(lines)

    debug.GetCellData().AddArray(cut_type)
    debug.GetCellData().AddArray(source_cell)
    debug.GetCellData().AddArray(cut_radius)
    debug.GetCellData().AddArray(attempt_number)
    debug.GetCellData().AddArray(planned_cut_index)

    write_polydata(debug, filename)


# ============================================================
# FIND THE ACTUAL TARGET PLANE/SURFACE INTERSECTION LOOP
#
# The plane may intersect several disconnected PA branches. vtkCutter
# creates ALL plane/surface intersections, then connectivity keeps ONLY
# the contour nearest the centerline cut point.
# ============================================================

class NoCutLoopFound(RuntimeError):
    pass


class SuspiciousCutLoop(RuntimeError):
    pass


def median_value(values):

    vals = sorted(float(v) for v in values)

    if len(vals) == 0:
        return None

    mid = len(vals) // 2

    if len(vals) % 2 == 1:
        return vals[mid]

    return 0.5 * (vals[mid - 1] + vals[mid])


def orientation_2d(a, b, c):

    return (
        (b[0] - a[0]) * (c[1] - a[1])
        - (b[1] - a[1]) * (c[0] - a[0])
    )


def point_distance_2d(a, b):

    dx = a[0] - b[0]
    dy = a[1] - b[1]
    return math.sqrt(dx * dx + dy * dy)


def cross_2d(a, b):
    return a[0] * b[1] - a[1] * b[0]


def point_to_segment_distance_2d(p, a, b):
    """Return distance, projection parameter, and closest point in XY."""
    vx, vy = b[0] - a[0], b[1] - a[1]
    denom = vx * vx + vy * vy
    if denom <= 1.0e-30:
        return point_distance_2d(p, a), 0.0, a
    t = ((p[0] - a[0])*vx + (p[1] - a[1])*vy) / denom
    tc = min(1.0, max(0.0, t))
    nearest = (a[0] + tc*vx, a[1] + tc*vy)
    return point_distance_2d(p, nearest), t, nearest


def segment_intersection_detail_2d(a, b, c, d, tol):
    """Detect actual crossing, on-segment touch, or collinear overlap.

    `tol` is a LENGTH. Unlike the previous test, a 2-D orientation
    determinant (L^2) is never directly compared to a length tolerance.
    Near-but-disjoint edges are not counted as physical intersections.
    """
    ab = (b[0] - a[0], b[1] - a[1])
    cd = (d[0] - c[0], d[1] - c[1])
    lab = math.hypot(*ab)
    lcd = math.hypot(*cd)
    if lab <= tol or lcd <= tol:
        return None

    if (max(a[0], b[0]) + tol < min(c[0], d[0])
        or max(c[0], d[0]) + tol < min(a[0], b[0])
        or max(a[1], b[1]) + tol < min(c[1], d[1])
        or max(c[1], d[1]) + tol < min(a[1], b[1])):
        return None

    # Signed perpendicular distances (L) of endpoints to each line.
    d_c = orientation_2d(a, b, c) / lab
    d_d = orientation_2d(a, b, d) / lab
    d_a = orientation_2d(c, d, a) / lcd
    d_b = orientation_2d(c, d, b) / lcd

    opposite_ab = (d_c > tol and d_d < -tol) or (d_c < -tol and d_d > tol)
    opposite_cd = (d_a > tol and d_b < -tol) or (d_a < -tol and d_b > tol)
    cp = (c[0] - a[0], c[1] - a[1])
    den = cross_2d(ab, cd)

    if opposite_ab and opposite_cd and den != 0.0:
        t = cross_2d(cp, cd) / den
        u = cross_2d(cp, ab) / den
        if 0.0 < t < 1.0 and 0.0 < u < 1.0:
            return {
                "kind": "crossing",
                "xy": (a[0] + t*ab[0], a[1] + t*ab[1])
            }

    # Exact or tolerance-level collinearity: check overlap of 1-D
    # projected intervals. Just being nearly parallel is not an overlap.
    if all(abs(s) <= tol for s in (d_c, d_d, d_a, d_b)):
        t0 = (cp[0]*ab[0] + cp[1]*ab[1]) / (lab*lab)
        t1 = t0 + (cd[0]*ab[0] + cd[1]*ab[1]) / (lab*lab)
        lo = max(0.0, min(t0, t1))
        hi = min(1.0, max(t0, t1))
        if (hi - lo)*lab > tol:
            tm = 0.5*(lo + hi)
            return {
                "kind": "overlap",
                "xy": (a[0]+tm*ab[0], a[1]+tm*ab[1])
            }

    # Non-adjacent endpoint lying on another segment. Explicitly require
    # proximity to the segment itself, not just an orientation sign.
    for p, p0, p1 in ((a, c, d), (b, c, d), (c, a, b), (d, a, b)):
        distance_to_seg, t, near = point_to_segment_distance_2d(p, p0, p1)
        if distance_to_seg <= tol and -tol/max(lab, lcd) <= t <= 1.0+tol/max(lab, lcd):
            return {"kind": "touch", "xy": near}

    return None


def segments_intersect_2d(a, b, c, d, tol):
    """Backward-compatible boolean used by earlier QC helpers."""
    return segment_intersection_detail_2d(a, b, c, d, tol) is not None


def inspect_nonadjacent_loop_intersections(xy, tol):
    """Return detailed confirmed intersections for debug output."""
    n = len(xy)
    details = []
    for i in range(n):
        jmin = i + 2
        for j in range(jmin, n):
            # Edges i and j must not share a topological loop vertex.
            if i == 0 and j == n - 1:
                continue
            detail = segment_intersection_detail_2d(
                xy[i], xy[(i + 1) % n],
                xy[j], xy[(j + 1) % n], tol
            )
            if detail is not None:
                details.append({"edge_i": i, "edge_j": j, **detail})
    return details


def count_nonadjacent_loop_intersections(xy, tol):
    return len(inspect_nonadjacent_loop_intersections(xy, tol))


def find_nonadjacent_near_duplicate_vertices(xy, tol):
    """Flag nearly coincident points only when far apart along the loop.

    In a finely sampled but valid contour, vertices only 2 or 3 edges apart
    can lie very close together. Such neighbors are NOT a folded/self-touching
    loop. Compare Euclidean distance against shortest along-loop arclength.
    """
    n = len(xy)
    if n < 4:
        return []
    cumulative = [0.0]
    for i in range(n):
        cumulative.append(cumulative[-1] +
                          point_distance_2d(xy[i], xy[(i+1) % n]))
    perimeter = cumulative[-1]
    hits = []
    for i in range(n):
        for j in range(i+2, n):
            if i == 0 and j == n-1:
                continue
            dist = point_distance_2d(xy[i], xy[j])
            if dist > tol:
                continue
            along = cumulative[j] - cumulative[i]
            shortest = min(along, perimeter - along)
            if shortest <= 5.0*tol:
                continue
            hits.append((i, j, ((xy[i][0]+xy[j][0])*0.5,
                                (xy[i][1]+xy[j][1])*0.5)))
    return hits


def count_nonadjacent_near_duplicate_vertices(xy, tol):
    return len(find_nonadjacent_near_duplicate_vertices(xy, tol))


def compute_loop_shape_metrics(loop, cut):

    # Join cutter line segments into an ordered closed polyline. A clean
    # vessel cross-section should normally collapse to one closed path.
    loop_clean = clean_polydata(
        loop,
        MERGE_TOLERANCE
    )

    # Total length of the raw cutter segments before vtkStripper joins them.
    # For a clean contour this should be essentially identical to the final
    # perimeter. Duplicate/overlapping triangles often create coincident
    # cutter segments and inflate this value.
    raw_segment_length = 0.0

    for raw_cell_id in range(loop_clean.GetNumberOfCells()):
        raw_cell = loop_clean.GetCell(raw_cell_id)

        if raw_cell is None or raw_cell.GetNumberOfPoints() < 2:
            continue

        for raw_i in range(raw_cell.GetNumberOfPoints() - 1):
            p0 = point_tuple(loop_clean, raw_cell.GetPointId(raw_i))
            p1 = point_tuple(loop_clean, raw_cell.GetPointId(raw_i + 1))
            raw_segment_length += distance(p0, p1)

    stripper = vtk.vtkStripper()
    stripper.SetInputData(loop_clean)
    stripper.JoinContiguousSegmentsOn()
    stripper.Update()

    stripped = vtk.vtkPolyData()
    stripped.DeepCopy(stripper.GetOutput())

    if stripped.GetNumberOfCells() != 1:
        raise SuspiciousCutLoop(
            "Selected contour did not reduce to one simple closed path "
            "after stripping ({} paths found)."
            .format(stripped.GetNumberOfCells())
        )

    cell = stripped.GetCell(0)

    if cell is None or cell.GetNumberOfPoints() < 3:
        raise SuspiciousCutLoop(
            "Selected contour has too few ordered points for shape QC."
        )

    ordered = []

    for i in range(cell.GetNumberOfPoints()):
        ordered.append(
            point_tuple(
                stripped,
                cell.GetPointId(i)
            )
        )

    # vtkStripper may repeat the first point at the end of a closed loop.
    if (
        len(ordered) >= 2
        and distance(ordered[0], ordered[-1])
            <= max(MERGE_TOLERANCE * 10.0, 1e-10)
    ):
        ordered = ordered[:-1]

    if len(ordered) < 3:
        raise SuspiciousCutLoop(
            "Selected contour collapsed to fewer than 3 unique points."
        )

    # Build an orthonormal basis spanning the cutting plane.
    normal = normalize(cut["normal"])

    if abs(normal[0]) < 0.9:
        reference = (1.0, 0.0, 0.0)
    else:
        reference = (0.0, 1.0, 0.0)

    u = normalize(cross(normal, reference))
    v = normalize(cross(normal, u))

    xy = []

    for p in ordered:
        rel = vec_sub(p, cut["center"])
        xy.append((dot(rel, u), dot(rel, v)))

    # Shoelace area and polyline perimeter.
    twice_area = 0.0
    perimeter = 0.0

    for i in range(len(ordered)):
        j = (i + 1) % len(ordered)

        x0, y0 = xy[i]
        x1, y1 = xy[j]

        twice_area += x0 * y1 - x1 * y0
        perimeter += distance(ordered[i], ordered[j])

    area = 0.5 * abs(twice_area)

    expected_area = (
        math.pi * float(cut["radius"]) * float(cut["radius"])
        if cut.get("radius") is not None and cut["radius"] > 0.0
        else 0.0
    )

    if expected_area > 0.0:
        area_ratio = area / expected_area
    else:
        area_ratio = float("inf")

    if perimeter > PLANE_EPS and area > 0.0:
        circularity = (
            4.0 * math.pi * area
            / (perimeter * perimeter)
        )
    else:
        circularity = 0.0

    # Assess one-sided protrusion about the loop's OWN center, not the
    # centerline point. This prevents simple centerline offset from being
    # mistaken for an asymmetric shape.
    loop_centroid = (
        sum(p[0] for p in ordered) / float(len(ordered)),
        sum(p[1] for p in ordered) / float(len(ordered)),
        sum(p[2] for p in ordered) / float(len(ordered))
    )

    shape_radii = []

    for p in ordered:
        rel = vec_sub(p, loop_centroid)
        axial = dot(rel, normal)
        radial_vec = vec_sub(
            rel,
            vec_scale(normal, axial)
        )
        shape_radii.append(magnitude(radial_vec))

    median_radius = median_value(shape_radii)
    max_shape_radius = max(shape_radii) if shape_radii else 0.0

    if median_radius is not None and median_radius > PLANE_EPS:
        radial_asymmetry = max_shape_radius / median_radius
    else:
        radial_asymmetry = float("inf")

    centroid_offset = distance(
        loop_centroid,
        cut["center"]
    )

    if cut.get("radius") is not None and cut["radius"] > 0.0:
        centroid_offset_ratio = centroid_offset / float(cut["radius"])
    else:
        centroid_offset_ratio = float("inf")

    effective_radius = (
        math.sqrt(area / math.pi)
        if area > 0.0
        else 0.0
    )

    reference_radius = cut.get(
        "reference_radius",
        cut.get("radius", None)
    )

    if (
        reference_radius is not None
        and math.isfinite(float(reference_radius))
        and float(reference_radius) > 0.0
    ):
        reference_radius = float(reference_radius)
        effective_radius_ratio = effective_radius / reference_radius
    else:
        reference_radius = None
        effective_radius_ratio = float("inf")

    if perimeter > PLANE_EPS:
        segment_redundancy_ratio = raw_segment_length / perimeter
    else:
        segment_redundancy_ratio = float("inf")

    overlap_scale = (
        reference_radius
        if reference_radius is not None
        else max(float(cut.get("radius", 0.0)), 1.0)
    )

    overlap_tol = max(
        overlap_scale * OVERLAP_VERTEX_TOLERANCE_FACTOR,
        MERGE_TOLERANCE * 100.0,
        1e-9
    )

    # A strict spatial tolerance is necessary for *actual* intersections.
    # The more generous overlap_tol remains dedicated to identifying nearly
    # coincident but not literally intersecting non-adjacent vertices.
    intersection_tol = max(
        overlap_scale * LOOP_INTERSECTION_TOLERANCE_RADIUS_FACTOR,
        LOOP_INTERSECTION_TOLERANCE_ABSOLUTE,
        MERGE_TOLERANCE * 0.1,
    )
    intersection_details = inspect_nonadjacent_loop_intersections(
        xy, intersection_tol
    )
    self_intersections = len(intersection_details)
    intersection_points_3d = []
    for hit in intersection_details:
        x, y = hit["xy"]
        p = vec_add(cut["center"],
                    vec_add(vec_scale(u, x), vec_scale(v, y)))
        intersection_points_3d.append((p, hit["kind"]))

    near_duplicate_details = find_nonadjacent_near_duplicate_vertices(
        xy, overlap_tol
    )
    near_duplicate_vertices = len(near_duplicate_details)
    near_duplicate_points_3d = []
    for _i, _j, (x, y) in near_duplicate_details:
        near_duplicate_points_3d.append(
            vec_add(cut["center"],
                    vec_add(vec_scale(u, x), vec_scale(v, y)))
        )

    return {
        "loop_area": area,
        "expected_area": expected_area,
        "area_ratio": area_ratio,
        "loop_perimeter": perimeter,
        "circularity": circularity,
        "radial_asymmetry": radial_asymmetry,
        "shape_centroid": loop_centroid,
        "shape_centroid_offset": centroid_offset,
        "shape_centroid_offset_ratio": centroid_offset_ratio,
        "effective_radius": effective_radius,
        "reference_radius": reference_radius,
        "effective_radius_ratio": effective_radius_ratio,
        "raw_segment_length": raw_segment_length,
        "segment_redundancy_ratio": segment_redundancy_ratio,
        "self_intersections": self_intersections,
        "intersection_points_3d": intersection_points_3d,
        "near_duplicate_vertices": near_duplicate_vertices,
        "near_duplicate_points_3d": near_duplicate_points_3d,
        "intersection_spatial_tolerance": intersection_tol,
    }


def find_target_cut_loop(surface, cut):

    plane = vtk.vtkPlane()
    plane.SetOrigin(cut["center"])
    plane.SetNormal(cut["normal"])

    cutter = vtk.vtkCutter()
    cutter.SetInputData(surface)
    cutter.SetCutFunction(plane)
    cutter.Update()

    cut_all = vtk.vtkPolyData()
    cut_all.DeepCopy(cutter.GetOutput())

    if (
        cut_all.GetNumberOfPoints() == 0
        or cut_all.GetNumberOfCells() == 0
    ):
        raise NoCutLoopFound(
            "Cut plane did not intersect the surface."
        )

    # Merge cutter segment endpoints so each physical contour becomes a
    # connected region for vtkPolyDataConnectivityFilter.
    cut_all = clean_polydata(
        cut_all,
        MERGE_TOLERANCE
    )

    connectivity = vtk.vtkPolyDataConnectivityFilter()
    connectivity.SetInputData(cut_all)
    connectivity.SetExtractionModeToClosestPointRegion()
    connectivity.SetClosestPoint(
        cut["center"][0],
        cut["center"][1],
        cut["center"][2]
    )
    connectivity.Update()

    loop = vtk.vtkPolyData()
    loop.DeepCopy(connectivity.GetOutput())

    if (
        loop.GetNumberOfPoints() < 3
        or loop.GetNumberOfCells() == 0
    ):
        raise NoCutLoopFound(
            "Could not isolate a valid cut contour near the centerline."
        )

    # --------------------------------------------------------
    # QC geometry for the selected loop.
    # --------------------------------------------------------
    max_radial = 0.0
    min_radial = None
    centroid = [0.0, 0.0, 0.0]

    for i in range(loop.GetNumberOfPoints()):

        p = point_tuple(loop, i)

        centroid[0] += p[0]
        centroid[1] += p[1]
        centroid[2] += p[2]

        radial = distance_from_axis(
            p,
            cut["center"],
            cut["normal"]
        )

        if radial > max_radial:
            max_radial = radial

        if min_radial is None or radial < min_radial:
            min_radial = radial

    npts = float(loop.GetNumberOfPoints())

    centroid = (
        centroid[0] / npts,
        centroid[1] / npts,
        centroid[2] / npts
    )

    centroid_offset = distance(
        centroid,
        cut["center"]
    )

    ratio = (
        max_radial / cut["cut_radius"]
        if cut["cut_radius"] > 0.0
        else float("inf")
    )

    # For true IsOutlet cuts, an already-pruned branch can leave no local
    # surface at the requested endpoint even though the infinite plane may
    # intersect some other distant vessel. Treat that situation as
    # "no local cross-section" so the outlet can be skipped safely.
    if cut.get("type") == "outlet":

        local_limit = (
            cut["cut_radius"]
            * OUTLET_LOCAL_SEARCH_FACTOR
        )

        if (
            min_radial is None
            or min_radial > local_limit
        ):
            raise NoCutLoopFound(
                "No local vessel cross-section was found near this "
                "IsOutlet point. The nearest plane contour was too far "
                "from the outlet centerline location."
            )

    retry_suspicious = (
        cut.get("type") == "inlet"
        or
        (
            cut.get("type") == "side_branch"
            and RETRY_SUSPICIOUS_SIDE_BRANCH_LOOPS
        )
        or
        (
            cut.get("type") == "terminal"
            and RETRY_SUSPICIOUS_TERMINAL_LOOPS
        )
        or
        (
            cut.get("type") == "outlet"
            and RETRY_SUSPICIOUS_OUTLET_LOOPS
        )
    )

    if ratio > MAX_LOOP_TO_EXPECTED_RADIUS_RATIO:

        message = (
            "Selected contour is larger than expected: "
            "max radial distance {:.6g}, expected QC radius {:.6g}, "
            "ratio {:.3f}."
            .format(
                max_radial,
                cut["cut_radius"],
                ratio
            )
        )

        if retry_suspicious:
            exc = SuspiciousCutLoop(message)
            exc.debug_loop = loop
            exc.debug_shape_stats = {}
            raise exc

        if FAIL_ON_SUSPICIOUS_LOOP:
            raise RuntimeError(message)
        else:
            print("    WARNING:", message)

    shape_stats = compute_loop_shape_metrics(
        loop,
        cut
    )

    if ENABLE_LOOP_SHAPE_QC:

        shape_failures = []

        if (
            shape_stats["area_ratio"]
            < MIN_LOOP_AREA_TO_EXPECTED_RATIO
        ):
            shape_failures.append(
                "area is too small (A/(pi*r^2)={:.3f} < {:.3f})"
                .format(
                    shape_stats["area_ratio"],
                    MIN_LOOP_AREA_TO_EXPECTED_RATIO
                )
            )

        if (
            shape_stats["area_ratio"]
            > MAX_LOOP_AREA_TO_EXPECTED_RATIO
        ):
            shape_failures.append(
                "area is too large (A/(pi*r^2)={:.3f} > {:.3f})"
                .format(
                    shape_stats["area_ratio"],
                    MAX_LOOP_AREA_TO_EXPECTED_RATIO
                )
            )

        if (
            shape_stats["circularity"]
            < MIN_LOOP_CIRCULARITY
        ):
            shape_failures.append(
                "shape circularity is low ({:.3f} < {:.3f})"
                .format(
                    shape_stats["circularity"],
                    MIN_LOOP_CIRCULARITY
                )
            )

        if (
            shape_stats["radial_asymmetry"]
            > MAX_LOOP_RADIAL_ASYMMETRY
        ):
            shape_failures.append(
                "one-sided radial asymmetry is high ({:.3f} > {:.3f})"
                .format(
                    shape_stats["radial_asymmetry"],
                    MAX_LOOP_RADIAL_ASYMMETRY
                )
            )

        if (
            shape_stats["shape_centroid_offset_ratio"]
            > MAX_LOOP_CENTROID_OFFSET_TO_RADIUS
        ):
            shape_failures.append(
                "loop centroid is too far from centerline "
                "({:.3f} local radii > {:.3f})"
                .format(
                    shape_stats["shape_centroid_offset_ratio"],
                    MAX_LOOP_CENTROID_OFFSET_TO_RADIUS
                )
            )

        if (
            shape_stats["effective_radius_ratio"]
            < MIN_OPENING_EFFECTIVE_RADIUS_TO_REFERENCE
        ):
            shape_failures.append(
                "opening is too small relative to its branch "
                "(effective radius/reference radius={:.3f} < {:.3f})"
                .format(
                    shape_stats["effective_radius_ratio"],
                    MIN_OPENING_EFFECTIVE_RADIUS_TO_REFERENCE
                )
            )

        if (
            shape_stats["segment_redundancy_ratio"]
            > MAX_LOOP_SEGMENT_REDUNDANCY_RATIO
        ):
            shape_failures.append(
                "overlapping/duplicate contour segments detected "
                "(raw length/perimeter={:.3f} > {:.3f})"
                .format(
                    shape_stats["segment_redundancy_ratio"],
                    MAX_LOOP_SEGMENT_REDUNDANCY_RATIO
                )
            )

        if (
            shape_stats["self_intersections"]
            > MAX_LOOP_SELF_INTERSECTIONS
        ):
            shape_failures.append(
                "self-intersecting contour detected "
                "({} non-adjacent edge intersection(s))"
                .format(shape_stats["self_intersections"])
            )

        if shape_stats["near_duplicate_vertices"] > 0:
            shape_failures.append(
                "self-touching/overlapping contour detected "
                "({} non-adjacent near-duplicate vertex pair(s))"
                .format(shape_stats["near_duplicate_vertices"])
            )

        if shape_failures:
            exc = SuspiciousCutLoop(
                "Suspicious cross-section shape: "
                + "; ".join(shape_failures)
            )
            # Keep actual rejected cutter polyline for inlet diagnostics.
            exc.debug_loop = loop
            exc.debug_shape_stats = shape_stats
            raise exc

    stats = {
        "num_loop_points": loop.GetNumberOfPoints(),
        "num_loop_cells": loop.GetNumberOfCells(),
        "min_radial": min_radial,
        "max_radial": max_radial,
        "centroid": centroid,
        "centroid_offset": centroid_offset
    }

    stats.update(shape_stats)

    return loop, plane, stats


# ============================================================
# TRIANGLE / PLANE INTERSECTION
#
# Used to map the selected vtkCutter contour back to the original surface
# triangles. This is deliberately geometric instead of relying on cutter
# cell-data propagation, which varies between VTK pipelines.
# ============================================================

def triangle_plane_intersection_points(surface, cell_id, plane):

    cell = surface.GetCell(cell_id)

    if cell.GetNumberOfPoints() != 3:
        return []

    pts = [
        point_tuple(surface, cell.GetPointId(i))
        for i in range(3)
    ]

    signed = [
        float(plane.EvaluateFunction(p))
        for p in pts
    ]

    # Completely coplanar triangles are extremely unlikely in a vessel
    # surface. Ignore them rather than ambiguously assigning them to a cut.
    if all(abs(s) <= PLANE_EPS for s in signed):
        return []

    raw = []

    for i, j in ((0, 1), (1, 2), (2, 0)):

        p0 = pts[i]
        p1 = pts[j]
        s0 = signed[i]
        s1 = signed[j]

        if abs(s0) <= PLANE_EPS:
            raw.append(p0)

        if (
            (s0 < -PLANE_EPS and s1 > PLANE_EPS)
            or
            (s0 > PLANE_EPS and s1 < -PLANE_EPS)
        ):

            denom = s0 - s1

            if abs(denom) > PLANE_EPS:

                t = s0 / denom

                p = (
                    p0[0] + t * (p1[0] - p0[0]),
                    p0[1] + t * (p1[1] - p0[1]),
                    p0[2] + t * (p1[2] - p0[2])
                )

                raw.append(p)

    # Deduplicate intersection points generated from shared/on-plane vertices.
    unique = []

    dedup_tol = max(
        MERGE_TOLERANCE * 10.0,
        1e-12
    )

    for p in raw:

        is_duplicate = False

        for q in unique:
            if distance(p, q) <= dedup_tol:
                is_duplicate = True
                break

        if not is_duplicate:
            unique.append(p)

    if len(unique) < 2:
        return []

    # A triangle/plane intersection is a line segment. If numerical edge
    # cases gave >2 unique points, take the farthest pair.
    best_pair = None
    best_dist = -1.0

    for i in range(len(unique)):
        for j in range(i + 1, len(unique)):

            d = distance(unique[i], unique[j])

            if d > best_dist:
                best_dist = d
                best_pair = (unique[i], unique[j])

    if best_pair is None or best_dist <= PLANE_EPS:
        return []

    return [best_pair[0], best_pair[1]]


# ============================================================
# MAP THE SELECTED CONTOUR BACK TO ITS SOURCE SURFACE TRIANGLES
#
# For candidate triangles near the loop, compute their exact plane
# intersection segment and ask whether its midpoint lies on the selected
# contour. Only those triangles are allowed to be clipped.
# ============================================================

def find_source_cells_for_loop(surface, loop, plane, cut):

    loop_locator = vtk.vtkCellLocator()
    loop_locator.SetDataSet(loop)
    loop_locator.BuildLocator()

    surface_locator = vtk.vtkStaticCellLocator()
    surface_locator.SetDataSet(surface)
    surface_locator.BuildLocator()

    bounds = list(loop.GetBounds())

    pad = max(
        cut["radius"] * 0.05,
        MERGE_TOLERANCE * 100.0,
        1e-10
    )

    expanded_bounds = (
        bounds[0] - pad,
        bounds[1] + pad,
        bounds[2] - pad,
        bounds[3] + pad,
        bounds[4] - pad,
        bounds[5] + pad
    )

    candidate_ids = vtk.vtkIdList()

    surface_locator.FindCellsWithinBounds(
        expanded_bounds,
        candidate_ids
    )

    match_tol = max(
        cut["radius"] * 1e-5,
        MERGE_TOLERANCE * 100.0,
        1e-9
    )

    crossing_cell_ids = set()

    for i in range(candidate_ids.GetNumberOfIds()):

        cell_id = int(candidate_ids.GetId(i))

        segment = triangle_plane_intersection_points(
            surface,
            cell_id,
            plane
        )

        if len(segment) != 2:
            continue

        midpoint = (
            0.5 * (segment[0][0] + segment[1][0]),
            0.5 * (segment[0][1] + segment[1][1]),
            0.5 * (segment[0][2] + segment[1][2])
        )

        closest = [0.0, 0.0, 0.0]
        closest_cell_id = vtk.reference(0)
        closest_sub_id = vtk.reference(0)
        dist2 = vtk.reference(0.0)

        loop_locator.FindClosestPoint(
            midpoint,
            closest,
            closest_cell_id,
            closest_sub_id,
            dist2
        )

        if float(dist2) <= match_tol * match_tol:
            crossing_cell_ids.add(cell_id)

    if len(crossing_cell_ids) == 0:
        raise RuntimeError(
            "Could not map the selected cut loop back to source triangles."
        )

    return crossing_cell_ids


# ============================================================
# FIND A SACRIFICIAL ONE-RING OF DOWNSTREAM NEIGHBOR CELLS
#
# The exact clip creates the desired planar boundary on the RETAINED side.
# We additionally discard one layer of triangles immediately downstream of
# the intersected ring. Those cells belong to the branch that will be
# removed anyway, and the tiny downstream gap guarantees that VTK
# connectivity cannot keep proximal and distal pieces joined through a
# plane vertex that happened to exist in the original mesh.
# ============================================================

def find_downstream_neighbor_cells(surface, crossing_ids, plane):

    surface.BuildLinks()

    candidate_neighbors = set()

    point_cells = vtk.vtkIdList()

    for cell_id in crossing_ids:

        cell = surface.GetCell(cell_id)

        for local_id in range(cell.GetNumberOfPoints()):

            point_id = cell.GetPointId(local_id)

            point_cells.Reset()
            surface.GetPointCells(
                point_id,
                point_cells
            )

            for i in range(point_cells.GetNumberOfIds()):

                neighbor_id = int(
                    point_cells.GetId(i)
                )

                if neighbor_id not in crossing_ids:
                    candidate_neighbors.add(neighbor_id)

    downstream_neighbors = set()

    for cell_id in candidate_neighbors:

        cell = surface.GetCell(cell_id)

        centroid = [0.0, 0.0, 0.0]
        npts = cell.GetNumberOfPoints()

        if npts == 0:
            continue

        for local_id in range(npts):

            p = point_tuple(
                surface,
                cell.GetPointId(local_id)
            )

            centroid[0] += p[0]
            centroid[1] += p[1]
            centroid[2] += p[2]

        centroid = (
            centroid[0] / float(npts),
            centroid[1] / float(npts),
            centroid[2] / float(npts)
        )

        signed_centroid = float(
            plane.EvaluateFunction(centroid)
        )

        # Normal points into the branch being discarded, so positive is
        # downstream. We remove only that distal one-ring.
        if signed_centroid > PLANE_EPS:
            downstream_neighbors.add(cell_id)

    return downstream_neighbors


# ============================================================
# EXTRACT A SUBSET OF POLYDATA CELLS BY A BINARY CELL FLAG
# ============================================================

def extract_cells_by_flag(surface, selected_ids, want_selected):
    """
    Extract selected or unselected cells using vtkExtractSelection.

    The older implementation deep-copied the complete surface, allocated a
    temporary flag array, and then looped over EVERY surface cell in Python for
    every cut. With ~600k triangles and ~1000 outlets that was one of the main
    runtime bottlenecks. vtkSelection performs the same subset/complement
    operation in compiled VTK code and only requires us to list the relatively
    small set of affected cell IDs.
    """

    selected_ids = sorted(int(cid) for cid in selected_ids)

    if len(selected_ids) == 0:
        if want_selected:
            return vtk.vtkPolyData()
        out = vtk.vtkPolyData()
        out.ShallowCopy(surface)
        return out

    id_array = vtk.vtkIdTypeArray()
    id_array.SetNumberOfComponents(1)
    id_array.SetNumberOfTuples(len(selected_ids))

    for i, cid in enumerate(selected_ids):
        id_array.SetValue(i, cid)

    node = vtk.vtkSelectionNode()
    node.SetFieldType(vtk.vtkSelectionNode.CELL)
    node.SetContentType(vtk.vtkSelectionNode.INDICES)
    node.SetSelectionList(id_array)

    if not want_selected:
        node.GetProperties().Set(vtk.vtkSelectionNode.INVERSE(), 1)

    selection = vtk.vtkSelection()
    selection.AddNode(node)

    extract = vtk.vtkExtractSelection()
    extract.SetInputData(0, surface)
    extract.SetInputData(1, selection)
    extract.Update()

    geometry = vtk.vtkGeometryFilter()
    geometry.SetInputConnection(extract.GetOutputPort())
    geometry.Update()

    out = vtk.vtkPolyData()
    out.ShallowCopy(geometry.GetOutput())

    # vtkExtractSelection may add bookkeeping arrays. They are harmless, but
    # removing them prevents repeated cuts from accumulating unused metadata.
    out.GetCellData().RemoveArray("vtkOriginalCellIds")
    out.GetPointData().RemoveArray("vtkOriginalPointIds")

    return out


# ============================================================
# BUILD A FLAT CAP FROM THE EXACT CUT CONTOUR
#
# vtkContourTriangulator uses the NEW interpolated contour points made by
# vtkCutter, so the cap boundary follows the true plane/surface oval rather
# than the pre-existing triangle vertices.
# ============================================================

def build_flat_cap(loop):

    loop_clean = clean_polydata(
        loop,
        MERGE_TOLERANCE
    )

    stripper = vtk.vtkStripper()
    stripper.SetInputData(loop_clean)
    stripper.JoinContiguousSegmentsOn()
    stripper.Update()

    triangulator = vtk.vtkContourTriangulator()
    triangulator.SetInputConnection(
        stripper.GetOutputPort()
    )
    triangulator.Update()

    cap = vtk.vtkPolyData()
    cap.DeepCopy(triangulator.GetOutput())

    cap = triangulate_polydata(cap)

    if cap.GetNumberOfCells() == 0:
        raise RuntimeError(
            "vtkContourTriangulator could not create a cap from the cut loop."
        )

    return cap


# ============================================================
# EXACT LOCAL PLANE CUT
#
# This is the refinement over the previous script.
#
# 1) vtkCutter creates exact interpolated plane/surface intersections.
# 2) Connectivity selects ONLY the contour nearest this centerline cut.
# 3) The source triangles producing that one contour are identified.
# 4) ONLY those triangles are run through vtkClipPolyData.
# 5) Their upstream portions are returned to the otherwise untouched model.
#
# This means the plane can intersect other distant PA branches without
# modifying them, while the target branch receives a true Meshmixer-like
# interpolated planar boundary.
# ============================================================

def apply_exact_local_plane_cut(surface, cut, cleanup_result=True):

    loop, plane, stats = (
        find_target_cut_loop(
            surface,
            cut
        )
    )

    crossing_ids = find_source_cells_for_loop(
        surface,
        loop,
        plane,
        cut
    )

    stats["num_crossing_cells"] = len(crossing_ids)

    downstream_neighbors = find_downstream_neighbor_cells(
        surface,
        crossing_ids,
        plane
    )

    stats["num_downstream_neighbor_cells_removed"] = len(
        downstream_neighbors
    )

    # Cells that generated ONLY the selected target contour.
    crossing_patch = extract_cells_by_flag(
        surface,
        crossing_ids,
        True
    )

    # Remove a sacrificial one-ring immediately downstream in addition to
    # replacing the actual crossing triangles. All other geometry remains
    # untouched, including distant branches crossed by the same plane.
    cells_removed_from_untouched = (
        crossing_ids | downstream_neighbors
    )

    untouched_surface = extract_cells_by_flag(
        surface,
        cells_removed_from_untouched,
        False
    )

    # The cut normal always points downstream / into the branch being
    # removed. Therefore the retained/proximal side has plane value <= 0.
    clipper = vtk.vtkClipPolyData()
    clipper.SetInputData(crossing_patch)
    clipper.SetClipFunction(plane)
    clipper.SetValue(0.0)
    clipper.InsideOutOn()
    clipper.GenerateClippedOutputOff()
    clipper.Update()

    proximal_fragments = vtk.vtkPolyData()
    proximal_fragments.DeepCopy(
        clipper.GetOutput()
    )

    proximal_fragments = triangulate_polydata(
        proximal_fragments
    )

    if proximal_fragments.GetNumberOfCells() == 0:
        raise RuntimeError(
            "Clipping the target ring produced no proximal surface fragments."
        )

    append = vtk.vtkAppendPolyData()
    append.AddInputData(untouched_surface)
    append.AddInputData(proximal_fragments)
    append.Update()

    result = vtk.vtkPolyData()
    result.ShallowCopy(append.GetOutput())

    # The two inputs are already triangles. Welding the entire pulmonary
    # surface is useful for long-term topology, but doing it after every single
    # outlet dominates runtime. The caller batches this operation instead.
    if cleanup_result:
        result = clean_polydata(
            result,
            MERGE_TOLERANCE
        )
        result = triangulate_polydata(result)

    result.GetCellData().RemoveArray("__TargetCutCell")
    result.GetCellData().RemoveArray("vtkOriginalCellIds")
    result.GetPointData().RemoveArray("vtkOriginalPointIds")

    # These per-cut caps are obsolete in the current workflow: final caps are
    # rebuilt from the ACTUAL final boundary loops after end-box and quality
    # cleanup. Avoid ~1000 unnecessary contour triangulations by default.
    if BUILD_LEGACY_PER_CUT_CAPS:
        cap = build_flat_cap(loop)
    else:
        cap = None

    return result, loop, cap, stats


# ============================================================
# INLET DEBUG OUTPUT + EARLY GEOMETRY PREFLIGHT
# ============================================================

def write_inlet_rejection_debug(records):
    """Write rejected contours and actual intersection locations.

    Each line segment has AttemptNumber, PathFraction, ValidationStage,
    and rejection count arrays for ParaView inspection. Input is never
    modified. Returns paths of files that were written.
    """
    if not records:
        return []

    points = vtk.vtkPoints()
    cells = vtk.vtkCellArray()
    count_arr = vtk.vtkIntArray()
    count_arr.SetName("SelfIntersections")
    duplicate_arr = vtk.vtkIntArray()
    duplicate_arr.SetName("NearDuplicateVertices")
    attempt_arr = vtk.vtkIntArray()
    attempt_arr.SetName("AttemptNumber")
    fraction_arr = vtk.vtkDoubleArray()
    fraction_arr.SetName("PathFraction")
    stage_arr = vtk.vtkIntArray()
    stage_arr.SetName("ValidationStage") # 0 preflight, 1 stage3

    ipoints = vtk.vtkPoints()
    iverts = vtk.vtkCellArray()
    itype = vtk.vtkIntArray()
    itype.SetName("IntersectionType") # 1 crossing, 2 touch, 3 overlap, 4 near duplicate
    iattempt = vtk.vtkIntArray()
    iattempt.SetName("AttemptNumber")
    ifraction = vtk.vtkDoubleArray()
    ifraction.SetName("PathFraction")
    istage = vtk.vtkIntArray()
    istage.SetName("ValidationStage")
    type_code = {"crossing": 1, "touch": 2, "overlap": 3}

    for candidate, error, stage in records:
        loop = getattr(error, "debug_loop", None)
        stats = getattr(error, "debug_shape_stats", {}) or {}
        attempt = int(candidate.get("attempt_number", 0))
        fraction = float(candidate.get("retry_path_fraction", INLET_CUT_PATH_FRACTION))
        stage_code = 0 if stage == "preflight" else 1
        if loop is not None:
            for cell_id in range(loop.GetNumberOfCells()):
                cell = loop.GetCell(cell_id)
                if cell is None or cell.GetNumberOfPoints() < 2:
                    continue
                line = vtk.vtkPolyLine()
                line.GetPointIds().SetNumberOfIds(cell.GetNumberOfPoints())
                for k in range(cell.GetNumberOfPoints()):
                    p = loop.GetPoint(cell.GetPointId(k))
                    pid = points.InsertNextPoint(p)
                    line.GetPointIds().SetId(k, pid)
                cells.InsertNextCell(line)
                count_arr.InsertNextValue(int(stats.get("self_intersections", 0)))
                duplicate_arr.InsertNextValue(int(stats.get("near_duplicate_vertices", 0)))
                attempt_arr.InsertNextValue(attempt)
                fraction_arr.InsertNextValue(fraction)
                stage_arr.InsertNextValue(stage_code)
        marker_positions = list(stats.get("intersection_points_3d", []))
        marker_positions.extend(
            (position, "near_duplicate")
            for position in stats.get("near_duplicate_points_3d", [])
        )
        for position, kind in marker_positions:
            pid = ipoints.InsertNextPoint(position)
            vert = vtk.vtkVertex()
            vert.GetPointIds().SetId(0, pid)
            iverts.InsertNextCell(vert)
            itype.InsertNextValue(
                4 if kind == "near_duplicate" else type_code.get(kind, 0)
            )
            iattempt.InsertNextValue(attempt)
            ifraction.InsertNextValue(fraction)
            istage.InsertNextValue(stage_code)

    written = []
    if cells.GetNumberOfCells() > 0:
        out = vtk.vtkPolyData()
        out.SetPoints(points)
        out.SetLines(cells)
        for arr in (count_arr, duplicate_arr, attempt_arr, fraction_arr, stage_arr):
            out.GetCellData().AddArray(arr)
        write_polydata(out, INLET_REJECTED_CONTOURS_DEBUG_FILE)
        written.append(INLET_REJECTED_CONTOURS_DEBUG_FILE)
    if iverts.GetNumberOfCells() > 0:
        out = vtk.vtkPolyData()
        out.SetPoints(ipoints)
        out.SetVerts(iverts)
        for arr in (itype, iattempt, ifraction, istage):
            out.GetCellData().AddArray(arr)
        write_polydata(out, INLET_INTERSECTION_POINTS_DEBUG_FILE)
        written.append(INLET_INTERSECTION_POINTS_DEBUG_FILE)
    for path in written:
        print("  Inlet diagnostic:", path)
    return written


def preflight_inlet_cut(surface, cut, collector):
    """Find first plausible inlet cross-section before the costly cut loop.

    Check the *original* triangulated surface and verify the selected cutter
    loop maps to some actual source surface triangles. The final Stage-3 cut
    still performs full QC against the modified post-outlet surface.
    """
    attempts = [dict(cut)] + [dict(c) for c in cut.get("retry_candidates", [])]
    print("\n============================================")
    print("EARLY INLET PREFLIGHT (before pruning)")
    print("============================================")
    for i, candidate in enumerate(attempts):
        candidate["attempt_number"] = i + 1
        try:
            loop, plane, stats = find_target_cut_loop(surface, candidate)
            crossing = find_source_cells_for_loop(surface, loop, plane, candidate)
            if len(crossing) == 0:
                raise RuntimeError("No intersected source triangles for inlet contour")
            fraction = candidate.get("retry_path_fraction", INLET_CUT_PATH_FRACTION)
            print("  Accepted preflight plane: {:.0f}% inlet trunk, {} crossing triangles".format(
                100*float(fraction), len(crossing)))
            return candidate
        except Exception as exc:
            collector.append((dict(candidate), exc, "preflight"))
            print("  Rejected inlet fraction {:.0f}%: {}".format(
                100 * float(candidate.get("retry_path_fraction", INLET_CUT_PATH_FRACTION)),
                exc))
    print("  WARNING: no inlet retry passed early validation.")
    return None


def prioritize_prevalidated_inlet_cut(cut, selected):
    """Try the prevalidated plane first; retain all remaining retries."""
    if selected is None:
        return
    variants = [dict(cut)] + [dict(c) for c in cut.get("retry_candidates", [])]
    selected_fraction = float(selected.get("retry_path_fraction", INLET_CUT_PATH_FRACTION))
    primary = next(c for c in variants if abs(float(c.get("retry_path_fraction", INLET_CUT_PATH_FRACTION)) - selected_fraction) < 1.e-12)
    rest = [c for c in variants if c is not primary]
    cut.update({key: value for key, value in primary.items() if key != "retry_candidates"})
    cut["retry_candidates"] = rest


class CutRetriesExhausted(RuntimeError):

    def __init__(self, cut, attempts, failure_messages):

        self.cut = cut
        self.attempts = attempts
        self.failure_messages = failure_messages

        super().__init__(
            "All {} candidate plane locations failed.\n{}"
            .format(
                len(attempts),
                "\n".join(failure_messages)
            )
        )


def apply_cut_with_retries(surface, cut, cleanup_result=True,
                           inlet_debug_collector=None):

    # Always try the originally planned plane first. Only if that attempt
    # fails do we walk to nearby positions on the same centerline cell.
    attempts = [dict(cut)]
    attempts.extend(
        dict(candidate)
        for candidate in cut.get("retry_candidates", [])
    )

    failure_messages = []
    failed_attempts = []

    for attempt_index, candidate in enumerate(attempts):

        candidate["attempt_number"] = attempt_index + 1
        candidate["planned_cut_index"] = cut.get(
            "planned_cut_index",
            -1
        )

        if attempt_index > 0 and VERBOSE_CUT_LOG:

            path_fraction = candidate.get("retry_path_fraction", None)
            offset = candidate.get("retry_radius_offset", None)

            if path_fraction is not None:
                offset_text = "{:.0f}% of branch length toward parent".format(
                    100.0 * path_fraction
                )
            elif offset is None:
                offset_text = "alternate centerline position"
            elif offset > 0.0:
                offset_text = "+{:.2f} local radii along cell".format(offset)
            else:
                offset_text = "{:.2f} local radii along cell".format(offset)

            print(
                "    Retrying {} cut (attempt {}/{}) at {}."
                .format(
                    candidate.get("type", "local"),
                    attempt_index + 1,
                    len(attempts),
                    offset_text
                )
            )

        try:

            result, loop, cap, stats = (
                apply_exact_local_plane_cut(
                    surface,
                    candidate,
                    cleanup_result=cleanup_result
                )
            )

            stats["retry_attempt"] = attempt_index
            stats["retry_distance_from_anchor"] = distance(
                candidate["center"],
                candidate["junction"]
            )

            return result, loop, cap, stats, candidate

        except Exception as exc:
            if (candidate.get("type") == "inlet"
                    and inlet_debug_collector is not None):
                inlet_debug_collector.append((dict(candidate), exc, "stage3"))

            # In this refinement/debugging version, geometric failures are
            # deliberately treated as retryable. This includes: no local
            # contour, oversized/compound contour, source-cell mapping
            # failure, empty proximal fragments, and cap triangulation.
            # The original input surface is not modified unless the function
            # returns successfully, so rejecting an attempt is safe.
            message = str(exc)

            failed_attempts.append(candidate)

            failure_messages.append(
                "attempt {}: {}".format(
                    attempt_index + 1,
                    message
                )
            )

            if VERBOSE_CUT_LOG:
                print(
                    "    Rejected this plane:",
                    message
                )

            continue

    raise CutRetriesExhausted(
        cut,
        failed_attempts,
        failure_messages
    )


# ============================================================
# ANNOTATE CUT LOOPS FOR DEBUG OUTPUT
# ============================================================

def annotate_loop(loop, cut_index, cut):

    out = vtk.vtkPolyData()
    out.DeepCopy(loop)

    cut_index_array = vtk.vtkIntArray()
    cut_index_array.SetName("CutIndex")
    cut_index_array.SetNumberOfTuples(out.GetNumberOfCells())

    cut_type_array = vtk.vtkIntArray()
    cut_type_array.SetName("CutType")
    cut_type_array.SetNumberOfTuples(out.GetNumberOfCells())

    source_cell_array = vtk.vtkIntArray()
    source_cell_array.SetName("SourceCellId")
    source_cell_array.SetNumberOfTuples(out.GetNumberOfCells())

    if cut["type"] == "side_branch":
        type_value = 1
    elif cut["type"] == "terminal":
        type_value = 2
    elif cut["type"] == "inlet":
        type_value = 5
    else:
        type_value = 3

    for i in range(out.GetNumberOfCells()):
        cut_index_array.SetTuple1(i, cut_index)
        cut_type_array.SetTuple1(i, type_value)
        source_cell_array.SetTuple1(
            i,
            int(cut["source_cell_id"])
        )

    out.GetCellData().AddArray(cut_index_array)
    out.GetCellData().AddArray(cut_type_array)
    out.GetCellData().AddArray(source_cell_array)

    return out


# ============================================================
# APPEND MANY POLYDATA OBJECTS
# ============================================================

def append_polydata_list(polydata_list, clean_result=False):

    append = vtk.vtkAppendPolyData()

    for poly in polydata_list:
        if poly is not None and poly.GetNumberOfCells() > 0:
            append.AddInputData(poly)

    append.Update()

    out = vtk.vtkPolyData()
    out.DeepCopy(append.GetOutput())

    if clean_result:
        out = clean_polydata(
            out,
            MERGE_TOLERANCE
        )

    return out


# ============================================================
# KEEP ONLY THE CONNECTED SURFACE NEAREST THE PA INLET
# ============================================================

def keep_inlet_connected_region(surface, inlet_xyz):

    connectivity = vtk.vtkPolyDataConnectivityFilter()
    connectivity.SetInputData(surface)
    connectivity.SetExtractionModeToClosestPointRegion()
    connectivity.SetClosestPoint(
        inlet_xyz[0],
        inlet_xyz[1],
        inlet_xyz[2]
    )
    connectivity.Update()

    out = vtk.vtkPolyData()
    out.DeepCopy(connectivity.GetOutput())

    return clean_polydata(
        out,
        MERGE_TOLERANCE
    )


# ============================================================
# CAP DATA-COMPATIBILITY HELPERS
#
# vtkAppendPolyData only keeps arrays that are available on every input.
# The cutter/triangulator naturally preserves most point data, but newly
# created cap triangles have no cell-data arrays. Add zero-valued versions
# of numeric arrays from the reference surface so useful metadata is not
# dropped simply because the cap is new geometry.
# ============================================================

def add_missing_numeric_arrays(target_attributes, reference_attributes, num_tuples):

    for i in range(reference_attributes.GetNumberOfArrays()):

        ref = reference_attributes.GetArray(i)

        # GetArray() returns vtkDataArray only, which is exactly what we
        # want here. String/other abstract arrays are intentionally ignored.
        if ref is None:
            continue

        name = ref.GetName()

        if name is None:
            continue

        existing = target_attributes.GetArray(name)

        if existing is not None:
            continue

        new_array = ref.NewInstance()
        new_array.SetName(name)
        new_array.SetNumberOfComponents(
            ref.GetNumberOfComponents()
        )
        new_array.SetNumberOfTuples(num_tuples)

        for comp in range(ref.GetNumberOfComponents()):
            new_array.FillComponent(comp, 0.0)

        target_attributes.AddArray(new_array)


def make_cap_compatible_with_surface(cap, reference_surface, cap_face_id=None):

    out = vtk.vtkPolyData()
    out.DeepCopy(cap)

    add_missing_numeric_arrays(
        out.GetPointData(),
        reference_surface.GetPointData(),
        out.GetNumberOfPoints()
    )

    add_missing_numeric_arrays(
        out.GetCellData(),
        reference_surface.GetCellData(),
        out.GetNumberOfCells()
    )

    if cap_face_id is not None:

        face_arr = out.GetCellData().GetArray(
            "ModelFaceID"
        )

        if face_arr is not None:
            for i in range(out.GetNumberOfCells()):
                face_arr.SetTuple1(i, cap_face_id)

    return out


def get_next_model_face_id(surface):

    arr = surface.GetCellData().GetArray(
        "ModelFaceID"
    )

    if arr is None or arr.GetNumberOfTuples() == 0:
        return None

    max_value = None

    for i in range(arr.GetNumberOfTuples()):
        value = int(round(arr.GetTuple1(i)))
        if max_value is None or value > max_value:
            max_value = value

    return max_value + 1


# ============================================================
# FINAL CAP + SIMVASCULAR FACENAME HELPERS
# ============================================================

def write_facenames_file(path, assignments):
    """Write the Tcl-style SimVascular .facenames sidecar file."""

    folder = os.path.dirname(path)
    if folder:
        os.makedirs(folder, exist_ok=True)

    with open(path, "w", newline="\n") as f:
        f.write("global gPolyDataFaceNames\n\n")
        for fid in sorted(assignments):
            f.write(
                "set gPolyDataFaceNames({}) {}\n".format(
                    int(fid), assignments[fid]
                )
            )


def set_uniform_model_face_id(surface, face_id):
    """Replace ModelFaceID on every cell of a PolyData surface."""

    out = vtk.vtkPolyData()
    out.DeepCopy(surface)

    out.GetCellData().RemoveArray("ModelFaceID")

    arr = vtk.vtkIntArray()
    arr.SetName("ModelFaceID")
    arr.SetNumberOfComponents(1)
    arr.SetNumberOfTuples(out.GetNumberOfCells())
    arr.FillComponent(0, int(face_id))
    out.GetCellData().AddArray(arr)

    return out


def extract_final_boundary_loops(surface):
    """
    Return one ordered vtkPolyData polyline per physical open boundary.

    The final model is already triangulated. vtkFeatureEdges extracts only
    boundary edges; vtkStripper joins each connected ring into one polyline.
    """

    edges = vtk.vtkFeatureEdges()
    edges.SetInputData(surface)
    edges.BoundaryEdgesOn()
    edges.FeatureEdgesOff()
    edges.ManifoldEdgesOff()
    edges.NonManifoldEdgesOff()
    edges.Update()

    boundary = vtk.vtkPolyData()
    boundary.DeepCopy(edges.GetOutput())
    boundary = clean_polydata(boundary, MERGE_TOLERANCE)

    stripper = vtk.vtkStripper()
    stripper.SetInputData(boundary)
    stripper.JoinContiguousSegmentsOn()
    stripper.Update()

    stripped = vtk.vtkPolyData()
    stripped.DeepCopy(stripper.GetOutput())

    loops = []

    for cell_id in range(stripped.GetNumberOfCells()):
        cell = stripped.GetCell(cell_id)
        if cell is None or cell.GetNumberOfPoints() < 3:
            continue

        ordered = []
        for j in range(cell.GetNumberOfPoints()):
            ordered.append(
                point_tuple(stripped, cell.GetPointId(j))
            )

        # vtkStripper normally repeats the first point for a closed ring.
        if (
            len(ordered) >= 2
            and distance(ordered[0], ordered[-1])
                <= max(10.0 * MERGE_TOLERANCE, 1.0e-10)
        ):
            ordered = ordered[:-1]

        if len(ordered) < 3:
            continue

        # Require a geometrically closed ring.  A non-closed boundary here is
        # more useful as an explicit error than as a silently malformed cap.
        closing_gap = distance(ordered[0], ordered[-1])
        edge_lengths = [
            distance(ordered[k], ordered[(k + 1) % len(ordered)])
            for k in range(len(ordered))
        ]
        positive_edges = [e for e in edge_lengths if e > PLANE_EPS]
        typical_edge = median_value(positive_edges) if positive_edges else 0.0

        # For an ordered closed polygon the last->first edge is expected and
        # can be much larger than MERGE_TOLERANCE because the repeated first
        # point was removed above.  The topology was already established by
        # vtkStripper, so we simply construct that final closing edge below.

        pts = vtk.vtkPoints()
        for p in ordered:
            pts.InsertNextPoint(p)

        polyline = vtk.vtkPolyLine()
        polyline.GetPointIds().SetNumberOfIds(len(ordered) + 1)
        for k in range(len(ordered)):
            polyline.GetPointIds().SetId(k, k)
        polyline.GetPointIds().SetId(len(ordered), 0)

        cells = vtk.vtkCellArray()
        cells.InsertNextCell(polyline)

        loop = vtk.vtkPolyData()
        loop.SetPoints(pts)
        loop.SetLines(cells)

        centroid = (
            sum(p[0] for p in ordered) / float(len(ordered)),
            sum(p[1] for p in ordered) / float(len(ordered)),
            sum(p[2] for p in ordered) / float(len(ordered)),
        )

        loops.append({
            "polydata": loop,
            "points": ordered,
            "centroid": centroid,
            "typical_edge": typical_edge,
        })

    return loops


def newell_polygon_normal(points):
    """Stable polygon normal for an ordered nearly planar boundary loop."""

    nx = ny = nz = 0.0
    n = len(points)
    for i in range(n):
        p = points[i]
        q = points[(i + 1) % n]
        nx += (p[1] - q[1]) * (p[2] + q[2])
        ny += (p[2] - q[2]) * (p[0] + q[0])
        nz += (p[0] - q[0]) * (p[1] + q[1])

    normal = (nx, ny, nz)
    if magnitude(normal) <= PLANE_EPS:
        return None
    return normalize(normal)


def build_radial_ring_cap(boundary_points, interior_rings=2):
    """
    Build a planar cap using concentric rings interpolated between the loop
    centroid and the exact wall boundary.

    This keeps every OUTER boundary point unchanged, so the cap welds exactly
    to the final open surface, while avoiding long ear-clipping diagonals that
    can generate one large triangle plus tiny sliver triangles.
    """

    ordered = list(boundary_points)
    n = len(ordered)
    if n < 3:
        raise RuntimeError("Cannot cap a boundary with fewer than 3 points.")

    normal = newell_polygon_normal(ordered)
    if normal is None:
        raise RuntimeError("Boundary loop is degenerate; cap normal is undefined.")

    center = (
        sum(p[0] for p in ordered) / float(n),
        sum(p[1] for p in ordered) / float(n),
        sum(p[2] for p in ordered) / float(n),
    )

    num_rings = max(0, int(interior_rings))

    pts = vtk.vtkPoints()
    center_id = pts.InsertNextPoint(center)

    ring_ids = []

    # Interior rings, from center outward.
    for ring_index in range(1, num_rings + 1):
        frac = float(ring_index) / float(num_rings + 1)
        ids = []
        for p in ordered:
            rp = vec_add(center, vec_scale(vec_sub(p, center), frac))
            ids.append(pts.InsertNextPoint(rp))
        ring_ids.append(ids)

    # Exact outer boundary is always the final ring.
    outer_ids = [pts.InsertNextPoint(p) for p in ordered]
    ring_ids.append(outer_ids)

    polys = vtk.vtkCellArray()

    # Center fan to the innermost ring.
    first_ring = ring_ids[0]
    for i in range(n):
        j = (i + 1) % n
        tri = vtk.vtkTriangle()
        tri.GetPointIds().SetId(0, center_id)
        tri.GetPointIds().SetId(1, first_ring[i])
        tri.GetPointIds().SetId(2, first_ring[j])
        polys.InsertNextCell(tri)

    # Two triangles per quad between successive rings.
    for r in range(len(ring_ids) - 1):
        inner = ring_ids[r]
        outer = ring_ids[r + 1]

        for i in range(n):
            j = (i + 1) % n

            tri1 = vtk.vtkTriangle()
            tri1.GetPointIds().SetId(0, inner[i])
            tri1.GetPointIds().SetId(1, outer[i])
            tri1.GetPointIds().SetId(2, outer[j])
            polys.InsertNextCell(tri1)

            tri2 = vtk.vtkTriangle()
            tri2.GetPointIds().SetId(0, inner[i])
            tri2.GetPointIds().SetId(1, outer[j])
            tri2.GetPointIds().SetId(2, inner[j])
            polys.InsertNextCell(tri2)

    cap = vtk.vtkPolyData()
    cap.SetPoints(pts)
    cap.SetPolys(polys)
    cap = triangulate_polydata(cap)

    return cap


def cap_edge_ratio_stats(cap):
    """Return maximum longest-edge/shortest-edge ratio over cap triangles."""

    worst = 1.0
    for cid in range(cap.GetNumberOfCells()):
        cell = cap.GetCell(cid)
        if cell is None or cell.GetNumberOfPoints() != 3:
            continue

        p0 = point_tuple(cap, cell.GetPointId(0))
        p1 = point_tuple(cap, cell.GetPointId(1))
        p2 = point_tuple(cap, cell.GetPointId(2))

        lengths = [
            distance(p0, p1),
            distance(p1, p2),
            distance(p2, p0),
        ]
        shortest = min(lengths)
        longest = max(lengths)
        if shortest <= PLANE_EPS:
            return float("inf")
        worst = max(worst, longest / shortest)

    return worst


def boundary_edge_count(surface):
    edges = vtk.vtkFeatureEdges()
    edges.SetInputData(surface)
    edges.BoundaryEdgesOn()
    edges.FeatureEdgesOff()
    edges.ManifoldEdgesOff()
    edges.NonManifoldEdgesOff()
    edges.Update()
    return edges.GetOutput().GetNumberOfCells()


def build_final_capped_face_model(
    open_surface,
    accepted_inlet_cut,
    accepted_cuts,
    trimmed_centerline
):
    """Cap actual final openings, name INLET and Left#/Right#/Other#.

    ModelFaceID=1 always Wall; =2 always INLET. All remaining faces are
    incrementally named within their individual centerline Sidedness classes.
    """
    loops = extract_final_boundary_loops(open_surface)
    if not loops:
        raise RuntimeError("Final open surface has no boundary loops to cap.")

    if accepted_inlet_cut is None:
        raise RuntimeError(
            "Cannot label final INLET: no accepted inlet cut exists."
        )

    inlet_cut_center = accepted_inlet_cut["center"]
    inlet_cut_radius = float(accepted_inlet_cut["radius"])
    inlet_loop_index = min(
        range(len(loops)),
        key=lambda i: distance(loops[i]["centroid"], inlet_cut_center)
    )
    inlet_distance = distance(
        loops[inlet_loop_index]["centroid"], inlet_cut_center
    )
    if inlet_distance > INLET_CAP_MAX_DISTANCE_RADII * inlet_cut_radius:
        raise RuntimeError(
            "No final boundary loop is close enough to the inlet cutting "
            "plane to label INLET. Nearest distance={:.6g}, local radius={:.6g}. "
            "Inspect the inlet and endpoint-box outputs."
            .format(inlet_distance, inlet_cut_radius)
        )

    # Stable ordering by corresponding accepted cuts, then loop center.
    cut_refs = sorted(
        (int(key), cut)
        for key, cut in accepted_cuts.items()
        if cut is not None and cut.get("center") is not None
        and cut.get("type") != "inlet"
    )
    remaining = [i for i in range(len(loops)) if i != inlet_loop_index]
    closest_cuts = {}
    for i in remaining:
        best = None
        for planned_index, cut in cut_refs:
            d = distance(loops[i]["centroid"], cut["center"])
            candidate = (d, planned_index)
            if best is None or candidate < best[0]:
                best = (candidate, cut)
        closest_cuts[i] = best

    def loop_sort_key(i):
        center = loops[i]["centroid"]
        matched = closest_cuts[i]
        if matched is None:
            return (float("inf"), float("inf"), *center)
        (distance_to_cut, planned_index), _cut = matched
        return (planned_index, distance_to_cut, *center)

    remaining.sort(key=loop_sort_key)
    ordered_indices = [inlet_loop_index] + remaining

    wall = set_uniform_model_face_id(
        triangulate_polydata(open_surface), FINAL_CAP_WALL_FACE_ID
    )
    assignments = {int(FINAL_CAP_WALL_FACE_ID): "Wall"}
    counters = {"Left": 0, "Right": 0, "Other": 0}
    sidedness_array = get_centerline_sidedness_array(trimmed_centerline)
    locator = build_sidedness_point_locator(trimmed_centerline)

    cap_surfaces = []
    debug_loops = []
    for ordinal, loop_index in enumerate(ordered_indices, start=2):
        info = loops[loop_index]
        face_id = int(ordinal)

        if loop_index == inlet_loop_index:
            face_name = FINAL_CAP_INLET_NAME
        else:
            matched = closest_cuts.get(loop_index)
            accepted_cut = None
            if matched is not None:
                (dist_to_cut, _planned_index), candidate_cut = matched
                # Don't attribute a distant remnant loop to an unrelated
                # outlet's source-cell Sidedness solely because it is the
                # nearest of the accepted cuts.
                radius = float(candidate_cut.get("radius", 0.0))
                if radius > 0.0 and dist_to_cut <= 3.0 * radius:
                    accepted_cut = candidate_cut

            group = classify_loop_sidedness(
                info["centroid"], trimmed_centerline,
                sidedness_array, locator, accepted_cut=accepted_cut
            )
            counters[group] += 1
            prefix = {
                "Left": LEFT_CAP_PREFIX,
                "Right": RIGHT_CAP_PREFIX,
                "Other": OTHER_CAP_PREFIX,
            }[group]
            face_name = "{}{}".format(prefix, counters[group])

        cap = build_radial_ring_cap(
            info["points"], FINAL_CAP_INTERIOR_RINGS
        )
        cap = make_cap_compatible_with_surface(cap, wall, face_id)
        cap_surfaces.append(cap)
        assignments[face_id] = face_name

        debug_loop = vtk.vtkPolyData()
        debug_loop.DeepCopy(info["polydata"])
        fid_arr = vtk.vtkIntArray()
        fid_arr.SetName("ModelFaceID")
        fid_arr.SetNumberOfTuples(debug_loop.GetNumberOfCells())
        fid_arr.FillComponent(0, face_id)
        debug_loop.GetCellData().AddArray(fid_arr)
        debug_loops.append(debug_loop)

        worst_ratio = cap_edge_ratio_stats(cap)
        if worst_ratio > FINAL_CAP_EDGE_RATIO_WARN:
            print(
                "  WARNING: final cap {} (fid={}) has max edge ratio {:.2f}"
                .format(face_name, face_id, worst_ratio)
            )

    capped = append_polydata_list(
        [wall] + cap_surfaces, clean_result=True
    )
    capped = triangulate_polydata(capped)

    remaining_boundary_edges = boundary_edge_count(capped)
    if remaining_boundary_edges != 0:
        raise RuntimeError(
            "Final capped model is not watertight: {} boundary edges remain."
            .format(remaining_boundary_edges)
        )

    normals = vtk.vtkPolyDataNormals()
    normals.SetInputData(capped)
    normals.ConsistencyOn()
    normals.AutoOrientNormalsOn()
    normals.SplittingOff()
    normals.ComputePointNormalsOn()
    normals.ComputeCellNormalsOn()
    normals.Update()

    capped_with_normals = vtk.vtkPolyData()
    capped_with_normals.DeepCopy(normals.GetOutput())
    if REGENERATE_GLOBAL_IDS:
        regenerate_global_ids(capped_with_normals)

    loop_debug = (
        append_polydata_list(debug_loops, clean_result=False)
        if debug_loops else vtk.vtkPolyData()
    )
    print(
        "Cap naming summary: INLET=1, Left={}, Right={}, Other={}".format(
            counters["Left"], counters["Right"], counters["Other"]
        )
    )
    return capped_with_normals, assignments, loop_debug


# ============================================================
# REGENERATE GLOBAL NODE / ELEMENT IDS AFTER GEOMETRY CHANGES
# ============================================================

def regenerate_global_ids(surface):

    pd = surface.GetPointData()
    cd = surface.GetCellData()

    if pd.GetArray("GlobalNodeID") is not None:
        pd.RemoveArray("GlobalNodeID")

    if cd.GetArray("GlobalElementID") is not None:
        cd.RemoveArray("GlobalElementID")

    node_ids = vtk.vtkIdTypeArray()
    node_ids.SetName("GlobalNodeID")
    node_ids.SetNumberOfComponents(1)
    node_ids.SetNumberOfTuples(surface.GetNumberOfPoints())

    for i in range(surface.GetNumberOfPoints()):
        node_ids.SetTuple1(i, i + 1)

    elem_ids = vtk.vtkIdTypeArray()
    elem_ids.SetName("GlobalElementID")
    elem_ids.SetNumberOfComponents(1)
    elem_ids.SetNumberOfTuples(surface.GetNumberOfCells())

    for i in range(surface.GetNumberOfCells()):
        elem_ids.SetTuple1(i, i + 1)

    pd.AddArray(node_ids)
    cd.AddArray(elem_ids)



# ============================================================
# FINAL MESH-QUALITY / SPIKE CLEANUP HELPERS
# ============================================================

def configure_mesh_quality_filter(mesh_quality):
    """Configure vtkMeshQuality to match the requested triangle measure."""

    name = str(MESH_QUALITY_TRIANGLE_MEASURE).strip().lower().replace(" ", "")

    setters = {
        "aspectratio": mesh_quality.SetTriangleQualityMeasureToAspectRatio,
        "edgeratio": mesh_quality.SetTriangleQualityMeasureToEdgeRatio,
        "radiusratio": mesh_quality.SetTriangleQualityMeasureToRadiusRatio,
        "aspectfrobenius": mesh_quality.SetTriangleQualityMeasureToAspectFrobenius,
        "condition": mesh_quality.SetTriangleQualityMeasureToCondition,
        "scaledjacobian": mesh_quality.SetTriangleQualityMeasureToScaledJacobian,
        "shape": mesh_quality.SetTriangleQualityMeasureToShape,
        "distortion": mesh_quality.SetTriangleQualityMeasureToDistortion,
        "minangle": mesh_quality.SetTriangleQualityMeasureToMinAngle,
        "maxangle": mesh_quality.SetTriangleQualityMeasureToMaxAngle,
        "equiangleskew": mesh_quality.SetTriangleQualityMeasureToEquiangleSkew,
    }

    if name not in setters:
        raise RuntimeError(
            "Unsupported MESH_QUALITY_TRIANGLE_MEASURE '{}'. "
            "Supported choices are: {}"
            .format(
                MESH_QUALITY_TRIANGLE_MEASURE,
                ", ".join(sorted(setters.keys()))
            )
        )

    setters[name]()


def compute_mesh_quality(surface):
    """Return a PolyData copy carrying vtkMeshQuality's cell-data Quality."""

    quality_filter = vtk.vtkMeshQuality()
    quality_filter.SetInputData(surface)
    configure_mesh_quality_filter(quality_filter)
    quality_filter.Update()

    out = vtk.vtkPolyData()
    out.DeepCopy(quality_filter.GetOutput())

    quality = out.GetCellData().GetArray("Quality")

    if quality is None:
        raise RuntimeError(
            "vtkMeshQuality did not create the expected cell-data array 'Quality'."
        )

    return out


def find_bad_quality_cell_ids(quality_surface):
    """Return cells with Quality > threshold, including non-finite values."""

    quality = quality_surface.GetCellData().GetArray("Quality")

    if quality is None:
        raise RuntimeError("Quality array is missing from quality surface.")

    bad = set()

    for cell_id in range(quality_surface.GetNumberOfCells()):

        value = float(quality.GetTuple1(cell_id))

        if (not math.isfinite(value)) or value > MESH_QUALITY_THRESHOLD:
            bad.add(cell_id)

    return bad


def cluster_bad_quality_cells(surface, bad_ids):
    """Group bad cells that touch through at least one mesh point."""

    if len(bad_ids) == 0:
        return []

    surface.BuildLinks()

    remaining = set(bad_ids)
    clusters = []
    point_cells = vtk.vtkIdList()

    while remaining:

        seed = remaining.pop()
        stack = [seed]
        cluster = {seed}

        while stack:

            cell_id = stack.pop()
            cell = surface.GetCell(cell_id)

            for local_id in range(cell.GetNumberOfPoints()):

                point_id = cell.GetPointId(local_id)
                point_cells.Reset()
                surface.GetPointCells(point_id, point_cells)

                for i in range(point_cells.GetNumberOfIds()):

                    neighbor_id = int(point_cells.GetId(i))

                    if neighbor_id in remaining:
                        remaining.remove(neighbor_id)
                        cluster.add(neighbor_id)
                        stack.append(neighbor_id)

        clusters.append(cluster)

    return clusters


def triangle_area_and_centroid(surface, cell_id):

    cell = surface.GetCell(cell_id)

    if cell.GetNumberOfPoints() != 3:
        return 0.0, (0.0, 0.0, 0.0)

    p0 = point_tuple(surface, cell.GetPointId(0))
    p1 = point_tuple(surface, cell.GetPointId(1))
    p2 = point_tuple(surface, cell.GetPointId(2))

    e1 = vec_sub(p1, p0)
    e2 = vec_sub(p2, p0)

    area = 0.5 * magnitude(cross(e1, e2))

    centroid = (
        (p0[0] + p1[0] + p2[0]) / 3.0,
        (p0[1] + p1[1] + p2[1]) / 3.0,
        (p0[2] + p1[2] + p2[2]) / 3.0
    )

    return area, centroid


def summarize_bad_quality_cluster(quality_surface, cluster_ids, cluster_index):
    """Area-weighted location and quality summary of one bad-cell cluster."""

    quality = quality_surface.GetCellData().GetArray("Quality")

    weighted = [0.0, 0.0, 0.0]
    total_area = 0.0
    max_quality = -float("inf")

    # Fallback if all bad triangles are numerically zero-area.
    fallback = [0.0, 0.0, 0.0]
    fallback_count = 0

    for cell_id in cluster_ids:

        q = float(quality.GetTuple1(cell_id))

        if not math.isfinite(q):
            max_quality = float("inf")
        elif q > max_quality:
            max_quality = q

        area, centroid = triangle_area_and_centroid(
            quality_surface,
            cell_id
        )

        fallback[0] += centroid[0]
        fallback[1] += centroid[1]
        fallback[2] += centroid[2]
        fallback_count += 1

        if area > 0.0 and math.isfinite(area):
            weighted[0] += area * centroid[0]
            weighted[1] += area * centroid[1]
            weighted[2] += area * centroid[2]
            total_area += area

    if total_area > 0.0:
        centroid = (
            weighted[0] / total_area,
            weighted[1] / total_area,
            weighted[2] / total_area
        )
    else:
        centroid = (
            fallback[0] / float(fallback_count),
            fallback[1] / float(fallback_count),
            fallback[2] / float(fallback_count)
        )

    return {
        "cluster_id": int(cluster_index),
        "cell_ids": set(cluster_ids),
        "num_cells": len(cluster_ids),
        "centroid": centroid,
        "surface_area": total_area,
        "max_quality": max_quality
    }


def build_centerline_cell_locator(centerline):

    locator = vtk.vtkStaticCellLocator()
    locator.SetDataSet(centerline)
    locator.BuildLocator()

    return locator


def closest_centerline_frame(centerline, locator, query_point):
    """
    Find the closest location on a centerline polyline and return a local
    tangent plus linearly interpolated Radius.
    """

    closest = [0.0, 0.0, 0.0]
    cell_id = vtk.reference(0)
    sub_id = vtk.reference(0)
    dist2 = vtk.reference(0.0)

    locator.FindClosestPoint(
        query_point,
        closest,
        cell_id,
        sub_id,
        dist2
    )

    cid = int(cell_id)
    sid = int(sub_id)

    if cid < 0 or cid >= centerline.GetNumberOfCells():
        raise RuntimeError("Could not locate a valid nearest centerline cell.")

    cell = centerline.GetCell(cid)
    npts = cell.GetNumberOfPoints()

    if npts < 2:
        raise RuntimeError(
            "Nearest centerline cell {} has fewer than two points."
            .format(cid)
        )

    # For vtkPolyLine cells, subId identifies the line segment. Clamp for
    # safety because different VTK versions can report an endpoint subId.
    if sid < 0:
        sid = 0
    if sid > npts - 2:
        sid = npts - 2

    id0 = cell.GetPointId(sid)
    id1 = cell.GetPointId(sid + 1)

    p0 = point_tuple(centerline, id0)
    p1 = point_tuple(centerline, id1)

    segment = vec_sub(p1, p0)
    seg_len2 = dot(segment, segment)

    if seg_len2 <= PLANE_EPS * PLANE_EPS:
        raise RuntimeError(
            "Nearest centerline segment in cell {} has zero length."
            .format(cid)
        )

    tangent = normalize(segment)

    # Interpolate Radius at the closest point projected onto this segment.
    rel = vec_sub(tuple(closest), p0)
    t = dot(rel, segment) / seg_len2

    if t < 0.0:
        t = 0.0
    elif t > 1.0:
        t = 1.0

    radius_array = centerline.GetPointData().GetArray(RADIUS_ARRAY)

    if radius_array is None:
        raise RuntimeError(
            "Trimmed centerline is missing required Radius array '{}'."
            .format(RADIUS_ARRAY)
        )

    r0 = float(radius_array.GetTuple1(id0))
    r1 = float(radius_array.GetTuple1(id1))

    valid0 = math.isfinite(r0) and r0 > 0.0
    valid1 = math.isfinite(r1) and r1 > 0.0

    if valid0 and valid1:
        radius = r0 + t * (r1 - r0)
    elif valid0:
        radius = r0
    elif valid1:
        radius = r1
    else:
        raise RuntimeError(
            "Nearest centerline segment has no valid Radius values."
        )

    return {
        "center": (
            float(closest[0]),
            float(closest[1]),
            float(closest[2])
        ),
        "tangent": tangent,
        "radius": float(radius),
        "centerline_cell_id": cid,
        "centerline_sub_id": sid,
        "distance_to_centerline": math.sqrt(float(dist2))
    }


def make_quality_box_plan(cluster_summary, centerline_frame):

    plan = dict(cluster_summary)
    plan.update(centerline_frame)

    plan["radial_half_width"] = (
        plan["radius"]
        * QUALITY_BOX_RADIAL_HALF_WIDTH_FACTOR
    )

    plan["axial_half_length"] = (
        plan["radius"]
        * QUALITY_BOX_AXIAL_HALF_LENGTH_FACTOR
    )

    return plan


def merge_nearby_quality_box_plans(plans):
    """Merge duplicate boxes proposed for the same local centerline region."""

    merged = []

    # Worst quality first, so the representative box for a merged group is
    # anchored to the most problematic cluster.
    plans_sorted = sorted(
        plans,
        key=lambda p: (
            float("inf") if not math.isfinite(p["max_quality"])
            else p["max_quality"]
        ),
        reverse=True
    )

    for plan in plans_sorted:

        found = None

        for existing in merged:

            merge_distance = (
                QUALITY_BOX_MERGE_DISTANCE_FACTOR
                * max(plan["radius"], existing["radius"])
            )

            if distance(plan["center"], existing["center"]) <= merge_distance:
                found = existing
                break

        if found is None:
            plan["merged_cluster_ids"] = [plan["cluster_id"]]
            plan["merged_bad_cell_count"] = plan["num_cells"]
            merged.append(plan)
        else:
            found["merged_cluster_ids"].append(plan["cluster_id"])
            found["merged_bad_cell_count"] += plan["num_cells"]

            if (
                not math.isfinite(plan["max_quality"])
                or plan["max_quality"] > found["max_quality"]
            ):
                found["max_quality"] = plan["max_quality"]

    return merged


def build_box_world_to_local_transform(center, tangent):
    """Create transform mapping world coordinates into the oriented box frame."""

    # Choose a reference vector not parallel to the tangent.
    if abs(tangent[0]) < 0.9:
        reference = (1.0, 0.0, 0.0)
    else:
        reference = (0.0, 1.0, 0.0)

    u = normalize(cross(tangent, reference))
    v = normalize(cross(tangent, u))

    matrix = vtk.vtkMatrix4x4()
    matrix.Identity()

    # Rows form the world -> local rotation. Translation moves the box center
    # to the local origin. vtkImplicitFunction applies this transform to each
    # input point before evaluating the axis-aligned vtkBox.
    for j in range(3):
        matrix.SetElement(0, j, u[j])
        matrix.SetElement(1, j, v[j])
        matrix.SetElement(2, j, tangent[j])

    matrix.SetElement(0, 3, -dot(u, center))
    matrix.SetElement(1, 3, -dot(v, center))
    matrix.SetElement(2, 3, -dot(tangent, center))

    transform = vtk.vtkTransform()
    transform.SetMatrix(matrix)

    return transform, u, v


def build_oriented_vtk_box(plan):

    box = vtk.vtkBox()

    radial = float(plan["radial_half_width"])
    axial = float(plan["axial_half_length"])

    box.SetBounds(
        -radial, radial,
        -radial, radial,
        -axial, axial
    )

    transform, u, v = build_box_world_to_local_transform(
        plan["center"],
        plan["tangent"]
    )

    box.SetTransform(transform)

    return box, u, v


def clip_out_quality_box(surface, plan, cleanup_result=True):
    """
    Remove the inside of one oriented box using vtkBoxClipDataSet.

    This is intentionally closer to ParaView's Box clipping behavior than an
    implicit vtkBox passed to vtkClipPolyData.  vtkBoxClipDataSet geometrically
    cuts cells against the six box planes, so a long skinny triangle cannot
    simply "jump across" the box because all of its original vertices happened
    to lie outside.
    """

    center = plan["center"]
    tangent = plan["tangent"]
    radial = float(plan["radial_half_width"])
    axial = float(plan["axial_half_length"])

    _, u, v = build_box_world_to_local_transform(
        center,
        tangent
    )

    # Six outward-facing planes of the oriented box.
    normals = [
        vec_scale(u, -1.0),
        u,
        vec_scale(v, -1.0),
        v,
        vec_scale(tangent, -1.0),
        tangent
    ]

    origins = [
        vec_add(center, vec_scale(u, -radial)),
        vec_add(center, vec_scale(u, radial)),
        vec_add(center, vec_scale(v, -radial)),
        vec_add(center, vec_scale(v, radial)),
        vec_add(center, vec_scale(tangent, -axial)),
        vec_add(center, vec_scale(tangent, axial))
    ]

    box_clipper = vtk.vtkBoxClipDataSet()
    box_clipper.SetInputData(surface)
    box_clipper.SetOrientation(1)
    box_clipper.SetBoxClip(
        normals[0], origins[0],
        normals[1], origins[1],
        normals[2], origins[2],
        normals[3], origins[3],
        normals[4], origins[4],
        normals[5], origins[5]
    )

    # Port 0 = inside the box. Port 1 = the portion clipped away/outside.
    # We want to REMOVE the box interior and retain port 1.
    box_clipper.GenerateClippedOutputOn()
    box_clipper.Update()

    outside_grid = box_clipper.GetOutput(1)

    geometry = vtk.vtkGeometryFilter()
    geometry.SetInputData(outside_grid)
    geometry.Update()

    out = vtk.vtkPolyData()
    out.ShallowCopy(geometry.GetOutput())

    # vtkBoxClipDataSet performs the geometric split itself. Global point
    # welding/triangulation can therefore be deferred safely during the large
    # endpoint sweep and done periodically + once at the end. Quality-cleanup
    # boxes still use cleanup_result=True by default.
    if cleanup_result:
        out = clean_polydata(out, MERGE_TOLERANCE)
        out = triangulate_polydata(out)

    return out


def make_debug_box_polydata(plan, box_index):
    """Create an oriented cube surface for visual QC in ParaView."""

    radial = float(plan["radial_half_width"])
    axial = float(plan["axial_half_length"])

    cube = vtk.vtkCubeSource()
    cube.SetBounds(
        -radial, radial,
        -radial, radial,
        -axial, axial
    )
    cube.Update()

    _, u, v = build_box_world_to_local_transform(
        plan["center"],
        plan["tangent"]
    )

    # Local -> world matrix (columns are local basis vectors).
    matrix = vtk.vtkMatrix4x4()
    matrix.Identity()

    for i in range(3):
        matrix.SetElement(i, 0, u[i])
        matrix.SetElement(i, 1, v[i])
        matrix.SetElement(i, 2, plan["tangent"][i])
        matrix.SetElement(i, 3, plan["center"][i])

    transform = vtk.vtkTransform()
    transform.SetMatrix(matrix)

    tf = vtk.vtkTransformPolyDataFilter()
    tf.SetInputConnection(cube.GetOutputPort())
    tf.SetTransform(transform)
    tf.Update()

    out = vtk.vtkPolyData()
    out.DeepCopy(tf.GetOutput())

    arrays = [
        ("QualityBoxIndex", int(box_index), vtk.vtkIntArray),
        ("BadClusterId", int(plan["cluster_id"]), vtk.vtkIntArray),
        ("CenterlineCellId", int(plan["centerline_cell_id"]), vtk.vtkIntArray),
        ("BadCellCount", int(plan.get("merged_bad_cell_count", plan["num_cells"])), vtk.vtkIntArray),
        ("LocalRadius", float(plan["radius"]), vtk.vtkDoubleArray),
        ("MaxBadQuality", float(plan["max_quality"]), vtk.vtkDoubleArray),
        ("DistanceToCenterline", float(plan["distance_to_centerline"]), vtk.vtkDoubleArray),
    ]

    for name, value, array_type in arrays:
        arr = array_type()
        arr.SetName(name)
        arr.SetNumberOfTuples(out.GetNumberOfCells())
        for i in range(out.GetNumberOfCells()):
            arr.SetTuple1(i, value)
        out.GetCellData().AddArray(arr)

    return out


def extract_bad_quality_debug_surface(quality_surface, bad_ids):

    if len(bad_ids) == 0:
        return vtk.vtkPolyData()

    return extract_cells_by_flag(
        quality_surface,
        bad_ids,
        True
    )


def build_optional_terminal_end_box_cuts(centerline):
    """Build box-only cut definitions from PointData TERMINAL_END_ARRAY."""

    pd = centerline.GetPointData()
    marker = pd.GetArray(TERMINAL_END_ARRAY)

    if marker is None:
        return []

    radius_array = pd.GetArray(RADIUS_ARRAY)
    if radius_array is None:
        return []

    is_outlet = pd.GetArray("IsOutlet")
    cuts = []

    for point_id in range(centerline.GetNumberOfPoints()):

        if int(round(marker.GetTuple1(point_id))) != 1:
            continue

        # True IsOutlet points are already covered by outlet_cuts.
        if (
            is_outlet is not None
            and int(round(is_outlet.GetTuple1(point_id))) == 1
        ):
            continue

        cell_ids = vtk.vtkIdList()
        centerline.GetPointCells(point_id, cell_ids)
        endpoint_candidates = []

        for j in range(cell_ids.GetNumberOfIds()):
            cell_id = int(cell_ids.GetId(j))
            cell = centerline.GetCell(cell_id)
            n = cell.GetNumberOfPoints()

            if n < 2:
                continue

            first_id = int(cell.GetPointId(0))
            last_id = int(cell.GetPointId(n - 1))

            if point_id == first_id:
                endpoint_candidates.append((cell_id, int(cell.GetPointId(1)), 0))
            elif point_id == last_id:
                endpoint_candidates.append((cell_id, int(cell.GetPointId(n - 2)), -1))

        if len(endpoint_candidates) == 0:
            continue

        source_cell_id, interior_id, endpoint_local_id = endpoint_candidates[0]
        endpoint = point_tuple(centerline, point_id)
        interior = point_tuple(centerline, interior_id)

        direction = vec_sub(endpoint, interior)
        if magnitude(direction) <= PLANE_EPS:
            continue

        normal = normalize(direction)

        r0 = float(radius_array.GetTuple1(point_id))
        r1 = float(radius_array.GetTuple1(interior_id))
        valid_radii = [
            r for r in (r0, r1)
            if math.isfinite(r) and r > 0.0
        ]

        if len(valid_radii) == 0:
            continue

        radius = max(valid_radii)
        reference_radius = robust_cell_reference_radius(
            centerline,
            source_cell_id,
            radius_array
        )

        if reference_radius is None:
            reference_radius = radius

        cuts.append({
            "type": "terminal_marker",
            "center": endpoint,
            "junction": endpoint,
            "normal": normal,
            "radius": radius,
            "reference_radius": reference_radius,
            "cut_radius": radius * CUT_RADIUS_FACTOR,
            "source_cell_id": source_cell_id,
            "source_endpoint_local_id": endpoint_local_id,
            "source_point_id": point_id,
            "removed_cell_ids": []
        })

    return cuts


def get_end_box_radius(cut):
    """Choose a robust physical radius for the final endpoint box."""

    values = []

    local_radius = cut.get("radius", None)
    if (
        local_radius is not None
        and math.isfinite(float(local_radius))
        and float(local_radius) > 0.0
    ):
        values.append(float(local_radius))

    reference_radius = cut.get("reference_radius", None)
    if (
        END_BOX_USE_REFERENCE_RADIUS
        and reference_radius is not None
        and math.isfinite(float(reference_radius))
        and float(reference_radius) > 0.0
    ):
        values.append(float(reference_radius))

    if len(values) == 0:
        return None

    return max(values)


def make_final_end_box_plan(cut, planned_cut_index):
    """
    Create a downstream-only box whose upstream face is the desired cut plane.
    """

    radius = get_end_box_radius(cut)
    if radius is None:
        return None

    tangent = normalize(cut["normal"])
    radial = radius * END_BOX_RADIAL_HALF_WIDTH_FACTOR
    axial = radius * END_BOX_AXIAL_HALF_LENGTH_FACTOR

    plane_center = cut["center"]

    # Shift the box downstream by one half-length. Therefore its upstream face
    # stays exactly at the accepted/planned cut location while the enlarged box
    # extends outward over the terminal junk we want removed.
    box_center = vec_add(
        plane_center,
        vec_scale(tangent, axial)
    )

    return {
        "center": box_center,
        "plane_center": plane_center,
        "tangent": tangent,
        "radius": radius,
        "radial_half_width": radial,
        "axial_half_length": axial,
        "cut_type": cut.get("type", "end"),
        "source_cell_id": int(cut.get("source_cell_id", -1)),
        "source_point_id": int(cut.get("source_point_id", -1)),
        "planned_cut_index": int(planned_cut_index),
    }


def make_end_box_debug_polydata(plan, box_index):
    """Create an oriented cube for visual QC of the final endpoint sweep."""

    radial = float(plan["radial_half_width"])
    axial = float(plan["axial_half_length"])

    cube = vtk.vtkCubeSource()
    cube.SetBounds(-radial, radial, -radial, radial, -axial, axial)
    cube.Update()

    _, u, v = build_box_world_to_local_transform(
        plan["center"],
        plan["tangent"]
    )

    matrix = vtk.vtkMatrix4x4()
    matrix.Identity()

    for i in range(3):
        matrix.SetElement(i, 0, u[i])
        matrix.SetElement(i, 1, v[i])
        matrix.SetElement(i, 2, plan["tangent"][i])
        matrix.SetElement(i, 3, plan["center"][i])

    transform = vtk.vtkTransform()
    transform.SetMatrix(matrix)

    tf = vtk.vtkTransformPolyDataFilter()
    tf.SetInputConnection(cube.GetOutputPort())
    tf.SetTransform(transform)
    tf.Update()

    out = vtk.vtkPolyData()
    out.DeepCopy(tf.GetOutput())

    type_map = {
        "side_branch": 1,
        "terminal": 2,
        "outlet": 3,
        "terminal_marker": 4,
        "inlet": 5,
    }

    arrays = [
        ("EndBoxIndex", int(box_index), vtk.vtkIntArray),
        ("PlannedCutIndex", int(plan["planned_cut_index"]), vtk.vtkIntArray),
        ("CutType", int(type_map.get(plan["cut_type"], 0)), vtk.vtkIntArray),
        ("SourceCellId", int(plan["source_cell_id"]), vtk.vtkIntArray),
        ("SourcePointId", int(plan["source_point_id"]), vtk.vtkIntArray),
        ("LocalRadius", float(plan["radius"]), vtk.vtkDoubleArray),
        ("RadialHalfWidth", float(plan["radial_half_width"]), vtk.vtkDoubleArray),
        ("AxialHalfLength", float(plan["axial_half_length"]), vtk.vtkDoubleArray),
    ]

    for name, value, array_type in arrays:
        arr = array_type()
        arr.SetName(name)
        arr.SetNumberOfTuples(out.GetNumberOfCells())
        for i in range(out.GetNumberOfCells()):
            arr.SetTuple1(i, value)
        out.GetCellData().AddArray(arr)

    return out


def run_final_end_box_sweep(
    surface,
    all_planned_cuts,
    accepted_cut_by_index,
    extra_terminal_cuts,
    inlet_xyz
):
    """
    Apply a large downstream-only oriented box to every intended end.

    Successful exact cuts use their ACTUAL accepted retry location. A skipped
    exact cut uses the safest/deepest available retry location, allowing the
    box pass to rescue difficult ends without requiring a valid oval contour.
    """

    working = vtk.vtkPolyData()
    working.DeepCopy(surface)

    locator_surface = vtk.vtkPolyData()
    locator_surface.DeepCopy(surface)

    locator = vtk.vtkStaticPointLocator()
    locator.SetDataSet(locator_surface)
    locator.BuildLocator()

    planned_entries = []

    for planned_index, planned_cut in enumerate(all_planned_cuts):

        if planned_index in accepted_cut_by_index:
            # Successful exact cut: place the final box at the ACTUAL accepted
            # retry location, not the original plan.
            chosen = accepted_cut_by_index[planned_index]

        else:
            # Exact cut was skipped. The box is our rescue mechanism, so use
            # the safest/deepest available retry location instead of blindly
            # returning to the original problematic plane.
            retry_candidates = planned_cut.get("retry_candidates", [])

            if len(retry_candidates) == 0:
                chosen = planned_cut

            elif planned_cut.get("type") == "side_branch":
                # Farther INTO the branch being removed is safer than moving
                # back toward its bifurcation with the retained tree.
                chosen = max(
                    retry_candidates,
                    key=lambda c: float(
                        c.get("retry_distance_from_anchor", 0.0)
                    )
                )

            else:
                # Terminal / IsOutlet retries are measured from the endpoint
                # back toward the parent. Choose the farthest parentward
                # available location to avoid tiny/overlapping terminal holes.
                chosen = max(
                    retry_candidates,
                    key=lambda c: float(
                        c.get("retry_distance_from_anchor", 0.0)
                    )
                )

        planned_entries.append((planned_index, chosen))

    base_count = len(all_planned_cuts)
    for i, cut in enumerate(extra_terminal_cuts):
        planned_entries.append((base_count + i, cut))

    debug_boxes = []
    applied = 0
    locality_skips = 0
    empty_skips = 0

    for local_index, (planned_index, cut) in enumerate(planned_entries):

        plan = make_final_end_box_plan(cut, planned_index)
        if plan is None:
            locality_skips += 1
            continue

        # Locality check uses the desired UPSTREAM FACE / cut plane location,
        # not the box center shifted downstream.
        nearest_id = locator.FindClosestPoint(plan["plane_center"])

        if nearest_id < 0:
            locality_skips += 1
            continue

        nearest = point_tuple(locator_surface, nearest_id)
        nearest_distance = distance(nearest, plan["plane_center"])

        if nearest_distance > END_BOX_LOCALITY_FACTOR * plan["radius"]:
            locality_skips += 1
            continue

        candidate = clip_out_quality_box(working, plan, cleanup_result=False)

        if (
            candidate.GetNumberOfPoints() == 0
            or candidate.GetNumberOfCells() == 0
        ):
            empty_skips += 1
            continue

        working = candidate

        debug_boxes.append(
            make_end_box_debug_polydata(plan, applied)
        )
        applied += 1

        if (
            END_BOX_GLOBAL_CLEAN_INTERVAL > 0
            and applied % END_BOX_GLOBAL_CLEAN_INTERVAL == 0
        ):
            working = clean_polydata(working, MERGE_TOLERANCE)
            working = triangulate_polydata(working)

        if (
            PROGRESS_PRINT_EVERY > 0
            and (local_index + 1) % PROGRESS_PRINT_EVERY == 0
        ):
            print(
                "Final end-box progress: {}/{} | applied={} skipped={}"
                .format(
                    local_index + 1,
                    len(planned_entries),
                    applied,
                    locality_skips + empty_skips
                )
            )

    # Discard all distal pieces severed by the endpoint boxes in one pass.
    working = keep_inlet_connected_region(working, inlet_xyz)
    working = clean_polydata(working, MERGE_TOLERANCE)
    working = triangulate_polydata(working)

    debug = append_polydata_list(debug_boxes, clean_result=False)

    return working, debug, applied, locality_skips, empty_skips


def run_mesh_quality_spike_cleanup(surface, centerline, inlet_xyz):
    """
    Detect Quality > threshold cells, group them, build oriented clipping boxes
    from the nearest centerline Radius/tangent, and keep the inlet-connected
    surface after each box.  Returns the cleaned surface plus debug objects.
    """

    working = vtk.vtkPolyData()
    working.DeepCopy(surface)

    centerline_locator = build_centerline_cell_locator(centerline)

    all_debug_boxes = []
    first_bad_debug = None

    for cleanup_pass in range(int(MESH_QUALITY_CLEANUP_PASSES)):

        quality_surface = compute_mesh_quality(working)
        bad_ids = find_bad_quality_cell_ids(quality_surface)

        print("")
        print("============================================")
        print(
            "MESH QUALITY CLEANUP PASS {}/{}"
            .format(
                cleanup_pass + 1,
                MESH_QUALITY_CLEANUP_PASSES
            )
        )
        print("============================================")
        print("Quality measure:", MESH_QUALITY_TRIANGLE_MEASURE)
        print("Quality threshold:", MESH_QUALITY_THRESHOLD)
        print("Bad cells found:", len(bad_ids))

        if cleanup_pass == 0:
            first_bad_debug = extract_bad_quality_debug_surface(
                quality_surface,
                bad_ids
            )

        if len(bad_ids) == 0:
            break

        clusters = cluster_bad_quality_cells(
            quality_surface,
            bad_ids
        )

        print("Connected bad-cell clusters:", len(clusters))

        plans = []

        for cluster_index, cluster_ids in enumerate(clusters):

            summary = summarize_bad_quality_cluster(
                quality_surface,
                cluster_ids,
                cluster_index
            )

            try:
                frame = closest_centerline_frame(
                    centerline,
                    centerline_locator,
                    summary["centroid"]
                )
            except Exception as exc:
                print(
                    "    SKIPPING bad cluster {}: could not find a valid "
                    "centerline frame: {}"
                    .format(cluster_index, exc)
                )
                continue

            plan = make_quality_box_plan(
                summary,
                frame
            )

            plans.append(plan)

            print("")
            print("Bad cluster {}".format(cluster_index))
            print("    Bad cells:", summary["num_cells"])
            print("    Max Quality:", summary["max_quality"])
            print("    Cluster centroid:", summary["centroid"])
            print("    Nearest centerline point:", frame["center"])
            print("    Distance to centerline:", frame["distance_to_centerline"])
            print("    Local Radius:", frame["radius"])
            print("    Centerline tangent:", frame["tangent"])

        plans = merge_nearby_quality_box_plans(plans)

        print("")
        print("Quality clip boxes after merging:", len(plans))

        for box_local_index, plan in enumerate(plans):

            before_cells = working.GetNumberOfCells()

            debug_box = make_debug_box_polydata(
                plan,
                len(all_debug_boxes)
            )
            all_debug_boxes.append(debug_box)

            print("")
            print(
                "Quality box {}/{}"
                .format(box_local_index + 1, len(plans))
            )
            print("    Center:", plan["center"])
            print("    Local Radius:", plan["radius"])
            print("    Radial half-width:", plan["radial_half_width"])
            print("    Axial half-length:", plan["axial_half_length"])
            print("    Merged bad clusters:", plan.get("merged_cluster_ids", []))

            clipped = clip_out_quality_box(
                working,
                plan
            )

            # This is the important topology step: a box spanning the complete
            # local vessel section disconnects the distal spike/branch. Keep
            # only what remains attached to the pulmonary-artery inlet.
            clipped = keep_inlet_connected_region(
                clipped,
                inlet_xyz
            )

            clipped = triangulate_polydata(clipped)

            working = clipped

            print("    Cells before box:", before_cells)
            print("    Cells after box/connectivity:", working.GetNumberOfCells())

    final_quality_surface = compute_mesh_quality(working)
    remaining_bad_ids = find_bad_quality_cell_ids(final_quality_surface)

    remaining_bad_debug = extract_bad_quality_debug_surface(
        final_quality_surface,
        remaining_bad_ids
    )

    if len(all_debug_boxes) > 0:
        boxes_debug = append_polydata_list(
            all_debug_boxes,
            clean_result=False
        )
    else:
        boxes_debug = vtk.vtkPolyData()

    print("")
    print("Remaining Quality > threshold cells:", len(remaining_bad_ids))

    return (
        working,
        first_bad_debug,
        boxes_debug,
        remaining_bad_debug,
        len(remaining_bad_ids)
    )


# ============================================================
# MAIN
# ============================================================

print("")
print("============================================")
print("Reading inputs")
print("============================================")

full_centerline = read_polydata(
    FULL_CENTERLINE_FILE
)

trimmed_centerline = read_polydata(
    TRIMMED_CENTERLINE_FILE
)

surface_original = read_polydata(
    SURFACE_FILE
)

print("Full centerline cells:", full_centerline.GetNumberOfCells())
print("Trimmed centerline cells:", trimmed_centerline.GetNumberOfCells())
print("Original surface points:", surface_original.GetNumberOfPoints())
print("Original surface cells:", surface_original.GetNumberOfCells())


# ============================================================
# OPTIONAL GENERATION INFORMATION
# ============================================================

generation_array = (
    trimmed_centerline
    .GetPointData()
    .GetArray("Generation")
)

if generation_array is not None:

    max_generation = None

    for i in range(generation_array.GetNumberOfTuples()):

        value = int(
            round(
                generation_array.GetTuple1(i)
            )
        )

        if max_generation is None or value > max_generation:
            max_generation = value

    print(
        "Maximum Generation in trimmed centerline:",
        max_generation
    )


# ============================================================
# RETAINED / REMOVED CENTERLINE CELLS
# ============================================================

retained_cell_ids = get_retained_cell_ids(
    trimmed_centerline,
    full_centerline
)

all_cell_ids = set(
    range(full_centerline.GetNumberOfCells())
)

removed_cell_ids = (
    all_cell_ids - retained_cell_ids
)

print("")
print(
    "Retained full-tree centerline cells:",
    len(retained_cell_ids)
)
print(
    "Removed full-tree centerline cells:",
    len(removed_cell_ids)
)


# ============================================================
# CUT DEFINITIONS
# ============================================================

# Stage 1: cuts inferred from branches removed between the full and trimmed
# centerline trees. These include side-branch pruning cuts and synthetic
# pruning-terminal cuts.
pruning_cuts, boundary_junction_count = build_cut_definitions(
    full_centerline,
    retained_cell_ids
)

side_cut_count = sum(
    1 for cut in pruning_cuts
    if cut["type"] == "side_branch"
)

terminal_cut_count = sum(
    1 for cut in pruning_cuts
    if cut["type"] == "terminal"
)

# Stage 2: true IsOutlet cuts from the trimmed centerline. These are applied
# AFTER the pruning stage has removed disconnected distal branches.
outlet_cuts, outlet_topology_skips = build_outlet_cut_definitions(
    trimmed_centerline
)

# Append the separate inlet cut LAST so all historical IsOutlet indices stay
# unchanged; the stage-3 inlet operation uses the next planned cut index.
if ENABLE_INLET_CUT:
    inlet_cut_definition, inlet_path = build_inlet_cut_definition(
        trimmed_centerline
    )
else:
    inlet_cut_definition, inlet_path = None, None

all_planned_cuts = pruning_cuts + outlet_cuts
if inlet_cut_definition is not None:
    all_planned_cuts.append(inlet_cut_definition)

print("")
print("============================================")
print("Cut planning summary")
print("============================================")
print("Pruning boundary junctions:", boundary_junction_count)
print("Side-branch pruning cuts:", side_cut_count)
print("Synthetic pruning-terminal cuts:", terminal_cut_count)
print("Trimmed-centerline IsOutlet cuts planned:", len(outlet_cuts))
print("IsOutlet points skipped for centerline topology/radius:", outlet_topology_skips)
if inlet_cut_definition is not None:
    print("Inlet cut planned at: {:.0f}% inlet-to-first-bifurcation arclength".format(
        100 * INLET_CUT_PATH_FRACTION
    ))
    print("Inlet trunk length:", inlet_path["total_length"])
    print("First bifurcation location:", inlet_path["branch_xyz"])
print("Total planned surface cuts:", len(all_planned_cuts))

if len(all_planned_cuts) == 0:
    raise RuntimeError(
        "No pruning-boundary cuts or IsOutlet cuts were found."
    )

if VERBOSE_CUT_LOG:

    for i, cut in enumerate(all_planned_cuts):

        print("")
        print("Planned cut {} / {}".format(i + 1, len(all_planned_cuts)))
        print("    Type:", cut["type"])
        print("    Center:", cut["center"])
        print("    Normal:", cut["normal"])
        print("    Centerline radius:", cut["radius"])
        print("    Expected QC radius:", cut["cut_radius"])
        print("    Source centerline cell:", cut["source_cell_id"])

        if cut["type"] in ("outlet", "inlet"):
            print("    Source terminal point ID:", cut["source_point_id"])
        else:
            print("    Removed cell(s):", cut["removed_cell_ids"])

write_cut_debug_file(
    all_planned_cuts,
    CUTS_DEBUG_FILE
)

print("")
print("Intended cut-circle debug file:")
print(CUTS_DEBUG_FILE)


# ============================================================
# FIND THE RETAINED CENTERLINE INLET
# ============================================================

is_inlet = (
    trimmed_centerline
    .GetPointData()
    .GetArray("IsInlet")
)

if is_inlet is None:
    raise RuntimeError(
        "Trimmed centerline is missing PointData array 'IsInlet'."
    )

inlet_ids = []

for i in range(trimmed_centerline.GetNumberOfPoints()):

    if int(round(is_inlet.GetTuple1(i))) == 1:
        inlet_ids.append(i)

if len(inlet_ids) != 1:
    raise RuntimeError(
        "Expected exactly one IsInlet point in the trimmed "
        "centerline, but found {}."
        .format(len(inlet_ids))
    )

inlet_xyz = point_tuple(
    trimmed_centerline,
    inlet_ids[0]
)

print("")
print("Trimmed-tree inlet coordinate:", inlet_xyz)


# ============================================================
# TRIANGULATE ORIGINAL SURFACE
# ============================================================

working_surface = triangulate_polydata(
    surface_original
)

working_surface = clean_polydata(
    working_surface,
    MERGE_TOLERANCE
)

print("")
print("Triangulated surface cells:", working_surface.GetNumberOfCells())

# Fail early if the inlet cross-section is already geometrically invalid.
# Crucially, the previous version did not attempt inlet QC until AFTER every
# pruning and outlet cut. Also retain *all* rejected inlet loops for ParaView.
inlet_rejected_debug_records = []
if inlet_cut_definition is not None and ENABLE_EARLY_INLET_PREFLIGHT:
    prevalidated_inlet_candidate = preflight_inlet_cut(
        working_surface, inlet_cut_definition, inlet_rejected_debug_records
    )
    if prevalidated_inlet_candidate is not None:
        prioritize_prevalidated_inlet_cut(
            inlet_cut_definition, prevalidated_inlet_candidate
        )
    else:
        write_inlet_rejection_debug(inlet_rejected_debug_records)
        if INLET_PREFLIGHT_FAIL_FAST and REQUIRE_SUCCESSFUL_INLET_CUT:
            raise RuntimeError(
                "Early inlet preflight rejected all nine candidate planes. "
                "Pruning was not started. Inspect the rejected-contour debug "
                "VTP and the point intersection diagnostics, if generated."
            )


# ============================================================
# APPLY EXACT LOCAL PLANE CUTS
#
# Stage 1: apply all full-vs-trimmed pruning cuts. These are required and
# remain fatal if they fail.
#
# Stage 2: first discard the distal pieces disconnected by Stage 1, then
# process every true IsOutlet point from the trimmed centerline. If an
# IsOutlet no longer has a LOCAL cross-section, it is skipped because the
# earlier pruning-terminal logic already removed that branch.
# ============================================================

selected_loops = []
successful_cut_indices = []
skipped_pruning_cut_indices = []
skipped_outlet_cut_indices = []
failed_cut_attempts = []
accepted_cut_by_index = {}
exact_cuts_since_global_clean = 0

print("")
print("============================================")
print("STAGE 1: Applying pruning-boundary VTK cuts")
print("============================================")

for cut_index, cut in enumerate(pruning_cuts):

    cut["planned_cut_index"] = cut_index

    if VERBOSE_CUT_LOG:
        print("")
        print(
            "Pruning cut {}/{} ({})"
            .format(
                cut_index + 1,
                len(pruning_cuts),
                cut["type"]
            )
        )
    elif (cut_index + 1) % PROGRESS_PRINT_EVERY == 0 or cut_index == 0:
        print(
            "Pruning progress: {}/{}".format(
                cut_index + 1, len(pruning_cuts)
            )
        )

    try:

        (
            working_surface,
            loop,
            _cap,
            stats,
            used_cut
        ) = apply_cut_with_retries(
            working_surface,
            cut,
            cleanup_result=False
        )

        selected_loops.append(
            annotate_loop(
                loop,
                cut_index,
                used_cut
            )
        )

        successful_cut_indices.append(
            cut_index
        )
        accepted_cut_by_index[cut_index] = dict(used_cut)

        exact_cuts_since_global_clean += 1
        if (
            EXACT_CUT_GLOBAL_CLEAN_INTERVAL > 0
            and exact_cuts_since_global_clean >= EXACT_CUT_GLOBAL_CLEAN_INTERVAL
        ):
            working_surface = clean_polydata(working_surface, MERGE_TOLERANCE)
            working_surface = triangulate_polydata(working_surface)
            exact_cuts_since_global_clean = 0

        if VERBOSE_CUT_LOG:
            print("    Target loop points:", stats["num_loop_points"])
            print("    Source triangles intersected:", stats["num_crossing_cells"])
            print(
                "    Downstream sacrificial one-ring triangles removed:",
                stats["num_downstream_neighbor_cells_removed"]
            )
            print("    Loop max radial distance:", stats["max_radial"])
            print("    Loop area:", stats.get("loop_area"))
            print(
                "    Area / expected local area:",
                stats.get("area_ratio")
            )
            print("    Loop circularity:", stats.get("circularity"))
            print(
                "    Loop radial asymmetry:",
                stats.get("radial_asymmetry")
            )
            print(
                "    Loop centroid offset from centerline:",
                stats["centroid_offset"]
            )
            if stats.get("retry_attempt", 0) > 0:
                print(
                    "    Accepted retry attempt:",
                    stats["retry_attempt"] + 1
                )
                print(
                    "    Accepted distance from anchor endpoint:",
                    stats["retry_distance_from_anchor"]
                )
            print(
                "    Working surface cells after cut:",
                working_surface.GetNumberOfCells()
            )

    except CutRetriesExhausted as exc:

        skipped_pruning_cut_indices.append(cut_index)
        failed_cut_attempts.extend(exc.attempts)

        if VERBOSE_CUT_LOG:
            print(
                "    SKIPPED ENTIRELY after {} failed plane location(s)."
                .format(len(exc.attempts))
            )

        if not SKIP_FAILED_PRUNING_CUTS:
            raise

    except Exception as exc:

        # Extra permissive safety net for this refinement phase. Missing
        # input arrays/top-level setup errors still fail before this loop,
        # but a single problematic anatomical cut does not kill the batch.
        failed_cut = dict(cut)
        failed_cut["attempt_number"] = 1
        failed_cut_attempts.append(failed_cut)
        skipped_pruning_cut_indices.append(cut_index)

        if VERBOSE_CUT_LOG:
            print("    SKIPPED ENTIRELY because of unexpected cut error:", exc)

        if not SKIP_FAILED_PRUNING_CUTS:
            raise


# ============================================================
# REMOVE DISTAL PIECES CREATED BY THE PRUNING CUTS
#
# This intermediate connectivity pass is deliberate. It ensures that an
# IsOutlet lying on a branch already removed in Stage 1 no longer has a local
# surface cross-section available when Stage 2 tries to process it.
# ============================================================

if len(pruning_cuts) > 0:

    # A batched exact cut may leave coincident seam points unwelded until the
    # next global clean. Weld once BEFORE connectivity so the retained proximal
    # cut fragments are recognized as part of the inlet-connected surface.
    if exact_cuts_since_global_clean > 0:
        working_surface = clean_polydata(working_surface, MERGE_TOLERANCE)
        working_surface = triangulate_polydata(working_surface)
        exact_cuts_since_global_clean = 0

    working_surface = keep_inlet_connected_region(
        working_surface,
        inlet_xyz
    )

    working_surface = triangulate_polydata(
        working_surface
    )
    exact_cuts_since_global_clean = 0

    print("")
    print(
        "Surface cells after Stage-1 inlet-connectivity cleanup:",
        working_surface.GetNumberOfCells()
    )


# ============================================================
# STAGE 2: TRUE IsOutlet CUTS
# ============================================================

print("")
print("============================================")
print("STAGE 2: Applying trimmed IsOutlet VTK cuts")
print("============================================")

outlet_success_count = 0
outlet_skip_count = 0

for outlet_local_index, cut in enumerate(outlet_cuts):

    cut_index = (
        len(pruning_cuts)
        + outlet_local_index
    )

    cut["planned_cut_index"] = cut_index

    if VERBOSE_CUT_LOG:
        print("")
        print(
            "IsOutlet cut {}/{} (point {})"
            .format(
                outlet_local_index + 1,
                len(outlet_cuts),
                cut["source_point_id"]
            )
        )
    elif (outlet_local_index + 1) % PROGRESS_PRINT_EVERY == 0 or outlet_local_index == 0:
        print(
            "IsOutlet progress: {}/{} | completed={} skipped={}"
            .format(
                outlet_local_index + 1,
                len(outlet_cuts),
                outlet_success_count,
                outlet_skip_count
            )
        )

    try:

        (
            working_surface,
            loop,
            _cap,
            stats,
            used_cut
        ) = apply_cut_with_retries(
            working_surface,
            cut,
            cleanup_result=False
        )

        selected_loops.append(
            annotate_loop(
                loop,
                cut_index,
                used_cut
            )
        )

        successful_cut_indices.append(
            cut_index
        )
        accepted_cut_by_index[cut_index] = dict(used_cut)

        outlet_success_count += 1
        exact_cuts_since_global_clean += 1
        if (
            EXACT_CUT_GLOBAL_CLEAN_INTERVAL > 0
            and exact_cuts_since_global_clean >= EXACT_CUT_GLOBAL_CLEAN_INTERVAL
        ):
            working_surface = clean_polydata(working_surface, MERGE_TOLERANCE)
            working_surface = triangulate_polydata(working_surface)
            exact_cuts_since_global_clean = 0

        if VERBOSE_CUT_LOG:
            print("    Target loop points:", stats["num_loop_points"])
            print("    Source triangles intersected:", stats["num_crossing_cells"])
            print(
                "    Downstream sacrificial one-ring triangles removed:",
                stats["num_downstream_neighbor_cells_removed"]
            )
            print("    Loop max radial distance:", stats["max_radial"])
            print("    Loop area:", stats.get("loop_area"))
            print(
                "    Area / expected local area:",
                stats.get("area_ratio")
            )
            print("    Loop circularity:", stats.get("circularity"))
            print(
                "    Loop radial asymmetry:",
                stats.get("radial_asymmetry")
            )
            print(
                "    Loop centroid offset from centerline:",
                stats["centroid_offset"]
            )
            if stats.get("retry_attempt", 0) > 0:
                print(
                    "    Accepted retry attempt:",
                    stats["retry_attempt"] + 1
                )
                print(
                    "    Accepted distance from outlet endpoint:",
                    stats["retry_distance_from_anchor"]
                )
            print(
                "    Working surface cells after cut:",
                working_surface.GetNumberOfCells()
            )

    except CutRetriesExhausted as exc:

        outlet_skip_count += 1
        skipped_outlet_cut_indices.append(cut_index)
        failed_cut_attempts.extend(exc.attempts)

        if VERBOSE_CUT_LOG:
            print(
                "    SKIPPED ENTIRELY after {} failed plane location(s)."
                .format(len(exc.attempts))
            )

        if not SKIP_FAILED_OUTLET_CUTS:
            raise

    except NoCutLoopFound as exc:

        outlet_skip_count += 1
        skipped_outlet_cut_indices.append(
            cut_index
        )

        if VERBOSE_CUT_LOG:
            print("    SKIPPED: no local cross-section was found:", exc)

    except Exception as exc:

        outlet_skip_count += 1
        skipped_outlet_cut_indices.append(cut_index)

        failed_cut = dict(cut)
        failed_cut["attempt_number"] = 1
        failed_cut_attempts.append(failed_cut)

        if VERBOSE_CUT_LOG:
            print(
                "    SKIPPED ENTIRELY because of unexpected cut error:",
                exc
            )

        if not SKIP_FAILED_OUTLET_CUTS:
            raise

print("")
# Final weld after batched exact cuts before connectivity/debug output.
if exact_cuts_since_global_clean > 0:
    working_surface = clean_polydata(working_surface, MERGE_TOLERANCE)
    working_surface = triangulate_polydata(working_surface)
    exact_cuts_since_global_clean = 0

print("Pruning cuts skipped:", len(skipped_pruning_cut_indices))
print("IsOutlet cuts completed:", outlet_success_count)
print("IsOutlet cuts skipped:", outlet_skip_count)

if VERBOSE_CUT_LOG:
    if len(skipped_pruning_cut_indices) > 0:
        print("Skipped pruning cut indices:", skipped_pruning_cut_indices)

    if len(skipped_outlet_cut_indices) > 0:
        print("Skipped IsOutlet cut indices:", skipped_outlet_cut_indices)


# ============================================================
# STAGE 3: CUT THE INLET HALFWAY TO THE FIRST BIFURCATION
# ============================================================
# This operation runs AFTER the existing outlet cuts. Recompute the inlet
# connectivity anchor on the RETAINED, branchward side before any later
# connectivity or endpoint-box operation.

accepted_inlet_cut = None
if inlet_cut_definition is not None:
    print("")
    print("============================================")
    print("STAGE 3: Cutting the inlet")
    print("============================================")
    inlet_cut_index = len(pruning_cuts) + len(outlet_cuts)
    inlet_cut_definition["planned_cut_index"] = inlet_cut_index

    try:
        (
            working_surface, loop, _cap, stats, used_cut
        ) = apply_cut_with_retries(
            working_surface, inlet_cut_definition, cleanup_result=True,
            inlet_debug_collector=inlet_rejected_debug_records,
        )
        selected_loops.append(
            annotate_loop(loop, inlet_cut_index, used_cut)
        )
        successful_cut_indices.append(inlet_cut_index)
        accepted_cut_by_index[inlet_cut_index] = dict(used_cut)
        accepted_inlet_cut = dict(used_cut)

        # The OLD IsInlet point is on the piece we just removed. The new
        # interior anchor remains safely downstream of every inlet retry.
        inlet_xyz = inlet_cut_definition["retained_anchor"]
        exact_cuts_since_global_clean = 0
        print("Inlet cut accepted at:", used_cut["center"])
        print("Inlet cut fraction:", used_cut.get(
            "retry_path_fraction", INLET_CUT_PATH_FRACTION
        ))
        print("New inlet connectivity anchor:", inlet_xyz)
    except Exception as exc:
        write_inlet_rejection_debug(inlet_rejected_debug_records)
        if REQUIRE_SUCCESSFUL_INLET_CUT:
            raise RuntimeError(
                "Required inlet cut failed after retries. Diagnostic contours "
                "were saved when a cutter loop was available. No capped "
                "model or .facenames file was written."
            ) from exc
        print("WARNING: Inlet cut failed:", exc)
    else:
        if inlet_rejected_debug_records:
            print("\nInlet retry contour diagnostics:")
            write_inlet_rejection_debug(inlet_rejected_debug_records)


# ============================================================
# WRITE FAILED/REJECTED CUT-ATTEMPT DEBUG FILE
# ============================================================

if len(failed_cut_attempts) > 0:

    write_cut_debug_file(
        failed_cut_attempts,
        FAILED_CUT_ATTEMPTS_DEBUG_FILE
    )

    print("")
    print("Failed/rejected cut-attempt debug file:")
    print(FAILED_CUT_ATTEMPTS_DEBUG_FILE)
    print(
        "Color by CutType, AttemptNumber, or PlannedCutIndex to "
        "inspect exactly where skipped cuts were attempted."
    )

else:

    print("")
    print("No failed cut attempts were recorded.")


# ============================================================
# WRITE ACTUAL CUT-LOOP DEBUG FILE
# ============================================================

loops_debug = append_polydata_list(
    selected_loops,
    clean_result=False
)

write_polydata(
    loops_debug,
    CUT_LOOPS_DEBUG_FILE
)

print("")
print("Actual exact cut-loop debug file:")
print(CUT_LOOPS_DEBUG_FILE)


# ============================================================
# WRITE PRE-FINAL-CONNECTIVITY EXACT-CUT SURFACE
#
# Stage-1 pruning distal pieces have already been removed before the IsOutlet
# stage. This debug surface can still contain the short distal stumps created
# by the new IsOutlet cuts, which the final connectivity pass removes below.
# The cut boundaries use interpolated plane-intersection points rather than
# whole missing triangles.
# ============================================================

write_polydata(
    working_surface,
    SEVERED_DEBUG_FILE
)

print("")
print("Exact severed debug surface:")
print(SEVERED_DEBUG_FILE)


# ============================================================
# COUNT CONNECTED REGIONS AFTER EXACT CUTTING
# ============================================================

connectivity_all = vtk.vtkPolyDataConnectivityFilter()
connectivity_all.SetInputData(working_surface)
connectivity_all.SetExtractionModeToAllRegions()
connectivity_all.ColorRegionsOn()
connectivity_all.Update()

print("")
print(
    "Connected regions after exact severing:",
    connectivity_all.GetNumberOfExtractedRegions()
)


# ============================================================
# KEEP ONLY THE PA TREE CONNECTED TO THE INLET
# ============================================================

final_open = keep_inlet_connected_region(
    working_surface,
    inlet_xyz
)

final_open = triangulate_polydata(
    final_open
)

if REGENERATE_GLOBAL_IDS:
    regenerate_global_ids(final_open)

write_polydata(
    final_open,
    OUTPUT_OPEN_FILE
)

print("")
print("============================================")
print("EXACT OPEN OUTPUT")
print("============================================")
print("Points:", final_open.GetNumberOfPoints())
print("Cells:", final_open.GetNumberOfCells())
print(OUTPUT_OPEN_FILE)


# ============================================================
# FINAL END-BOX SAFETY SWEEP
# ============================================================

optional_terminal_end_cuts = build_optional_terminal_end_box_cuts(
    trimmed_centerline
)

end_box_open = final_open
end_boxes_debug = vtk.vtkPolyData()
end_box_applied_count = 0
end_box_locality_skips = 0
end_box_empty_skips = 0

if ENABLE_FINAL_END_BOX_SWEEP:

    print("")
    print("============================================")
    print("FINAL END-BOX SAFETY SWEEP")
    print("============================================")
    print(
        "Endpoint box half-widths: radial={:.2f} R, axial={:.2f} R"
        .format(
            END_BOX_RADIAL_HALF_WIDTH_FACTOR,
            END_BOX_AXIAL_HALF_LENGTH_FACTOR
        )
    )

    (
        end_box_open,
        end_boxes_debug,
        end_box_applied_count,
        end_box_locality_skips,
        end_box_empty_skips
    ) = run_final_end_box_sweep(
        final_open,
        all_planned_cuts,
        accepted_cut_by_index,
        optional_terminal_end_cuts,
        inlet_xyz
    )

    if REGENERATE_GLOBAL_IDS:
        regenerate_global_ids(end_box_open)

    if end_boxes_debug.GetNumberOfCells() > 0:
        write_polydata(end_boxes_debug, END_BOXES_DEBUG_FILE)

    write_polydata(end_box_open, END_BOX_OPEN_OUTPUT_FILE)

    print("Endpoint boxes applied:", end_box_applied_count)
    print("Endpoint boxes skipped by locality:", end_box_locality_skips)
    print("Endpoint boxes skipped because result was empty:", end_box_empty_skips)
    print("Additional IsTerminalEnd markers:", len(optional_terminal_end_cuts))
    print("Endpoint-box debug file:", END_BOXES_DEBUG_FILE)
    print("Endpoint-box open surface:", END_BOX_OPEN_OUTPUT_FILE)


# ============================================================
# FINAL CAPPING IS DEFERRED UNTIL AFTER ALL CLEANUP
# ============================================================
#
# The older version capped the exact-cut loops here.  That is intentionally
# disabled: the final endpoint-box sweep and mesh-quality cleanup can move an
# opening after those loops were recorded.  The final capping pass below
# re-detects the ACTUAL boundary loops of the finished open surface, creates
# better-quality planar caps, stamps unique ModelFaceID values, and writes the
# matching SimVascular .facenames file.

# ============================================================
# FINAL MESH-QUALITY / SPIKE CLEANUP
# ============================================================

# The endpoint-box sweep is now the direct input to mesh-quality cleanup.
# This variable was accidentally left undefined in the previous revision.
quality_cleanup_input = end_box_open

if ENABLE_MESH_QUALITY_CLEANUP:

    (
        final_quality_cleaned,
        bad_quality_debug,
        quality_boxes_debug,
        remaining_bad_quality_debug,
        remaining_bad_quality_count
    ) = run_mesh_quality_spike_cleanup(
        quality_cleanup_input,
        trimmed_centerline,
        inlet_xyz
    )

    # Recompute normals after box clipping.  The quality-cleaned result is
    # intentionally left OPEN where a box severed a bad terminal/spike region.
    quality_normals = vtk.vtkPolyDataNormals()
    quality_normals.SetInputData(final_quality_cleaned)
    quality_normals.ConsistencyOn()
    quality_normals.SplittingOff()
    quality_normals.ComputePointNormalsOn()
    quality_normals.ComputeCellNormalsOn()
    quality_normals.Update()

    final_quality_cleaned_with_normals = vtk.vtkPolyData()
    final_quality_cleaned_with_normals.DeepCopy(
        quality_normals.GetOutput()
    )

    if REGENERATE_GLOBAL_IDS:
        regenerate_global_ids(
            final_quality_cleaned_with_normals
        )

    write_polydata(
        final_quality_cleaned_with_normals,
        OUTPUT_QUALITY_CLEANED_FILE
    )

    if (
        bad_quality_debug is not None
        and bad_quality_debug.GetNumberOfCells() > 0
    ):
        write_polydata(
            bad_quality_debug,
            BAD_QUALITY_CELLS_DEBUG_FILE
        )
        print("")
        print("Initial bad-quality cells debug file:")
        print(BAD_QUALITY_CELLS_DEBUG_FILE)

    if quality_boxes_debug.GetNumberOfCells() > 0:
        write_polydata(
            quality_boxes_debug,
            QUALITY_BOXES_DEBUG_FILE
        )
        print("")
        print("Quality cleanup boxes debug file:")
        print(QUALITY_BOXES_DEBUG_FILE)

    if remaining_bad_quality_debug.GetNumberOfCells() > 0:
        write_polydata(
            remaining_bad_quality_debug,
            REMAINING_BAD_QUALITY_CELLS_DEBUG_FILE
        )
        print("")
        print("Remaining bad-quality cells debug file:")
        print(REMAINING_BAD_QUALITY_CELLS_DEBUG_FILE)

    print("")
    print("============================================")
    print("QUALITY-CLEANED OPEN OUTPUT")
    print("============================================")
    print(
        "Points:",
        final_quality_cleaned_with_normals.GetNumberOfPoints()
    )
    print(
        "Cells:",
        final_quality_cleaned_with_normals.GetNumberOfCells()
    )
    print(
        "Remaining cells above quality threshold:",
        remaining_bad_quality_count
    )
    print(OUTPUT_QUALITY_CLEANED_FILE)



# ============================================================
# FINAL CAP + MODELFACEID + .FACENAMES PASS
# ============================================================
# Cap the ACTUAL final open model, after endpoint-box and quality cleanup.
# This avoids using stale pre-cleanup cap loops and gives every cap a unique
# SimVascular ModelFaceID plus a matching .facenames entry.

if ENABLE_MESH_QUALITY_CLEANUP:
    final_open_for_capping = final_quality_cleaned_with_normals
else:
    final_open_for_capping = end_box_open

if WRITE_CAPPED_OUTPUT:
    print("")
    print("============================================")
    print("FINAL CAPS + SIMVASCULAR FACE IDS")
    print("============================================")

    (
        final_capped_with_normals,
        final_face_assignments,
        final_cap_loops_debug,
    ) = build_final_capped_face_model(
        final_open_for_capping,
        accepted_inlet_cut,
        accepted_cut_by_index,
        trimmed_centerline,
    )

    write_polydata(
        final_capped_with_normals,
        OUTPUT_CAPPED_FILE
    )

    write_facenames_file(
        OUTPUT_CAPPED_FACENAMES,
        final_face_assignments
    )

    if final_cap_loops_debug.GetNumberOfCells() > 0:
        write_polydata(
            final_cap_loops_debug,
            FINAL_CAP_LOOPS_DEBUG_FILE
        )

    print("Final caps created:", len(final_face_assignments) - 1)
    print("Wall ModelFaceID:", FINAL_CAP_WALL_FACE_ID)
    print("Boundary edges after capping:", boundary_edge_count(final_capped_with_normals))
    print("Capped face model:", OUTPUT_CAPPED_FILE)
    print("Facenames file:", OUTPUT_CAPPED_FACENAMES)
    print("Cap-loop debug:", FINAL_CAP_LOOPS_DEBUG_FILE)


print("")
print("============================================")
print("CUT SUMMARY")
print("============================================")
print("Pruning cuts completed:", len(pruning_cuts))
print("IsOutlet cuts planned:", len(outlet_cuts))
print("IsOutlet cuts completed:", outlet_success_count)
print("IsOutlet cuts skipped:", outlet_skip_count)
print("Inlet cut completed:", accepted_inlet_cut is not None)
print("Final endpoint boxes applied:", end_box_applied_count)
print("Final endpoint-box locality skips:", end_box_locality_skips)

print("")
print("============================================")
print("FINISHED")
print("============================================")
print("")
if ENABLE_MESH_QUALITY_CLEANUP:
    print("Primary final OPEN output:", OUTPUT_QUALITY_CLEANED_FILE)
else:
    print("Primary final OPEN output:", END_BOX_OPEN_OUTPUT_FILE)
if WRITE_CAPPED_OUTPUT:
    print("Primary final CAPPED face model:", OUTPUT_CAPPED_FILE)
    print("SimVascular facenames:", OUTPUT_CAPPED_FACENAMES)
print("Failed-cut debug output:", FAILED_CUT_ATTEMPTS_DEBUG_FILE)
