import vtk
import os
import math
from collections import defaultdict


# ============================================================
# USER SETTINGS
# ============================================================

# Full, untrimmed centerline tree.
FULL_CENTERLINE_FILE = r"C:\path\to\Network_model_V2.vtp"

# Centerline tree after trimming to the extent you want to keep.
TRIMMED_CENTERLINE_FILE = r"C:\path\to\JL_Centerline_mini-cut.vtp"

# Original pulmonary artery surface to prune.
SURFACE_FILE = r"C:\path\to\Pulmonary_Artery_Surface.vtp"

# Exact-clipped PA tree with the new outlet rings left OPEN.
# This is the most useful QC output because you can inspect the new
# plane-cut boundary before any cap is added.
OUTPUT_OPEN_FILE = r"C:\path\to\Pulmonary_Artery_Surface_Pruned_Exact_Open.vtp"

# Exact-clipped PA tree with flat triangulated caps added at the new cuts.
OUTPUT_CAPPED_FILE = r"C:\path\to\Pulmonary_Artery_Surface_Pruned_Exact_Capped.vtp"

# Debug output: the intended circular search/QC disk for every cut.
CUTS_DEBUG_FILE = r"C:\path\to\Pulmonary_Artery_Cuts_Debug.vtp"

# Debug output: the ACTUAL plane/surface intersection loop selected for
# every cut. These should look like clean oval/circular outlet rings.
CUT_LOOPS_DEBUG_FILE = r"C:\path\to\Pulmonary_Artery_Exact_Cut_Loops_Debug.vtp"

# Debug output: after all exact cuts, but BEFORE connectivity removes the
# detached distal vessel pieces.
SEVERED_DEBUG_FILE = r"C:\path\to\Pulmonary_Artery_Exact_Severed_Debug.vtp"


# ============================================================
# CENTERLINE / CUT SETTINGS
# ============================================================

# In the supplied trimmed centerline file, CellIds map back to the
# original full centerline cells.
CELL_ID_ARRAY = "CellIds"

# Radius array on the full centerline.
RADIUS_ARRAY = "Radius"

# Used as a QC radius for deciding whether the plane intersection loop
# selected by VTK is plausibly the vessel surrounding this centerline.
# The exact cut itself is made from the true plane/surface intersection;
# it is NOT approximated by this circle.
CUT_RADIUS_FACTOR = 1.35

# When a side branch disappears from the trimmed centerline, move the
# plane slightly INTO the removed branch so we are not cutting directly
# through the bifurcation crotch.
#
# 0.0 = exactly at the branch junction
# 0.35 = 35% of the way from the junction to the next point in the
#        removed centerline branch
SIDE_BRANCH_OFFSET_FRACTION = 0.35

# If the retained centerline terminates at a junction, make one transverse
# cut using the last retained branch direction. Positive values move the
# plane slightly upstream into the retained parent branch.
TERMINAL_OFFSET_FRACTION = 0.0

# Coordinate tolerance used to group coincident centerline endpoints that
# use different VTK point IDs at the same anatomical bifurcation.
COORDINATE_TOLERANCE = 1e-5

# Small geometric tolerance.
PLANE_EPS = 1e-10

# Debug circle resolution only. This does NOT control the actual cut.
DEBUG_CIRCLE_RESOLUTION = 64

# If the selected plane/surface intersection loop extends much farther
# from the centerline than expected, print a warning. The cut still
# proceeds unless FAIL_ON_SUSPICIOUS_LOOP is True.
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

            cuts.append(
                {
                    "type": "terminal",
                    "center": cut_center,
                    "junction": junction,
                    "normal": normal,
                    "radius": radius,
                    "cut_radius": radius * CUT_RADIUS_FACTOR,
                    "source_cell_id": kept["cell_id"],
                    "removed_cell_ids": [
                        e["cell_id"]
                        for e in removed_entries
                    ]
                }
            )

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

                cuts.append(
                    {
                        "type": "side_branch",
                        "center": cut_center,
                        "junction": junction,
                        "normal": normal,
                        "radius": radius,
                        "cut_radius": radius * CUT_RADIUS_FACTOR,
                        "source_cell_id": removed["cell_id"],
                        "removed_cell_ids": [
                            removed["cell_id"]
                        ]
                    }
                )

    return cuts, boundary_junction_count


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
        else:
            cut_type.InsertNextValue(2)

        source_cell.InsertNextValue(
            int(cut["source_cell_id"])
        )

        cut_radius.InsertNextValue(
            float(cut["cut_radius"])
        )

    debug = vtk.vtkPolyData()
    debug.SetPoints(points)
    debug.SetLines(lines)

    debug.GetCellData().AddArray(cut_type)
    debug.GetCellData().AddArray(source_cell)
    debug.GetCellData().AddArray(cut_radius)

    write_polydata(debug, filename)


# ============================================================
# FIND THE ACTUAL TARGET PLANE/SURFACE INTERSECTION LOOP
#
# The plane may intersect several disconnected PA branches. vtkCutter
# creates ALL plane/surface intersections, then connectivity keeps ONLY
# the contour nearest the centerline cut point.
# ============================================================

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
        raise RuntimeError(
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
        raise RuntimeError(
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

        if FAIL_ON_SUSPICIOUS_LOOP:
            raise RuntimeError(message)
        else:
            print("    WARNING:", message)

    stats = {
        "num_loop_points": loop.GetNumberOfPoints(),
        "num_loop_cells": loop.GetNumberOfCells(),
        "min_radial": min_radial,
        "max_radial": max_radial,
        "centroid": centroid,
        "centroid_offset": centroid_offset
    }

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

    flagged = vtk.vtkPolyData()
    flagged.DeepCopy(surface)

    arr = vtk.vtkIntArray()
    arr.SetName("__TargetCutCell")
    arr.SetNumberOfComponents(1)
    arr.SetNumberOfTuples(flagged.GetNumberOfCells())

    for cell_id in range(flagged.GetNumberOfCells()):

        is_selected = cell_id in selected_ids

        arr.SetTuple1(
            cell_id,
            1 if is_selected else 0
        )

    flagged.GetCellData().AddArray(arr)

    threshold = vtk.vtkThreshold()
    threshold.SetInputData(flagged)
    threshold.SetInputArrayToProcess(
        0,
        0,
        0,
        vtk.vtkDataObject.FIELD_ASSOCIATION_CELLS,
        "__TargetCutCell"
    )

    target_value = 1.0 if want_selected else 0.0

    threshold.SetLowerThreshold(target_value - 0.1)
    threshold.SetUpperThreshold(target_value + 0.1)
    threshold.SetThresholdFunction(
        vtk.vtkThreshold.THRESHOLD_BETWEEN
    )
    threshold.Update()

    geometry = vtk.vtkGeometryFilter()
    geometry.SetInputConnection(
        threshold.GetOutputPort()
    )
    geometry.Update()

    out = vtk.vtkPolyData()
    out.DeepCopy(geometry.GetOutput())

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

def apply_exact_local_plane_cut(surface, cut):

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
    result.DeepCopy(append.GetOutput())

    result = clean_polydata(
        result,
        MERGE_TOLERANCE
    )

    result = triangulate_polydata(result)

    # Remove the internal selection flag before the next cut iteration.
    result.GetCellData().RemoveArray("__TargetCutCell")

    cap = build_flat_cap(loop)

    return result, loop, cap, stats


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

    type_value = 1 if cut["type"] == "side_branch" else 2

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

cuts, boundary_junction_count = build_cut_definitions(
    full_centerline,
    retained_cell_ids
)

side_cut_count = sum(
    1 for cut in cuts
    if cut["type"] == "side_branch"
)

terminal_cut_count = sum(
    1 for cut in cuts
    if cut["type"] == "terminal"
)

print("")
print("============================================")
print("Pruning boundary summary")
print("============================================")
print("Boundary junctions:", boundary_junction_count)
print("Side-branch cuts:", side_cut_count)
print("Terminal cuts:", terminal_cut_count)
print("Total surface cuts:", len(cuts))

if len(cuts) == 0:
    raise RuntimeError(
        "No retained/removed centerline boundaries were found."
    )

for i, cut in enumerate(cuts):

    print("")
    print("Cut {} / {}".format(i + 1, len(cuts)))
    print("    Type:", cut["type"])
    print("    Center:", cut["center"])
    print("    Normal:", cut["normal"])
    print("    Centerline radius:", cut["radius"])
    print("    Expected QC radius:", cut["cut_radius"])
    print("    Source full-tree cell:", cut["source_cell_id"])
    print("    Removed cell(s):", cut["removed_cell_ids"])

write_cut_debug_file(
    cuts,
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


# ============================================================
# APPLY EACH EXACT LOCAL PLANE CUT
# ============================================================

selected_loops = []
flat_caps = []
successful_cut_indices = []

print("")
print("============================================")
print("Applying exact VTK plane cuts")
print("============================================")

for cut_index, cut in enumerate(cuts):

    print("")
    print(
        "Exact cut {}/{} ({})"
        .format(
            cut_index + 1,
            len(cuts),
            cut["type"]
        )
    )

    try:

        (
            working_surface,
            loop,
            cap,
            stats
        ) = apply_exact_local_plane_cut(
            working_surface,
            cut
        )

        selected_loops.append(
            annotate_loop(
                loop,
                cut_index,
                cut
            )
        )

        flat_caps.append(
            (cut_index, cut, cap)
        )

        successful_cut_indices.append(
            cut_index
        )

        print(
            "    Target loop points:",
            stats["num_loop_points"]
        )
        print(
            "    Source triangles intersected:",
            stats["num_crossing_cells"]
        )
        print(
            "    Downstream sacrificial one-ring triangles removed:",
            stats["num_downstream_neighbor_cells_removed"]
        )
        print(
            "    Loop max radial distance:",
            stats["max_radial"]
        )
        print(
            "    Loop centroid offset from centerline:",
            stats["centroid_offset"]
        )
        print(
            "    Working surface cells after cut:",
            working_surface.GetNumberOfCells()
        )

    except Exception as exc:

        # For this refinement stage, stop rather than silently producing a
        # partially pruned model. The debug cut circles tell you which plane
        # needs adjustment.
        raise RuntimeError(
            "Exact cut {} ({}) failed:\n{}"
            .format(
                cut_index,
                cut["type"],
                exc
            )
        )


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
# WRITE PRE-CONNECTIVITY EXACT-CUT SURFACE
#
# This still contains the detached distal vessel pieces. The difference
# from the old script is that the cut boundaries now contain interpolated
# plane-intersection points rather than whole missing triangles.
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
# ADD FLAT CAPS
#
# The caps are created from the exact vtkCutter loops. Therefore their
# boundary uses the same interpolated plane-intersection geometry as the
# refined open surface.
# ============================================================

if WRITE_CAPPED_OUTPUT:

    next_face_id = get_next_model_face_id(
        final_open
    )

    compatible_caps = []

    for local_cap_index, item in enumerate(flat_caps):

        cut_index, cut, cap = item

        if next_face_id is None:
            cap_face_id = None
        else:
            cap_face_id = (
                next_face_id + local_cap_index
            )

        cap_compatible = make_cap_compatible_with_surface(
            cap,
            final_open,
            cap_face_id
        )

        compatible_caps.append(
            cap_compatible
        )

    capped_inputs = [final_open] + compatible_caps

    final_capped = append_polydata_list(
        capped_inputs,
        clean_result=True
    )

    final_capped = triangulate_polydata(
        final_capped
    )

    # A final connectivity pass removes any cap that for some unexpected
    # reason was not attached to the retained inlet-connected PA tree.
    final_capped = keep_inlet_connected_region(
        final_capped,
        inlet_xyz
    )

    final_capped = triangulate_polydata(
        final_capped
    )

    if REGENERATE_GLOBAL_IDS:
        regenerate_global_ids(final_capped)

    # Recompute consistent normals. Splitting is disabled so no geometry is
    # duplicated along the new flat cap boundary.
    normals = vtk.vtkPolyDataNormals()
    normals.SetInputData(final_capped)
    normals.ConsistencyOn()
    normals.SplittingOff()
    normals.ComputePointNormalsOn()
    normals.ComputeCellNormalsOn()
    normals.Update()

    final_capped_with_normals = vtk.vtkPolyData()
    final_capped_with_normals.DeepCopy(
        normals.GetOutput()
    )

    write_polydata(
        final_capped_with_normals,
        OUTPUT_CAPPED_FILE
    )

    print("")
    print("============================================")
    print("EXACT CAPPED OUTPUT")
    print("============================================")
    print(
        "Points:",
        final_capped_with_normals.GetNumberOfPoints()
    )
    print(
        "Cells:",
        final_capped_with_normals.GetNumberOfCells()
    )
    print(OUTPUT_CAPPED_FILE)


print("")
print("============================================")
print("FINISHED")
print("============================================")
print("")
print(
    "Recommended QC order: CUTS_DEBUG_FILE -> "
    "CUT_LOOPS_DEBUG_FILE -> OUTPUT_OPEN_FILE -> OUTPUT_CAPPED_FILE"
)
