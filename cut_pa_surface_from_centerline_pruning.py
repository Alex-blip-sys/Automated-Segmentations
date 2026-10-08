import vtk
import os
import math
from collections import defaultdict


# ============================================================
# USER SETTINGS
# ============================================================

# Full, untrimmed centerline tree.
FULL_CENTERLINE_FILE = r"C:\path\to\Network_model_V2.vtp"

# Centerline tree after trimming to the generation you want to keep.
TRIMMED_CENTERLINE_FILE = r"C:\path\to\JL_Centerline_mini-cut.vtp"

# Original pulmonary artery surface to prune.
SURFACE_FILE = r"C:\path\to\Pulmonary_Artery_Surface.vtp"

# Final pruned surface.
OUTPUT_FILE = r"C:\path\to\Pulmonary_Artery_Surface_Pruned.vtp"

# Debug output: surface after the local cuts but before connectivity
# removes the distal pieces.
SEVERED_DEBUG_FILE = r"C:\path\to\Pulmonary_Artery_Surface_Severed_Debug.vtp"

# Debug output: circles showing every local cut.
CUTS_DEBUG_FILE = r"C:\path\to\Pulmonary_Artery_Cuts_Debug.vtp"


# ============================================================
# SETTINGS
# ============================================================

# The trimmed file supplied for this workflow contains CellIds
# that correspond to cell indices in the full centerline file.
CELL_ID_ARRAY = "CellIds"

# Radius array on the full centerline.
RADIUS_ARRAY = "Radius"

# Radius of each local cutting disk:
#
#       cut radius = centerline radius * CUT_RADIUS_FACTOR
#
# Increase slightly if a cut does not fully sever a vessel.
# Decrease if a nearby vessel is accidentally touched.
CUT_RADIUS_FACTOR = 1.35

# When a whole side branch has disappeared from the trimmed
# centerline, move the cut slightly INTO that removed branch.
#
# 0.0 = exactly at the bifurcation
# 0.35 = 35% of the way from the bifurcation to the next
#        centerline point in the removed branch.
SIDE_BRANCH_OFFSET_FRACTION = 0.35

# If only one retained centerline cell reaches a pruning junction,
# the retained tree terminates there. In that case one transverse
# cut is made using the retained branch direction instead of making
# separate cuts into every removed daughter.
#
# 0.0 places the plane exactly at the terminal junction / IsOutlet
# location. A value such as 0.20 would move it slightly upstream
# into the retained parent branch.
TERMINAL_OFFSET_FRACTION = 0.0

# Coordinate tolerance used only to group coincident centerline
# endpoints belonging to different VTK polyline cells.
COORDINATE_TOLERANCE = 1e-5

# Plane intersection tolerance.
PLANE_EPS = 1e-10

# Debug circle resolution.
DEBUG_CIRCLE_RESOLUTION = 64


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


# ============================================================
# COORDINATE KEY
#
# Full centerline branches are separate vtkPolyLine cells.
# At bifurcations, separate point IDs can occupy the same
# physical coordinate. This key lets us recognize them as the
# same anatomical junction.
# ============================================================

def coordinate_key(point):

    tol = COORDINATE_TOLERANCE

    return (
        int(round(point[0] / tol)),
        int(round(point[1] / tol)),
        int(round(point[2] / tol))
    )


# ============================================================
# CELL ENDPOINT INFORMATION
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
#
# In the supplied trimmed file, CellIds are inherited from the
# original full centerline. This gives us a very clean way to
# distinguish:
#
#       retained full centerline cells
#       removed full centerline cells
#
# without relying on nearest-neighbor geometry.
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
#
# Each dictionary entry represents one physical endpoint
# coordinate. Several vtkPolyLine cells may meet there.
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

            endpoint_map[key].append(
                entry
            )

    return endpoint_map


# ============================================================
# BUILD CUT DEFINITIONS
#
# At every full-tree junction containing BOTH retained and
# removed centerline cells:
#
# CASE 1:
#   >= 2 retained incident cells
#
#   The retained tree continues through the junction and only
#   one or more SIDE branches disappeared. Cut each removed
#   branch slightly downstream of the junction.
#
# CASE 2:
#   exactly 1 retained incident cell
#
#   The retained tree terminates at this junction. Make ONE
#   transverse cut using the final retained branch direction.
#   This avoids leaving a Y-shaped stump when both daughter
#   branches were removed.
# ============================================================

def build_cut_definitions(
    full_centerline,
    retained_cell_ids
):

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

        # Not a pruning boundary.
        if (
            len(kept_entries) == 0
            or len(removed_entries) == 0
        ):
            continue

        boundary_junction_count += 1

        # All entries at this junction occupy the same physical
        # location, so use the first entry to get the coordinate.
        junction = point_tuple(
            full_centerline,
            entries[0]["point_id"]
        )

        # ----------------------------------------------------
        # CASE 1:
        # Retained tree TERMINATES here.
        #
        # One retained branch reaches the junction and all
        # other incident branches were removed.
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

            # Direction from retained parent toward its terminal
            # junction.
            normal = normalize(
                vec_sub(
                    junction,
                    interior
                )
            )

            # Default offset = 0, so this is exactly at the
            # retained terminal point. Positive offset moves
            # upstream into the retained branch.
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
                    "cut_radius": (
                        radius
                        * CUT_RADIUS_FACTOR
                    ),
                    "source_cell_id": kept["cell_id"],
                    "removed_cell_ids": [
                        e["cell_id"]
                        for e in removed_entries
                    ]
                }
            )

        # ----------------------------------------------------
        # CASE 2:
        # Retained tree continues through the junction.
        #
        # Each missing branch gets its own local cut.
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

                # Move slightly into the removed branch so the
                # cutting disk is not centered directly on the
                # bifurcation.
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

                # Use the smaller of the two valid values when
                # possible because this is specifically the
                # daughter branch being removed, not the parent.
                valid_radii = [
                    r for r in (r0, r1)
                    if math.isfinite(r) and r > 0.0
                ]

                if len(valid_radii) == 0:
                    raise RuntimeError(
                        "Invalid radius for removed full "
                        "centerline cell {}."
                        .format(
                            removed["cell_id"]
                        )
                    )

                radius = min(valid_radii)

                cuts.append(
                    {
                        "type": "side_branch",
                        "center": cut_center,
                        "junction": junction,
                        "normal": normal,
                        "radius": radius,
                        "cut_radius": (
                            radius
                            * CUT_RADIUS_FACTOR
                        ),
                        "source_cell_id": removed["cell_id"],
                        "removed_cell_ids": [
                            removed["cell_id"]
                        ]
                    }
                )

    return cuts, boundary_junction_count


# ============================================================
# CUT GEOMETRY TEST
# ============================================================

def distance_from_axis(
    point,
    center,
    axis
):

    rel = vec_sub(
        point,
        center
    )

    axial = dot(
        rel,
        axis
    )

    radial_vector = vec_sub(
        rel,
        vec_scale(axis, axial)
    )

    return magnitude(
        radial_vector
    )


def cell_intersects_local_cut(
    surface,
    cell_id,
    center,
    normal,
    cut_radius
):

    cell = surface.GetCell(cell_id)

    npts = cell.GetNumberOfPoints()

    if npts < 2:
        return False

    pts = []
    signed = []

    for j in range(npts):

        pid = cell.GetPointId(j)

        p = point_tuple(
            surface,
            pid
        )

        pts.append(p)

        signed.append(
            dot(
                vec_sub(p, center),
                normal
            )
        )

    # If every vertex is strictly on the same side of the
    # plane, the cell cannot cross it.
    all_positive = True
    all_negative = True

    for s in signed:

        if s <= PLANE_EPS:
            all_positive = False

        if s >= -PLANE_EPS:
            all_negative = False

    if all_positive or all_negative:
        return False

    # A vertex lying directly on the plane.
    for j in range(npts):

        if abs(signed[j]) <= PLANE_EPS:

            if (
                distance_from_axis(
                    pts[j],
                    center,
                    normal
                )
                <= cut_radius
            ):
                return True

    # Check edge-plane intersections.
    for j in range(npts):

        k = (j + 1) % npts

        s0 = signed[j]
        s1 = signed[k]

        crosses = (
            (
                s0 < -PLANE_EPS
                and s1 > PLANE_EPS
            )
            or
            (
                s0 > PLANE_EPS
                and s1 < -PLANE_EPS
            )
        )

        if not crosses:
            continue

        denom = s0 - s1

        if abs(denom) <= PLANE_EPS:
            continue

        t = s0 / denom

        p0 = pts[j]
        p1 = pts[k]

        intersection = (
            p0[0]
            + t * (p1[0] - p0[0]),

            p0[1]
            + t * (p1[1] - p0[1]),

            p0[2]
            + t * (p1[2] - p0[2])
        )

        if (
            distance_from_axis(
                intersection,
                center,
                normal
            )
            <= cut_radius
        ):
            return True

    return False


# ============================================================
# DEBUG CUT CIRCLES
# ============================================================

def add_debug_circle(
    center,
    normal,
    radius,
    points,
    lines
):

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

    for i in range(
        DEBUG_CIRCLE_RESOLUTION
    ):

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

    for i in range(
        DEBUG_CIRCLE_RESOLUTION
    ):

        polyline.GetPointIds().SetId(
            i,
            ids[i]
        )

    polyline.GetPointIds().SetId(
        DEBUG_CIRCLE_RESOLUTION,
        ids[0]
    )

    lines.InsertNextCell(
        polyline
    )


def write_cut_debug_file(
    cuts,
    filename
):

    points = vtk.vtkPoints()
    lines = vtk.vtkCellArray()

    cut_type = vtk.vtkIntArray()
    cut_type.SetName("CutType")

    source_cell = vtk.vtkIntArray()
    source_cell.SetName("SourceCellId")

    cut_radius = vtk.vtkDoubleArray()
    cut_radius.SetName("CutRadius")

    for cut in cuts:

        add_debug_circle(
            cut["center"],
            cut["normal"],
            cut["cut_radius"],
            points,
            lines
        )

        # 1 = branch-specific side cut
        # 2 = retained-tree terminal cut
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

    debug.GetCellData().AddArray(
        cut_type
    )

    debug.GetCellData().AddArray(
        source_cell
    )

    debug.GetCellData().AddArray(
        cut_radius
    )

    write_polydata(
        debug,
        filename
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


print(
    "Full centerline cells:",
    full_centerline.GetNumberOfCells()
)

print(
    "Trimmed centerline cells:",
    trimmed_centerline.GetNumberOfCells()
)

print(
    "Original surface points:",
    surface_original.GetNumberOfPoints()
)

print(
    "Original surface cells:",
    surface_original.GetNumberOfCells()
)


# Optional informational check.
generation_array = (
    trimmed_centerline
    .GetPointData()
    .GetArray("Generation")
)

if generation_array is not None:

    max_generation = None

    for i in range(
        generation_array.GetNumberOfTuples()
    ):

        value = int(
            round(
                generation_array.GetTuple1(i)
            )
        )

        if (
            max_generation is None
            or value > max_generation
        ):
            max_generation = value

    print(
        "Maximum Generation in trimmed centerline:",
        max_generation
    )


# ============================================================
# FIND RETAINED / REMOVED CENTERLINE CELLS
# ============================================================

retained_cell_ids = get_retained_cell_ids(
    trimmed_centerline,
    full_centerline
)

all_cell_ids = set(
    range(
        full_centerline.GetNumberOfCells()
    )
)

removed_cell_ids = (
    all_cell_ids
    - retained_cell_ids
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
# CREATE CUT DEFINITIONS
# ============================================================

cuts, boundary_junction_count = (
    build_cut_definitions(
        full_centerline,
        retained_cell_ids
    )
)


side_cut_count = sum(
    1
    for cut in cuts
    if cut["type"] == "side_branch"
)

terminal_cut_count = sum(
    1
    for cut in cuts
    if cut["type"] == "terminal"
)


print("")
print("============================================")
print("Pruning boundary summary")
print("============================================")

print(
    "Boundary junctions:",
    boundary_junction_count
)

print(
    "Side-branch cuts:",
    side_cut_count
)

print(
    "Terminal cuts:",
    terminal_cut_count
)

print(
    "Total surface cuts:",
    len(cuts)
)


if len(cuts) == 0:

    raise RuntimeError(
        "No retained/removed centerline boundaries were found."
    )


for i, cut in enumerate(cuts):

    print("")
    print(
        "Cut {} / {}".format(
            i + 1,
            len(cuts)
        )
    )

    print(
        "    Type:",
        cut["type"]
    )

    print(
        "    Center:",
        cut["center"]
    )

    print(
        "    Normal:",
        cut["normal"]
    )

    print(
        "    Centerline radius:",
        cut["radius"]
    )

    print(
        "    Cutting radius:",
        cut["cut_radius"]
    )

    print(
        "    Source full-tree cell:",
        cut["source_cell_id"]
    )

    print(
        "    Removed cell(s):",
        cut["removed_cell_ids"]
    )


# ============================================================
# DEBUG CUT CIRCLES
# ============================================================

write_cut_debug_file(
    cuts,
    CUTS_DEBUG_FILE
)


print("")
print(
    "Cut-circle debug file written to:"
)

print(
    CUTS_DEBUG_FILE
)


# ============================================================
# FIND TRIMMED CENTERLINE INLET
# ============================================================

is_inlet = (
    trimmed_centerline
    .GetPointData()
    .GetArray("IsInlet")
)

if is_inlet is None:

    raise RuntimeError(
        "Trimmed centerline is missing PointData "
        "array 'IsInlet'."
    )


inlet_ids = []

for i in range(
    trimmed_centerline.GetNumberOfPoints()
):

    if int(
        round(
            is_inlet.GetTuple1(i)
        )
    ) == 1:

        inlet_ids.append(i)


if len(inlet_ids) != 1:

    raise RuntimeError(
        "Expected exactly one IsInlet point in the "
        "trimmed centerline, but found {}."
        .format(len(inlet_ids))
    )


inlet_xyz = point_tuple(
    trimmed_centerline,
    inlet_ids[0]
)


print("")
print(
    "Trimmed-tree inlet coordinate:",
    inlet_xyz
)


# ============================================================
# TRIANGULATE SURFACE
# ============================================================

triangle_filter = vtk.vtkTriangleFilter()

triangle_filter.SetInputData(
    surface_original
)

triangle_filter.PassLinesOff()
triangle_filter.PassVertsOff()

triangle_filter.Update()


surface = vtk.vtkPolyData()

surface.DeepCopy(
    triangle_filter.GetOutput()
)


print("")
print(
    "Triangulated surface cells:",
    surface.GetNumberOfCells()
)


# ============================================================
# BUILD A CELL LOCATOR
#
# Instead of testing every surface triangle against every cut,
# use a local bounding box around each cut first.
# ============================================================

cell_locator = vtk.vtkStaticCellLocator()

cell_locator.SetDataSet(
    surface
)

cell_locator.BuildLocator()


# ============================================================
# IDENTIFY SURFACE CELLS TO REMOVE
# ============================================================

cells_to_remove = set()

print("")
print("============================================")
print("Finding local surface cutting rings")
print("============================================")


for i, cut in enumerate(cuts):

    center = cut["center"]
    radius = cut["cut_radius"]

    # A cube enclosing the circular cutting disk.
    search_radius = (
        radius * 1.10
    )

    bounds = (
        center[0] - search_radius,
        center[0] + search_radius,
        center[1] - search_radius,
        center[1] + search_radius,
        center[2] - search_radius,
        center[2] + search_radius
    )

    candidate_ids = vtk.vtkIdList()

    cell_locator.FindCellsWithinBounds(
        bounds,
        candidate_ids
    )

    cut_cells = set()

    for j in range(
        candidate_ids.GetNumberOfIds()
    ):

        cell_id = int(
            candidate_ids.GetId(j)
        )

        if cell_intersects_local_cut(
            surface,
            cell_id,
            cut["center"],
            cut["normal"],
            cut["cut_radius"]
        ):

            cut_cells.add(
                cell_id
            )

    cells_to_remove.update(
        cut_cells
    )


    print(
        "Cut {}/{} ({}) removed {} surface triangles."
        .format(
            i + 1,
            len(cuts),
            cut["type"],
            len(cut_cells)
        )
    )


    if len(cut_cells) == 0:

        print(
            "    WARNING: this cut did not intersect the "
            "surface. Inspect the debug circle."
        )


if len(cells_to_remove) == 0:

    raise RuntimeError(
        "None of the proposed cuts intersected the surface."
    )


print("")
print(
    "Unique surface triangles removed:",
    len(cells_to_remove)
)


# ============================================================
# REMOVE CUTTING-RING TRIANGLES
# ============================================================

keep_array = vtk.vtkIntArray()

keep_array.SetName(
    "Keep_After_Centerline_Prune"
)

keep_array.SetNumberOfComponents(1)

keep_array.SetNumberOfTuples(
    surface.GetNumberOfCells()
)


for cell_id in range(
    surface.GetNumberOfCells()
):

    if cell_id in cells_to_remove:

        keep_array.SetTuple1(
            cell_id,
            0
        )

    else:

        keep_array.SetTuple1(
            cell_id,
            1
        )


surface_flagged = vtk.vtkPolyData()

surface_flagged.DeepCopy(
    surface
)

surface_flagged.GetCellData().AddArray(
    keep_array
)


threshold = vtk.vtkThreshold()

threshold.SetInputData(
    surface_flagged
)

threshold.SetInputArrayToProcess(
    0,
    0,
    0,
    vtk.vtkDataObject.FIELD_ASSOCIATION_CELLS,
    "Keep_After_Centerline_Prune"
)

threshold.SetLowerThreshold(0.5)
threshold.SetUpperThreshold(1.5)

threshold.SetThresholdFunction(
    vtk.vtkThreshold.THRESHOLD_BETWEEN
)

threshold.Update()


geometry = vtk.vtkGeometryFilter()

geometry.SetInputConnection(
    threshold.GetOutputPort()
)

geometry.Update()


severed_surface = vtk.vtkPolyData()

severed_surface.DeepCopy(
    geometry.GetOutput()
)


clean_severed = vtk.vtkCleanPolyData()

clean_severed.SetInputData(
    severed_surface
)

clean_severed.PointMergingOn()

clean_severed.Update()


severed_surface = vtk.vtkPolyData()

severed_surface.DeepCopy(
    clean_severed.GetOutput()
)


# ============================================================
# SAVE THE PRE-CONNECTIVITY SURFACE
#
# This contains BOTH the retained tree and all distal pieces,
# separated by the local cutting rings.
# ============================================================

write_polydata(
    severed_surface,
    SEVERED_DEBUG_FILE
)


print("")
print(
    "Severed debug surface written to:"
)

print(
    SEVERED_DEBUG_FILE
)


# ============================================================
# COUNT CONNECTED REGIONS AFTER CUTTING
# ============================================================

connectivity_all = (
    vtk.vtkPolyDataConnectivityFilter()
)

connectivity_all.SetInputData(
    severed_surface
)

connectivity_all.SetExtractionModeToAllRegions()

connectivity_all.ColorRegionsOn()

connectivity_all.Update()


print("")
print(
    "Connected regions after severing:",
    connectivity_all.GetNumberOfExtractedRegions()
)


# ============================================================
# KEEP ONLY THE SURFACE CONNECTED TO THE PA INLET
# ============================================================

connectivity_keep = (
    vtk.vtkPolyDataConnectivityFilter()
)

connectivity_keep.SetInputData(
    severed_surface
)

connectivity_keep.SetExtractionModeToClosestPointRegion()

connectivity_keep.SetClosestPoint(
    inlet_xyz[0],
    inlet_xyz[1],
    inlet_xyz[2]
)

connectivity_keep.Update()


final_surface = vtk.vtkPolyData()

final_surface.DeepCopy(
    connectivity_keep.GetOutput()
)


# ============================================================
# FINAL CLEAN
# ============================================================

final_clean = vtk.vtkCleanPolyData()

final_clean.SetInputData(
    final_surface
)

final_clean.PointMergingOn()

final_clean.Update()


final_surface = vtk.vtkPolyData()

final_surface.DeepCopy(
    final_clean.GetOutput()
)


# ============================================================
# WRITE FINAL RESULT
# ============================================================

write_polydata(
    final_surface,
    OUTPUT_FILE
)


print("")
print("============================================")
print("FINISHED")
print("============================================")

print("")
print(
    "Final surface points:",
    final_surface.GetNumberOfPoints()
)

print(
    "Final surface cells:",
    final_surface.GetNumberOfCells()
)

print("")
print(
    "Final pruned surface:"
)

print(
    OUTPUT_FILE
)

print("")
print(
    "The new pruning boundaries are intentionally OPEN."
)

print(
    "Inspect CUTS_DEBUG_FILE and SEVERED_DEBUG_FILE "
    "before doing any final capping/smoothing."
)
