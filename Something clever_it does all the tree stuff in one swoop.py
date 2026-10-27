from paraview.simple import *
import os

centerlines_path = r"C:\Users\Admin\Documents\3D_Slicer_Files\JL_Testing\Network_model_V4.vtp"

# ENDPOINT SELECTION: choose "Generations" or "Radius".
ENDPOINT_METHOD = "Radius"
GENERATION_THRESHOLD = 8
RADIUS_THRESHOLD = 1  # Same length units as the Radius array (usually cm).

method = ENDPOINT_METHOD.strip().lower()
if method not in ("generations", "radius"):
    raise ValueError('ENDPOINT_METHOD must be "Generations" or "Radius".')
if method == "radius" and RADIUS_THRESHOLD <= 0:
    raise ValueError("RADIUS_THRESHOLD must be positive.")

print("Endpoint method: {}".format(ENDPOINT_METHOD))
print("Endpoint cutoff: {}".format(
    GENERATION_THRESHOLD if method == "generations" else RADIUS_THRESHOLD))

# -------------------------
# 1) Load data sources
# -------------------------

if not os.path.exists(centerlines_path):
    raise RuntimeError(f"Missing file: {centerlines_path}")

centerlines_src = XMLPolyDataReader(FileName=[centerlines_path])
centerlines_src.UpdatePipeline()

# -------------------------
# 2) Clean filter to connect all centerline segments
# -------------------------

clean_pf1 = Clean(Input=centerlines_src)

# Merge points within this distance
clean_pf1.Tolerance = 0.00001

# Recommended defaults (shown for clarity)
clean_pf1.ConvertLinesToPoints = 0   # leave topology intact
clean_pf1.ConvertPolysToLines = 0
clean_pf1.ConvertStripsToPolys = 0

# If you want to guarantee point merging:
clean_pf1.PointMerging = 1           # default = On

clean_pf1.UpdatePipeline()

# -------------------------
# 3) Generates IDs so there are references for the points
# -------------------------

ID1 = PointAndCellIds(Input=clean_pf1)
ID1.UpdatePipeline()

# -------------------------
# 4) Programmable Filter 1 
# Assigns an inlet, outlets, and lengths to the centerlines
# -------------------------

pf1 = ProgrammableFilter(Input=[ID1])

pf1.Script = r"""
import vtk
from vtkmodules.util import numpy_support
import numpy as np
from collections import deque

# Get input and output
input = self.GetInputDataObject(0, 0)
output = self.GetOutput()

# Deep copy input to output to preserve geometry and existing arrays
output.DeepCopy(input)

# Get 'CenterlineIds' from CellData
id_array = input.GetPointData().GetArray("PointIds")
if id_array is None:
    raise RuntimeError("PointIds array not found in CellData.")

ids = numpy_support.vtk_to_numpy(id_array)


# Build connectivity map from cells
connectivity = {i: [] for i in range(input.GetNumberOfPoints())}
for cell_id in range(input.GetNumberOfCells()):
    cell = input.GetCell(cell_id)
    num_pts = cell.GetNumberOfPoints()
    for i in range(num_pts - 1):
        p1 = cell.GetPointId(i)
        p2 = cell.GetPointId(i + 1)
        connectivity[p1].append(p2)
        connectivity[p2].append(p1)

# Map cell IDs to point IDs
point_ids = np.full(input.GetNumberOfPoints(), -1)
for cell_id in range(input.GetNumberOfCells()):
    cell = input.GetCell(cell_id)
    for i in range(cell.GetNumberOfPoints()):
        pid = cell.GetPointId(i)
        point_ids[pid] = ids[cell_id]

# Find inlet point: first point from a cell with GroupId == 0
inlet_index = None
for cell_id in range(input.GetNumberOfCells()):
    if ids[cell_id] == 0:
        cell = input.GetCell(cell_id)
        inlet_index = cell.GetPointId(0)
        break

if inlet_index is None:
    raise RuntimeError("No inlet point found with GroupId == 0.")

# Identify outlet points (one neighbor, not inlet)
outlet_indices = [i for i in range(input.GetNumberOfPoints()) if len(connectivity[i]) == 1 and i != inlet_index]

# Mark inlet and outlet points
inlet_array = np.zeros(input.GetNumberOfPoints())
outlet_array = np.zeros(input.GetNumberOfPoints())

inlet_array[inlet_index] = 1
for oid in outlet_indices:
    outlet_array[oid] = 1

inlet_vtk = numpy_support.numpy_to_vtk(inlet_array)
inlet_vtk.SetName("IsInlet")
output.GetPointData().AddArray(inlet_vtk)

outlet_vtk = numpy_support.numpy_to_vtk(outlet_array)
outlet_vtk.SetName("IsOutlet")
output.GetPointData().AddArray(outlet_vtk)

# New BFS: compute distance from inlet to all reachable points
def bfs_all_distances(start, connectivity, input):
    visited = set()
    queue = deque([(start, 0.0)])  # (point_id, cumulative_distance)
    distances = np.full(input.GetNumberOfPoints(), -1.0)
    distances[start] = 0.0

    while queue:
        current, dist = queue.popleft()
        visited.add(current)

        for neighbor in connectivity[current]:
            if neighbor not in visited:
                p1 = np.array(input.GetPoint(current))
                p2 = np.array(input.GetPoint(neighbor))
                edge_length = np.linalg.norm(p2 - p1)
                new_dist = dist + edge_length

                if distances[neighbor] < 0 or new_dist < distances[neighbor]:
                    distances[neighbor] = new_dist
                    queue.append((neighbor, new_dist))

    return distances

# Compute segment distances from inlet
segment_distance = bfs_all_distances(inlet_index, connectivity, input)

# Compute straight-line distance from inlet
straight_line_distance = np.zeros(input.GetNumberOfPoints())
inlet_coord = np.array(input.GetPoint(inlet_index))

for i in range(input.GetNumberOfPoints()):
    pt = np.array(input.GetPoint(i))
    straight_line_distance[i] = np.linalg.norm(pt - inlet_coord)

linear_vtk = numpy_support.numpy_to_vtk(straight_line_distance)
linear_vtk.SetName("LinearDistanceFromInlet")
output.GetPointData().AddArray(linear_vtk)

# Compute tortuosity only for reachable points
tortuosity = np.full(input.GetNumberOfPoints(), -1.0)
for i in range(input.GetNumberOfPoints()):
    seg_dist = segment_distance[i]
    straight_dist = straight_line_distance[i]
    if seg_dist >= 0 and straight_dist > 0:
        tortuosity[i] = seg_dist / straight_dist
tortuosity_vtk = numpy_support.numpy_to_vtk(tortuosity)
tortuosity_vtk.SetName("Tortuosity")
output.GetPointData().AddArray(tortuosity_vtk)

segment_vtk = numpy_support.numpy_to_vtk(segment_distance)
segment_vtk.SetName("CenterlineLength")
output.GetPointData().AddArray(segment_vtk)

# Summary
print(f"Inlet point index (GroupId == 0): {inlet_index}")
print(f"Number of outlet points: {len(outlet_indices)}")
print("CenterlineLength and Tortuosity arrays added.")
"""
pf1.RequestInformationScript = ''
pf1.RequestUpdateExtentScript = ''
pf1.PythonPath = ''
pf1.UpdatePipeline()

# -------------------------
# 5) Programmable Filter 2 
# Counts the generations of the tree
# -------------------------

pf2 = ProgrammableFilter(Input=[pf1])

pf2.Script = r"""
from vtkmodules.util import numpy_support
import numpy as np
from collections import deque

# Get input
input = self.GetInput()

# Get arrays
points = input.GetPoints()
num_points = input.GetNumberOfPoints()

is_inlet_vtk = input.GetPointData().GetArray("IsInlet")
if is_inlet_vtk is None:
    raise RuntimeError("Required array 'IsInlet' not found.")

is_inlet = numpy_support.vtk_to_numpy(is_inlet_vtk)

# Build connectivity graph
connectivity = {i: [] for i in range(num_points)}
for cell_id in range(input.GetNumberOfCells()):
    cell = input.GetCell(cell_id)
    for j in range(cell.GetNumberOfPoints() - 1):
        p1 = cell.GetPointId(j)
        p2 = cell.GetPointId(j + 1)
        connectivity[p1].append(p2)
        connectivity[p2].append(p1)

# Compute number of connections per point
num_connections = np.array([len(connectivity[i]) for i in range(num_points)])

# Initialize Generation array
generation_array = np.full(num_points, -1, dtype=np.int32)  # -1 means unvisited

# BFS traversal from inlet points
inlet_indices = np.where(is_inlet == 1)[0]
queue = deque()

for inlet in inlet_indices:
    generation_array[inlet] = 0
    queue.append((inlet, 0, True))  # (point_id, generation, is_first_step)

while queue:
    current, gen, is_first_step = queue.popleft()

    for neighbor in connectivity[current]:

        if generation_array[neighbor] == -1:  # not visited

            # Do not increment immediately after the inlet
            if is_first_step:
                next_gen = gen

            else:
                # Increment ONLY when entering a bifurcation point.
                # A bifurcation has more than two connected neighbors.
                if num_connections[neighbor] > 2:
                    next_gen = gen + 1
                else:
                    next_gen = gen

            generation_array[neighbor] = next_gen
            queue.append((neighbor, next_gen, False))

# Attach Generation array
generation_vtk = numpy_support.numpy_to_vtk(generation_array)
generation_vtk.SetName("Generation")
input.GetPointData().AddArray(generation_vtk)

# Output
output.ShallowCopy(input)
"""
pf2.RequestInformationScript = ''
pf2.RequestUpdateExtentScript = ''
pf2.PythonPath = ''
pf2.UpdatePipeline()

# -------------------------
# 6) Trim downstream branches by generation OR first radius crossing
# -------------------------
# In Generations mode, preserve the original ParaView threshold behavior.
# In Radius mode, traverse the connected centerline graph from the inlet.
# Keep a point when first encountered with Radius <= RADIUS_THRESHOLD, but
# do not visit its descendants -- even if their radius later becomes larger.
# This MUST be a traversal, not a Threshold(Radius) filter: a pointwise
# threshold could reintroduce distal vessels beyond a stenotic region.

if method == "generations":
    thresh_filter = Threshold(Input=pf2)
    thresh_filter.Scalars = ['POINTS', 'Generation']
    thresh_filter.ThresholdMethod = 'Between'
    thresh_filter.LowerThreshold = -1
    thresh_filter.UpperThreshold = GENERATION_THRESHOLD
    thresh_filter.AllScalars = 1
    thresh_filter.UpdatePipeline()
    trimming_result = thresh_filter

else:
    radius_trim = ProgrammableFilter(Input=[pf2])
    radius_trim.OutputDataSetType = 'vtkPolyData'
    radius_trim.Script = r"""
import vtk
import numpy as np
import builtins as py_builtins  # Avoid ParaView shadowing Python sum/min/max
from collections import deque

RADIUS_THRESHOLD = __RADIUS_THRESHOLD__

src = self.GetInputDataObject(0, 0)
dst = self.GetOutputDataObject(0)
if not isinstance(src, vtk.vtkPolyData):
    raise RuntimeError("Radius trim expects centerlines as vtkPolyData.")

n = src.GetNumberOfPoints()
if n < 2:
    raise RuntimeError("Radius trim requires at least two centerline points.")

pd = src.GetPointData()
ra = pd.GetArray("Radius")
if ra is None:
    raise RuntimeError("Radius trim: required POINT array 'Radius' not found.")
if ra.GetNumberOfComponents() != 1:
    raise RuntimeError("Radius trim: 'Radius' must be a scalar point array.")

inlet = pd.GetArray("IsInlet")
if inlet is None:
    raise RuntimeError("Radius trim: IsInlet missing. Check upstream filters.")
inlets = [pid for pid in range(n) if inlet.GetTuple1(pid) >= 0.5]
if len(inlets) != 1:
    raise RuntimeError("Radius trim expects one inlet; found {}".format(len(inlets)))
root = inlets[0]

# Deduplicated graph; map every topological edge to a source cell for copying
# its cell arrays, even if the input is one or more long vtkPolyLine cells.
adj = [set() for _ in range(n)]
edge_cell = {}
for cid in range(src.GetNumberOfCells()):
    cell = src.GetCell(cid)
    for j in range(cell.GetNumberOfPoints() - 1):
        a = int(cell.GetPointId(j))
        b = int(cell.GetPointId(j + 1))
        if a == b:
            continue
        adj[a].add(b)
        adj[b].add(a)
        edge_cell.setdefault((a, b) if a < b else (b, a), cid)

# BFS builds a rooted tree, stops at the FIRST qualifying radius on each
# inlet-to-terminal path, and never resumes after that terminal.
parent = {root: -1}
queue = deque([root])
cut_points = set()
missing_radius = []
while queue:
    a = queue.popleft()
    value = float(ra.GetTuple1(a))
    if not np.isfinite(value):
        missing_radius.append(a)
        continue
    # The inlet itself is not treated as a branch termination. In particular,
    # this prevents a missing/zero inlet radius from deleting the whole tree.
    if a != root and value <= RADIUS_THRESHOLD:
        cut_points.add(a)
        continue
    for b in sorted(adj[a]):
        if b not in parent:
            parent[b] = a
            queue.append(b)

if missing_radius:
    raise RuntimeError(
        "Radius trim: non-finite Radius at traversed centerline IDs {}. "
        "Fix or interpolate Radius before trimming.".format(missing_radius[:12])
    )

if len(parent) < 2:
    raise RuntimeError("Radius trim retained no inlet-to-vessel edge.")

# Make a compact output with ONLY retained points, not orphaned vertices;
# each retained original edge becomes a vtkLine. This also handles source
# polylines that cross a cutoff in the middle of a cell.
keep_ids = sorted(parent)
new_id = {old: i for i, old in enumerate(keep_ids)}
pts = vtk.vtkPoints()
pts.SetDataType(src.GetPoints().GetDataType())
pts.SetNumberOfPoints(len(keep_ids))
for old, new in new_id.items():
    pts.SetPoint(new, src.GetPoint(old))

lines = vtk.vtkCellArray()
old_cell_ids = []
for old in keep_ids:
    p = parent[old]
    if p < 0:
        continue
    edge = (p, old) if p < old else (old, p)
    if edge not in edge_cell:
        raise RuntimeError("Cannot locate source cell for edge {}".format(edge))
    line = vtk.vtkLine()
    line.GetPointIds().SetId(0, new_id[p])
    line.GetPointIds().SetId(1, new_id[old])
    lines.InsertNextCell(line)
    old_cell_ids.append(edge_cell[edge])

result = vtk.vtkPolyData()
result.SetPoints(pts)
result.SetLines(lines)

# Copy ALL centerline point/cell arrays (numeric and string), including
# Generation, Radius, original PointIds, and VMTK labels.
result.GetPointData().CopyAllocate(src.GetPointData(), len(keep_ids))
for old, new in new_id.items():
    result.GetPointData().CopyData(src.GetPointData(), old, new)

result.GetCellData().CopyAllocate(src.GetCellData(), len(old_cell_ids))
for new_cid, old_cid in enumerate(old_cell_ids):
    result.GetCellData().CopyData(src.GetCellData(), old_cid, new_cid)

# Ensure the source data arrays are maintained through the compacted tree.
# PointIds are optional provenance; downstream filters do not require them.
if result.GetPointData().GetArray("Radius") is None:
    raise RuntimeError("Radius trim lost its Radius array while copying PointData.")
if result.GetPointData().GetArray("Generation") is None:
    raise RuntimeError("Radius trim lost its Generation array while copying PointData.")

# Refresh the outlets BEFORE PF3, so downstream code can immediately inspect
# new radius-cut terminal points. PF3 recomputes the flags again as normal.
deg = [0] * len(keep_ids)
for old in keep_ids:
    p = parent[old]
    if p >= 0:
        deg[new_id[p]] += 1
        deg[new_id[old]] += 1

out_flags = vtk.vtkIntArray()
out_flags.SetName("IsOutlet")
out_flags.SetNumberOfTuples(len(keep_ids))
cut_flags = vtk.vtkIntArray()
cut_flags.SetName("RadiusCutTerminal")
cut_flags.SetNumberOfTuples(len(keep_ids))
for old in keep_ids:
    new = new_id[old]
    out_flags.SetValue(new, int(old != root and deg[new] == 1))
    cut_flags.SetValue(new, int(old in cut_points))
result.GetPointData().RemoveArray("IsOutlet")
result.GetPointData().AddArray(out_flags)
result.GetPointData().RemoveArray("RadiusCutTerminal")
result.GetPointData().AddArray(cut_flags)

print("Radius trim: threshold <= {}".format(RADIUS_THRESHOLD))
print("Radius trim: kept {}/{} points, {} edges".format(
    len(keep_ids), n, len(old_cell_ids)))
print("Radius trim: first-radius cutoff terminals {}, total outlets {}".format(
    len(cut_points), py_builtins.sum(deg[new_id[old]] == 1 for old in keep_ids if old != root)))
print("Radius trim: original connected graph edges {}, unvisited/cut descendants {}".format(
    len(edge_cell), n-len(keep_ids)))

# Fail closed if no connected tubes were kept or if the graph is not tree-like.
# Duplicated overlapping centerline paths are normal (deduped above), but
# extra unique graph edges would make the chosen route ambiguous.
# The intended input is a tree. If it has a cycle, warn because the parent
# traversal deliberately selects one spanning-tree route rather than using
# ambiguous alternative paths.
seen_all = {root}
q_all = deque([root])
while q_all:
    a = q_all.popleft()
    for b in adj[a]:
        if b not in seen_all:
            seen_all.add(b)
            q_all.append(b)
reachable_edges = py_builtins.sum(len(adj[a] & seen_all) for a in seen_all) // 2
if reachable_edges >= len(seen_all):
    print("WARNING: centerline contains cycle(s); pruning uses an inlet-rooted BFS spanning tree.")

dst.ShallowCopy(result)
""".replace("__RADIUS_THRESHOLD__", repr(float(RADIUS_THRESHOLD)))
    radius_trim.RequestInformationScript = ''
    radius_trim.RequestUpdateExtentScript = ''
    radius_trim.PythonPath = ''
    radius_trim.UpdatePipeline()
    trimming_result = radius_trim

# -------------------------
# -------------------------
# 7) Programmable Filter 3 
# 
# -------------------------

pf3 = ProgrammableFilter(Input=[trimming_result])

pf3.Script = r"""

import vtk
from vtkmodules.util import numpy_support
import numpy as np
from collections import deque
import math
import builtins as py_builtins  # ParaView may override min/max/sum/any


# Get input and output
input = self.GetInputDataObject(0, 0)
output = self.GetOutput()

# Deep copy input to output to preserve geometry and existing arrays
output.DeepCopy(input)

# The downstream filter operates on current centerline point IDs.
# The original PointIds array is preserved when present, but it is not
# required for connectivity, length, outlet detection, or sidedness.
# In particular, no original IDs should be indexed by current cell IDs.

# Build connectivity map from cells
connectivity = {i: [] for i in range(input.GetNumberOfPoints())}
for cell_id in range(input.GetNumberOfCells()):
    cell = input.GetCell(cell_id)
    num_pts = cell.GetNumberOfPoints()
    for i in range(num_pts - 1):
        p1 = cell.GetPointId(i)
        p2 = cell.GetPointId(i + 1)
        connectivity[p1].append(p2)
        connectivity[p2].append(p1)

# Find inlet from the preserved flag, rather than assuming that the
# original PointIds still match current cell IDs after trimming.
inlet_index = None
existing_inlet = input.GetPointData().GetArray("IsInlet")
if existing_inlet is not None:
    candidates = [pid for pid in range(input.GetNumberOfPoints())
                  if existing_inlet.GetTuple1(pid) >= 0.5]
    if len(candidates) == 1:
        inlet_index = candidates[0]
if inlet_index is None:
    raise RuntimeError("Expected one IsInlet point after trimming.")

# Identify outlet points (one neighbor, not inlet)
outlet_indices = [i for i in range(input.GetNumberOfPoints()) if len(connectivity[i]) == 1 and i != inlet_index]

# Mark inlet and outlet points
inlet_array = np.zeros(input.GetNumberOfPoints())
outlet_array = np.zeros(input.GetNumberOfPoints())

inlet_array[inlet_index] = 1
for oid in outlet_indices:
    outlet_array[oid] = 1

inlet_vtk = numpy_support.numpy_to_vtk(inlet_array)
inlet_vtk.SetName("IsInlet")
output.GetPointData().AddArray(inlet_vtk)

outlet_vtk = numpy_support.numpy_to_vtk(outlet_array)
outlet_vtk.SetName("IsOutlet")
output.GetPointData().AddArray(outlet_vtk)

# New BFS: compute distance from inlet to all reachable points
def bfs_all_distances(start, connectivity, input):
    visited = set()
    queue = deque([(start, 0.0)])  # (point_id, cumulative_distance)
    distances = np.full(input.GetNumberOfPoints(), -1.0)
    distances[start] = 0.0

    while queue:
        current, dist = queue.popleft()
        visited.add(current)

        for neighbor in connectivity[current]:
            if neighbor not in visited:
                p1 = np.array(input.GetPoint(current))
                p2 = np.array(input.GetPoint(neighbor))
                edge_length = np.linalg.norm(p2 - p1)
                new_dist = dist + edge_length

                if distances[neighbor] < 0 or new_dist < distances[neighbor]:
                    distances[neighbor] = new_dist
                    queue.append((neighbor, new_dist))

    return distances

# --------------------------------
# Connectivity and segment utilities (no NumPy)
# --------------------------------
def build_connectivity(polydata):
    
    npts = polydata.GetNumberOfPoints()
    conn = [[] for _ in range(npts)]
    for cid in range(polydata.GetNumberOfCells()):
        cell = polydata.GetCell(cid)
        m = cell.GetNumberOfPoints()
        for j in range(m - 1):
            a = cell.GetPointId(j)
            b = cell.GetPointId(j + 1)
            if b not in conn[a]:
                conn[a].append(b)
            if a not in conn[b]:
                conn[b].append(a)
    return conn

def compute_degrees(connectivity):
    return [len(nbs) for nbs in connectivity]

def find_segments(centerline, connectivity):
    
    n = centerline.GetNumberOfPoints()
    deg = compute_degrees(connectivity)
    junctions = [i for i in range(n) if deg[i] != 2]

    visited_edge = set()  # undirected edges as (min(a,b), max(a,b))
    segments = []

    def mark_edge(a, b):
        visited_edge.add((a, b) if a < b else (b, a))

    def edge_seen(a, b):
        return ((a, b) if a < b else (b, a)) in visited_edge

    # Walk outward from each junction along each neighbor, stopping at the next junction.
    for j in junctions:
        for nb in connectivity[j]:
            if edge_seen(j, nb):
                continue

            seg = [j]
            prev = j
            curr = nb
            mark_edge(prev, curr)

            while True:
                seg.append(curr)
                if deg[curr] != 2:
                    break
                nbs = connectivity[curr]
                nxt = nbs[0] if nbs[1] == prev else nbs[1]
                if edge_seen(curr, nxt):
                    break
                prev, curr = curr, nxt
                mark_edge(prev, curr)

            segments.append(seg)

    # Optional: detect/handle pure loops with all degree==2.

    return segments

# --- DFS restricted to one side ---
def dfs_with_tag(start, tag_label, temp_tags):
    stack = [start]
    visited = set()

    while stack:
        current = stack.pop()
        if current in visited:
            continue
        visited.add(current)

        if generation[current] > 0.9 and temp_tags[current] is None:
            temp_tags[current] = tag_label

        for neighbor in connectivity[current]:
            # Do not cross trunk
            if generation[neighbor] < 0.1:
                continue
            if neighbor not in visited:
                stack.append(neighbor)

# Compute segment distances from inlet
segment_distance = bfs_all_distances(inlet_index, connectivity, input)

# Compute straight-line distance from inlet
straight_line_distance = np.zeros(input.GetNumberOfPoints())
inlet_coord = np.array(input.GetPoint(inlet_index))

for i in range(input.GetNumberOfPoints()):
    pt = np.array(input.GetPoint(i))
    straight_line_distance[i] = np.linalg.norm(pt - inlet_coord)

linear_vtk = numpy_support.numpy_to_vtk(straight_line_distance)
linear_vtk.SetName("LinearDistanceFromInlet")
output.GetPointData().AddArray(linear_vtk)

# Compute tortuosity only for reachable points
tortuosity = np.full(input.GetNumberOfPoints(), -1.0)
for i in range(input.GetNumberOfPoints()):
    seg_dist = segment_distance[i]
    straight_dist = straight_line_distance[i]
    if seg_dist >= 0 and straight_dist > 0:
        tortuosity[i] = seg_dist / straight_dist
tortuosity_vtk = numpy_support.numpy_to_vtk(tortuosity)
tortuosity_vtk.SetName("Tortuosity")
output.GetPointData().AddArray(tortuosity_vtk)

segment_vtk = numpy_support.numpy_to_vtk(segment_distance)
segment_vtk.SetName("CenterlineLength")
output.GetPointData().AddArray(segment_vtk)

# Summary
print(f"Number of outlet points: {len(outlet_indices)}")
print("CenterlineLength and Tortuosity arrays added.")

# --------------------------------
# Inputs
# --------------------------------
# Accept either a PolyData directly, or a MultiBlock with centerline in block 0

num_points = input.GetNumberOfPoints()

# Arrays
pd = input.GetPointData()


# --------------------------------
# Build segments
# --------------------------------
connectivity = build_connectivity(input)
segments = find_segments(input, connectivity)


# --------------------------------
# Prepare output arrays (per-point, segment-wise constant values)
# --------------------------------
seg_id_arr = vtk.vtkIntArray()
seg_id_arr.SetName("segment_ID")
seg_id_arr.SetNumberOfComponents(1)
seg_id_arr.SetNumberOfTuples(num_points)
seg_id_arr.FillComponent(0, -1)

# --- Get required array ---
generation_vtk = input.GetPointData().GetArray("Generation")

if generation_vtk is None:
    raise RuntimeError("Missing required array: 'Generation'.")

generation = numpy_support.vtk_to_numpy(generation_vtk)


# ============================================================
# ROBUST SIDEDNESS FROM GENERATION-1 BRANCH TOPOLOGY
# ============================================================
# A single anatomical branch can contain many Generation==1 point IDs, all
# immediately adjacent to Generation 0. Do NOT count point IDs as branches.
# Instead, find connected downstream trees, or split their shared root.
#
# Label assignment is deterministic (smallest point ID first), not anatomical.
# Flip the switch below after checking the left/right orientation in ParaView.
SWAP_LEFT_RIGHT = False
GENERATION_1_TOL = 0.10

# Build one *unique* adjacency set per point. The same edge may appear in
# multiple VMTK centerline cells, and a list with duplicates inflates degree.
connectivity = {pid: set() for pid in range(num_points)}
for cell_id in range(input.GetNumberOfCells()):
    cell = input.GetCell(cell_id)
    for j in range(cell.GetNumberOfPoints() - 1):
        a = int(cell.GetPointId(j))
        b = int(cell.GetPointId(j + 1))
        if a != b:
            connectivity[a].add(b)
            connectivity[b].add(a)

trunk_points = {pid for pid in range(num_points) if generation[pid] < 0.1}
nontrunk_points = {pid for pid in range(num_points) if generation[pid] >= 0.9}
gen1_points = {
    pid for pid in range(num_points)
    if abs(float(generation[pid]) - 1.0) < GENERATION_1_TOL
}


def find_components(allowed_nodes, blocked_nodes=()):
    allowed = set(allowed_nodes) - set(blocked_nodes)
    components = []

    while allowed:
        first = py_builtins.min(allowed)
        stack = [first]
        group = set()
        allowed.remove(first)

        while stack:
            pid = stack.pop()
            group.add(pid)
            for nb in connectivity[pid]:
                if nb in allowed:
                    allowed.remove(nb)
                    stack.append(nb)

        components.append(group)

    return sorted(components, key=lambda group: py_builtins.min(group))


def describe_components(components, label):
    print("{}: {} connected groups".format(label, len(components)))
    print("  group sizes:", [len(c) for c in components[:15]])
    if len(components) > 15:
        print("  (only first 15 group sizes shown)")


def find_shared_root_split(whole_component):

    candidate_points = sorted(
        (pid for pid in whole_component if pid in gen1_points),
        key=lambda pid: (
            float(segment_distance[pid]) if segment_distance[pid] >= 0.0
            else float('inf'),
            pid,
        )
    )

    for junction in candidate_points:
        neighbors = connectivity[junction] & whole_component
        if len(neighbors) < 2:
            continue

        # Deleting this point must disconnect downstream limbs. We find the
        # resulting regions, then identify the ones whose direct junction
        # neighbor is farther from the inlet than the junction itself.
        fragments = find_components(whole_component, {junction})
        if len(fragments) < 2:
            continue

        downstream_fragments = []
        for fragment in fragments:
            attached = neighbors & fragment
            if not attached:
                continue

            if py_builtins.any(
                segment_distance[nb] > segment_distance[junction] + 1e-9
                for nb in attached
            ):
                downstream_fragments.append(fragment)

        if len(downstream_fragments) == 2:
            left_right = sorted(
                downstream_fragments, key=lambda group: py_builtins.min(group)
            )
            shared = whole_component - left_right[0] - left_right[1]
            print("  Shared Generation-1 junction:", junction)
            print("  Shared/unlabelled points:", len(shared))
            print("  Downstream component sizes:",
                  len(left_right[0]), len(left_right[1]))
            return left_right, shared

    return None, None


# Primary strategy: remove the Generation-0 trunk and count the downstream
# connected trees. Several Gen-1 points can touch the trunk within each tree.
nontrunk_components = find_components(nontrunk_points)

print("Generation-1 point IDs:", sorted(gen1_points))
print("Gen-1 points immediately adjacent to trunk:", sorted(
    pid for pid in gen1_points
    if py_builtins.any(nb in trunk_points for nb in connectivity[pid])
))
describe_components(nontrunk_components, "All downstream points")

shared_points = set()
if len(nontrunk_components) == 2:
    branch_groups = nontrunk_components
    print("  Using the two downstream connected trees.")

elif len(nontrunk_components) == 1:
    branch_groups, shared_points = find_shared_root_split(
        nontrunk_components[0]
    )
    if branch_groups is None:
        raise RuntimeError(
            "All downstream points form one connected region, but no "
            "unambiguous two-arm Generation-1 junction was found. "
            "Inspect centerline connectivity and duplicate point geometry."
        )

else:
    # A disconnected higher-generation tip can create extra components. If
    # Generation 1 still forms exactly two groups, assign their connected
    # downstream trees and leave detached regions unassigned for inspection.
    gen1_components = find_components(gen1_points)
    describe_components(gen1_components, "Generation-1-only points")

    if len(gen1_components) == 2:
        selected = []
        for gen1_component in gen1_components:
            matching = [
                group for group in nontrunk_components
                if bool(group & gen1_component)
            ]
            if len(matching) != 1:
                raise RuntimeError(
                    "Generation-1 points span multiple disconnected "
                    "regions. Cannot identify two unique branches."
                )
            selected.append(matching[0])
        if selected[0] is selected[1]:
            raise RuntimeError("Both Generation-1 groups share the same region.")
        branch_groups = selected
        print("  WARNING: Some higher-generation components are detached "
              "and will retain Sidedness='null'.")
    else:
        raise RuntimeError(
            "Found {} downstream components, including {} Generation-1 "
            "components. The graph does not unambiguously resolve to "
            "two primary branches; refusing to assign arbitrary labels. "
            "A centerline VTP would be needed to resolve the topology."
            .format(len(nontrunk_components), len(gen1_components))
        )

# The two graph components are initially ordered by point ID. Depending on
# point numbering, Left/Right may need to be swapped based on anatomy.
if SWAP_LEFT_RIGHT:
    branch_groups = [branch_groups[1], branch_groups[0]]

branch1 = py_builtins.min(branch_groups[0])
branch2 = py_builtins.min(branch_groups[1])
temp_tags = [None] * num_points

for pid in branch_groups[0]:
    temp_tags[pid] = 'Right'
for pid in branch_groups[1]:
    if temp_tags[pid] is not None:
        raise RuntimeError("Left/Right groups overlap at point {}.".format(pid))
    temp_tags[pid] = 'Left'

# Trunk and any shared Generation-1 bifurcation point stay 'null'; a single
# vtkStringArray value cannot represent Left and Right at the same point ID.
for pid in trunk_points | shared_points:
    temp_tags[pid] = None

sidedness_vtk = vtk.vtkStringArray()
sidedness_vtk.SetName('Sidedness')
sidedness_vtk.SetNumberOfValues(num_points)
for pid, tag in enumerate(temp_tags):
    sidedness_vtk.SetValue(pid, tag if tag is not None else 'null')

# Avoid accumulating a stale Sidedness array when the filter is rerun.
output.GetPointData().RemoveArray('Sidedness')
output.GetPointData().AddArray(sidedness_vtk)

print('Left points:', temp_tags.count('Left'))
print('Right points:', temp_tags.count('Right'))
print('Unassigned/shared points:', temp_tags.count(None))
print('Unassigned non-trunk points:', py_builtins.sum(
    temp_tags[pid] is None for pid in nontrunk_points
))
print('Sidedness array added.')


# --------------------------------
# Attach arrays & output
# --------------------------------

# seg_id_arr is already a vtkIntArray; do not pass it through numpy_to_vtk.
output.GetPointData().RemoveArray('segment_ID')
output.GetPointData().AddArray(seg_id_arr)

"""
pf3.RequestInformationScript = ''
pf3.RequestUpdateExtentScript = ''
pf3.PythonPath = ''
pf3.UpdatePipeline()

# -------------------------
# 8) Extract Surface
# 
# -------------------------

Extract_for_vtp = ExtractSurface(Input=[pf3])
Extract_for_vtp.UpdatePipeline()
