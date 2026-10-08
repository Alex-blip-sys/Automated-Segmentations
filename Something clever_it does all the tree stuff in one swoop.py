from paraview.simple import *
import os

centerlines_path = r"C:\Users\Admin\Documents\3D_Slicer_Files\JL_Testing\Network_model_V2.vtp"
thres = 15

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

ID1 = GenerateIds(Input=clean_pf1)
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
# 6) Threshold to chop off a certain number of branches 
# 
# -------------------------

# 2. Create the Threshold filter
# Input can also be the active source: Threshold(Input=GetActiveSource())
thresh_filter = Threshold(Input=pf2)

# 3. Configure the Threshold parameters
# Format: ['ATTRIBUTE_MODE', 'ARRAY_NAME']
# ATTRIBUTE_MODE can be 'POINTS' (Point Data) or 'CELLS' (Cell Data)
thresh_filter.Scalars = ['POINTS', 'Generation']

# Set the method and range (ParaView 5.X / 6.X standard)
thresh_filter.ThresholdMethod = 'Between'
thresh_filter.LowerThreshold = -1
thresh_filter.UpperThreshold = thres

# Optional: If using Point Data, decide whether ALL points of a cell 
# must meet the criteria, or just ANY point. (1 = True, 0 = False)
thresh_filter.AllScalars = 1 

# 4. Update the pipeline to calculate the threshold
thresh_filter.UpdatePipeline()

# -------------------------
# 7) Programmable Filter 3 
# 
# -------------------------

pf3 = ProgrammableFilter(Input=[thresh_filter])

pf3.Script = r"""

import vtk
from vtkmodules.util import numpy_support
import numpy as np
from collections import deque
import math


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


# --------------------------------
# Compute per-segment values and paint them on all points in the segment
# --------------------------------
for seg_idx, seg in enumerate(segments):

    # Paint values to all points in the segment
    for pid in seg:
        seg_id_arr.SetTuple1(pid, int(seg_idx))

# --------------------------------
# Attach arrays & output
# --------------------------------

segment_vtk = numpy_support.numpy_to_vtk(seg_id_arr)
segment_vtk.SetName("segment_ID")
output.GetPointData().AddArray(segment_vtk)

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