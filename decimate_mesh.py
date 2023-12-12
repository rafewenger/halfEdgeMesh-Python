## \file decimate_mesh.py
#  Some simple mesh decimation routines.
#  - Collapse, join, split edges and cells.
#  - Uses data structure HALF_EDGE_MESH_EDIT_BASE.
#  @version 0.4.0

import math
from math import sqrt
from math import acos
import sys
from time import time, localtime, strftime

import half_edge_mesh_edit
import half_edge_mesh_IO

from half_edge_mesh_elements\
    import HMESH_VERTEX_BASE, HMESH_HALF_EDGE_BASE, HMESH_CELL_BASE
from half_edge_mesh_edit import HALF_EDGE_MESH_EDIT_BASE
from half_edge_mesh_compute import print_coordP
from half_edge_mesh_compute\
    import compute_min_max_edge_length_squared,\
        compute_cell_min_max_edge_length_squared,\
        compute_min_cell_edge_length_ratio_squared,\
        compute_cos_min_max_cell_angles,\
        compute_angle_info
from half_edge_mesh_angle_info import HMESH_ANGLE_INFO


def main(argv):

    global input_filename, output_filename
    global flag_silent, flag_terse, flag_no_warn, flag_time
    global flag_collapse_short_edges, flag_split_long_edges
    global flag_split_cells, flag_join_cells
    global flag_triangulate_cells
    global flag_allow_non_manifold, flag_fail_on_non_manifold
    global flag_reduce_checks

    begin_time = time()

    # Number of cells in a "large" data set.
    LARGE_DATA_NUM_CELLS = 1000

    # Initialize
    input_filename = None
    output_filename = None
    InitFlags()

    parse_command_line(sys.argv)

    mesh = HALF_EDGE_MESH_EDIT_BASE\
        (HMESH_VERTEX_BASE,HMESH_HALF_EDGE_BASE,HMESH_CELL_BASE)

    half_edge_mesh_IO.open_and_read_off_file(input_filename, mesh)

    check_input_mesh(mesh, input_filename, flag_silent or flag_no_warn)

    time2 = time()

    try:
        if (mesh.NumCells() == 0):
            print("Mesh in file " + input_filename + " has no cells.")
            print("  Nothing to edit.")
            print("  Exiting...")
            return 0

        if not(flag_reduce_checks):
            flag_reduce_checks =\
                reduce_checks_on_large_datasets\
                    (mesh, flag_no_warn, LARGE_DATA_NUM_CELLS)

        if (flag_collapse_short_edges):
            collapse_shortest_edge_in_each_cell\
                (mesh, flag_terse, flag_no_warn)

        if (flag_split_long_edges):
            split_longest_edge_in_each_cell(mesh, flag_terse, flag_no_warn)

        if (flag_split_cells):
            split_each_cell(mesh, flag_terse, flag_no_warn)

        if (flag_join_cells):
            join_each_cell(mesh, flag_terse, flag_no_warn)

        if (flag_triangulate_cells):
            triangulate_each_cell(mesh, flag_terse, flag_no_warn)

        passed_check = check_mesh(mesh, flag_silent or flag_no_warn)

        if not(flag_silent):
            if (passed_check):
                print("Mesh data structure passed check.")

        if not(flag_silent):
            print()
            print_mesh_info(mesh)

    except Exception as e:
        print(e)
        sys.stderr.write("Exiting.")
        exit(-1)

    time3 = time()

    if (output_filename is None):
        output_filename = "out.off"
    if (output_filename == input_filename):
        output_filename = "out2.off"

    if not(flag_silent):
        print()
        print("Writing file: " + output_filename + ".")

    half_edge_mesh_IO.open_and_write_file(output_filename, mesh)

    end_time = time()

    if (flag_time):
        print_time("Time to read file:    ", (time2-begin_time))
        print_time("Time to process mesh: ", (time3-time2))
        print_time("Time to write file:   ", (end_time-time3))
        print_time("Total time:           ", (end_time-begin_time))

    return 0


# *** Collapse edge routines ***

## Collapse edge
def collapse_edge(mesh, ihalf_edge, flag_terse, flag_no_warn, flag_check):
    global flag_allow_non_manifold

    half_edge = mesh.HalfEdge(ihalf_edge)

    if mesh.IsIllegalEdgeCollapseH(ihalf_edge):
        if not(flag_no_warn):
            print_illegal_edge_collapse(half_edge)
        return

    if mesh.DoesEdgeCollapseChangeMeshTopologyH(ihalf_edge):
        if not(flag_no_warn):
            print_edge_collapse_topology_change(mesh, half_edge)
        if not(flag_allow_non_manifold):
            if not(flag_no_warn):
                print_skipped_edge_collapse(half_edge)
            return

    if not(flag_terse):
        print("Collapsing edge (" + half_edge.EndpointsStr(",") + ").")

    vnew = mesh.CollapseEdge(ihalf_edge)
    if (vnew is None):
        print("Skipped illegal collapse of edge (" +\
                half_edge.EndpointsStr(",") + ").")

    if flag_check:
        check_mesh(mesh, flag_no_warn)


## Collapse shortest cell edge.
def collapse_shortest_cell_edge\
    (mesh, icell, flag_terse, flag_no_warn, flag_check):

    cell = mesh.Cell(icell)
    if (cell is None):
        return

    minL, maxL, ihalf_edge_min, ihalf_edge_max =\
        compute_cell_min_max_edge_length_squared(cell)
    collapse_edge(mesh, ihalf_edge_min, flag_terse, flag_no_warn, flag_check)


## Collapse shortest edge in each cell.
def collapse_shortest_edge_in_each_cell(mesh, flag_terse, flag_no_warn):
    global flag_reduce_checks

    n = mesh.NumCells()
    flag_check = not(flag_reduce_checks)

    # Create a list of the cell indices.
    cell_indices_list = list(mesh.CellIndices());
    cell_indices_list.sort()

    # DO NOT iterate over CellIndices() directly, since collapse/split/join
    #   may delete cells.
    kount = 0
    for icell in cell_indices_list:

        # Check if cell index is valid.
        #   cell may be none if it has been deleted from the cell dictionary
        cell = mesh.Cell(icell)
        if (cell is None):
            continue

        collapse_shortest_cell_edge\
            (mesh, icell, flag_terse, flag_no_warn, flag_check)
        kount = kount+1

        if (flag_reduce_checks):
            if (kount == n/2):
                check_mesh(mesh, flag_no_warn)


## Print illegal edge collapse message.
def print_illegal_edge_collapse(half_edge):
    print("Collapse of edge  (" + half_edge.EndpointsStr(",") +\
            ") is illegal.")
    print("  Some cell contains vertices " +\
            half_edge.EndpointsStr(" and ") + " but not edge (" +\
            half_edge.EndpointsStr(",") + ").")


## Print information about any topology changes caused by edge collapse.
#  @param mesh Half edge mesh.
#  @param half_edge Half edge to be collapsed.
#    - Note: Data structure, not half edge index.
def print_edge_collapse_topology_change(mesh, half_edge):
    ihalf_edge = half_edge.Index()
    icell = half_edge.CellIndex()
    flag, ivC = mesh.FindTriangleHole(ihalf_edge)
    if (flag):
        print("Collapsing edge (" + half_edge.EndpointsStr(",") +\
                ") will change the mesh topology.")
        print("  Vertices (" + half_edge.EndpointsStr(", ") +\
                ", " + str(ivC) + ") form a triangle hole.")

    if (mesh.IsInteriorEdgeWithBoundaryVertices(ihalf_edge)):
        print("Collapsing edge (" + half_edge.EndpointsStr(",") +\
                ") merges two non-adjacent boundary vertices.")

    if mesh.IsIsolatedTriangle(icell):
        print("Collapsing edge (" + half_edge.EndpointsStr(",") +\
                f") will delete isolated triangle {icell}.")

    if mesh.IsInTetrahedron(icell):
        print("Collapsing edge (" + half_edge.EndpointsStr(",") +\
                ") will collapse a tetrahedron.")


## Print skipped edge collapse.
def print_skipped_edge_collapse(half_edge):
    print("Skipped collapse of edge (" +\
            half_edge.EndpointsStr(",") + ").")


# *** Split edge routines. ***

## Split edge.
def split_edge(mesh, half_edge, flag_terse, flag_no_warn, flag_check):

    if not(flag_terse):
        print("Splitting edge (" + half_edge.EndpointsStr(",") + ").")

    vnew = mesh.SplitEdge(half_edge.Index())
    if (vnew is None):
        print("Split of edge (" + half_edge.EndpointsStr(",") + ") failed.")

    if (flag_check):
        check_mesh(mesh, flag_no_warn)


## Split longest cell edge.
def split_longest_cell_edge\
    (mesh, icell, flag_terse, flag_no_warn, flag_check):
    cell = mesh.Cell(icell)
    if (cell is None):
        return

    minL, maxL, ihalf_edge_min, ihalf_edge_max =\
        compute_cell_min_max_edge_length_squared(cell)

    half_edge_max = mesh.HalfEdge(ihalf_edge_max)

    split_edge(mesh, half_edge_max, flag_terse, flag_no_warn, flag_check)


## Split longest edge in each cell.
def split_longest_edge_in_each_cell(mesh, flag_terse, flag_no_warn):

    n = mesh.NumCells()
    flag_check = not(flag_reduce_checks)

    # Create a list of the cell indices.
    cell_indices_list = list(mesh.CellIndices());
    cell_indices_list.sort()

    # DO NOT iterate over CellIndices() directly, since collapse/split/join
    #   may delete cells.
    kount = 0
    for icell in cell_indices_list:

        # Check if cell index is valid.
        #   cell may be none if it has been deleted from the cell dictionary
        cell = mesh.Cell(icell)
        if (cell is None):
            continue

        split_longest_cell_edge\
            (mesh, icell, flag_terse, flag_no_warn, flag_check)
        kount = kount + 1

        if (flag_reduce_checks):
            # Check mesh halfway through.
            if (kount == n/2):
                check_mesh(mesh, flag_no_warn)


# *** Split cell routines. ***

## Split cell with diagonal (half_edgeA.FromVertex(), half_edgeB.FromVertex())
#  - Returns split edge.
#  - Returns None if split fails.
def split_cell(mesh, ihalf_edgeA, ihalf_edgeB,\
                flag_terse, flag_no_warn, flag_check):
    half_edgeA = mesh.HalfEdge(ihalf_edgeA)
    half_edgeB = mesh.HalfEdge(ihalf_edgeB)
    ivA = half_edgeA.FromVertexIndex()
    ivB = half_edgeB.FromVertexIndex()
    icell = half_edgeA.CellIndex()

    flag = check_split_cell(mesh, ihalf_edgeA, ihalf_edgeB, flag_no_warn)

    if (mesh.IsIllegalSplitCell(ihalf_edgeA, ihalf_edgeB)):
        return None

    if (flag or flag_allow_non_manifold):
        if not(flag_terse):
            print(f"Splitting cell {icell} at diagonal ({ivA},{ivB}).")

        split_edge = mesh.SplitCell(ihalf_edgeA, ihalf_edgeB)
        if (split_edge is None):
            print(f"Split of cell {icell} at diagonal ({ivA},{ivB}) failed.")

        if (flag_check):
            check_mesh(mesh, flag_no_warn)

        return split_edge
    else:
        if not(flag_terse):
            print(f"Skipping split of cell {icell} at diagonal ({ivA},{ivB}).")

        return None


## Split cell at largest angle.
#  - Split cell at vertex forming the largest angle.
def split_cell_at_largest_angle\
        (mesh, icell, flag_terse, flag_no_warn, flag_check):
    cell = mesh.Cell(icell)
    if (cell.IsTriangle()):
        return

    cos_minA, cos_maxA, ihalf_edge_min, ihalf_edge_max, flag_zero =\
        compute_cos_min_max_cell_angles(cell)

    half_edgeA = mesh.HalfEdge(ihalf_edge_max)
    numv = cell.NumVertices();
    half_numv = numv//2
    half_edgeB = half_edgeA
    for i in range(half_numv):
        half_edgeB = half_edgeB.NextHalfEdgeInCell()

    ihalf_edgeA = half_edgeA.Index()
    ihalf_edgeB = half_edgeB.Index()
    split_edge = split_cell(mesh, ihalf_edgeA, ihalf_edgeB,\
                                flag_terse, flag_no_warn, flag_check)
    if (split_edge is None):
        # Cannot split cell at largest angle.
        return

    if (flag_check):
        check_mesh(mesh, flag_no_warn)


## Split each cell.
def split_each_cell(mesh, flag_terse, flag_no_warn):
    global flag_reduce_checks

    n = mesh.MaxCellIndex()
    flag_check = not(flag_reduce_checks)

    # Create a list of the cell indices.
    cell_indices_list = list(mesh.CellIndices());
    cell_indices_list.sort()

    # DO NOT iterate over CellIndices() directly, since collapse/split/join
    #   may delete cells.
    kount = 0
    for icell in cell_indices_list:

        # Check if cell index is valid.
        #   cell may be none if it has been deleted from the cell dictionary
        cell = mesh.Cell(icell)
        if (cell is None):
            continue

        split_cell_at_largest_angle\
            (mesh, icell, flag_terse, flag_no_warn, flag_check)
        kount = kount+1

        if (flag_reduce_checks):
            if (kount == n/2):
                check_mesh(mesh, flag_no_warn)


## Triangulate cell from vertex half_edge0.FromVertex()
#  - Return true if triangulation succeeds.
def triangulate_cell_from_vertex\
        (mesh, ihalf_edge, flag_no_warn, flag_check):
    half_edge = mesh.HalfEdge(ihalf_edge)
    iv = half_edge.FromVertexIndex()
    icell = half_edge.CellIndex()
    cell = half_edge.Cell()
    if (cell.IsTriangle()):
        # Cell is already a triangle.
        return False

    if mesh.DoesTriangulateCellFromVertexChangeTopology(ihalf_edge):
        if not(flag_no_warn):
            print(f"Triangulating cell {icell} from vertex {iv} changes mesh topology.")
            print(f"  Skipping triangulation of cell {icell} from vertex {iv}.")
        return False

    if not(flag_terse):
        print(f"Triangulating cell {icell} from vertex {iv}.")

    mesh.TriangulateCellFromVertex(ihalf_edge)

    if (flag_check):
        check_mesh(mesh, flag_no_warn)

    return True


## Triangulate cell containing half edge.
#  - Attempt to triangulate from half_edgeA.FromVertex().
#  - If triangulation fails, try other vertices.
def triangulate_cell\
        (mesh, ihalf_edgeA, flag_terse, flag_no_warn, flag_check):
    half_edge = mesh.HalfEdge(ihalf_edgeA)
    cellA = half_edge.Cell()
    for i in range(0,cellA.NumVertices()):
        ihalf_edge = half_edge.Index()
        if triangulate_cell_from_vertex\
            (mesh, ihalf_edge, flag_no_warn, flag_check):
            return True
        half_edge = half_edge.NextHalfEdgeInCell()

    # Unable to find acceptable triangulation.
    return False


## Triangulate cell from vertex with largest angle.
#  - Attempt to triangulate from vertex with largest angle.
#  - If triangulation fails, try other vertices.
def triangulate_cell_from_vertex_with_largest_angle\
        (mesh, icell, flag_terse, flag_no_warn, flag_check):
    cell = mesh.Cell(icell)
    cos_minA, cos_maxA, ihalf_edge_min, ihalf_edge_max, flag_zero =\
        compute_cos_min_max_cell_angles(cell)

    return triangulate_cell\
        (mesh, ihalf_edge_max, flag_terse, flag_no_warn, flag_check)


## Triangulate each cell.
def triangulate_each_cell(mesh, flag_terse, flag_no_warn):
    global flag_reduce_checks

    NUM_VERT_PER_TRIANGLE = 3;
    n = mesh.MaxCellIndex()
    flag_check = not(flag_reduce_checks)

    # Create a list of the cell indices.
    cell_indices_list = list(mesh.CellIndices());
    cell_indices_list.sort()

    # DO NOT iterate over CellIndices() directly, since collapse/split/join
    #   may delete cells.
    kount = 0
    for icell in cell_indices_list:

        # Check if cell index is valid.
        #   cell may be none if it has been deleted from the cell dictionary
        cell = mesh.Cell(icell)
        if (cell is None):
            continue

        if (cell.NumVertices() <= NUM_VERT_PER_TRIANGLE):
            continue

        flag = triangulate_cell_from_vertex_with_largest_angle\
            (mesh, icell, flag_terse, flag_no_warn, flag_check)

        kount = kount+1

        if (flag_reduce_checks):
            if (kount == n/2):
                check_mesh(mesh, flag_no_warn)


# *** Join cell routines ***

def join_two_cells(mesh, ihalf_edge, flag_terse, flag_no_warn, flag_check):
    half_edge = mesh.HalfEdge(ihalf_edge)
    half_edgeX = half_edge.NextHalfEdgeAroundEdge()
    icell = half_edge.CellIndex()
    icellX = half_edgeX.CellIndex()

    flag = check_join_cell(mesh, ihalf_edge, flag_no_warn)

    if (mesh.IsIllegalJoinTwoCells(ihalf_edge)):
        return

    if flag:
        if not(flag_terse):
            print(f"Joining cell {icell} to cell {icellX} by deleting edge (" +\
                    half_edge.EndpointsStr(",") + ").")

        new_cell = mesh.JoinTwoCells(half_edge.Index())
        if (new_cell is None):
            print(f"Join of cell {icell} to cell {icellX} failed.")
        else:
            if (flag_check):
                check_mesh(mesh, flag_no_warn)

        return
    else:
        if not(flag_terse):
            print(f"Skipping join of cell {icell} with cell {icellX}.");


## Attempt to join each cell by deleting longest edge.
def join_each_cell(mesh, flag_terse, flag_no_warn):
    # Don't join cells with MAX_NUMV or more vertices.
    MAX_NUMV = 6
    n = mesh.NumCells()

    flag_check = not(flag_reduce_checks)

    # Create a list of the cell indices.
    cell_indices_list = list(mesh.CellIndices());
    cell_indices_list.sort()

    # DO NOT iterate over CellIndices() directly, since collapse/split/join
    #   may delete cells.
    kount = 0
    for icell in cell_indices_list:

        # Check if cell index is valid.
        #   cell may be none if it has been deleted from the cell dictionary
        cell = mesh.Cell(icell)
        if (cell is None):
            continue

        if (cell.NumVertices() >= MAX_NUMV):
            # Don't let the cell get to large.
            continue

        Lmin, Lmax, ihalf_edge_min, ihalf_edge_max =\
            compute_cell_min_max_edge_length_squared(cell)

        half_edge = mesh.HalfEdge(ihalf_edge_max)

        half_edgeX = half_edge.NextHalfEdgeAroundEdge()
        if (half_edgeX.Cell().NumVertices() >= MAX_NUMV):
            # Don't let the cell get too large.
            continue

        ihalf_edge = half_edge.Index()
        join_two_cells\
            (mesh, ihalf_edge, flag_terse, flag_no_warn, flag_check)
        kount = kount+1

        if (flag_reduce_checks):
            if (kount == n/2):
                check_mesh(mesh, flag_no_warn)


# *** Check routines ***

def exit_on_non_manifold():
    sys.stderr.write\
        ("Detected non-manifold or inconsistent orientations.\n")
    sys.stderr.write("  Exiting...")

    exit(-1)


def check_oriented_manifold(mesh, flag_no_warn):

    flag_manifold_vertex, flag_manifold_edge, iv, ihalf_edgeA =\
        mesh.CheckManifold()

    if not(flag_manifold_edge):
        if not(flag_no_warn):
            mesh.PrintNonManifoldEdge(sys.stderr, "Warning:", ihalf_edgeA)

        # Non-manifold edge automatically implies inconsistent orientations.
        return False

    flag_orientation, ihalf_edgeB = mesh.CheckOrientation()

    if flag_orientation:
        if not(flag_manifold_vertex):
            if not(flag_no_warn):
                mesh.PrintNonManifoldVertex(sys.stderr, "Warning:", iv)

            return False
    else:
        if not(flag_no_warn):
            if not(flag_manifold_vertex):
                sys.stderr.write\
                    (f"Warning: Non-manifold vertex or inconsistent orientations in cells incident on vertex {iv}.\n")
            else:
                mesh.PrintNotOriented(sys.stderr, "Warning:", ihalf_edgeB)

        return False

    return True


def check_mesh(mesh, flag_no_warn):
    global flag_fail_on_non_manifold

    flag, error_msg = mesh.CheckAll()
    if not(flag):
        sys.stderr.write("Error detected in mesh data structure.\n")
        if not(error_msg is None):
            sys.stderr.write(error_msg + "\n")

        exit(-1)

    if (not(flag_no_warn) or flag_fail_on_non_manifold):
        flag_oriented_manifold = check_oriented_manifold(mesh, flag_no_warn)

        if (flag_fail_on_non_manifold and not(flag_oriented_manifold)):
            exit_on_non_manifold()

        return flag_oriented_manifold

    else:
        return True


## Check if input is an oriented manifold
def check_input_oriented_manifold(mesh, input_filename, flag_no_warn):

    flag_manifold_vertex, flag_manifold_edge, iv, ihalf_edgeA =\
        mesh.CheckManifold()

    warning_prefix = "Warning: Mesh input from file " + input_filename
    if not(flag_manifold_edge):
        if not(flag_no_warn):
            sys.stderr.write\
                (warning_prefix + " has non-manifold edges.\n")
            mesh.PrintNonManifoldEdge(sys.stderr, "  ", ihalf_edgeA)

        # Non-manifold edge automatically implies inconsistent orientations.
        return False

    flag_orientation, ihalf_edgeB = mesh.CheckOrientation()

    if flag_orientation:
        if not(flag_manifold_vertex):
            if not(flag_no_warn):
                sys.stderr.write\
                    (warning_prefix + " has non-manifold vertices.\n")
                mesh.PrintNonManifoldVertex(sys.stderr, "  ", iv)

            return False
    else:
        if not(flag_no_warn):
            if not(flag_manifold_vertex):
                sys.stderr.write\
                    (warning_prefix + " has non-manifold vertex or inconsistent cell orientations.\n")
                sys.stderr.write\
                    (f"Warning: Non-manifold vertex or inconsistent orientations in cells incident on vertex {iv}.\n")
            else:
                sys.stderr.write\
                    (warning_prefix + " has inconsistent cell orientations.\n")
                mesh.PrintNotOriented(sys.stderr, "  ", ihalf_edgeB)

        return False

    return True


## Check mesh read from input_filename
def check_input_mesh(mesh, input_filename, flag_no_warn):
    global flag_fail_on_non_manifold

    flag, error_msg = mesh.CheckAll()
    if not(flag):
        sys.stderr.write("Error detected in mesh data structure created\n")
        sys.stderr.write("  from input file " + input_filename + ".")
        if not(error_msg is None):
            sys.stderr.write(error_msg + "\n")

        exit(-1)

    if (not(flag_no_warn) or flag_fail_on_non_manifold):
        flag_oriented_manifold = \
            check_input_oriented_manifold(mesh, input_filename, flag_no_warn)
        if not(flag_oriented_manifold):
            sys.stderr.write("\n")

        if (flag_fail_on_non_manifold and not(flag_oriented_manifold)):
            exit_on_non_manifold()

        return flag_oriented_manifold

    else:
        return True


## Print a warning message if collapsing half_edge is illegal or
#    will change mesh topology.
#  - Return True if collapse is not illegal and does not change
#    mesh topology.
def check_edge_collapse(mesh, ihalf_edge, flag_no_warn):
    half_edge = mesh.HalfEdge(ihalf_edge)
    icell = half_edge.CellIndex()
    return_flag = True

    if (mesh.IsIllegalEdgeCollapseH(ihalf_edge)):
        if not(flag_no_warn):
            print("Collapse of edge  (" + half_edge.EndpointsStr(",") +\
                    ") is illegal.")
            print("  Some cell contains vertices " +\
                    half_edge.EndpointsStr(" and ") + " but not edge (" +\
                    half_edge.EndpointsStr(",") + ").")

        return_flag = False

    flag, ivC = mesh.FindTriangleHole(ihalf_edge)
    if (flag):
        if not(flag_no_warn):
            print("Collapsing edge (" + half_edge.EndpointsStr(",") +\
                    ") will change the mesh topology.")
            print("  Vertices (" + half_edge.EndpointsStr(", ") +\
                    ", " + str(ivC) + ") form a triangle hole.")

        return_flag = False

    if (mesh.IsInteriorEdgeWithBoundaryVertices(ihalf_edge)):
        if not(flag_no_warn):
            print("Collapsing edge (" + half_edge.EndpointsStr(",") +\
                    ") merges two non-adjacent boundary vertices.")

        return_flag = False;

    if mesh.IsIsolatedTriangle(icell):
        if not(flag_no_warn):
            print("Collapsing edge(" + half_edge.EndpointsStr(",") +\
                    f") will delete isolated cell {icell}.")

        return_flag = False

    if mesh.IsInTetrahedron(icell):
        if not(flag_no_warn):
            print("Collapsing edge (" + half_edge.EndpointsStr(",") +\
                    ") will collapse a tetrahedron.")

        return_flag = False

    return return_flag

## Print a warning message if splitting cell at diagonal
#    (half_edgeA.FromVertex(), half_edgeB.FromVertex())
#    will change the mesh topology.
#  - Return True if split does not change manifold topology.
def check_split_cell(mesh, ihalf_edgeA, ihalf_edgeB, flag_no_warn):
    half_edgeA = mesh.HalfEdge(ihalf_edgeA)
    half_edgeB = mesh.HalfEdge(ihalf_edgeB)
    vA = half_edgeA.FromVertex()
    vB = half_edgeB.FromVertex()
    ivA = vA.Index()
    ivB = vB.Index()
    icell = half_edgeA.CellIndex()
    half_edgeC = mesh.FindEdge(vA, vB)

    flag_cell_edge = False
    return_flag = True

    if (mesh.IsIllegalSplitCell(ihalf_edgeA, ihalf_edgeB)):
        if (vA is half_edgeB.ToVertex()) or (vB is half_edgeA.FromVertex()):
            flag_cell_edge = True

        if not(flag_no_warn):
            if flag_cell_edge:
                print(f"({ivA},{ivB}) is a cell edge, not a cell diagonal.")
            else:
                print(f"Illegal split of cell {icell} with diagonal ({ivA}{ivB}).")
        return_flag = False;

    if not(half_edgeC is None) and not(flag_cell_edge):
        if not(flag_no_warn):
            sys.stdout.write\
                (f"Splitting cell {icell} with diagonal ({ivA},{ivB})");
            sys.stdout.write\
                (" creates an edge incident on three or more cells.\n")
        return_flag = False

    return return_flag


## Print a warning if joining cells separated by half_edge is illegal.
#  - Return true if join is legal.
def check_join_cell(mesh, ihalf_edge, flag_no_warn):
    half_edge = mesh.HalfEdge(ihalf_edge)
    TWO = 2
    return_flag = True

    if (mesh.IsIllegalJoinTwoCells(ihalf_edge)):
        half_edgeX = half_edge.NextHalfEdgeAroundEdge()

        if not(flag_no_warn):
            ivfrom = half_edge.FromVertexIndex()
            ivto = half_edge.ToVertexIndex()
            if (half_edge.IsBoundary()):
                print("Only one cell contains edge (" +\
                    half_edge.EndpointsStr(",") + ").")
            elif not(mesh.IsVertexIncidentOnMoreThanTwoEdges(ivfrom)):
                print(f"Half edge endpoint {ivfrom} is incident on only two edges.")
            elif not(mesh.IsVertexIncidentOnMoreThanTwoEdges(ivto)):
                print(f"Half edge endpoint {ivto} is incident on only two edges.")
            elif not(half_edge is half_edgeX.NextHalfEdgeAroundEdge()):
                print("More than two cells are incident on edge (" +\
                        half_edge.EndpointsStr(",") + ").")
            else:
                icell = half_edge.CellIndex()
                icellX = half_edgeX.CellIndex()
                num_shared_vertices =\
                    mesh.CountNumVerticesSharedByTwoCells(icell, icellX)
                print("Join of two cells " + str(icell) + " and "\
                        + str(icellX) + " incident on edge (" +\
                        half_edge.EndpointsStr(",") + ") is illegal.")
                if (num_shared_vertices > TWO):
                    print(f"  Cells {icell} and {icellX} share {num_shared_vertices} vertices.")

        return_flag = False

    return return_flag


## Return True and print warning message if data set is large
def reduce_checks_on_large_datasets\
    (mesh, flag_no_warn, large_data_num_cells):

    num_cells = mesh.NumCells()
    if (num_cells >= large_data_num_cells):
        if not(flag_no_warn):
            print(f"Warning: Large data set with {num_cells} cells.")
            print("  Reducing checks (using -flag_reduce_checks.)")

        return True
    else:
        return False


# *** Init/parse/print functions. ***

## Initialize global flags.
def InitFlags():
    global input_filename, output_filename
    global flag_silent, flag_terse, flag_no_warn, flag_time
    global flag_collapse_short_edges, flag_split_long_edges
    global flag_split_cells, flag_join_cells
    global flag_triangulate_cells
    global flag_allow_non_manifold, flag_fail_on_non_manifold
    global flag_reduce_checks

    # Initialize
    input_filename = None
    output_filename = None
    flag_silent = False
    flag_terse = False
    flag_no_warn = False
    flag_time = False
    flag_collapse_short_edges = False
    flag_split_cells = False
    flag_join_cells = False
    flag_split_long_edges = False
    flag_triangulate_cells = False
    flag_allow_non_manifold = False
    flag_fail_on_non_manifold = False
    flag_reduce_checks = False


def parse_command_line(argv):
    global input_filename, output_filename
    global flag_silent, flag_terse, flag_no_warn, flag_time
    global flag_collapse_short_edges, flag_split_long_edges
    global flag_split_cells, flag_join_cells
    global flag_triangulate_cells
    global flag_allow_non_manifold, flag_fail_on_non_manifold
    global flag_reduce_checks

    iarg = 1
    while (iarg < len(argv) and argv[iarg][0] == '-'):
        s = argv[iarg]
        if (s == "-collapse_short_edges"):
            flag_collapse_short_edges = True
        elif (s == "-split_long_edges"):
            flag_split_long_edges = True
        elif (s == "-split_cells"):
            flag_split_cells = True
        elif (s == "-join_cells"):
            flag_join_cells = True
        elif ((s == "-triangulate_cells") or (s == "-triangulate")):
            flag_triangulate_cells = True
        elif (s == "-allow_non_manifold"):
            flag_allow_non_manifold = True
        elif (s == "-fail_on_non_manifold"):
            flag_fail_on_non_manifold = True
        elif (s == "-reduce_checks"):
            flag_reduce_checks = True
        elif (s == "-s"):
            flag_silent = True
            flag_terse = True
        elif (s == "-terse"):
            flag_terse = True
        elif (s == "-no_warn"):
            flag_no_warn = True
        elif (s == "-time"):
            flag_time = True
        elif (s == "-h"):
            help()
        else:
            sys.stderr.write("Usage error. Option " + s + " is undefined.\n")
            usage_error();

        iarg = iarg+1

    if (iarg >= len(argv) or iarg+2 < len(argv)):
        usage_error()

    input_filename = argv[iarg]

    if (iarg+1 < len(argv)):
        output_filename = argv[iarg+1]

    if not(flag_collapse_short_edges or flag_split_long_edges or\
            flag_split_cells or flag_join_cells or\
            flag_triangulate_cells):
        print("No edit operation specified.")
        print("Specify -collapse_short_edges, -split_long_edges, -split_cells")
        print("  -join_cells or -triangulate_cells.")
        usage_error()

## Print timeX in seconds.
def print_time(label, timeX):
    print(label + "{:.4f}".format(timeX) + " seconds")

## Print mesh information, such as number of vertices, edges, cells
#    min and max edge lengths and cell angles.
def print_mesh_info(mesh):
    FIVE = 5
    num_vertices = mesh.NumVertices()
    num_edges = mesh.CountNumEdges()
    num_boundary_edges = mesh.CountNumBoundaryEdges()
    num_cells = mesh.NumCells()
    num_triangles = mesh.CountNumTriangles()
    num_quads = mesh.CountNumQuads()
    num_large_cells = mesh.CountNumCellsOfSizeGE(FIVE)

    minL_squared, maxL_squared, ihmin, ihmax =\
        compute_min_max_edge_length_squared(mesh)
    min_ratio_squared, icmin, Lmin, Lmax, ihmin, ihmax =\
        compute_min_cell_edge_length_ratio_squared(mesh)
    angle_info = compute_angle_info(mesh, [1, 5, 10], [175, 170])

    flag_manifoldV, flag_manifoldE, iv, ie = mesh.CheckManifold()
    is_oriented, iv = mesh.CheckOrientation()

    print("Number of vertices: ", num_vertices)
    print("Number of mesh edges: ", num_edges)
    print("Number of boundary edges: ", num_boundary_edges)
    print("Number of mesh cells: ", num_cells)
    print("  Number of mesh triangles: ", num_triangles)
    print("  Number of mesh quadrilaterals: ", num_quads)
    print("  Number of mesh cells with > 4 vertices: ", num_large_cells)
    print(f"Min edge length: {sqrt(minL_squared):.4f}")
    print(f"Max edge length: {sqrt(maxL_squared):.4f}")
    print(f"Min cell edge length ratio: {sqrt(min_ratio_squared):.4f}")
    min_angle = math.degrees(acos(angle_info.cos_min_angle))
    print(f"Minimum cell angle: {min_angle:.4f}")
    for j in range(len(angle_info.small_angle_bounds)):
        A = angle_info.small_angle_bounds[j]
        num_cells = angle_info.num_cells_with_angle_le_small[j]
        print(f"  Num cells with angles <= {A}: {num_cells}.")

    max_angle = math.degrees(acos(angle_info.cos_max_angle))
    print(f"Maximum cell angle: {max_angle:.4f}")
    for j in range(len(angle_info.large_angle_bounds)):
        A = angle_info.large_angle_bounds[j]
        num_cells = angle_info.num_cells_with_angle_ge_large[j]
        print(f"  Num cells with angles >= {A}: {num_cells}.")

    if (flag_manifoldV and flag_manifoldE and is_oriented):
        print("Mesh is an oriented manifold.")
    else:
        print("Mesh is non-manifold or has inconsistent cell orientations.")


def usage_msg(out):
    out.write("Usage: python decimate_mesh.py [OPTIONS] <input filename> [<output_filename>]\n")
    out.write("Options:\n")
    out.write("  [-collapse_short_edges] [-split_long_edges]\n")
    out.write("  [-split_cells] [-join_cells]\n")
    out.write("  [-triangulate_cells | -triangulate]\n")
    out.write("  [-allow_non_manifold] [-fail_on_non_manifold]\n")
    out.write("  [-s | -terse] [-no_warn] [-reduce_checks] [-time] [-h]\n")

def usage_error():
    usage_msg(sys.stderr)
    sys.stderr.flush()
    exit(-1)

def help():
    usage_msg(sys.stdout)
    print()
    print("decimate_mesh.py -- Decimate_mesh.")
    print("    Collapse/split/join mesh edges or cells.")
    print()
    print("Options:")
    print("-collapse_short_edges: Attempt to collapse short edge in each cell.")
    print("-split_long_edges: Split longest edge in each cell.")
    print("-split_cells:      Split cells at vertices with largest angles.")
    print("-join_cells:       Join cells sharing long edges.")
    print("-triangulate_cells | -triangulate:")
    print("      Triangulate cells from vertex with largest angle.")
    print("-terse:   Terse output. Suppress messages output after each")
    print("      collapse/join/split iteration.")
    print("  Does not suppress warning messages at each iteration.")
    print("  Does not suppress final mesh information.")
    print("-s:       Silent. Output only warnings and error messages.")
    print("-no_warn: Do not ouptut non-manifold or inconsistent orientation warnings.")
    print("-time:    Report run time.")
    print("-h:       Output this help message and exit.")

    exit(0)

if __name__ == '__main__':
    main(sys.argv)
