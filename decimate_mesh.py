## \file decimate_mesh.py
#  Some simple mesh decimation routines.
#  Use data structure HALF_EDGE_MESH_DCMT_BASE (DCMT = decimate).

import math
from math import sqrt
from math import acos
import sys
from time import time, localtime, strftime

import half_edge_mesh
import half_edge_mesh_DCMT
import half_edge_mesh_IO
from half_edge_mesh_DCMT import HALF_EDGE_MESH_DCMT_BASE
from half_edge_mesh_DCMT\
    import VERTEX_DCMT_BASE, HALF_EDGE_DCMT_BASE, CELL_DCMT_BASE


def main(argv):

    global input_filename, output_filename
    global flag_silent, flag_terse, flag_no_warn, flag_time
    global flag_collapse_edges, flag_collapse_short_edges
    global flag_split_cells, flag_split_all_cells
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

    mesh = HALF_EDGE_MESH_DCMT_BASE\
        (VERTEX_DCMT_BASE,HALF_EDGE_DCMT_BASE,CELL_DCMT_BASE)

    half_edge_mesh_IO.open_and_read_off_file(input_filename, mesh)

    time2 = time()

    try:
        if not(flag_reduce_checks):
            flag_reduce_checks =\
                reduce_checks_on_large_datasets\
                    (mesh, flag_no_warn, LARGE_DATA_NUM_CELLS)

        if (flag_collapse_edges):
            prompt_and_collapse_edges(mesh, flag_terse, flag_no_warn)

        if (flag_split_cells):
            prompt_and_split_cells(mesh, flag_terse, flag_no_warn)

        if (flag_collapse_short_edges):
            collapse_shortest_edge_in_each_cell\
                (mesh, flag_terse, flag_no_warn)

        if (flag_split_all_cells):
            split_all_cells(mesh, flag_terse, flag_no_warn)


        passed_check = check_mesh(mesh, flag_silent and flag_no_warn)

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
def collapse_edge(mesh, half_edge, flag_terse, flag_no_warn, flag_check):
    global flag_allow_non_manifold

    flag = check_edge_collapse(mesh, half_edge, flag_no_warn)

    if mesh.IsIllegalEdgeCollapseH(half_edge):
        return

    if flag or flag_allow_non_manifold:

        if not(flag_terse):
            print("Collapsing edge (" + half_edge.EndpointsStr(",") + ").")

        vnew = mesh.CollapseEdge(half_edge.Index())
        if (vnew is None):
            print("Skipped illegal collapse of edge (" +\
                    half_edge.EndpointsStr(",") + ").")

        if flag_check:
            check_mesh(mesh, flag_no_warn)
    else:
        if not(flag_terse):
            print("Skipped collapse of edge (" +\
                    half_edge.EndpointsStr(",") + ").")


## Prompt and collapse edges.
def prompt_and_collapse_edges(mesh, flag_terse, flag_no_warn):

    while (True):
        half_edge0 = prompt_for_mesh_edge(mesh, False)

        if (half_edge0 is None):
            # End.
            print()
            return

        collapse_edge(mesh, half_edge0, flag_terse, flag_no_warn, True)

        print()

## Collapse shortest cell edge.
def collapse_shortest_cell_edge\
    (mesh, icell, flag_terse, flag_no_warn, flag_check):

    cell = mesh.Cell(icell)
    if (cell is None):
        return

    minL, maxL, ihalf_edge_min, ihalf_edge_max =\
        cell.ComputeMinMaxEdgeLengthSquared()

    half_edge_min = mesh.HalfEdge(ihalf_edge_min)
    collapse_edge(mesh, half_edge_min, flag_terse, flag_no_warn, flag_check)


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


# *** Split cell routines.

## Split cell with diagonal (half_edgeA.FromVertex(), half_edgeB.FromVertex())
#  - Returns split edge.
#  - Returns None if split fails.
def split_cell(mesh, half_edgeA, half_edgeB,\
                flag_terse, flag_no_warn, flag_check):
    ihalf_edgeA = half_edgeA.Index()
    ihalf_edgeB = half_edgeB.Index()
    vA = half_edgeA.FromVertex()
    vB = half_edgeB.FromVertex()
    ivA = vA.Index()
    ivB = vB.Index()
    icell = half_edgeA.CellIndex()

    flag = check_split_cell(mesh, half_edgeA, half_edgeB, flag_no_warn)

    if (mesh.IsIllegalSplitCell(half_edgeA, half_edgeB)):
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


## Get at most max_num cells with more than three vertices.
def get_cells_with_more_than_three_vertices(mesh, max_num):
    THREE = 3

    cell_list = []
    if (max_num < 1):
        return cell_list

    for icell in range(0,mesh.MaxCellIndex()):
        cell = mesh.Cell(icell)
        if (icell is None):
            continue

        if (cell.NumVertices() > THREE):
            cell_list.append(icell)

        if (len(cell_list) >= max_num):
            # Stop getting more cells.
            return cell_list

    return cell_list


## Prompt and split cells.
def prompt_and_split_cells(mesh, flag_terse, flag_no_warn):
    THREE = 3
    MAX_NUM = 10

    cell_list = get_cells_with_more_than_three_vertices(mesh, MAX_NUM)

    if (len(cell_list) == 0):
        if not(flag_no_warn):
            print("All cells are triangle. No cells can be split.")
        return

    print_cells_with_more_than_three_vertices(MAX_NUM, cell_list)
    print()

    while(True):

        icell = int(input("Enter cell (-1 to end): "))

        if (icell < 0):
            return

        if (icell > mesh.MaxCellIndex()):
            print(f"No cell has index {icell}.")
            print(f"  Maximum cell index: {mesh.MaxCellIndex()}")
            continue

        cell = mesh.Cell(icell)
        if (cell is None):
            print(f"No cell has index {icell}.")
            print()
            continue

        if (cell.NumVertices() <= THREE):
            print(f"Cell {icell} has fewer than four vertices and cannot be split.")
            print()
            continue

        half_edge = cell.HalfEdge()
        sys.stdout.write(f"Vertices in cell {icell}:")
        for k in range(0,cell.NumVertices()):
            sys.stdout.write(f"  {half_edge.FromVertexIndex()}")
            half_edge = half_edge.NextHalfEdgeInCell()
        sys.stdout.write("\n")

        while True:
            input_list = input("Enter two distinct vertex indices (-1 to end): ").split()
            if (len(input_list) >= 2):
                ivA = int(input_list[0])
                ivB = int(input_list[1])
                break
            elif (len(input_list) == 1):
                ivA = int(input_list[0])
                if (ivA < 0):
                    return
                ivB = int(input("Enter a second vertex index (-1 to end):"))
                break

        if (ivA < 0) or (ivB < 0):
            return

        if (ivA == ivB):
            print()
            print("Vertices are not distinct. Start again.")
            print()
            continue

        half_edgeA = None
        half_edgeB = None
        half_edge = cell.HalfEdge()
        for k in range(0, cell.NumVertices()):
            if (half_edge.FromVertexIndex() == ivA):
                half_edgeA = half_edge
            if (half_edge.FromVertexIndex() == ivB):
                half_edgeB = half_edge

            half_edge = half_edge.NextHalfEdgeInCell()

        if (half_edgeA is None) or (half_edgeB is None):
            print()
            print(f"Vertices are not in cell {icell}.")
            print("Start again.")
            print()
            continue

        if (half_edgeA.ToVertexIndex() == ivB) or\
            (half_edgeB.ToVertexIndex() == ivA):
            print ()
            print(f"({ivA},{ivB}) is a cell edge, not a cell diagaonal.")
            print("  Vertices must not be adjacent.")
            print("Start again.")
            print()
            continue

        split_cell\
            (mesh, half_edgeA, half_edgeB, flag_terse, flag_no_warn, True)

        print()


## Split cell at largest angle.
#  - Split cell at vertex forming the largest angle.
#  - Split as many times as necessary to triangulate.
def split_cell_at_largest_angle\
        (mesh, cell, flag_terse, flag_no_warn, flag_check):

    cos_minA, cos_maxA, ihalf_edge_min, ihalf_edge_max =\
        cell.ComputeCosMinMaxAngle()

    half_edgeA = mesh.HalfEdge(ihalf_edge_max)

    while not(half_edgeA.Cell().IsTriangle()):
        half_edgeB = half_edgeA.PrevHalfEdgeInCell().PrevHalfEdgeInCell()
        vA = half_edgeA.FromVertex()
        vB = half_edgeB.FromVertex()
        ivA = vA.Index()
        ivB = vB.Index()

        split_edge = split_cell(mesh, half_edgeA, half_edgeB,\
                                flag_terse, flag_no_warn, flag_check)
        if (split_edge is None):
            # Cannot split cell at largest angle.
            return

        if (flag_check):
            check_mesh(mesh, flag_no_warn)

        # Get largest angle in remaining cell.
        cellA = split_edge.Cell()
        cos_minA, cos_maxA, ihalf_edge_min, ihalf_edge_max =\
            cell.ComputeCosMinMaxAngle()
        half_edgeA = mesh.HalfEdge(ihalf_edge_max)


## Split all cells.
def split_all_cells(mesh, flag_terse, flag_no_warn):
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
            (mesh, cell, flag_terse, flag_no_warn, flag_check)
        kount = kount+1

        if (flag_reduce_checks):
            if (kount == n/2):
                check_mesh(mesh, flag_no_warn)


# *** Check routines ***

def check_oriented_manifold(mesh, flag_no_warn):
    flag_non_manifold_vertex, flag_non_manifold_edge, iv, ihalf_edgeA =\
        mesh.CheckManifold()


    if flag_non_manifold_edge:
        if not(flag_no_warn):
            half_edgeA = mesh.HalfEdge(ihalf_edgeA)
            sys.stderr.write("Warning: Non-manifold edge (" +\
                                half_edgeA.EndpointsStr(",") + ").\n")

        # Non-manifold edge automatically implies inconsistent orientations.
        return False

    flag_orientation, ihalf_edgeB = mesh.CheckOrientation()

    if flag_orientation:
        if flag_non_manifold_vertex:
            if not(flag_no_warn):
                sys.stderr.write(f"Warning: Non-manifold vertex {iv}.")

            return False
    else:
        if flag_non_manifold_vertex:
            if not(flag_no_warn):
                sys.stderr.write\
                    (f"Warning: Non-manifold vertex or inconsistent orientations in cells incident on vertex {iv}.\n")
        else:
            sys.stderr.write\
                ("Warning: Inconsistent orientation of cells incident on edge (" +\
                half_edgeB.EndpointsStr(",") + ").")

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
            if not(flag_no_warn):
                sys.stderr.write\
                    ("Detected non-manifold or inconsistent orientations")
                sys.stderr.write("Exiting.")

            exit(-1)

        return flag_oriented_manifold

    else:
        return True


## Print a warning message if collapsing half_edge is illegal or
#    will change mesh topology.
#  - Return true if collapse is not illegal and does not change
#    mesh topology.
def check_edge_collapse(mesh, half_edge, flag_no_warn):
    icell = half_edge.CellIndex()
    return_flag = True

    if (mesh.IsIllegalEdgeCollapseH(half_edge)):
        if not(flag_no_warn):
            print("Collapse of edge  (" + half_edge.EndpointsStr(",") +\
                    ") is illegal.")
            print("  Some cell contains vertices " +\
                    half_edge.EndpointsStr(" and ") + " but not edge (" +\
                    half_edge.EndpointsStr(",") + ").")

        return_flag = False

    flag, ivC = mesh.FindTriangleHole(half_edge)
    if (flag):
        if not(flag_no_warn):
            print("Collapsing edge (" + half_edge.EndpointsStr(",") +\
                    ") will change the mesh topology.")
            print("  Vertices (" + half_edge.EndpointsStr(", ") +\
                    ", " + str(ivC) + ") form a triangle hole.")

        return_flag = False

    if not(half_edge.IsBoundary()):
        if half_edge.FromVertex().IsBoundary() and\
            half_edge.ToVertex().IsBoundary():

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
#  - Return true if split does not change manifold topology.
def check_split_cell(mesh, half_edgeA, half_edgeB, flag_no_warn):
    vA = half_edgeA.FromVertex()
    vB = half_edgeB.FromVertex()
    ivA = vA.Index()
    ivB = vB.Index()
    icell = half_edgeA.CellIndex()
    half_edgeC = mesh.FindEdge(vA, vB)

    flag_cell_edge = False
    return_flag = True

    if (mesh.IsIllegalSplitCell(half_edgeA, half_edgeB)):
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
                (" creates an edge incident on three or nmore cells.\n")
        return_flag = False

    return return_flag


## Return true and print warning message if data set is large
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


# *** Init/parse/print/prompt functions. ***

## Initialize global flags.
def InitFlags():
        global input_filename, output_filename
        global flag_silent, flag_terse, flag_no_warn, flag_time
        global flag_collapse_edges, flag_collapse_short_edges
        global flag_split_cells, flag_split_all_cells
        global flag_allow_non_manifold, flag_fail_on_non_manifold
        global flag_reduce_checks

        # Initialize
        input_filename = None
        output_filename = None
        flag_silent = False
        flag_terse = False
        flag_no_warn = False
        flag_time = False
        flag_collapse_edges = False
        flag_collapse_short_edges = False
        flag_split_cells = False
        flag_split_all_cells = False
        flag_allow_non_manifold = False
        flag_fail_on_non_manifold = False
        flag_reduce_checks = False


def parse_command_line(argv):
    global input_filename, output_filename
    global flag_silent, flag_terse, flag_no_warn, flag_time
    global flag_collapse_edges, flag_collapse_short_edges
    global flag_split_cells, flag_split_all_cells
    global flag_allow_non_manifold, flag_fail_on_non_manifold
    global flag_reduce_checks

    iarg = 1
    while (iarg < len(argv) and argv[iarg][0] == '-'):
        s = argv[iarg]
        if (s == "-collapse_edges"):
            flag_collapse_edges = True
        elif (s == "-collapse_short_edges"):
            flag_collapse_short_edges = True
        elif (s == "-split_cells"):
            flag_split_cells = True
        elif (s == "-split_all_cells"):
            flag_split_all_cells = True
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


## Prompt for mesh edge.
#  - Return None if user enters -1.
def prompt_for_mesh_edge(mesh, flag_only_internal):

    while (True):
        iv0 = int(input("Enter vertex (-1 to end): "))
        if (iv0 < 0):
            return None

        flag, error_msg = mesh.CheckVertexIndex(iv0)
        if not(flag):
            if not(error_msg is None):
                print(error_msg)
            continue;

        v0 = mesh.Vertex(iv0)

        if (v0.NumHalfEdgesFrom() == 0):
            print(f"Vertex {iv0} is not incident on any cell.")
            continue

        if (flag_only_internal):
            print("NOT YET IMPLEMENTED")
        else:
            sys.stdout.write(f"Half edges from {iv0}:")
            for k in range(0,v0.NumHalfEdgesFrom()):
                half_edge = v0.KthHalfEdgeFrom(k)
                sys.stdout.write("  (" + half_edge.EndpointsStr(",") + ")")
            print()

        iv1 = int(input(f"Enter vertex adjacent to {iv0} (-1 to end): "))

        if (iv1 < 0):
            return None

        flag, error_msg = mesh.CheckVertexIndex(iv1)
        if not(flag):
            if not(error_msg is None):
                print(error_msg)
            continue;

        half_edge0 = v0.FindIncidentHalfEdge(iv1)
        if (half_edge0 is None):
            print(f"Mesh does not have a half edge ({iv0},{iv1}).")
            continue;

        if (flag_only_internal and half_edge0.IsBoundary()):
            print("Half edge (" + half_edge0.EndpointsStr(",") +\
                    ") is a boundary half edge.")
            print("Start again.")
            print()
            continue

        return half_edge0


## Print cells with more than three vertices.
def print_cells_with_more_than_three_vertices(max_num, cell_list):
    sys.stdout.write("Cells with more than three vertices")
    if (len(cell_list) >= max_num):
        sys.stdout.write(" (partial list)")
    sys.stdout.write(": ")
    for i in range(0,len(cell_list)):
        sys.stdout.write(f"  {cell_list[i]}")
    sys.stdout.write("\n")


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
        mesh.ComputeMinMaxEdgeLengthSquared()
    min_ratio_squared, icmin, Lmin, Lmax, ihmin, ihmax =\
        mesh.ComputeMinCellEdgeLengthRatioSquared()
    cos_minA, cos_maxA, ihmin, ihmax = mesh.ComputeCosMinMaxAngle()

    flag_non_manifoldV, flag_non_manifoldE, iv, ie = mesh.CheckManifold()
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
    print(f"Minimum cell angle: {math.degrees(acos(cos_minA)):.4f}")
    print(f"Minimum cell angle: {math.degrees(acos(cos_maxA)):.4f}")

    if (not(flag_non_manifoldV) and not(flag_non_manifoldE) and is_oriented):
        print("Mesh is an oriented manifold.")
    else:
        print("Mesh is non-manifold or has inconsistent cell orientations.")


def usage_msg(out):
    out.write("Usage: python3 decimate_mesh.py [OPTIONS] <input filename> [<output_filename>]\n")
    out.write("Options:\n")
    out.write("  [-collapse_edges] [-collapse_short_edges]\n")
    out.write("  [-split_cells] [-split_all_cells]\n")
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
    print("-collapse_edges:  Prompt and collapse edges.")
    print("-collapse_short_edges:  Attempt to collapse short edge in each cell.")
    print("-split_cells:     Prompt and split cells across diagonals.")
    print("-split_all_cells: Attempt to split all cells.")
    print("-terse:   Terse output. Suppress messages output after each")
    print("    collapse/join/split iteration.")
    print("  Does not suppress warning messages at each iteration.")
    print("  Does not suppress final mesh information.")
    print("-s:       Silent. Output only wanrings and error messages.")
    print("-no_warn: Do not ouptut non-manifold or inconsistent orientation warnings.")
    print("-time:    Report run time.")
    print("-h:       Output this help message and exit.")

    exit(0)

if __name__ == '__main__':
    main(sys.argv)
