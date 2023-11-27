## \file meshinfo.py
#  Print mesh information, including number of vertices, edges, cells,
#  and minimum and maximum edge lengths and cell angles.

import math
from math import sqrt
from math import acos
import sys
import half_edge_mesh
import half_edge_mesh_IO
import half_edge_mesh_compute
from half_edge_mesh import HALF_EDGE_MESH_BASE
from half_edge_mesh_elements\
    import HMESH_VERTEX_BASE, HMESH_HALF_EDGE_BASE, HMESH_CELL_BASE
from half_edge_mesh_compute import print_coordP
from half_edge_mesh_compute\
    import compute_min_max_edge_length_squared,\
        compute_min_cell_edge_length_ratio_squared,\
        compute_angle_info

def main(argv):

    global input_filename, flag_more_info


    # Initialize
    input_filename = None
    flag_more_info = False

    parse_command_line(sys.argv)

    mesh = HALF_EDGE_MESH_BASE\
        (HMESH_VERTEX_BASE, HMESH_HALF_EDGE_BASE, HMESH_CELL_BASE)

    half_edge_mesh_IO.open_and_read_off_file(input_filename, mesh)

    try:
        print_mesh_size(mesh, flag_more_info)
        print_min_max_edge_lengths(mesh, flag_more_info)
        print_min_cell_edge_length_ratio(mesh, flag_more_info)
        print_min_max_angles(mesh, flag_more_info)
        print_manifold_info(mesh, flag_more_info)

    except Exception as e:
        print(e)
        sys.stderr.write("Exiting.")
        exit(-1)


# *** Print Mesh Information ***

## Print number of vertices, number of edges, number of cells, etc.
def print_mesh_size(mesh, flag_more_info):
    num_vertices = mesh.NumVertices()
    num_isolated_vertices = mesh.CountNumIsolatedVertices()
    num_edges = mesh.CountNumEdges()
    num_boundary_edges = mesh.CountNumBoundaryEdges()
    num_cells = mesh.NumCells()

    print("Number of mesh vertices: ", num_vertices-num_isolated_vertices)
    if (num_isolated_vertices > 0):
        print("Total number of vertices in input file: ", num_vertices)

    print("Number of mesh edges: ", num_edges)
    print("Number of boundary mesh edges: ", num_boundary_edges)
    print("Number of mesh cells: ", num_cells)

    if (flag_more_info):
        num_triangles = mesh.CountNumTriangles()
        num_quads = mesh.CountNumQuads()
        num_pentagons = mesh.CountNumPentagons()
        num_large_cells = mesh.CountNumCellsOfSizeGE(6)

        print("  Number of mesh triangles: ", num_triangles)
        print("  Number of mesh quadrilaterals: ", num_quads)
        if (num_pentagons > 0):
            print("  Number of mesh pentagons: ", num_pentagons)
            print("  Number of cells with > 5 vertices: ", num_large_cells)
        else:
            print("  Number of cells with > 4 vertices: ", num_large_cells)


## Print minimum and maximum edge lengths.
def print_min_max_edge_lengths(mesh, flag_more_info):

    min_Lsquared, max_Lsquared, ihalf_edge_min, ihalf_edge_max =\
        compute_min_max_edge_length_squared(mesh)

    print(f"Min edge length: {sqrt(min_Lsquared):.4f}")
    if (flag_more_info):
        half_edge_min = mesh.HalfEdge(ihalf_edge_min)
        print("  Min length = length of edge (" +\
                half_edge_min.EndpointsStr(",") +\
                ") in cell: " + str(half_edge_min.CellIndex()) + ".")

    print(f"Max edge length: {sqrt(max_Lsquared):.4f}")
    if (flag_more_info):
        half_edge_max = mesh.HalfEdge(ihalf_edge_max)
        print("  Max length = length of edge (" +\
                half_edge_max.EndpointsStr(",") +\
                ") in cell: " + str(half_edge_max.CellIndex()) + ".")


## Print minimum ratio of shortest to longest edge in a cell.
def print_min_cell_edge_length_ratio(mesh, flag_more_info):

    min_ratio_squared, icell, Lmin_sq, Lmax_sq, ihalfE_min, ihalfE_max = \
        compute_min_cell_edge_length_ratio_squared(mesh)

    min_ratio = sqrt(min_ratio_squared)
    print(f"Min cell edge length ratio: {min_ratio:.4f}");
    if (flag_more_info):
        half_edge_min = mesh.HalfEdge(ihalfE_min)
        half_edge_max = mesh.HalfEdge(ihalfE_max)
        Lmin = sqrt(Lmin_sq)
        Lmax = sqrt(Lmax_sq)
        print(f"  In cell: {icell}.")
        print(f"  Min cell edge length: {Lmin:.4f}.  Edge: (" +\
                half_edge_min.EndpointsStr(",") + ").")
        print(f"  Max cell edge length: {Lmax:.4f}.  Edge: (" +\
                half_edge_max.EndpointsStr(",") + ").")


def print_min_max_angles(mesh, flag_more_info):

    small_angle_bounds = [1, 5, 10]
    large_angle_bounds = [175, 170]
    angle_info = compute_angle_info\
                    (mesh, small_angle_bounds, large_angle_bounds)

    print(f"Min angle: {math.degrees(acos(angle_info.cos_min_angle)):.4f}")
    if (flag_more_info):
        half_edge_minA = mesh.HalfEdge(angle_info.imin)
        iv = half_edge_minA.FromVertexIndex()
        icell = half_edge_minA.CellIndex()
        print(f"  At vertex {iv} in cell {icell}.")

    for j in range(len(angle_info.small_angle_bounds)):
        A = angle_info.small_angle_bounds[j]
        num_cells = angle_info.num_cells_with_angle_le_small[j]
        print(f"Num cells with angles <= {A}: {num_cells}.")

    print(f"Max angle: {math.degrees(acos(angle_info.cos_max_angle)):.4f}")
    if (flag_more_info):
        half_edge_maxA = mesh.HalfEdge(angle_info.imax)
        iv = half_edge_maxA.FromVertexIndex()
        icell = half_edge_maxA.CellIndex()
        print(f"  At vertex {iv} in cell {icell}.")

    for j in range(len(angle_info.large_angle_bounds)):
        A = angle_info.large_angle_bounds[j]
        num_cells = angle_info.num_cells_with_angle_ge_large[j]
        print(F"Num cells with angles >= {A}: {num_cells}.")


## Print manifold and orientation information.
def print_manifold_info(mesh, flag_more_info):

    flag_manifoldV, flag_manifoldE, ivA, ihalf_edgeA =\
        mesh.CheckManifold()
    flag_orientation, ihalf_edgeB = mesh.CheckOrientation()

    if (flag_manifoldE and flag_manifoldV and flag_orientation):
        print("Mesh is an oriented manifold.")
    elif not(flag_manifoldE):
        half_edgeA = mesh.HalfEdge(ihalf_edgeA)
        print("Mesh has a non-manifold edge (" +\
                half_edgeA.EndpointsStr(",") + ").")
    elif not(flag_manifoldV) and flag_orientation:
        print(f"Mesh has a non-manifold vertex {ivA}.")
    elif not(flag_manifoldV):
        print(f"Non-manifold or inconsistent orientation at vertex {ivA}.")
    elif not(flag_orientation):
        half_edgeB = mesh.HalfEdge(ihalf_edgeB)
        half_edgeBX = half_edgeB.NextHalfEdgeAroundEdge()
        icellB = half_edgeB.CellIndex();
        icellBX = half_edgeBX.CellIndex();
        print("Mesh is a manifold.");
        print(f"Inconsistent orientation of cells {icellB} and {icellBX}.")


# *** Parse/output functions ***

def parse_command_line(argv):
    global input_filename, flag_more_info, flag_midpoints, flag_centroids

    iarg = 1
    while (iarg < len(argv) and argv[iarg][0] == '-'):
        s = argv[iarg]
        if (s == "-more"):
            flag_more_info = True
        elif ((s == "-h") or (s == "-help")):
            help()
        else:
            sys.stderr.write("Usage error. Option " + s + " is undefined.\n")
            usage_error();

        iarg = iarg+1

    if (iarg >= len(argv) or iarg+1 < len(argv)):
        usage_error()

    input_filename = argv[iarg]


def usage_msg(out):
    out.write("Usage: python3 meshinfo [-more] [-h] <input filename>")

def usage_error():
    usage_msg(sys.stderr)
    sys.stderr.flush()
    exit(-1)

def help():
    usage_msg(sys.stdout)
    print("\n\nmeshinfo.py -- Print mesh information.")
    print("")
    print("Options:")
    print("-more:      Print additional information.")
    print("-h:         Output this help message and exit.")

    exit(0)

if __name__ == '__main__':
    main(sys.argv)
