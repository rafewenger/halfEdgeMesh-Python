## \file test_hmesh_split_edges.py
#    Test splitting edges in half edge mesh data structure.
#
#  - Reads a .off file into HALF_EDGE_MESH_BASE.
#  - Checks mesh, manifold and orientations properties of the mesh.
#  - Splits each mesh edge.
#  - Rechecks mesh, manifold and orientation properties of the mesh.
#  - Writes the mesh to a .off file.


from time import time, localtime, strftime
import sys
import half_edge_mesh
import half_edge_mesh_IO
from half_edge_mesh_elements\
    import HMESH_VERTEX_BASE, HMESH_HALF_EDGE_BASE, HMESH_CELL_BASE
from half_edge_mesh import HALF_EDGE_MESH_BASE

def main(argv):

    global input_filename, output_filename, flag_no_warn
    global flag_silent, flag_verbose
    global flag_time

    begin_time = time()

    # Initialize
    input_filename = None
    output_filename = None
    flag_silent = False
    flag_verbose = False
    flag_no_warn = False
    flag_time = False

    parse_command_line(sys.argv)

    mesh = HALF_EDGE_MESH_BASE\
        (HMESH_VERTEX_BASE,HMESH_HALF_EDGE_BASE,HMESH_CELL_BASE)

    half_edge_mesh_IO.open_and_read_off_file(input_filename, mesh)

    time2 = time()

    flag, error_msg = mesh.CheckAll()
    if not flag:
        mesh.PrintErrorMessage(sys.stderr, error_msg)
        sys.stderr.write("  Exiting...\n")
        exit(-1)

    # Check manifold and orientation.
    flag_manifold_vertices, flag_manifold_edges, iv, ihalf_edgeM =\
        mesh.CheckManifold()
    flag_manifold = flag_manifold_vertices and flag_manifold_edges
    flag_oriented, ihalf_edgeO = mesh.CheckOrientation()

    if (not flag_no_warn):
        if not(flag_manifold_vertices):
            mesh.PrintNonManifoldVertex(sys.stderr, "Warning:", iv)
        if not(flag_manifold_edges):
            mesh.PrintNonManifoldEdge(sys.stderr, "Warning:", ihalf_edgeM)
        if flag_manifold_edges and not(flag_oriented):
            mesh.PrintNotOriented(sys.stderr, "Warning:", ihalf_edgeO)

    if (not flag_silent):
        if (flag_manifold and flag_oriented):
            print("Mesh data structure passed check.")

    num_boundary_edges = mesh.CountNumBoundaryEdges()

    # Split edges
    edge_list = mesh.GetEdgeList()

    for half_edge in edge_list:
        ihalf_edge = half_edge.Index()
        if (flag_verbose):
            s = half_edge.IndexAndEndpointsStr(",")
            print(f"Splitting half edge {s}.")
        mesh.SplitEdge(ihalf_edge)

    if (not flag_silent):
        print("Split completed.")

    # Check mesh after splitting vertices.
    flag, error_msg = mesh.CheckAll()
    if not flag:
        sys.stderr.write\
            ("Error detected in mesh data structure after splitting edges.\n")
        if not(error_msg is None):
            sys.stderr.write(error_msg + "\n")
        exit(-1)

    # Check mesh after splitting edges.
    check_mesh_after_split\
        (mesh, flag_manifold, flag_oriented, num_boundary_edges)

    if (output_filename is None):
        output_filename = "out.off"
    if (output_filename == input_filename):
        output_filename = "out2.off"

    time3 = time()

    if (not flag_silent):
        print("Writing file: " + output_filename)

    half_edge_mesh_IO.open_and_write_file(output_filename, mesh)

    end_time = time()

    if (flag_time):
        print_time("Time to read file:  ", (time2-begin_time))
        print_time("Time to check mesh: ", (time3-time2))
        print_time("Time to write file: ", (end_time-time3))
        print_time("Total time:         ", (end_time-begin_time))


# ****** SUBROUTINES ******

## Check mesh after split of edges.
#  @param flag_manifold True if original mesh is manifold.
#  @param flag_oriented True if original mesh is oriented.
#  @param num_boundary_edges Number of boundary edges
#       in original mesh.
def check_mesh_after_split\
    (mesh, flag_manifold, flag_oriented, num_boundary_edges):

    flag, error_msg = mesh.CheckAll()
    if not flag:
        sys.stderr.write\
            ("Error detected in mesh data structure after splitting edges.\n")
        if not(error_msg is None):
            sys.stderr.write(error_msg + "\n")
        sys.stderr.write("  Exiting...")
        exit(-1)

    flag_new_manifold_vertices, flag_new_manifold_edges, iv, ihalf_edgeM =\
        mesh.CheckManifold()
    flag_new_manifold = flag_new_manifold_vertices and flag_new_manifold_edges
    flag_new_oriented, ihalf_edgeO = mesh.CheckOrientation()

    if flag_manifold:
        if not(flag_new_manifold):
            sys.stderr.write\
                ("Error. Input mesh is a manifold but mesh after split is not.\n")
        if not(flag_new_manifold_vertices):
            mesh.PrintNonManifoldVertex(sys.stderr, "Error:", iv)
        if not(flag_new_manifold_edges):
            mesh.PrintNonManifoldEdge(sys.stderr, "Error:", ihalf_edgeM)
    if flag_manifold and flag_oriented:
        if flag_new_manifold_edges and not(flag_new_oriented):
            sys.stderr.write\
                ("Error. Input mesh is oriented but mesh after split is not.\n")
            mesh.PrintNotOriented(sys.stderr, "Error:", ihalf_edgeO)

    # Check num boundary edges.
    numh = mesh.CountNumBoundaryEdges()
    if (numh != 2*num_boundary_edges):
        sys.stderr.write\
            ("Error detected in boundary after splitting cells.\n")
        sys.stderr.write\
            ("  Original mesh has " + str(num_boundary_edges) +\
                " boundary edges.\n")
        sys.stderr.write\
            ("  New mesh has " + str(numh) + " boundary edges.\n")
        sys.stderr.write\
            ("  New mesh should have " + str(2*num_boundary_edges) +\
                " boundary edges.\n")

# Print timeX in seconds.
def print_time(label, timeX):
    print(label + "{:.4f}".format(timeX) + " seconds")


def parse_command_line(argv):
    global input_filename, output_filename, flag_no_warn
    global flag_silent, flag_verbose
    global flag_time

    iarg = 1
    while (iarg < len(argv) and argv[iarg][0] == '-'):
        s = argv[iarg]
        if (s == "-s"):
            flag_silent = True
        elif (s == "-verbose"):
            flag_verbose = True
        elif (s == "-no_warn"):
            flag_no_warn = True
        elif (s == "-time"):
            flag_time = True
        elif (s == "-h"):
            help()
        else:
            sys.stderr.write("Usage error. Option " + s + " is undefined.\n")
            usage_error()

        iarg += 1

    if (iarg >= len(argv) or iarg+2 < len(argv)):
        usage_error()

    input_filename = argv[iarg]

    if (iarg+1 < len(argv)):
        output_filename = argv[iarg+1]

    return


def usage_msg(out):
    out.write("Usage: python3 test_hmesh_split_edges.py [-s] [-no_warn] [-time] [-h] <input filename> [<output_filename>]")

def usage_error():
    usage_msg(sys.stderr)
    sys.stderr.flush()
    exit(-1)

def help():
    usage_msg(sys.stdout)
    print("\n\ntest_hmesh_split_edges.py -- Test splitting edgs in HALF_EDGE_MESH.")
    print("  Read mesh from a .off file, split each edge,")
    print("  and write the mesh to a .off file.")
    print("\n\nOptions:")
    print("-s:        Silent. Output only warnings and error messages.")
    print("-no_warn:  Do not output non-manifold or inconsistent orientation warnings.")
    print("-time:     Report run time.")
    print("-h:        Output this help message and exit.")
    exit(0)

if __name__ == '__main__':
    main(sys.argv)
