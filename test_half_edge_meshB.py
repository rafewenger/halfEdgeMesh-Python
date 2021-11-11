## \file test_half_edge_meshB.py
#    Test using classes derived from VERTEX_BASE, HALF_EDGE_BASE,
#    CELL_BASE and HALF_EDGE_MESH_BASE.
#
#  - Reads in a mesh.
#  - Stores "new values" in the new mesh vertices, half edges
#    and cells.
#  - Prints some sample new mesh information.
#  - Creates a new mesh whose vertex and cell identifiers equal
#    the new values stored in the old mesh.
#  - Writes the new mesh to a .off file.
#  - The new .off file will have 5 vertices with coordinates (0C0)
#    since the new mesh has no vertices equal to 0, 1, 2, 3, or 4.

from time import time, localtime, strftime
import sys
import half_edge_mesh
import half_edge_mesh_IO
from half_edge_mesh import HALF_EDGE_MESH_BASE


def main(argv):

    global input_filename, output_filename, flag_silent, flag_no_warn
    global flag_time

    begin_time = time()

    # Initialize
    input_filename = None
    output_filename = None
    flag_silent = False
    flag_no_warn = False
    flag_time = False

    parse_command_line(sys.argv)

    mesh = HALF_EDGE_MESH_BASE(VERTEX_B,HALF_EDGE_B,CELL_B)

    half_edge_mesh_IO.open_and_read_off_file(input_filename, mesh)

    time2 = time()

    check_mesh(mesh, flag_no_warn)

    mesh = test_new_classes(mesh)
    check_mesh(mesh, flag_no_warn)
    if (not flag_silent):
        write_mesh_info(mesh)

    new_mesh = create_new_mesh(mesh)

    flag_passed_check = check_mesh(new_mesh, flag_no_warn)

    if (not flag_silent):
        if flag_passed_check:
            print("Mesh data structure passed check.")

    if (output_filename is None):
        output_filename = "out.off"
    if (output_filename == input_filename):
        output_filename = "out2.off"

    time3 = time()

    if (not flag_silent):
        print("Writing file: " + output_filename)

    half_edge_mesh_IO.open_and_write_file(output_filename, new_mesh)

    end_time = time()

    if (flag_time):
        print_time("Time to read file:  ", (time2-begin_time))
        print_time("Time to check mesh: ", (time3-time2))
        print_time("Time to write file: ", (end_time-time3))
        print_time("Total time:         ", (end_time-begin_time))


# ****** New classes ******

class VERTEX_B(half_edge_mesh.VERTEX_BASE):
    pass

    def __init__(self):
        # Initialize VERTEX_BASE.
        half_edge_mesh.VERTEX_BASE.__init__(self)

        # New field.
        self.new_val=0;

class HALF_EDGE_B(half_edge_mesh.HALF_EDGE_BASE):
    pass

    def __init__(self):

        # Initialize HALF_EDGE_BASE.
        half_edge_mesh.HALF_EDGE_BASE.__init__(self)

        # New field.
        self.new_val=0;

class CELL_B(half_edge_mesh.CELL_BASE):
    pass

    def __init__(self):

        # Initialize CELL_BASE.
        half_edge_mesh.CELL_BASE.__init__(self)

        # New field.
        self.new_val=0;


# ****** Test new classes ******

def test_new_classes(mesh):

    for iv in mesh.VertexIndices():
        v = mesh.Vertex(iv)
        v.new_val = v.Index()+5;

    for ihalf_edge in mesh.HalfEdgeIndices():
        half_edge = mesh.HalfEdge(ihalf_edge)
        half_edge.new_val = half_edge.Index()+5;

    for icell in mesh.CellIndices():
        cell = mesh.Cell(icell);
        cell.new_val = cell.Index()+5;

    return mesh;

def write_mesh_info(mesh):
    iv = int(mesh.NumVertices()/2)
    v = mesh.Vertex(iv)
    if (v is not(None)):
        print("Vertex(" + str(iv) + ").new_val = " + str(v.new_val))

    ihalf_edge = int(mesh.NumHalfEdges()/2)
    half_edge = mesh.HalfEdge(ihalf_edge)
    if (half_edge is not(None)):
        print("HalfEdge(" + str(ihalf_edge) + ").new_val = " + str(half_edge.new_val))

    icell = int(mesh.NumCells()/2)
    cell = mesh.Cell(icell)
    if (cell is not(None)):
        print("Cell(" + str(icell) + ").new_val = " + str(cell.new_val))


## Create a new mesh using indices at new_val.
def create_new_mesh(mesh):

    # Create new mesh.
    new_mesh = HALF_EDGE_MESH_BASE(VERTEX_B,HALF_EDGE_B,CELL_B)

    # Create vertices in the new mesh.
    for iv in mesh.VertexIndices():
        v = mesh.Vertex(iv)
        ivnew = v.new_val
        new_mesh.AddVertex(ivnew)
        new_mesh.SetCoord(ivnew, v.coord)

    # Create cells in the new mesh.
    for icell in mesh.CellIndices():
        cell = mesh.Cell(icell)
        half_edge = cell.HalfEdge()
        cell_vlist = []
        for j in range(0,cell.NumVertices()):
            iv = half_edge.FromVertexIndex()
            v = mesh.Vertex(iv)
            iv_new = v.new_val;
            cell_vlist.append(iv_new)
            half_edge = half_edge.NextHalfEdgeInCell()
        new_mesh.AddCell(cell.new_val, cell_vlist)

    return new_mesh


# ****** SUBROUTINES ******

## Check mesh and mesh manifold and orientation properties.
#  - Return true if passed check.
def check_mesh(mesh, flag_no_warn):

    flag, error_msg = mesh.CheckAll()
    if not flag:
        sys.stderr.write("Error detected in mesh data structure.\n")
        if not(error_msg is None):
            sys.stderr.write(error_msg + "\n")
        exit(-1)

    flag_non_manifold = False
    flag_not_oriented = False
    if (not flag_no_warn):
        flag_non_manifold, flag_not_oriented = \
            warn_non_manifold_or_not_oriented(mesh)

    return (flag and not(flag_non_manifold) and not(flag_not_oriented))


## Print a warning message if the mesh is not a manifold or not oriented.
#  - Returns flag_non_manifold, flag_not_oriented.
def warn_non_manifold_or_not_oriented(mesh):

    flag_non_manifold_vertex, flag_non_manifold_edge, iv, ihalf_edge =\
        mesh.CheckManifold()

    if flag_non_manifold_vertex:
        sys.stderr.write("Warning: Non-manifold vertex " + str(iv) + ".\n")

    if flag_non_manifold_edge:
        sys.stderr.write("Warning: Non-manifold edge " + str(ihalf_edge) + ".\n")

    flag_non_manifold = flag_non_manifold_vertex or flag_non_manifold_edge

    # Initialize
    flag_not_oriented = False;

    if (flag_non_manifold_edge):
        flag_not_oriented = True;
    else:
        flag, ihalf_edge = mesh.CheckOrientation()
        flag_not_oriented = not(flag);
        if (flag_not_oriented):
            flag_oriented = False
            sys.stderr.write("Warning: Inconsistent orientation of cells incident on edge (" +\
                mesh.HalfEdge(ihalf_edge).EndpointsStr(",") + ").\n")

    return flag_non_manifold, flag_not_oriented


# Print timeX in seconds.
def print_time(label, timeX):
    print(label + "{:.4f}".format(timeX) + " seconds")


def parse_command_line(argv):
    global input_filename, output_filename, flag_silent, flag_no_warn
    global flag_time

    iarg = 1
    while (iarg < len(argv) and argv[iarg][0] == '-'):
        s = argv[iarg]
        if (s == "-s"):
            flag_silent = True
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
    out.write("Usage: python3 test_half_edge_meshB [-s] [-no_warn] [-time] <input filename> [<output_filename>]")

def usage_error():
    usage_msg(sys.stderr)
    sys.stderr.flush()
    exit(-1)

def help():
    usage_msg(sys.stdout)
    print("\n\ntest_half_edge_meshB.py -- Test derived classes in HALF_EDGE_MESH.")
    exit(0)

if __name__ == '__main__':
    main(sys.argv)
