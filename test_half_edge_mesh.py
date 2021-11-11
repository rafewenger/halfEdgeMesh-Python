## \file test_half_edge_mesh.py
#    Test program for half_edge_mesh.py.
#
#  - Reads a .off file into HALF_EDGE_MESH_BASE.
#  - Checks mesh, manifold and orientations properties of the mesh.
#  - Writes the mesh to a .off file.


from time import time, localtime, strftime
import sys
import half_edge_mesh
import half_edge_mesh_IO


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

    mesh = half_edge_mesh.HALF_EDGE_MESH_BASE()

    half_edge_mesh_IO.open_and_read_off_file(input_filename, mesh)

    time2 = time()

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

    if (not flag_silent):
        if not(flag_non_manifold or flag_not_oriented):
            print("Mesh data structure passed check.")

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
    out.write("Usage: python3 test_half_edge_mesh [-s] [-no_warn] [-time] <input filename> [<output_filename>]")

def usage_error():
    usage_msg(sys.stderr)
    sys.stderr.flush()
    exit(-1)

def help():
    usage_msg(sys.stdout)
    print("\n\ntest_half_edge_mesh.py -- Test the HALF_EDGE_MESH and associated classes")
    print("  and I/O routines by reading a .off file to the mesh,")
    print("  running check mesh, check manifold and check orientation routines")
    print("  and then writing the mesh to a .off file.")
    print("\n\nOptions:")
    print("-s:        Silent. Output only warnings and error messages.")
    print("-no_warn:  Do not output non-manifold or inconsistent orientation warnings.")
    print("-time:     Report run time.")
    print("-h:        Output this help message and exit.")
    exit(0)

if __name__ == '__main__':
    main(sys.argv)
