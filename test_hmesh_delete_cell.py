## \file test_hmesh_delete_cell.py
#    Test program for deleting cells from half_edge_mesh.py.
#
#  - Reads a .off file into HALF_EDGE_MESH_BASE.
#  - Checks mesh data structure consistency.
#  - Deletes cells in order, checking mesh data structure consistency
#    at each iteration.


from time import time, localtime, strftime
import sys
import half_edge_mesh
import half_edge_mesh_IO
from half_edge_mesh_elements\
    import HMESH_VERTEX_BASE, HMESH_HALF_EDGE_BASE, HMESH_CELL_BASE
from half_edge_mesh import HALF_EDGE_MESH_BASE

def main(argv):

    global input_filename, output_filename, flag_no_warn
    global flag_verbose, flag_silent
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
        mesh.PrintErrorMessge(sys.stderr, error_msg)
        sys.stderr("  Exiting...\n")
        exit(-1)

    if (not flag_silent):
        print("Mesh data structure passed initial check.")

    cell_list = list(mesh.CellIndices())

    for icell in cell_list:
        if (flag_verbose):
            print(f"Deleting cell {icell}.")

        mesh.DeleteCell(icell)
        flag, error_msg = mesh.CheckAll()
        if not flag:
            sys.stderr.write\
                ("Error detected in mesh data structure after deleting cell " + str(icell) + ".\n")
            if not(error_msg is None):
                sys.stderr.write(error_msg + "\n")
            exit(-1)

    if (mesh.NumCells() != 0):
        sty.stderr.write\
            ("Error. All cells deleted but mesh.NumCells() = " + \
                str(mesh.NumCells()) + ".\n")

    if (not flag_silent):
        print("Mesh data structure passed all checkes after all deletions.")

    end_time = time()

    if (flag_time):
        print_time("Time to read file:  ", (time2-begin_time))
        print_time("Time to check delete cells: ", (time3-time2))
        print_time("Time to write file: ", (end_time-time3))
        print_time("Total time:         ", (end_time-begin_time))


# ****** SUBROUTINES ******

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
    out.write("Usage: python3 test_hmesh_delete_cell [-s] [-no_warn] [-time] <input filename>")

def usage_error():
    usage_msg(sys.stderr)
    sys.stderr.flush()
    exit(-1)

def help():
    usage_msg(sys.stdout)
    print("\n\ntest_hmesh_delete_cell.py -- Test the HALF_EDGE_MESH delete cell operations")
    print("  by reading a .off file to the mesh, repeatedly deleting cells,")
    print("  and running check mesh routines after each deletion.")
    print("  (Does not check manifold or orientation properties.)")
    print("\n\nOptions:")
    print("-s:        Silent. Output only warnings and error messages.")
    print("-time:     Report run time.")
    print("-h:        Output this help message and exit.")
    exit(0)

if __name__ == '__main__':
    main(sys.argv)
