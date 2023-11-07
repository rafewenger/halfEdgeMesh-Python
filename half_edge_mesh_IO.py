## \file half_edge_mesh_IO.py
# Functions for 2D half edge mesh read/write.
# - Requires python3.
# - Version 0.0.1

#  Copyright (C) 2021 Rephael Wenger
#
#  This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# (LGPL) as published by the Free Software Foundation; either
#  version 2.1 of the License, or any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

import half_edge_mesh

## Get first line that is not blank or a comment.
# - First non-blank character in a comment line is '#'.
def get_next_non_comment_line(infile):

    line = infile.readline()
    while line:
        line = line.lstrip()
        line = line.rstrip()
        if (len(line) > 0):
            if (line[0] != '#'): return(line)
        line = infile.readline()
    return(line)


## Read off file.
#  - mesh is derived from HALF_EDGE_MESH_BASE.
#  @pre Dimension of vertices is 3.
def read_off_file(infile, mesh):

    line = infile.readline()
    line = line.rstrip()

    if (line != "OFF"):
        raise Exception("Read error. File does not being with OFF.")

    line = get_next_non_comment_line(infile)
    if (not line):
        raise Exception("Read error. File does not contain line with number of vertices and polygons.")

    listX = line.split()
    if (len(listX) == 0):
        raise Exception("Incorrect file format. File missing number of vertices and polygons.")
    elif (len(listX) == 1):
        raise Exception("Incorrect file format. File missing number of polygons.")

    numv = int(listX[0])
    numpoly = int(listX[1])

    mesh.AddVertices(numv)
    for iv in range(0,numv):
        line = get_next_non_comment_line(infile)
        if (not line):
            raise Exception("Read error. File is missing some vertex coordinates.")

        listX = line.split()
        if (len(listX) < 3):
            raise Exception("Read error.  Error reading vertex coordinates.")

        coord = [ float(listX[0]), float(listX[1]), float(listX[2])]
        mesh.SetCoord(iv, coord)

    for ipoly in range(0,numpoly):
        line = get_next_non_comment_line(infile)
        if (not line):
            raise Exception("Read error. File is missing some polygon vertices.")

        listX = line.split()
        if (len(listX) < 1):
            raise Exception("Read error. Error reading polygon vertices.")

        num_poly_vert = int(listX[0])
        if (len(listX) < num_poly_vert+1):
            raise Exception("Read error. Error reading polygon vertices.")

        cell_vlist = []
        for k in range(1,num_poly_vert+1):
            cell_vlist.append(int(listX[k]))

        mesh.AddCell(ipoly, cell_vlist)



## Open and read off file into mesh.
#  - mesh is derived from HALF_EDGE_MESH_BASE.
#  @pre Dimension of vertices is 3.
def open_and_read_off_file(input_filename, mesh):

    try:
        with open(input_filename, 'r') as infile:
            read_off_file(infile, mesh)

    except Exception as e:
        print("Error reading file " + input_filename + ".")
        print(e)
        exit(-1)


## Write off file.
#  - mesh is derived from HALF_EDGE_MESH_BASE.
#  @pre Dimension of vertices is 3.
def write_off_file(outfile, mesh):

    # Write OFF label.
    outfile.write("OFF\n")

    # Write number of vertices and polygons.
    outfile.write(str(mesh.MaxVertexIndex()+1) + " " + str(mesh.NumCells()) + " 0\n")
    outfile.write("\n")

    # Write vertex coordinates.
    # Points
    buffer = ""
    for iv in range(0,mesh.MaxVertexIndex()+1):
        v = mesh.Vertex(iv)
        if (v is None):
            buffer += "0 0 0\n"
        else:
            buffer += v.CoordStr() + "\n"

    buffer += "\n"
    outfile.write(buffer)

    # Write polygon vertices.
    num_poly = 0
    poly_indices = list(mesh.CellIndices())
    poly_indices.sort()
    buffer = ""
    for ipoly in poly_indices:
        poly = mesh.Cell(ipoly)
        if (poly is None):
            continue
        else:

            if (num_poly >= mesh.NumCells()):
                # Error. Number of cells in mesh.cell_dict does not equal NumCells()?!?
                raise Exception("Error in write_off_file. Incorrect mesh.NumCells().")

            s = str(poly.NumVertices()) + " ";
            half_edge = poly.HalfEdge()
            for k in range (0, poly.NumVertices()):
                s += " " + str(half_edge.FromVertexIndex())
                half_edge = half_edge.NextHalfEdgeInCell()
            buffer += s + "\n"
            num_poly += 1

    outfile.write(buffer)

#    # Write polygon vertices.
#    num_poly = 0
#    poly_indices = list(mesh.CellIndices())
#    poly_indices.sort()
#    for ipoly in poly_indices:
#        poly = mesh.Cell(ipoly)
#        if (poly is None):
#            continue
#        else:
#
#            if (num_poly >= mesh.NumCells()):
#                # Error. Number of cells in mesh.cell_dict does not equal NumCells()?!?
#                raise Exception("Error in write_off_file. Incorrect mesh.NumCells().")
#
##            half_edge = poly.HalfEdge()
#            for k in range (0, poly.NumVertices()):
#                s += " " + str(half_edge.FromVertexIndex())
#                half_edge = half_edge.NextHalfEdgeInCell()
#            outfile.write(s + "\n")
#            num_poly += 1

    return


## Open and write mesh into off file.
#  - mesh is derived from HALF_EDGE_MESH_BASE.
#  @pre Dimension of vertices is 3.
def open_and_write_file(output_filename, mesh):

    try:
        with open(output_filename, 'w') as outfile:
            write_off_file(outfile, mesh)

    except Exception as e:
        print("Error write file " + output_filename + ".")
        print(e)
        exit(-1)
