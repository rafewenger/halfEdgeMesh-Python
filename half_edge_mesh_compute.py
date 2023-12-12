## \file half_edge_mesh_compute.py
#  Functions for computations
#  (midpoints, edge lengths, angles, etc.)
#  - Uses numpy
#  - For simple calculations, works directly on lists.

#  Copyright (C) 2021-2023 Rephael Wenger
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

import numpy as np
import math
from math import sqrt
from math import cos
from math import radians

from half_edge_mesh_angle_info import HMESH_ANGLE_INFO


# *** MIDPOINT COMPUTATION ***

## Return the midpoint of coord0[] and coord1[].
#  - Does not use numpy.
#  @param coord0[] Python list (not a numpy array).
#  @param coord1[] Python list (not a numpy array).
#  @pre len(coord1) == len(coord0).
def compute_midpoint(coord0, coord1, coord2):
    for ic in range(0,len(coord0)):
        coord2[ic] = (coord0[ic] + coord1[ic])/2.0

    return


# *** DISTANCE COMPUTATIONS ***

## Return squared distance between two coordinates.
#  @param coord0[] Python list (not a numpy array).
#  @param coord1[] Python list (not a numpy array).
#  @pre len(coord1) == len(coord0).
def compute_squared_distance(coord0, coord1):
    c0 = np.array(coord0)
    c1 = np.array(coord1)
    c2 = c0-c1
    return np.inner(c2,c2)


## Return edge length squared.
#  @param half_edge Half edge.
def compute_edge_length_squared(half_edge):
    L = compute_squared_distance\
        (half_edge.FromVertex().coord, half_edge.ToVertex().coord)
    return L


## Return square of min and max edge lengths and indices of half edges with min and max lengths.
def compute_min_max_edge_length_squared(mesh):
    elist = mesh.GetEdgeList()
    if (len(elist) == 0):
        return 0.0, 0.0, 0, 0

    # Initialize
    min_Lsquared = compute_edge_length_squared(elist[0])
    max_Lsquared = min_Lsquared
    half_edge_min = elist[0]
    half_edge_max = half_edge_min

    for i in range(len(elist)):
        if (i > 0):
            half_edge = elist[i]
            L = compute_edge_length_squared(half_edge)
            if (L < min_Lsquared):
                min_Lsquared = L
                half_edge_min = half_edge
            elif (L > max_Lsquared):
                max_Lsquared = L
                half_edge_max = half_edge

    ihalf_edge_min = half_edge_min.Index()
    ihalf_edge_max = half_edge_max.Index()

    return min_Lsquared, max_Lsquared, ihalf_edge_min, ihalf_edge_max


## Return square of min and max edge lengths in given CELL and indices of half edges with min and max lengths.
#  @param cell Cell. Compute square of cell edge lengths.
def compute_cell_min_max_edge_length_squared(cell):

    if (cell.NumVertices() < 1):
        return 0.0, 0.0, 0, 0

    half_edge = cell.HalfEdge()
    min_Lsquared = compute_edge_length_squared(half_edge)
    max_Lsquared = min_Lsquared
    half_edge_min = half_edge
    half_edge_max = half_edge
    for i in range(cell.NumVertices()):
        if (i > 0):
            L = compute_edge_length_squared(half_edge)
            if (L < min_Lsquared):
                min_Lsquared = L
                half_edge_min = half_edge
            elif (L > max_Lsquared):
                max_Lsquared = L
                half_edge_max = half_edge
        half_edge = half_edge.NextHalfEdgeInCell()

    imin_half_edge = half_edge_min.Index()
    imax_half_edge = half_edge_max.Index()

    return min_Lsquared, max_Lsquared, imin_half_edge, imax_half_edge


## Return ratio of the min squared edge length to max squared edge length in cell.
#  - Returns ratio, min edge length squared, max edge length squared,
#    indices of half edges with min max edge lengths.
#  - Returns ratio 1.0 if all edges have length 0.
#  @param cell Cell.
def compute_cell_edge_length_ratio_squared(cell):

    min_Lsquared, max_Lsquared, imin_half_edge, imax_half_edge =\
        compute_cell_min_max_edge_length_squared(cell)

    ratio = 1.0
    if (max_Lsquared > 0.0):
        ratio = min_Lsquared/max_Lsquared

    return ratio, min_Lsquared, max_Lsquared, imin_half_edge, imax_half_edge


## Return min squared ratio of the min to max edge in any cell.
#  - Ignores cells with all edge lengths 0.
#  - Return also cell index and length and indices
#      of shortest and longest half edges in the cell.
#  - Returns 1.0 if there are no cells or all edges are length 0.
def compute_min_cell_edge_length_ratio_squared(mesh):

    # Initialize.
    min_edge_length_ratio_squared = 1.0
    icell_min_ratio = 0
    min_edge_length_squared = 0.0
    max_edge_length_squared = 0.0
    ihalf_edge_min = 0
    ihalf_edge_max = 0

    flag_is_set = False
    for icell in mesh.CellIndices():
        cell = mesh.Cell(icell);
        if (cell is None):
            # Shouldn't happen but just in case.
            continue

        ratio, min_Lsquared, max_Lsquared, ihalf_min, ihalf_max =\
            compute_cell_edge_length_ratio_squared(cell)

        if (max_Lsquared == 0.0 or cell.NumVertices() == 0):
            continue

        if not(flag_is_set) or (ratio < min_edge_length_ratio_squared):
            min_edge_length_ratio_squared = ratio
            icell_min_ratio = icell
            min_edge_length_squared = min_Lsquared
            max_edge_length_squared = max_Lsquared
            ihalf_edge_min = ihalf_min
            ihalf_edge_max = ihalf_max
            flag_is_set = True

    return min_edge_length_ratio_squared, icell_min_ratio,\
            min_edge_length_squared, max_edge_length_squared,\
            ihalf_edge_min, ihalf_edge_max


# *** ANGLE COMPUTATIONS ***


## Return normalized vector and magnitude of the original vector.
#  - If the vector has magnitude 0, returns vector (1,0,0,...) and 0.
#  @param vect[] numpy array.
def normalize_vector(vectA):
    magnitudeA = sqrt(np.inner(vectA,vectA))

    if (abs(magnitudeA) == 0.0):
        vectB = np.zeros(len(vectA))
        vectB[0] = 1.0
        return vectB, 0.0

    vectC = np.abs(vectA)
    imax = np.argmax(vectC)
    maxc = vectC[imax]
    vectC = vectA/maxc
    if (vectC[imax] < 0):
        vectC[imax] = -1
    else:
        vectC[imax] = 1

    # Since vectC[imax] is 1 or -1, magnitudeC >= 1.
    magnitudeC = sqrt(np.inner(vectC,vectC))

    # Divide by magnitudeC.
    vectC = vectC/magnitudeC

    return vectC, magnitudeC


## Returns cosine of angle between two vectors.
#  - Returns also flag_zero, true if either vector is 0
#    (or very, very near 0).
#  - If either vector is zero, returns 0 as cosine.
#  @param vect0[] numpy array.  Function modifies vect0[].
#  @param vect1[] numpy array.  Function modifies vect1[].
def compute_cos_angle(vect0, vect1):

    vect0, magnitude0 = normalize_vector(vect0)
    vect1, magnitude1 = normalize_vector(vect1)

    if (magnitude0 == 0.0) or (magnitude1 == 0.0):
        return 0, True

    cos_angle = np.inner(vect0,vect1)

    # Clamp to [-1,1] to handle numerical error.
    if (cos_angle < -1):
        return -1, False
    if (cos_angle > 1):
        return 1, False

    return cos_angle, False


## Compute cos of triangle angle at coord1[]
#    in triangle (coord0[],coord1[],coord2[]).
#  - Returns also flag_zero, true if either coord0[] or
#    coord2[] are at coord1[] (or very, very close to coord1[].)
#  - If coord0[] == coord1[] (or is very, very close,) or
#    coord2[] == coord1[] (or is very, very close,)
#    returns 0 as cosine.
#  @param coord0[] Python list (not a numpy array).
#  @param coord1[] Python list (not a numpy array).
#  @pre len(coord1) == len(coord0).
#  @param coord2[] Python list (not a numpy array).
#  @pre len(coord2) == len(coord0).
def compute_cos_triangle_angle(coord0, coord1, coord2):
    c0 = np.array(coord0)
    c1 = np.array(coord1)
    c2 = np.array(coord2)

    return compute_cos_angle(c0-c1, c2-c1)


## Compute cos of angle at half_edge1.FromVertex() in cell half_edge1.Cell().
#  - Returns cosine of angle and flag_zero.
#  - Returns also flag_zero, true if either coord0[] or
#    coord2[] are at coord1[] (or very, very close to coord1[].)
def compute_cos_vertex_angle(half_edge1):
    half_edge0 = half_edge1.PrevHalfEdgeInCell()
    half_edge2 = half_edge1.NextHalfEdgeInCell()
    v0 = half_edge0.FromVertex()
    v1 = half_edge1.FromVertex()
    v2 = half_edge2.FromVertex()

    return compute_cos_triangle_angle(v0.coord, v1.coord, v2.coord)


## Compute cos of min and max angles in cell.
#  - Returns cos_min, cos_max and half edges whose from vertices
#    have the min and max angles.
#  - Returns also flag_zero, true if some edge has zero length.
#  - Note: min angle has maximum cosine and max angle has minimum cosine.
#  - Returns 0.0, 0.0 if cell has no vertices.
#  - Ignore vertices incident on 0 length edges.
def compute_cos_min_max_cell_angles(cell):
    if (cell.NumVertices() < 1):
        return 0.0, 0.0, 0, 0, False

    half_edge = cell.HalfEdge()
    flag_found = False

    # Initialize
    cos_minA, flag = compute_cos_vertex_angle(half_edge)
    flag_zero = flag
    cos_maxA = cos_minA
    half_edge_minA = half_edge
    half_edge_maxA = half_edge

    flag_found = False
    if not(flag):
        flag_found = True

    for i in range(cell.NumVertices()):
        cosA, flag = compute_cos_vertex_angle(half_edge)
        flag_zero = (flag_zero or flag)

        if not(flag):
            if not(flag_found):
                # First vertex incident on edges with non-zero lengths.
                cos_minA = cosA
                cos_maxA = cosA
                half_edge_minA = half_edge
                half_edge_maxA = half_edge
            elif (cosA > cos_minA):
                # Note: if cos angle is large, then angle is small.
                cos_minA = cosA
                half_edge_minA = half_edge
            elif (cosA < cos_maxA):
                # Note: if cos angle is small, then angle is large.
                cos_maxA = cosA
                half_edge_maxA = half_edge
        half_edge = half_edge.NextHalfEdgeInCell()

    iminA_half_edge = half_edge_minA.Index()
    imaxA_half_edge = half_edge_maxA.Index()

    return cos_minA, cos_maxA, iminA_half_edge, imaxA_half_edge, flag_zero


## Compute cos of min and max angles in mesh.
#  - Returns cos_min_max_angle_info.
#  - Ignore vertices incident on 0 length edges.
#  @param small_angle_bounds[] Array of small angle bounds.
#    - Computes number of cells with angles less than
#      or equal to small_angle_bounds[i] for each i.
#  @param large_angle_bounds[] Array of large angle bounds.
#    - Computes number of cells with angles greater than
#      or equal to large_angle_bounds[i] for each i.
def compute_angle_info\
    (mesh, small_angle_bounds, large_angle_bounds):
    angle_info = HMESH_ANGLE_INFO()

    angle_info.Initialize()
    angle_info.SetSmallAngleBounds(small_angle_bounds)
    angle_info.SetLargeAngleBounds(large_angle_bounds)

    cos_small_angle_bounds = []
    for j in range(len(angle_info.small_angle_bounds)):
        A = angle_info.small_angle_bounds[j]
        cos_small_angle_bounds.append(cos(radians(A)))

    cos_large_angle_bounds = []
    for j in range(len(angle_info.large_angle_bounds)):
        A = angle_info.large_angle_bounds[j]
        cos_large_angle_bounds.append(cos(radians(A)))

    for icell in mesh.CellIndices():
        cell = mesh.Cell(icell);
        if (cell is None):
            # Shouldn't happen but just in case.
            continue

        cos_minA, cos_maxA, ihalf_edge_minA, ihalf_edge_maxA, flag =\
            compute_cos_min_max_cell_angles(cell)

        angle_info.flag_zero = (angle_info.flag_zero or flag)

        if (cos_minA > angle_info.cos_min_angle):
            # Note: if cos angle is large, then angle is small.
            angle_info.SetMinAngle(cos_minA, ihalf_edge_minA)

        if (cos_maxA < angle_info.cos_max_angle):
            # Note: if cos angle is small, then angle is large.
            angle_info.SetMaxAngle(cos_maxA, ihalf_edge_maxA)

        for j in range (len(angle_info.small_angle_bounds)):
            if (cos_minA >= cos_small_angle_bounds[j]):
                # Note: If cos_minA >= cos_small_angle_bounds[j],
                #   then minA <= small_angle_bounds[j].
                angle_info.IncrementNumCellsWithAngleLESmall(j)

        for j in range (len(angle_info.large_angle_bounds)):
            if (cos_maxA  <= cos_large_angle_bounds[j]):
                # Note: If cos_maxA <= cos_large_angle_bounds[j],
                #   then maxA >= Large_angle_bounds[j].
                angle_info.IncrementNumCellsWithAngleGELarge(j)

    return angle_info


# *** PRINT FUNCTIONS ***

## Print coordinate array.
#  - Does not print brackets around array.
#  @param out Output stream.
#  @param coord Array of floating point coordinates.
#  @param separator Separator between coordinates.
#  @param precision Number of digits after decimal point.
#    - Non-negative integer.
def print_coord(out, coord, separator, precision):
    output_format = "{:." + str(precision) + "f}"

    for i in range(len(coord)):
        if (i > 0):
            out.write(separator)
        s = output_format.format(coord[i])
        if ("." in s):
            s = s.rstrip('0')
            if (s[-1] == '.'):
                s = s + "0"
        out.write(s)


## Print coordinate array surrounded by parentheses.
def print_coordP(out, coord, separator, precision):
    out.write("(")
    print_coord(out, coord, separator, precision)
    out.write(")")
