## \file half_edge_mesh_coord.py
#  Functions for coordinate computations
#  (midpoints, edge lengths, angles, etc.)
#  - Uses numpy
#  - For simple calculations, works directly on lists.

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

import numpy as np
import math
from math import sqrt


## Return the midpoint of coord0[] and coord1[].
#  - Does not use numpy.
#  @param coord0[] Python list (not a numpy array).
#  @param coord1[] Python list (not a numpy array).
#  @pre len(coord1) == len(coord0).
def compute_midpoint(coord0, coord1, coord2):
    for ic in range(0,len(coord0)):
        coord2[ic] = (coord0[ic] + coord1[ic])/2.0

    return


## Return squared distance between two coordinates.
#  @param coord0[] Python list (not a numpy array).
#  @param coord1[] Python list (not a numpy array).
#  @pre len(coord1) == len(coord0).
def compute_squared_distance(coord0, coord1):
    c0 = np.array(coord0)
    c1 = np.array(coord1)
    c2 = c0-c1
    return np.inner(c2,c2)

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
#  - Returns also flag_zero, indicating if either vector is 0
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
#  - Returns also flag_zero, indicating if either coord0[] or
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
