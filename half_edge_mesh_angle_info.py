## \file half_edge_mesh_angle_info.py
#  Class for returning angle information.

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


## Class storing information about min/max angles.
class HMESH_MIN_MAX_ANGLE_INFO:

    ## Cosine of min angle.
    #  - Note: Min angle has maximum cosine.
    cos_min_angle = 0.0

    ## Cosine of max angle.
    #  - Note: Max angle has minimum cosine.
    cos_max_angle = 0.0

    ## Index of half edge whose from vertex forms min angle.
    imin = 0

    ## Index of half edge whose from vertex forms max angle.
    imax = 0

    ## True if computation encounters a zero length edge.
    flag_zero = False

    ## Initialize structure.
    #  - Set cos_min_angle to -1.
    #  - Set cos_max_angle to 1.
    def Initialize(self):
        self.cos_min_angle = -1
        self.cos_max_angle = 1
        self.imin = 0
        self.imax = 0
        self.flag_zero = False

    ## Set cos_min_angle and half edge whose from vertex is min angle.
    def SetMinAngle(self, cos_angle, ihalf_edge):
        self.cos_min_angle = cos_angle
        self.imin = ihalf_edge

    ## Set cos_max_angle and half edge whose from vertex is max angle.
    def SetMaxAngle(self, cos_angle, ihalf_edge):
        self.cos_max_angle = cos_angle
        self.imax = ihalf_edge

    ## Set both min and max angle to cos_angle.
    def SetMinMaxAngle(self, cos_angle, ihalf_edge):
        self.SetMinAngle(cos_angle, ihalf_edge)
        self.SetMaxAngle(cos_angle, ihalf_edge)

    ## Copy min angle, max angle, imin, and imax.
    def Copy(self, info):
        self.SetMinAngle(info.cos_min_angle, info.imin)
        self.SetMaxAngle(info.cos_max_angle, info.imax)


## Class storing information about min/max angles
#    and number of cells with angles less than or equal to
#    or greater than or equal to some bound.
class HMESH_ANGLE_INFO(HMESH_MIN_MAX_ANGLE_INFO):

    ## List of small angle bounds.
    small_angle_bounds = []

    ## List of large angle bounds.
    large_angle_bounds = []

    ## num_cells_with_angles_le[i] =
    #    Number of cells with angles less than or equal to
    #    small_angle_bounds[i].
    num_cells_with_angle_le_small = []

    ## num_cells_with_angles_ge[i] =
    #    Number of cells with angles greater than or equal to
    #    large_angle_bounds[i].
    num_cells_with_angle_ge_large = []

    ## Clear all lists.
    def ClearLists(self):
        self.small_angle_bounds = []
        self.large_angle_bounds = []
        self.num_cells_with_angle_le_small = []
        self.num_cells_with_angle_ge_large = []

    ## Initialize structure.
    #  - Set cos_min_angle to -1.
    #  - Set cos_max_angle to 1.
    def Initialize(self):
        super().Initialize()
        self.ClearLists()

    ## Set list small_angle_bounds to _small_angle_bounds.
    def SetSmallAngleBounds(self, _small_angle_bounds):
        self.small_angle_bounds.clear()
        self.small_angle_bounds.extend(_small_angle_bounds)
        self.num_cells_with_angle_le_small =\
            [0] * len(self.small_angle_bounds)

    ## Set list large_angle_bounds to _large_angle_bounds.
    def SetLargeAngleBounds(self, _large_angle_bounds):
        self.large_angle_bounds.clear()
        self.large_angle_bounds.extend(_large_angle_bounds)
        self.num_cells_with_angle_ge_large =\
            [0] * len(self.large_angle_bounds)


    ## Increment num_cells_with_angle_le_small[i].
    def IncrementNumCellsWithAngleLESmall(self, i):
        self.num_cells_with_angle_le_small[i] =\
            self.num_cells_with_angle_le_small[i]+1

    ## Increment num_cells_with_angle_ge_large[i].
    def IncrementNumCellsWithAngleGELarge(self, i):
        self.num_cells_with_angle_ge_large[i] =\
            self.num_cells_with_angle_ge_large[i]+1
