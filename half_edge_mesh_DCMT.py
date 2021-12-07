## \file half_edge_mesh_DCMT.py
# Extension of half edge mesh supporting decimation (DCMT) operations.

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
import half_edge_mesh_coord
from half_edge_mesh_coord import compute_squared_distance
from half_edge_mesh_coord import compute_cos_triangle_angle
from half_edge_mesh_coord import compute_midpoint


## Vertex class supporting mesh decimation.
class VERTEX_DCMT_BASE(half_edge_mesh.VERTEX_BASE):

    ## Initialize
    def __init__(self):

        super(VERTEX_DCMT_BASE,self).__init__()

        ## Internal flag used for detecting vertex adjacencies.
        self._visited_flag = False

    ## Return visited_flag.
    def IsVisited(self):
        return self._visited_flag

    ## SetVisitedFlag to False.
    def ClearVisitedFlag(self):
        self._visited_flag = False

    ## Set visited_flag to flag in all neighbors of self.
    def SetVisitedFlagsInAdjacentVertices(self, flag):

        for k in range(0,self.NumHalfEdgesFrom()):
            half_edgeA = self.KthHalfEdgeFrom(k)
            half_edgeA.ToVertex()._visited_flag = flag

            half_edgeB = half_edgeA.PrevHalfEdgeInCell()

            # Set half_edgeB.FromVertex()._visited_fag in case
            #   of boundary edges or cells with arbitrary orientations.
            half_edgeB.FromVertex()._visited_flag = flag


    ## Set visited_flag to False in all neighbors of self.
    def ClearVisitedFlagsInAdjacentVertices(self):
        self.SetVisitedFlagsInAdjacentVertices(False)


    ## Compare first and last half edges in half_edges_from[].
    #  - Swap half edges if last half edge is a boundary edge
    #    and first half edge is internal or if both are internal,
    #    but half_edge_from[-1].PrevHalfEdgeInCell() is boundary,
    #    while half_edge_from[0].PrevHalfEdgeInCell()) is interior.
    def _ProcessFirstLastHalfEdgesFrom(self):
        if (self.NumHalfEdgesFrom() < 2):
            # No swap
            return

        if (self.KthHalfEdgeFrom(0).IsBoundary()):
            # No swap
            return

        if (self.KthHalfEdgeFrom(-1).IsBoundary()):
            self._SwapHalfEdgesInHalfEdgeFromList(0,-1)
            return

        if (self.KthHalfEdgeFrom(0).PrevHalfEdgeInCell().IsBoundary()):
            # No swap.
            return

        if (self.KthHalfEdgeFrom(-1).PrevHalfEdgeInCell().IsBoundary()):
            self._SwapHalfEdgesInHalfEdgeFromList(0,-1)
            return


    # *** Public ***

    ## Return True if vertex is incident on more than edges.
    def IsIncidentOnMoreThanTwoEdges(self):
        TWO = 2
        num_half_edges_from = self.NumHalfEdgesFrom()

        if (num_half_edges_from > TWO):
            return True

        if not(self.IsBoundary()):
            return False

        if (num_half_edges_from == TWO):
            # Boundary vertex in two cells must have at least
            #   three incident edges.
            return True
        else:
            # Boundary vertex is in just one cell, and has exactly
            #   two incident edges
            return False


## Half edge class supporting mesh decimation.
class HALF_EDGE_DCMT_BASE(half_edge_mesh.HALF_EDGE_BASE):

    ## Compute square of edge length.
    def ComputeLengthSquared(self):
        return compute_squared_distance\
                (self.FromVertex().coord, self.ToVertex().coord)

    ## Returns cosine of angle at FromVertex() in triangle formed by
    #  PrevHalfEdge().FromVertex(), FromVertex(), ToVertex().
    #  - Returns also flag_zero indicating if middle vertex has same
    #    coordinates as one of the other two (or middle vertex is
    #    very, very close to one of the other two.)
    def ComputeCosAngleAtFromVertex(self):
        prev_half_edge = self.PrevHalfEdgeInCell()
        v0 = prev_half_edge.FromVertex()
        v1 = self.FromVertex()
        v2 = self.ToVertex()

        return compute_cos_triangle_angle\
                (v0.coord, v1.coord, v2.coord)

## Cell class supporting mesh decimation.
class CELL_DCMT_BASE(half_edge_mesh.CELL_BASE):

    ## Set visited_flag to flag in all cell vertices.
    def SetVisitedFlagsInAllVertices(self, flag):
        half_edge = self.HalfEdge()
        for k in range(0,self.NumVertices()):
            half_edge.FromVertex()._visited_flag = flag
            half_edge = half_edge.NextHalfEdgeInCell()

    ## Set visited_flag to False in all cell vertices.
    def ClearVisitedFlagsInAllVertices(self):
        self.SetVisitedFlagsInAllVertices(False)


    ## Return min and max squared edge lengths in a cell.
    #  - Return also half edges with min and max edge lengths.
    def ComputeMinMaxEdgeLengthSquared(self):
        if (self.NumVertices() == 0):
            # Empty cell.
            return 0.0, 0.0, 0, 0

        half_edge = self.HalfEdge()
        min_edge_length_squared = half_edge.ComputeLengthSquared()
        max_edge_length_squared = min_edge_length_squared
        ihalf_edge_min = half_edge.Index();
        ihalf_edge_max = ihalf_edge_min

        for i in range(1,self.NumVertices()):
            half_edge = half_edge.NextHalfEdgeInCell()

            length_squared = half_edge.ComputeLengthSquared()
            if (length_squared < min_edge_length_squared):
                min_edge_length_squared = length_squared
                ihalf_edge_min = half_edge.Index()

                # Note: min_edge_length_squared <= max_edge_length_squared, so
                #   if length_squared < min_edge_length_squared, then
                #   (length_squared > max_edge_length_squared) is False.
            elif (length_squared > max_edge_length_squared):
                max_edge_length_squared = length_squared
                ihalf_edge_max = half_edge.Index()

        return min_edge_length_squared, max_edge_length_squared,\
                ihalf_edge_min, ihalf_edge_max


    ## Return cosine of min and max angles between consecutive cell edges.
    #  - Return also half edges whose from vertices are incident
    #      on the two cell edges forming the min and max angles.
    #  - Note: The smallest angle has the largest cosine and
    #      the largest angle has the smallest cosine.
    def ComputeCosMinMaxAngle(self):
        if (self.NumVertices() == 0):
            # Empty cell.
            return 0.0, 0.0, 0, 0

        # Initialize
        cos_min_angle = 0.0
        cos_max_angle = 0.0
        ihalf_edge_min = 0
        ihalf_edge_max = 0

        half_edge = self.HalfEdge()

        flag_found = False
        for i in range(0,self.NumVertices()):

            cos_angle, flag_zero = half_edge.ComputeCosAngleAtFromVertex()

            if (not(flag_zero)):
                if (not(flag_found)):
                    cos_min_angle = cos_angle
                    cos_max_angle = cos_angle
                    ihalf_edge_min = half_edge.Index()
                    ihalf_edge_max = ihalf_edge_min
                    flag_found = True
                elif (cos_angle > cos_min_angle):
                    # Remember: Small angles have large cos values.
                    cos_min_angle = cos_angle
                    ihalf_edge_min = half_edge.Index()

                    # Note: cos_min_angle >= cos_max_angle, so
                    #   if cos_angle > cos_min_angle, then
                    #   (cos_angle < cos_max_angle) is False.
                elif (cos_angle < cos_max_angle):
                    cos_max_angle = cos_angle
                    ihalf_edge_max = half_edge.Index()

            half_edge = half_edge.NextHalfEdgeInCell()

        return cos_min_angle, cos_max_angle,\
                ihalf_edge_min, ihalf_edge_max


# Half edge mesh decimation class.
#  - Initialized with vertex, half edge and cell classes.
#  - These classes should be derived from VERTEX_DCMT_BASE,
#      HALF_EDGE_DCMT_BASE and CELL_DCMT_BASE.
class HALF_EDGE_MESH_DCMT_BASE(half_edge_mesh.HALF_EDGE_MESH_BASE):

    # Public find edge function.

    ## Return half edge (v0,v1) or (v1,v0) if it exists
    #  - Return None if no edge found.
    def FindEdge(self, v0, v1):
        half_edge = v0.FindIncidentHalfEdge(v1.Index())

        if not(half_edge is None):
            return half_edge

        half_edge = v1.FindIncidentHalfEdge(v0.Index())

        return half_edge


    # Private member functions.

    def __init__(self, classV, classHE, classC):
        super(HALF_EDGE_MESH_DCMT_BASE,self).__init__(classV, classHE, classC)

    ## Remove half edge from the half_edge_from list of its from_vertex.
    #  - Does not ensure that first element is a boundary half edge.
    #  - Call _MoveBoundaryHalfEdgeToIncidentHalfEdgeFrom0() to ensure
    #    that first half edge in vertex list is a boundary edge.
    def _RemoveHalfEdgeFromVertexList(self, half_edge0):
        v0 = half_edge0.FromVertex()
        list_length = v0.NumHalfEdgesFrom()
        ilast = list_length-1

        for k in range(0,list_length):
            half_edge = v0.KthHalfEdgeFrom(k)
            if (half_edge0 is half_edge):
                if (k != ilast):
                    # Replace half_edge0 with last entry.
                    v0.half_edge_from[k] = v0.half_edge_from[ilast]

                v0.half_edge_from.pop()
                return


    ## Move half edges in vA.half_edge_from[] to vB.half_edge_from[].
    #  - Clear vA.half_edge_from[]
    def _MoveVertexHalfEdgeFromList(self, vA, vB):
        for k in range(0,vA.NumHalfEdgesFrom()):
            half_edge = vA.KthHalfEdgeFrom(k)
            half_edge.from_vertex = vB

        # Add vA.half_edge_from[] to vB.half_edge_from[].
        vB.half_edge_from.extend(vA.half_edge_from)

        vA.half_edge_from.clear()


    ## Swap next_half_edge_around_edge.
    def _SwapNextHalfEdgeAroundEdge(self, half_edgeA, half_edgeB):
        tempA = half_edgeA.NextHalfEdgeAroundEdge()
        tempB = half_edgeB.NextHalfEdgeAroundEdge()

        half_edgeA.next_half_edge_around_edge = tempB
        half_edgeB.next_half_edge_around_edge = tempA

    ## Return previous half edge around edge.
    def _FindPrevHalfEdeAroundEdge(self, half_edge0):
        max_numh = half_edge0.FromVertex().NumHalfEdgesFrom() +\
                    half_edge0.ToVertex().NumHalfEdgesFrom()

        half_edge = half_edge0.NextHalfEdgeAroundEdge()
        for k in range(0,max_numh):
            if (half_edge0 is half_edge.NextHalfEdgeAroundEdge()):
                return half_edge

            half_edge = half_edge.NextHalfEdgeAroundEdge()

        # Should never reach here. Data structure inconsistency.
        raise Exception\
            ("Programming error. Unable to find previous half edge around edge.")


    ## Find some half edge (v0,v1) or (v1,v0) and link with half_edgeA
    #    in half edge around edge cycle.
    #  - If half edge not found, then do nothing.
    def _FindAndLinkHalfEdgeAroundEdge(self, v0, v1, half_edgeA):
        half_edgeB = self.FindEdge(v0, v1)
        if (half_edgeB is None):
            return

        self._SwapNextHalfEdgeAroundEdge(half_edgeA, half_edgeB)

        # Reorder v1.half_edge_from in case half_edgeA
        #   is no longer a boundary edge.
        v1.MoveBoundaryHalfEdgeToIncidentHalfEdge0()


    ## Link half edges merged by merging v0 and v1.
    #  - Search for all possible merged edges, not just triangles
    #    containing v0 and v1, to handle non-manifolds.
    #  - Running time: O(v0.NumHalfEdgesFrom() + v1.NumHalfEdgesFrom())
    def _LinkHalfEdgesAroundMergedEdges(self,v0,v1):
        v1.ClearVisitedFlagsInAdjacentVertices()
        v0.SetVisitedFlagsInAdjacentVertices(True)

        for k in range(0, v1.NumHalfEdgesFrom()):
            half_edgeA = v1.KthHalfEdgeFrom(k)
            vtoA = half_edgeA.ToVertex()

            if (vtoA.IsVisited()):
                # vtoA is a neighbor of v0 and v1.

                self._FindAndLinkHalfEdgeAroundEdge(v0, vtoA, half_edgeA)

                # Set vtoA.visited_flag to False so that vtoA
                #   will not be processed twice.
                vtoA.ClearVisitedFlag()

            half_edgeB = half_edgeA.PrevHalfEdgeInCell()
            vfromB = half_edgeB.FromVertex()

            # Check vfromB to handle boundary edges and/or cells
            #   with arbitrary orientations.
            if (vfromB.IsVisited()):
                # vfromB is a neighbor of v0 and v1.

                self._FindAndLinkHalfEdgeAroundEdge(v0, vfromB, half_edgeB)

                # Set vfromB.visited_flag to False so that vfromB
                #   will not be processed twice.
                vfromB.ClearVisitedFlag()


    ## Reliink half edges in cell.
    #  - Overwrites previous links.
    def _RelinkHalfEdgesInCell(self, hprev, hnext):
        hprev.next_half_edge_in_cell = hnext
        hnext.prev_half_edge_in_cell = hprev


    ## Delete vertex.
    def _DeleteVertex(self, v):
        if (v is None):
            # Can't delete None.
            return

        iv = v.Index()
        self._vertex_dict.pop(iv,0)


    ## Delete half edge.
    def _DeleteHalfEdge(self, half_edge0):

        if not(half_edge0.IsBoundary()):
            next_half_edge_around_edge = half_edge0.NextHalfEdgeAroundEdge()

            prev_half_edge_around_edge =\
                self._FindPrevHalfEdeAroundEdge(half_edge0)

            prev_half_edge_around_edge.next_half_edge_around_edge =\
                next_half_edge_around_edge

        self._RemoveHalfEdgeFromVertexList(half_edge0)

        ihalf_edge0 = half_edge0.Index()
        self._half_edge_dict.pop(ihalf_edge0,0)


    ## Delete half edges around edge.
    #  - max_numh is an upper bound on the number of half edges
    #      around the edge containing half_edge0.
    def _DeleteHalfEdgesAroundEdge(self, half_edge0, max_numh):
        half_edge = half_edge0

        for k in range(0, max_numh):
            next_half_edge_around_edge = half_edge.NextHalfEdgeAroundEdge()

            if (next_half_edge_around_edge is half_edge):
                # Delete half edge.
                ihalf_edge = half_edge.Index()
                self._RemoveHalfEdgeFromVertexList(half_edge)
                self._half_edge_dict.pop(ihalf_edge)
                return
            else:
                # Delete next_half_edge_around_edge.
                half_edge.next_half_edge_around_edge =\
                    next_half_edge_around_edge.NextHalfEdgeAroundEdge()
                inext_half_edge = next_half_edge_around_edge.Index()
                self._RemoveHalfEdgeFromVertexList(next_half_edge_around_edge)
                self._half_edge_dict.pop(inext_half_edge,0)


    ## Delete cell
    def _DeleteCell(self, cell):
        icell = cell.Index()
        self._cell_dict.pop(icell,0)


    # *** Internal split functions.

    def _SplitInternalEdge(self, half_edgeA):
        if (half_edgeA is None):
            raise Exception\
                ("Programming error. Argument to _SplitInternalEdge is None.")

        half_edgeB = half_edgeA.NextHalfEdgeAroundEdge()
        if not(half_edgeB.NextHalfEdgeAroundEdge() is half_edgeA):
            raise Exception\
                ("Programming error. Half edge passed to _SplitInternalEdge is in an edge shared by three or more cells.")

        if (half_edgeB is half_edgeA):
            raise Exception\
                ("Programming error. Half edge passed to _SplitInternalEdge is a boundary edge. Call _SplitBoundaryEdge().")

        vA = half_edgeA.FromVertex()
        vB = half_edgeB.FromVertex()
        cellA = half_edgeA.Cell()
        cellB = half_edgeB.Cell()
        numvA = cellA.NumVertices()
        numvB = cellB.NumVertices()
        nextA = half_edgeA.NextHalfEdgeInCell()
        nextB = half_edgeB.NextHalfEdgeInCell()

        # Create a new vertex.
        ivnew = self.MaxVertexIndex()+1
        newv = self.AddVertex(ivnew)

        if (newv is None):
            raise Exception("Error creating new vertex. Out of memory?")

        # Set newv to midpoint of (vA,vB).
        compute_midpoint(vA.coord, vB.coord, newv.coord)

        # Create two new half edges.
        inew_half_edgeA = self.MaxHalfEdgeIndex()+1

        # _AddHalfEdge increments cellA.num_vertices.
        new_half_edgeA = self._AddHalfEdge(inew_half_edgeA, cellA, newv)

        inew_half_edgeB = self.MaxHalfEdgeIndex()+1

        # _AddHalfEdge increments cellB.num_vertices.
        new_half_edgeB = self._AddHalfEdge(inew_half_edgeB, cellB, newv)

        newv.half_edge_from.append(new_half_edgeA)
        newv.half_edge_from.append(new_half_edgeB)

        # Relink half edges in cell.
        self._RelinkHalfEdgesInCell(half_edgeA, new_half_edgeA)
        self._RelinkHalfEdgesInCell(half_edgeB, new_half_edgeB)
        self._RelinkHalfEdgesInCell(new_half_edgeA, nextA)
        self._RelinkHalfEdgesInCell(new_half_edgeB, nextB)

        # Unlink half_edgeA.next_half_edge_around_edge and
        #   half_edgeB.next_half_edge_around_edge.
        half_edgeA.next_half_edge_around_edge = half_edgeA
        half_edgeB.next_half_edge_around_edge = half_edgeB

        # Link half edges around edge.
        self._LinkHalfEdgesAroundEdge(half_edgeA, new_half_edgeB)
        self._LinkHalfEdgesAroundEdge(half_edgeB, new_half_edgeA)

        # half_edgeA and half_edgeB are not boundary edges,
        #   but the previous edges in the cell might be boundary edges.
        vA.MoveBoundaryHalfEdgeToIncidentHalfEdge0()
        vB.MoveBoundaryHalfEdgeToIncidentHalfEdge0()

        return newv


    ## Split a boundary edge.
    #  - Returns new vertex.
    #  @pre half_edgeA is a boundary edge.
    def _SplitBoundaryEdge(self, half_edgeA):
        if (half_edgeA is None):
            raise Exception\
                ("Programming error. Argument to _SplitBoundaryEdge is None.")

        if not(half_edgeA.IsBoundary()):
            raise Exception\
                ("Programming error. Half edge passed to _SplitBoundaryEdges is not a boundary edge. Call _SplitInternalEdge().")

        vA = half_edgeA.FromVertex()
        vB = half_edgeA.ToVertex()
        cellA = half_edgeA.Cell()
        numvA = cellA.NumVertices()
        nextA = half_edgeA.NextHalfEdgeInCell()

        # Create a new vertex.
        ivnew = self.MaxVertexIndex()+1
        newv = self.AddVertex(ivnew)

        if (newv is None):
            raise Exception("Error creating new vertex. Out of memory?")

        # Set newv to midpoint of (vA,vB).
        compute_midpoint(vA.coord, vB.coord, newv.coord)

        # Create new half edge.
        inew_half_edgeA = self.MaxHalfEdgeIndex()+1

        # _AddHalfEdge() increments cellA.num_vertices.
        new_half_edgeA = self._AddHalfEdge(inew_half_edgeA,cellA,newv)

        newv.half_edge_from.append(new_half_edgeA)

        # Relink half edges in cell.
        self._RelinkHalfEdgesInCell(half_edgeA, new_half_edgeA)
        self._RelinkHalfEdgesInCell(new_half_edgeA, nextA)

        # No need to move edges in half_edge_from[] lists.

        return newv


    # Public member functions.

    # *** Collapse/join/split functions ***

    ## Collapse edge, mapping two vertices to a single vertex.
    def CollapseEdge(self, ihalf_edge0):
        half_edge0 = self.HalfEdge(ihalf_edge0);
        if (half_edge0 is None):
            # Can't collapse a half edge that doesn't exist.
            return None;

        # Don't collapse half_edge0 if its two endoints (vA,vB) are
        #   in some cell, but edge (vA,vB) is not in the cell.
        if (self.IsIllegalEdgeCollapseH(half_edge0)):
            return None

        max_num_half_edges_around_edge =\
            half_edge0.FromVertex().NumHalfEdgesFrom() +\
            half_edge0.ToVertex().NumHalfEdgesFrom()

        v0 = half_edge0.FromVertex()
        v1 = half_edge0.ToVertex()

        # Update *.next_half_edge_around_edge links.
        self._LinkHalfEdgesAroundMergedEdges(v0,v1)

        # Update *.prev_half_edge_in_cell and *.next_half_edge_in_cell links.
        # Delete triangles cells containing edge (v0,v1).
        half_edge = half_edge0
        k = 0
        while True:
            cell = half_edge.Cell()
            prev_half_edge = half_edge.PrevHalfEdgeInCell()
            next_half_edge = half_edge.NextHalfEdgeInCell()

            if (cell.HalfEdge() == half_edge):
                cell.half_edge = next_half_edge

            if (cell.IsTriangle()):
                # Cell is a triangle.
                # Collapsing half edge removes cell.
                v2 = prev_half_edge.FromVertex()
                self._DeleteHalfEdge(prev_half_edge)
                self._DeleteHalfEdge(next_half_edge)
                self._DeleteCell(cell)

                v2.MoveBoundaryHalfEdgeToIncidentHalfEdge0()

            else:
                # Link previous and next half edge.
                self._RelinkHalfEdgesInCell(prev_half_edge, next_half_edge)

                cell.num_vertices = cell.num_vertices-1

            half_edge = half_edge.NextHalfEdgeAroundEdge()
            k = k+1

            if (k >= max_num_half_edges_around_edge) or\
                (half_edge is half_edge0):
                break

        # Compute an upper bound on the number of half edges
        #  with endpoints (v0,v1).
        max_numh0 = v0.NumHalfEdgesFrom() + v1.NumHalfEdgesFrom()

        # Move all half edges from v0 to be from v1.
        # - Moves v0.half_edge_from to v1.half_edge_from
        self._MoveVertexHalfEdgeFromList(v0, v1)

        # Set v1 to midpoint of (v0,v1).
        compute_midpoint(v0.coord, v1.coord, v1.coord)

        # Delete half edges around edge (v0,v1).
        self._DeleteHalfEdgesAroundEdge(half_edge0, max_numh0)

        self._DeleteVertex(v0)

        v1.MoveBoundaryHalfEdgeToIncidentHalfEdge0()

        return v1


    ## Split cell with diagonal connecting the two from vertices.
    #  - Diagonal (half_edgeA.FromVertex(), half_edgeB.FromVertex()).
    def SplitCell(self, ihalf_edgeA, ihalf_edgeB):
        half_edgeA = self.HalfEdge(ihalf_edgeA)
        half_edgeB = self.HalfEdge(ihalf_edgeB)

        if (half_edgeA is None) or (half_edgeB is None):
            raise Exception("Programming error. Arguments to SplitCell are not half edge indices.")

        if self.IsIllegalSplitCell(half_edgeA, half_edgeB):
            return None

        vA = half_edgeA.FromVertex()
        vB = half_edgeB.FromVertex()
        half_edgeC = self.FindEdge(vA, vB)

        cellA = half_edgeA.Cell()
        numvA = cellA.NumVertices()    # Store before adding diagA to mesh
        icellB = self.MaxCellIndex()+1
        cellB = self._AddCell(icellB)
        idiagA = self.MaxHalfEdgeIndex()+1
        diagA = self._AddHalfEdge(idiagA, cellA, vB)
        idiagB = self.MaxHalfEdgeIndex()+1
        diagB = self._AddHalfEdge(idiagB, cellB, vA)

        # Link diagA and diagB around edge.
        diagA.next_half_edge_around_edge = diagB
        diagB.next_half_edge_around_edge = diagA

        if not(half_edgeC is None):
            # Link half_edge_around_edge cycle of half_edgeC and diagA/diagB
            self._SwapNextHalfEdgeAroundEdge(half_edgeC, diagA)

        # Add diagA and diagB to vertex half_edge_from[] lists.
        diagA.from_vertex.half_edge_from.append(diagA)
        diagB.from_vertex.half_edge_from.append(diagB)

        # Change cell of half edges from half_edgeB to half_edgeA.
        half_edge = half_edgeB
        k = 0
        while (k < numvA) and not(half_edge is half_edgeA):
            half_edge.cell = cellB
            half_edge = half_edge.NextHalfEdgeInCell()
            k = k+1

        # Set num_vertices in cellA and cellB
        cellB.num_vertices = k+1
        cellA.num_vertices = numvA+1-k

        # Set cellB.half_edge.
        cellB.half_edge = half_edgeB

        # Change cellA.half_edge, if necessary.
        if not(cellA.HalfEdge().Cell() is cellA):
            cellA.half_edge = half_edgeA

        hprevA = half_edgeA.PrevHalfEdgeInCell()
        hprevB = half_edgeB.PrevHalfEdgeInCell()

        # Link half edges in cell.
        self._RelinkHalfEdgesInCell(hprevB, diagA)
        self._RelinkHalfEdgesInCell(diagA, half_edgeA)
        self._RelinkHalfEdgesInCell(hprevA, diagB)
        self._RelinkHalfEdgesInCell(diagB, half_edgeB)

        # Swap first and last edges in half_edge_list[], if necessary.
        # diagA and diagB are not boundary edges, but
        #  diagA.PrevEdgeInCell() or diagB.PrevEdgeInCell() could
        #  be boundary edges.
        diagA.FromVertex()._ProcessFirstLastHalfEdgesFrom()
        diagB.FromVertex()._ProcessFirstLastHalfEdgesFrom()

        return diagA


    ## Joint two cells sharing an edge.
    #  - Returns edge incident on the joined cell.
    def JoinTwoCells(self, ihalf_edgeA):
        half_edgeA = self.HalfEdge(ihalf_edgeA)
        if (half_edgeA is None):
            raise Exception("Programming error. Argument to JoinTwoCells is not a cell index.")

        if (self.IsIllegalJoinCells(half_edgeA)):
            return None

        half_edgeB = half_edgeA.NextHalfEdgeAroundEdge()
        if (half_edgeB.NextHalfEdgeAroundEdge() != half_edgeA):
            raise Exception("Programming error. Half edge passed to JoinToCells is in an edge shared by three or more cells.")

        vA = half_edgeA.FromVertex()
        vB = half_edgeB.FromVertex()
        cellA = half_edgeA.Cell()
        cellB = half_edgeB.Cell()
        numvA = cellA.NumVertices()
        numvB = cellB.NumVertices()
        prevA = half_edgeA.PrevHalfEdgeInCell()
        prevB = half_edgeB.PrevHalfEdgeInCell()
        nextA = half_edgeA.NextHalfEdgeInCell()
        nextB = half_edgeB.NextHalfEdgeInCell()

        if not(vA.IsIncidentOnMoreThanTwoEdges()) or\
            not(vB.IsIncidentOnMoreThanTwoEdges()):
            # Can't remove an edge if some endpoint only has degree 2.
            return None

        # Change cellA.HalfEdge() if necessary.
        if (cellA.HalfEdge() is half_edgeA):
            cellA.half_edge = half_edgeA.NextHalfEdgeInCell()

        # Change edges in cellB to be in cellA.
        half_edge = half_edgeB.NextHalfEdgeInCell()
        for k in range(0, numvB-1):
            half_edge.cell = cellA
            half_edge = half_edge.NextHalfEdgeInCell()

        # Set number of vertices in cell.
        cellA.num_vertices = numvA + numvB - 2

        # Relink half edges in cell.
        self._RelinkHalfEdgesInCell(prevA, nextB)
        self._RelinkHalfEdgesInCell(prevB, nextA)

        # Delete cellB and half_edgeA and half_edgeB.
        self._DeleteHalfEdge(half_edgeA)
        self._DeleteHalfEdge(half_edgeB)
        self._DeleteCell(cellB)

        # half_edgeA and half_edgeB are not boundary edges,
        #   but the previous edges in the cell might be boundary edges.
        vA.MoveBoundaryHalfEdgeToIncidentHalfEdge0()
        vB.MoveBoundaryHalfEdgeToIncidentHalfEdge0()

        return nextA


    ## Split edge at midpoint.
    #  - Returns new vertex.
    def SplitEdge(self, ihalf_edgeA):
        half_edgeA = self.HalfEdge(ihalf_edgeA)
        if (half_edgeA is None):
            raise Exception\
                ("Programming error. Argument to SplitEdge is not a half edge index.")

        if (half_edgeA.IsBoundary()):
            return self._SplitBoundaryEdge(half_edgeA)
        else:
            return self._SplitInternalEdge(half_edgeA)


    # *** Functions to check potential edge collapses. ***

    ## Return True if half edge endpoints and v are in a mesh triangle.
    def IsInTriangle(self, half_edge0, v):
        # Cannot have more than max_numh half edges around an edge.
        max_numh = half_edge0.FromVertex().NumHalfEdgesFrom() +\
                    half_edge0.ToVertex().NumHalfEdgesFrom()

        half_edge = half_edge0
        k = 0
        while (True):
            if (half_edge.Cell().IsTriangle()):
                prev_half_edge = half_edge.PrevHalfEdgeInCell()

                if (prev_half_edge.FromVertex() is v):
                    return True

            half_edge = half_edge.NextHalfEdgeAroundEdge()
            k = k+1

            if (k >= max_numh) or (half_edge is half_edge0):
                break

        return False

    ## Return True if both endpoints (vfrom,vto) of half_edge
    #    are neighbors of some vertex vC, but (vfrom, vto, vC)
    #    is not a mesh triangle.
    #  - Returns also ivC, the index of the third vertex vC.
    def FindTriangleHole(self, half_edge):
        # Initialize.
        ivC = 0

        vfrom = half_edge.FromVertex()
        vto = half_edge.ToVertex()
        vto.ClearVisitedFlagsInAdjacentVertices()
        vfrom.SetVisitedFlagsInAdjacentVertices(True)

        for k in range(0, vto.NumHalfEdgesFrom()):
            half_edgeA = vto.KthHalfEdgeFrom(k)
            vtoA = half_edgeA.ToVertex()

            if (vtoA.IsVisited()):
                # vtoA is a neighbor of vfrom and vto.
                if not(self.IsInTriangle(half_edge, half_edgeA.ToVertex())):
                    ivC = vtoA.Index()
                    return True, ivC

            half_edgeB = half_edgeA.PrevHalfEdgeInCell()
            vfromB = half_edgeB.FromVertex()

            ## Check vfromB to handle boundary edges and/or cells
            #    with arbitrary orientations.
            if vfromB.IsVisited():
                # vfromB is a neighbor of vfrom and vto.
                if not(self.IsInTriangle(half_edge, vfromB)):
                    ivC = vfromB.Index()
                    return True, ivC

        return False, 0


    ## Return True if cell icell is a triangle whose 3 edges
    #    are boundary edges.
    def IsIsolatedTriangle(self, icell):
        THREE = 3

        cell = self.Cell(icell)
        if (cell is None):
            return False

        if not(cell.IsTriangle()):
            return False

        half_edge = cell.HalfEdge()
        for i in range(0,THREE):
            if not(half_edge.IsBoundary()):
                return False

        # Cell has three vertices (and three edges) and all edges
        #   are boundary edges.
        return True;


    ## Return True if cell icell is in the boundary of a tetrahedron.
    def IsInTetrahedron(self, icell):
        cell0 = self.Cell(icell)
        if (cell0 is None):
            return False

        if not(cell0.IsTriangle()):
            return False

        half_edge0 = cell0.HalfEdge()
        v2 = half_edge0.PrevHalfEdgeInCell().FromVertex()

        # Cannot have more than max_numh half edges around an edge.
        max_numh = half_edge0.FromVertex().NumHalfEdgesFrom() +\
                    half_edge0.ToVertex().NumHalfEdgesFrom()

        half_edge = half_edge0.NextHalfEdgeAroundEdge()
        k = 0
        while (k < max_numh and not(half_edge is half_edge0)):
            cell = half_edge.Cell()
            if cell.IsTriangle():
                prev_half_edge = half_edge.PrevHalfEdgeInCell()
                next_half_edge = half_edge.NextHalfEdgeInCell()

                if self.IsInTriangle(prev_half_edge, v2) and\
                    self.IsInTriangle(next_half_edge, v2):
                    # cell0, cell, and two triangles form a tetrahedron.
                    return True

            k = k+1
            half_edge = half_edge.NextHalfEdgeAroundEdge()

        return False


    ## Count number of vertices shared by two cells.
    def CountNumVerticesSharedByTwoCells(self, cellA, cellB):
        num_shared_vertices = 0

        cellB.ClearVisitedFlagsInAllVertices()
        cellA.SetVisitedFlagsInAllVertices(True)

        half_edgeB = cellB.HalfEdge()
        for k in range(0, cellB.NumVertices()):
            v = half_edgeB.FromVertex()
            if (v.IsVisited()):
                num_shared_vertices = num_shared_vertices+1

            half_edgeB = half_edgeB.NextHalfEdgeInCell()

        return num_shared_vertices


    ## Return True if edge collapse is illegal.
    # - Edge collapse (vA,vB) is illegal if some cell contains
    #   both vA and vB but not edge (vA,vB).
    # - Version that takes two vertices.
    # - NOTE: Function suffix is 'V' (for vertex arguments).
    def IsIllegalEdgeCollapseV(self, vA, vB):
        if (vA.NumHalfEdgesFrom() > vB.NumHalfEdgesFrom()):
            # Swap vA and vB to reduce number of cells processed.
            return self.IsIllegalEdgeCollapseV(vB, vA)
        else:
            for k in range(0, vA.NumHalfEdgesFrom()):
                half_edge0 = vA.KthHalfEdgeFrom(k)
                cell = half_edge0.Cell()
                if (cell.NumVertices() < 4):
                    # All pairs of cell vertices form an edge.
                    continue

                half_edge =\
                    (half_edge0.NextHalfEdgeInCell()).NextHalfEdgeInCell()

                for i in range(2,cell.NumVertices()-1):
                    if (half_edge.FromVertex() == vB):
                        return True

            return False


    ## Return True if edge collapse is illegal.
    # - NOTE: Function suffix is 'H' (for half_edge argument).
    def IsIllegalEdgeCollapseH(self, half_edge):
        return self.IsIllegalEdgeCollapseV\
                    (half_edge.FromVertex(), half_edge.ToVertex())


    ## Return True if split cell is illegal.
    # - Split cell is illegal
    #     if half_edgeA and half_edgeB are in different cells or
    #     if half_edgeA.FromVertex() and half_edgeB.FromVertex()
    #       are adjacent vertices.
    def IsIllegalSplitCell(self, half_edgeA, half_edgeB):
        if not(half_edgeA.Cell() is half_edgeB.Cell()):
            return True

        if (half_edgeA is half_edgeB):
            return True

        if (half_edgeA.FromVertex() is half_edgeB.ToVertex()):
            return True

        if (half_edgeA.ToVertex() is half_edgeB.FromVertex()):
            return True

        return False


    ## Return True if join cells is illegal.
    #  - Join cells is illegal if half_edge is a boundary half edge
    #    or more than two cells are incident on the edge
    #    or some endpoint of half edge has degree 2.
    def IsIllegalJoinCells(self, half_edge):
        TWO = 2

        if (half_edge.IsBoundary()):
            return True

        if not(half_edge.FromVertex().IsIncidentOnMoreThanTwoEdges()):
            return True

        if not(half_edge.ToVertex().IsIncidentOnMoreThanTwoEdges()):
            return True

        half_edgeX = half_edge.NextHalfEdgeAroundEdge()
        if not(half_edge is half_edgeX.NextHalfEdgeAroundEdge()):
            # More than two cells are incident on edge
            #  (half_edge.FromVertex(), half_edge.ToVertex()).
            return True

        if (self.CountNumVerticesSharedByTwoCells\
                (half_edge.Cell(), half_edgeX.Cell()) > TWO):
                # Cells share more than two vertices.
                return True

        # Join is LEGAL
        return False


    # *** Compute mesh information. ***

    ## Return min and max squared edge lengths over all mesh edges.
    #  - Return also half edges with min and max edge lengths.
    def ComputeMinMaxEdgeLengthSquared(self):

        flag_found = False;

        # Initialize
        min_edge_length_squared = 0.0;
        max_edge_length_squared = 0.0;
        ihalf_edge_min = 0;
        ihalf_edge_max = 0;

        for ihalf_edge in self.HalfEdgeIndices():
            half_edge = self.HalfEdge(ihalf_edge)
            if (half_edge is None):
                # Shouldn't happen but just in case.
                continue;

            length_squared = half_edge.ComputeLengthSquared()
            if (not(flag_found) or (length_squared < min_edge_length_squared)):
                min_edge_length_squared = length_squared
                ihalf_edge_min = half_edge.Index()

            if (not(flag_found) or length_squared > max_edge_length_squared):
                max_edge_length_squared = length_squared
                ihalf_edge_max = half_edge.Index()

            flag_found = True

        return min_edge_length_squared, max_edge_length_squared,\
                ihalf_edge_min, ihalf_edge_max


    ## Return min squared ratio of the min to max edge in any cell.
    #  - Ignores cells with all edge lengths 0.
    #  - Return also cell index and length and indices
    #      of shortest and longest half edges in the cell.
    #  - Returns 1.0 if there are no cells or all edges are length 0.
    def ComputeMinCellEdgeLengthRatioSquared(self):

        # Initialize.
        min_edge_length_ratio_squared = 1.0
        icell_min_ratio = 0
        min_edge_length_squared = 0.0
        max_edge_length_squared = 0.0
        ihalf_edge_min = 0
        ihalf_edge_max = 0

        for icell in self.CellIndices():
            cell = self.Cell(icell);
            if (cell is None):
                # Shouldn't happen but just in case.
                continue

            min_Lsquared, max_Lsquared, ihalf_min, ihalf_max =\
                cell.ComputeMinMaxEdgeLengthSquared()

            if (max_Lsquared == 0 or cell.NumVertices() == 0):
                continue

            ratio = min_Lsquared/max_Lsquared;
            if (ratio < min_edge_length_ratio_squared):
                min_edge_length_ratio_squared = ratio
                icell_min_ratio = icell
                min_edge_length_squared = min_Lsquared
                max_edge_length_squared = max_Lsquared
                ihalf_edge_min = ihalf_min
                ihalf_edge_max = ihalf_max

        return min_edge_length_ratio_squared, icell_min_ratio,\
                min_edge_length_squared, max_edge_length_squared,\
                ihalf_edge_min, ihalf_edge_max


    ## Compute cosine min and max cell angles over all mesh cells.
    #  - Return also half edges whose from vertices are incident
    #      on the two cell edges forming the min and max angles.
    #  - Note: The smallest angle has the largest cosine and
    #      the largest angle has the smallest cosine.
    def ComputeCosMinMaxAngle(self):

        # Initialize
        cos_min_angle = 0.0;
        cos_max_angle = 0.0;
        ihalf_edge_min = 0;
        ihalf_edge_max = 0;

        flag_found = False;
        for icell in self.CellIndices():
            cell = self.Cell(icell);
            if (cell is None):
                # Shouldn't happen but just in case.
                continue

            cos_minA, cos_maxA, ihalf_min, ihalf_max =\
                cell.ComputeCosMinMaxAngle()

            if not(flag_found) or (cos_minA > cos_min_angle):
                cos_min_angle = cos_minA
                ihalf_edge_min = ihalf_min

            if not(flag_found) or (cos_maxA < cos_max_angle):
                cos_max_angle = cos_maxA
                ihalf_edge_max = ihalf_max

            flag_found = True

        return cos_min_angle, cos_max_angle, ihalf_edge_min, ihalf_edge_max
