## \file half_edge_mesh_elements.py
# Mesh elements (vertices, edges, cells) for half edge mesh data structure.
# Classes for vertices, edges, cells of half edge mesh.
#
# \mainpage Half Edge Mesh Elements (Python implementation):
# - Class HMESH_VERTEX_BASE: Half edge mesh vertices.
# - Class HMESH_HALF_EDGE_BASE: Mesh half edges.
# - Class HMESH_CELL_BASE: Half edge mesh cells (2-cells).
# The mesh is stored in HALF_EDGE_MESH_BASE, including lists
# - All allocations of vertices, half edges and cells should be done
#   in HALF_EDGE_MESH_BASE or a subclass of HALF_EDGE_MESH_BASE.
# - Each vertex, half edge and cell can be identified by a pointer
#   to the object containing the vertex, half edge or cell, or
#   by an integer index (identifier) of the vertex, half edge or cell.
# - This is NOT a very efficient/compact implementation of half edges.
# - This implementation is meant to be simple and (hopefully) robust
#   for use in OSU CSE homeworks and prototypes.
# - Note: Many of the simpler get functions do not check their arguments,
#   e.g. Objects that are None or indices in range.
#    Such checks would be too time consuming for large meshes.
#    The calling function is responsible to ensure that objects
#    not None and indices are in a specified range.
# - Version 0.0.3

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

class HMESH_VERTEX_BASE:

    ## Dimension = Number of coordinates.
    # - All vertices have the same dimension.
    _DIM = 3;

    def Dimension(self) :
        return self._DIM

    def KthCoord(self, k) :
        return self.coord[k]

    def Index(self) :
        return self.index

    ## Return number of half edges with from vertex
    #    equal to this.
    def NumHalfEdgesFrom(self):
        return len(self.half_edge_from)

    ## Return k'th incident half edge with from vertex equal to this.
    #  @pre k < NumHalfEdgesFrom() and k >= 0.
    def KthHalfEdgeFrom(self, k):
        return self.half_edge_from[k]


    ## Return true if vertex is on the boundary or vertex is
    #    is not incident on any cells.
    def IsBoundary(self):
        if (self.NumHalfEdgesFrom() == 0):
            # Vertex is not incident on any cells.
            return True

        # Check PrevHalfEdgeInCell() in case cells are inconsistently oriented
        #  and both boundary half edges point to vertex.
        return (self.KthHalfEdgeFrom(0).IsBoundary() or
                self.KthHalfEdgeFrom(0).PrevHalfEdgeInCell().IsBoundary());


    ## Find half edge with from vertex self and
    #     ToVertexIndex() iv.
    #  - Return None if no half edge found.
    def FindHalfEdgeTo(self, iv):
        for k in range(0, self.NumHalfEdgesFrom()):
            half_edge = self.half_edge_from[k]
            if (half_edge.ToVertex()).Index() == iv:
                return(half_edge)

        return(None)

    ## Count number of half edges with from vertex self
    #    and ToVertexIndex() iv.
    def CountNumIncidentHalfEdges(self, iv):
        num = 0
        for half_edge in self.half_edge_from:
            if (half_edge.ToVertexIndex() == iv):
                num = num + 1

        return num

    ## Swap half edges in list half_edge_from[].
    #  - Internal (protected) routine.
    def _SwapHalfEdgesInHalfEdgeFromList(self, k0, k1):
        temp = self.half_edge_from[k0]
        self.half_edge_from[k0] = self.half_edge_from[k1]
        self.half_edge_from[k1] = temp


    ## Move boundary half edge to half_edge_from[0].
    #  - If there are no boundary half edges in half_edge_from[],
    #    but half_edge_from[k].PrevHalfEdgeInCell() is a boundary
    #    half edge, move half_edge_from[k] to half_edge_from[0].
    #  - If cells are inconsistenly oriented and there are no boundary
    #    half edges in half_edge_from[] but
    #    half_edge_from[k].PrevHalfEdgeInCell() is a boundary half edge
    #    move half_edge_from[k] to half_edge_from[0].
    #  - Does nothing if half_edge_from[0] is a boundary half edge
    #    or if vertex is not incident on any boundary half edges (from or to).
    def MoveBoundaryHalfEdgeToHalfEdgeFrom0(self):

        if (self.NumHalfEdgesFrom() < 1):
            return

        if (self.half_edge_from[0]).IsBoundary():
            # half_edge_from[0] is already a boundary half edge.
            # Do nothing.
            return

        for k in range(1,self.NumHalfEdgesFrom()):
            half_edge = self.half_edge_from[k]

            if (half_edge.IsBoundary()):
                self._SwapHalfEdgesInHalfEdgeFromList(0,k)
                return

        # No boundary half edges found.

        # Extra processing in case cells are inconsistently oriented.
        # Check if half_edge_from[k].PrevHalfEdgeInCell()
        #   is a boundary half edge for some k.
        prev_half_edge0 = self.half_edge_from[0]

        if (prev_half_edge0.IsBoundary()):
            # prev_half_edge0 is already a boundary half edge.
            # Do nothing.
            return

        for k in range(1,self.NumHalfEdgesFrom()):
            half_edge = self.half_edge_from[k]
            if (half_edge.PrevHalfEdgeInCell().IsBoundary):
                self._SwapHalfEdgesInHalfEdgeFromList(0,k)
                return


    ## Return string containing coordinates separated by spaces.
    def CoordStr(self):
        if (self.Dimension() < 1):
            return ""

        s = str(self.KthCoord(0))
        for ic in range(1,self.Dimension()):
            s = s + " " + str(self.KthCoord(ic))

        return s

    ## Initialize
    def __init__(self):

        ## Unique non-negative integer identifying the vertex.
        self.index = None

        ## List of all half edges originating at vertex in case mesh
        #  is not a manifold and cells incident on vertex do not
        #  form a fan.
        self.half_edge_from = []

        ## Vertex coordinates.
        self.coord = []

        for ic in range(0,self.Dimension()):
            self.coord.append(0)


class HMESH_HALF_EDGE_BASE:

    def Index(self):
        return self.index

    ## Return previous half edge in cell.
    def PrevHalfEdgeInCell(self):
        return self.prev_half_edge_in_cell

    ## Return next half edge in cell.
    def NextHalfEdgeInCell(self):
        return self.next_half_edge_in_cell

    ## Return next half edge around edge.
    def NextHalfEdgeAroundEdge(self):
        return self.next_half_edge_around_edge

    ## Return previous half edge around edge.
    #  Takes time proportional to number of half edges around edge.
    def PrevHalfEdgeAroundEdge(self):

        #  - max_numh is an upper bound on the number of half edges
        #      around the edge containing the half edge.
        max_numh = self.FromVertex().NumHalfEdgesFrom() +\
                    self.ToVertex().NumHalfEdgesFrom()

        half_edge = self.NextHalfEdgeAroundEdge()
        for k in range(0,max_numh):
            if (self is half_edge.NextHalfEdgeAroundEdge()):
                return half_edge

            half_edge = half_edge.NextHalfEdgeAroundEdge()

        if (self is not half_edge.NextHalfEdgeAroundEdge()):
            raise Exception \
            ("Half edge mesh data structure error.\nError in linked list of half edges around edge.")


    ## Return from vertex.
    def FromVertex(self):
        return self.from_vertex

    ## Return to vertex.
    def ToVertex(self):
        return (self.next_half_edge_in_cell).FromVertex()

    ## Return index of from vertex.
    def FromVertexIndex(self):
        return (self.FromVertex()).Index()

    ## Return index of to vertex.
    def ToVertexIndex(self):
        return (self.ToVertex()).Index()

    ## Return true if half edge is boundary half edge.
    def IsBoundary(self):
        return (self == self.next_half_edge_around_edge)

    ## Return cell containing half edge.
    def Cell(self):
        return self.cell

    ## Return index of cell.
    #  - Added: 12-02-2021 - RW
    def CellIndex(self):
        return self.Cell().Index()


    ## Count number of half edges around edge.
    def CountNumHalfEdgesAroundEdge(self):
        # Cannot have more than max_num half edges around an edge.
        max_num = self.FromVertex().NumHalfEdgesFrom() +\
            self.ToVertex().NumHalfEdgesFrom()

        num = 1
        half_edge = self.NextHalfEdgeAroundEdge()
        while (half_edge is not self) and (num <= max_num):
            half_edge = half_edge.NextHalfEdgeAroundEdge()
            num = num + 1

        return num

    ## Return half edge with minimum index in cycle of half edges around edge.
    #  - Useful in getting a unique half edge representing an edge.
    #  - Could return self.
    def MinIndexHalfEdgeAroundEdge(self):

        # Cannot have more than max_numh half edges around this edge.
        max_numh = self.FromVertex().NumHalfEdgesFrom() +\
                    self.ToVertex().NumHalfEdgesFrom()

        # Initialize
        min_index_half_edge = self;
        min_index = self.Index()

        half_edge = self.NextHalfEdgeAroundEdge()
        k = 0;
        while (half_edge != self) and (k < max_numh):
            if (half_edge.Index() < min_index):
                min_index_half_edge = half_edge;
                min_index = half_edge.Index()

            half_edge = half_edge.NextHalfEdgeAroundEdge()
            k = k+1

        return min_index_half_edge


    ## Return true if half_edgeB has same endpoints as current half_edge (self).
    def SameEndpoints(self, half_edgeB):
        if (self.FromVertex() is half_edgeB.ToVertex() and
            self.ToVertex() is half_edgeB.FromVertex()):
            return True

        if (self.FromVertex() is half_edgeB.FromVertex() and
            self.ToVertex() is half_edgeB.ToVertex()):
            return True

        return False

    ## Return previous half edge around from vertex.
    # - Returns PrevHalfEdgeInCell() if PrevHalfEdgeInCell()
    #   is a boundary half edge.
    # - NextHalfEdgeAroundFromVertex() is not defined, since
    #   PrevHalfEdgeAroundFromVertex() should always be used in moving
    #   around a vertex.
    def PrevHalfEdgeAroundFromVertex(self):
        return self.PrevHalfEdgeInCell().NextHalfEdgeAroundEdge()

    ## Returns previous half edge around vertex iv.
    #  - Returns PrevHalfEdgeInCell() if PrevHalfEdgeInCell()
    #    is a boundary half edge.
    #  - Note: If iv == ToVertexIndex(), returns
    #    NextHalfEdgeInCell().NextHalfEdgeAroundEdge() so that
    #    repeated calls to PrevHalfEdgeAroundVertex() move in a
    #    consistent direction.
    def PrevHalfEdgeAroundVertex(self, iv):
        if (self.FromVertexIndex() == iv):
            return self.PrevHalfEdgeAroundFromVertex()
        else:
            return self.NextHalfEdgeInCell().NextHalfEdgeAroundEdge()


    ## Compute coordinates of edge midpoint.
    def ComputeMidpointCoord(self):
        ivfrom = self.FromVertex()
        ivto = self.ToVertex()
        midpoint_coord = []
        for ic in range(len(ivfrom.coord)):
            midc = (ivfrom.coord[ic]+ivto.coord[ic])/2.0
            midpoint_coord.append(midc)

        return midpoint_coord


    ## Return string of half edge endpoints.
    def EndpointsStr(self, separator):
        s = str(self.FromVertexIndex()) + separator
        if (self.NextHalfEdgeInCell() is None):
            # Cannot determine ToVertex(). Replace vertex index with "*".
            s = s + "*"
        else:
            s = s + str(self.ToVertexIndex())
        return s

    ## Return string of half edge index and endpoints.
    def IndexAndEndpointsStr(self, separator):
        s = str(self.Index()) + " (" + self.EndpointsStr(separator) + ")"
        return s

    def __init__(self):

        ## Unique non-negative integer identifying the half-edge.
        self.index = None;

        ## Next half edge in cell.
        self.next_half_edge_in_cell = None

        ## Previous half edge in cell.
        self.prev_half_edge_in_cell = None

        ## Next half edge around edge.
        self.next_half_edge_around_edge = self

        ## From vertex.
        self.from_vertex = None

        ## Cell containing half edge.
        self.cell = None


class HMESH_CELL_BASE:

    ## Return cell index.
    def Index(self):
        return(self.index)

    ## Return number of cell vertices.
    def NumVertices(self):
        return(self.num_vertices)

    ## Return one cell half edge.
    def HalfEdge(self):
        return(self.half_edge)

    ## Return true if cell has exactly 3 vertices.
    def IsTriangle(self):
        return (self.NumVertices() == 3)


    ## Return list of cell vertex indices.
    def GetVertexIndices(self):
        ivlist = []
        half_edge = self.HalfEdge()
        for i in range(self.NumVertices()):
            ivlist.append(half_edge.FromVertexIndex())
            half_edge = half_edge.NextHalfEdgeInCell()

        return ivlist


    ## Compute centroid of cell.
    def ComputeCentroid(self):
        if (self.NumVertices() < 1):
            raise Exception("Cannot compute centroid of cell with 0 vertices.")

        half_edge = self.HalfEdge()
        dimension = half_edge.FromVertex().Dimension()
        centroid_coord = [ 0.0 ] * dimension

        for i in range(self.NumVertices()):
            for ic in range(dimension):
                from_vertex = half_edge.FromVertex()
                centroid_coord[ic] = centroid_coord[ic] + from_vertex.coord[ic]
            half_edge = half_edge.NextHalfEdgeInCell()

        for ic in range(dimension):
            centroid_coord[ic] = centroid_coord[ic]/self.NumVertices()

        return centroid_coord


    ## Return string of cell vertices.
    def VerticesStr(self, separator):
        s = ""
        half_edge = self.HalfEdge()
        for i in range(self.NumVertices()):
            if (i > 0):
                s = s + ", "

            s = s + str(half_edge.FromVertexIndex())
            half_edge = half_edge.NextHalfEdgeInCell()

        return s


    ## Initialize
    def __init__(self):

        ## Unique non-negative integer identifying the cell.
        self.index = None

        ## Some half edge in the cell.
        self.half_edge = None

        ## Number of cell vertices.
        self.num_vertices = 0
