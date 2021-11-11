## \file half_edge_mesh.py
# Half edge mesh data structure.
# Classes for half edge mesh.
# This is an implementation of a half edge mesh.
#
# \mainpage Half Edge Mesh (Python implementation):
# The mesh is stored in HALF_EDGE_MESH_BASE, including lists
# of all the vertices, half edges and cells in the mesh.
# - Each vertex, half edge and cell is in its own class.
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

class VERTEX_BASE:

    ## Dimension = Number of coordinates.
    # - All vertices have the same dimension.
    _DIM = 3;

    def Dimension(self) :
        return self._DIM

    def Coord(self, ic) :
        return self.coord[ic]

    def Index(self) :
        return self.index

    def NumHalfEdgesFrom(self):
        return len(self.half_edge_from)

    def KthHalfEdgeFrom(self, k):
        return self.half_edge_from[k]

    ## Find half edge with from vertex self and
    #     ToVertexIndex() iv.
    def FindIncidentHalfEdge(self, iv):
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

    def MoveBoundaryHalfEdgeToIncidentHalfEdge0(self):
        if (self.NumHalfEdgesFrom() < 1):
            return

        if (self.half_edge_from[0]).IsBoundary():
            return

        for k in range(1,self.NumHalfEdgesFrom()):
            half_edge = self.half_edge_from[k]
            if (half_edge.IsBoundary()):
                temp = self.half_edge_from[0]
                self.half_edge_from[0] = half_edge
                self.half_edge_from[k] = temp
                return

    def CoordStr(self):
        if (self.Dimension() < 1):
            return ""

        s = str(self.Coord(0))
        for ic in range(1,self.Dimension()):
            s = s + " " + str(self.Coord(ic))

        return s

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


class HALF_EDGE_BASE:

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
        return(self == self.next_half_edge_around_edge)

    ## Return cell containing half edge.
    def Cell(self):
        return(self.cell)

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


class CELL_BASE:

    def Index(self):
        return(self.index)

    def NumVertices(self):
        return(self.num_vertices)

    def HalfEdge(self):
        return(self.half_edge)

    def __init__(self):

        ## Unique non-negative integer identifying the cell.
        self.index = None

        ## Some half edge in the cell.
        self.half_edge = None

        ## Number of cell vertices.
        self.num_vertices = 0


# Pass vertex, half edge and cell classes to HALF_EDGE_MESH_BASE.
class HALF_EDGE_MESH_BASE:

    # Private member functions.

    def __init__(self, classV, classHE, classC):

        ## Set class types of vertices, half edges and cells.
        self.VERTEX_TYPE = classV;
        self.HALF_EDGE_TYPE = classHE;
        self.CELL_TYPE = classC;

        ## Dictionary of vertices.
        self._vertex_dict = dict()

        ## Dictionary of half edges.
        self._half_edge_dict = dict()

        ## Dictionary of cells.
        self._cell_dict = dict()

        ## Upper bound on the vertex index.
        #  - Could be greater than the maximum index if some vertices are deleted.
        self._max_vertex_index = -1

        ## Upper bound on the half edge index.
        # - Could be greater than the maximum index if some half edges are deleted.
        self._max_half_edge_index = -1

        ## Upper bound on the cell index.
        # - Could be greater than the maximum index if some cells are deleted.
        self._max_cell_index = -1


    ## Create vertex with index iv, if vertex iv does not exist.
    # - Return vertex.
    # = Returns vertex, even if vertex already exists.
    def _CreateVertex(self,iv):
        if (iv < 0):
            raise Exception("Illegal argument to _CreateVertex().  Vertex index must be non-negative.")

        self._max_vertex_index = max(self._max_vertex_index, iv)

        if (iv in self._vertex_dict):
            return self._vertex_dict[iv]

        v = self.VERTEX_TYPE()
        v.index = iv
        self._vertex_dict[iv] = v

        return v


    ## Add half edge with index ihalf_edge.
    #  - Raises an exception if some half edge already has index ihalf_edge.
    #  - Use MaxHalfEdgeIndex() to find an unused half edge index.
    def _AddHalfEdge(self, ihalf_edge, cell, vfrom):

        if (ihalf_edge in self.HalfEdgeIndices()):
            raise Exception("Illegal argument to AddHalfEdge(). Half edge with index " +\
                str(ihalf_edge) + " already exists.")

        half_edge = self.HALF_EDGE_TYPE()
        half_edge.index = ihalf_edge
        self._half_edge_dict[ihalf_edge] = half_edge
        half_edge.cell = cell
        cell.num_vertices = cell.num_vertices + 1
        half_edge.from_vertex = vfrom
        self._max_half_edge_index = max(self._max_half_edge_index, ihalf_edge)

        return half_edge


    def _AddAndLinkHalfEdge(self, ihalf_edge, cell, vfrom, vto, hprev):

        half_edge = self._AddHalfEdge(ihalf_edge, cell, vfrom)

        half_edge.prev_half_edge_in_cell = hprev
        if (hprev is not None):
            hprev.next_half_edge_in_cell = half_edge

        half_edgeB = vto.FindIncidentHalfEdge(vfrom.Index())
        if (half_edgeB is None):
            # Check whether there is a half edge around the edges
            #   with same orientation as half_edge.
            half_edgeC = vfrom.FindIncidentHalfEdge(vto.Index())
            if (half_edgeC is None):
                half_edge.next_half_edge_around_edge = half_edge
            else:
                self._LinkHalfEdgesAroundEdge(half_edgeC, half_edge)
        else:
            self._LinkHalfEdgesAroundEdge(half_edgeB, half_edge)

        vfrom.half_edge_from.append(half_edge)
        vfrom.MoveBoundaryHalfEdgeToIncidentHalfEdge0()
        vto.MoveBoundaryHalfEdgeToIncidentHalfEdge0()

        return half_edge


    def _LinkHalfEdgesInCell(self, hprev, hnext):
        if (hprev.Cell() != hnext.Cell()):
            raise Exception("Link of half edges failed in _LinkHalfEdgesInCell.  Half edges are in different cells.")

        hprev.next_half_edge_in_cell = hnext
        hnext.prev_half_edge_in_cell = hprev

    def _LinkHalfEdgesAroundEdge(self, half_edgeA, half_edgeB):
        half_edgeC = half_edgeA.next_half_edge_around_edge;
        half_edgeA.next_half_edge_around_edge = half_edgeB;
        half_edgeB.next_half_edge_around_edge = half_edgeC;


    ## Add cell with index icell.
    #  - Raises an exception if some cell already has index icell.
    #  - Use MaxCellIndex() to find an unused cell index.
    def _AddCell(self, icell):
        cell = self.CELL_TYPE()
        cell.index = icell
        self._cell_dict[icell] = cell
        self._max_cell_index = max(self._max_cell_index, icell)

        return cell


    # Public member Functions.

    # Get functions.

    ## Return vertex with index iv.
    # - Returns None if no vertex has index iv.
    def Vertex(self, iv):
        v = self._vertex_dict.get(iv)
        return v

    ## Return list of vertex indices (_vertex_dict keys).
    #  - Note: This returns a view object, not a python list.
    #  - Use list(obj.VertexIndices()) to get a python list.
    def VertexIndices(self):
        return self._vertex_dict.keys()

    ## Return half edge with index ihalf_edge.
    # - Returns None if no half edge has index ihalf_edge.
    def HalfEdge(self, ihalf_edge):
        half_edge = self._half_edge_dict.get(ihalf_edge)
        return half_edge

    ## Return list of half edge indices (_half_edge_dict keys).
    #  - Note: This returns a view object, not a python list.
    #  - Use list(obj.HalfEdgeIndices()) to get a python list.
    def HalfEdgeIndices(self):
        return self._half_edge_dict.keys()

    ## Return cell with index icell.
    # - Returns None if no cell has index icell.
    def Cell(self, icell):
        return self._cell_dict.get(icell)

    ## Return list of cell indices (_cell_dict keys).
    #  - Note: This returns a view object, not a python list.
    #  - Use list(obj.CellIndices()) to get a python list.
    def CellIndices(self):
        return self._cell_dict.keys()


    ## Return number of vertices.
    # - Note: If vertices are deleted, MaxVertexIndex() may be greater than NumVertices().
    def NumVertices(self):
        return len(self._vertex_dict)

    def NumHalfEdges(self):
        return len(self._half_edge_dict)

    def NumCells(self):
        return len(self._cell_dict)

    ## Return max index (key) of vertices in _vertex_dict.
    #  - Return -1 if _vertex_dict is empty.
    def MaxVertexIndex(self):
        return self._max_vertex_index

    ## Return max index (key) of half edges in _half_edge_dict.
    #  - Return -1 if _half_edge_dict is empty.
    def MaxHalfEdgeIndex(self):
        return self._max_half_edge_index

    ## Return max index (key) of cells in _cell_dict.
    #  - Return -1 if _cell_dict is empty.
    def MaxCellIndex(self):
        return self._max_cell_index


    # Functions to add vertices or cells.

    ## Add vertex with index iv.
    #  @param iv Vertex index. Should not already be in VertexIndices().
    #  - Returns vertex.
    def AddVertex(self, iv):
        if (iv in self._vertex_dict):
            raise Exception("Illegal argument to AddVertex(). Vertex " + str(iv) + " already exists.")

        return self._CreateVertex(iv)


    ## Add numv vertices ith indices from 0 to numv-1.
    #  @pre Half edge mesh has no vertices.
    def AddVertices(self, numv):
        if (len(self._vertex_dict) > 0):
            raise Exception("Error. Cannot call AddVertices() if mesh already has some vertices.")

        for iv in range(0,numv):
            self._CreateVertex(iv)


    ## Add cell with index icell.
    # param icell Cell index. Should not already be in CellIndices().
    def AddCell(self, icell, cell_vertex):

        if (len(cell_vertex) < 3):
            raise Exception("Illegal argument to AddCell(). List cell_vertex[] must have 3 or more vertices.")

        if (icell in self.CellIndices()):
            raise Exception("Illegal argument to AddCell(). Cell with index " +\
                str(icell) + " already exists.")

        cell = self._AddCell(icell)

        iv0 = cell_vertex[0]
        iv1 = cell_vertex[1]
        v0 = self._CreateVertex(iv0)
        v1 = self._CreateVertex(iv1)

        ihalf_edge0 = self.MaxHalfEdgeIndex() + 1
        half_edge0 = self._AddAndLinkHalfEdge(ihalf_edge0, cell, v0, v1, None)
        cell.half_edge = half_edge0

        hprev = half_edge0
        ihalf_edge = ihalf_edge0
        for i0 in range(1,len(cell_vertex)):
            i1 = (i0+1)%len(cell_vertex)
            iv0 = cell_vertex[i0]
            iv1 = cell_vertex[i1]
            v0 = self._vertex_dict[iv0]
            v1 = self._CreateVertex(iv1)
            ihalf_edge += 1
            half_edge = self._AddAndLinkHalfEdge(ihalf_edge, cell, v0, v1, hprev)
            hprev = half_edge

        # Lin last half edge (hprev) and first half edge (half_edge0)
        self._LinkHalfEdgesInCell(hprev, half_edge0)

        if (len(cell_vertex) != cell.NumVertices()):
            raise Exception("Error in AddCell().  Incorrect number of vertices in cell.")


    ## Set coordinates of iv to coord[].
    def SetCoord(self, iv, coord):
        v = self._vertex_dict[iv]

        if (v == None):
            raise Exception("Error in SetCoord(). Vertex " + str(iv) + " does not exist.")

        if (len(coord) < v.Dimension()):
            raise Exception("Error in SetCoord(). Argument coord[] has too few coordinates.")

        for ic in range(0,v.Dimension()):
            v.coord[ic] = coord[ic]


    # Check routines.

    ## Check data structure vertices.
    # - Return true if no problems found.
    # - Return index of problem vertex and error messge.
    def CheckVertices(self):

        # Initialize
        error_msg = None
        iv = 0

        ivmax = max(self.VertexIndices(),default=-1)
        if (self.MaxVertexIndex() < ivmax):
            error_msg = "Incorrect value (" + str(self.MaxVertexIndex()) +\
                ") of _max_vertex_index.  Max vertex is " + str(ivmax) + "."
            return False, ivmax, error_msg

        for jv in self._vertex_dict:
            iv = jv
            v = self.Vertex(iv)

            if (v.Index() != jv):
                error_msg = "Incorrect vertex index for vertex " + str(iv) + "."
                return False, iv, error_msg

            flag_boundary = False
            boundary_half_edge = None
            for k in range(0,v.NumHalfEdgesFrom()):
                half_edge = v.KthHalfEdgeFrom(k)
                if (half_edge is None):
                    error_msg = "Vertex " + str(iv) + " list half_edge_from[" + str(k) + "] = None."
                    return False, iv, error_msg

                if (half_edge.FromVertex() != v):
                    error_msg = "Error in list half_edge_from[] for vertex " + str(iv)\
                        + ". List incorrectly includes half edge " + int(half_edge.Index()) + "."
                    return False, iv, error_msg

                if (half_edge.IsBoundary()):
                    flag_boundary = True
                    boundary_half_edge = half_edge

            if (flag_boundary):
                # Note: v.half_edge_from[] must contain at least
                #   one half edge for flag_boundary to be true.
                half_edge = v.half_edge_from[0]
                if not(half_edge.IsBoundary()):
                    error_msg = "Vertex " + str(iv) + " is on a boundary half edge " + \
                        boundary_half_edge.IndexAndEndpointsStr(",") + \
                        " but first incident half edge " + \
                        " is not a boundary half edge."
                    return False, iv, error_msg

        return True, 0, None


    ## Check data structure half edges.
    # - Return true if no problems found.
    # - Return index of problem half edge and error message.
    def CheckHalfEdges(self):

        # Initialize
        error_msg = None
        ihalf_edge = 0

        ihmax = max(self.HalfEdgeIndices(),default=-1)
        if (self.MaxHalfEdgeIndex() < ihmax):
            error_msg = "Incorrect value (" + str(self.MaxHalfEdgeIndex()) +\
                ") of _max_half_edge_index.  Max half edge is " + str(ihmax) + "."
            return False, ihmax, error_msg

        for j in self._half_edge_dict:
            ihalf_edge = j

            half_edge = self.HalfEdge(ihalf_edge)

            if (half_edge.Index() != ihalf_edge):
                error_msg = "Incorrect half edge index for half edge " +\
                    str(ihalf_edge) + "."
                return False, ihalf_edge, error_msg

            v = half_edge.FromVertex()
            if (v is None):
                error_msg = "Missing (None) from vertex in half edge " +\
                    str(ihalf_edge) + "."
                return False, ihalf_edge, error_msg

            num_match = v.half_edge_from.count(half_edge)
            if (num_match < 1):
                error_msg = "Half edge " + half_edge.IndexAndEndpointsStr(",") +\
                    " does not appear in half_edge_from[] list for vertex " +\
                    str(v.Index()) + "."
                return False, ihalf_edge, error_msg
            elif (num_match > 1):
                error_msg = "Half edge appears more than once in half_edge_from[] list for vertex " +\
                    str(v.Index()) + "."
                return False, ihalf_edge, error_msg

            half_edge0 = v.KthHalfEdgeFrom(0)
            if half_edge.IsBoundary() and not(half_edge0.IsBoundary()):
                error_msg = "Half edge " + str(ihalf_edge) + " is a boundary half edge " +\
                    "but v.KthHalfEdgeFrom(0) is not a boundary half edge."
                return False, ihalf_edge, error_msg

            cell = half_edge.Cell()
            prev_half_edge = half_edge.PrevHalfEdgeInCell()
            next_half_edge = half_edge.NextHalfEdgeInCell()

            if (cell is None):
                error_msg = "Half edge " + str(ihalf_edge) +\
                    " missing cell containing half edge."
                return False, ihalf_edge, error_msg

            if (prev_half_edge is None):
                error_msg = "Half edge " + str(ihalf_edge) +\
                    " missing previous half edge in cell."
                return False, ihalf_edge, error_msg

            if (next_half_edge is None):
                error_msg = "Half edge " + str(ihalf_edge) +\
                    " missing next half edge in cell."
                return False, ihalf_edge, error_msg

            if (prev_half_edge.Cell() != cell):
                error_msg = "Half edge " + str(ihalf_edge) +\
                    " and previous half edge " + str(previous_half_edge.Index()) +\
                    " are in different cells."
                return False, ihalf_edge, error_msg

            if (next_half_edge.Cell() != cell):
                error_msg = "Half edge " + str(ihalf_edge) +\
                    " and next half edge " + str(next_half_edge.Index()) +\
                    " are in different cells."
                return False, ihalf_edge, error_msg

            half_edgeX = half_edge.NextHalfEdgeAroundEdge()

            if (half_edgeX is None):
                error_msg = "Half edge " + half_edge.IndexAndEndpointsStr(",") +\
                    " missing next half edge around edge."
                return False, ihalf_edge, error_msg

            if (half_edgeX != half_edge):
                if (not half_edge.SameEndpoints(half_edgeX)):
                    error_msg = "Error. Two half edges around edge have different endpoints." +\
                        " Half edge:" + half_edge.IndexAndEndpointsStr(",") + " "\
                        + half_edgeX.IndexAndEndpointsStr(",")
                    return False, ihalf_edge, error_msg

                if (half_edgeX.Cell() == cell):
                    error_msg = "Error. Two half edges around edge are in the same cell. " +\
                        " Half edges: " + half_edge.IndexAndEndpointsStr(",") + " "\
                        + half_edgeX.IndexAndEndpointsStr(",")
                    return False, ihalf_edge, error_msg

        # Dictionary tracking visited half edges.
        # Avoid revisiting (reprocessing) visited half edges.
        is_visited = dict()
        for j in self._half_edge_dict:
            ihalf_edge = j

            if ihalf_edge in is_visited:
                continue

            half_edge = self.HalfEdge(ihalf_edge)
            numh = half_edge.CountNumHalfEdgesAroundEdge()
            vfrom = half_edge.FromVertex()
            vto = half_edge.ToVertex()
            ivfrom = vfrom.Index()
            ivto = vto.Index()
            numh2 = vto.CountNumIncidentHalfEdges(ivfrom) +\
                    vfrom.CountNumIncidentHalfEdges(ivto)

            if (numh != numh2):
                error_msg = "Inconsistency between half edges around edge " +\
                    "and vertex incident lists for edge (" +\
                    half_edge.EndpointsStr(",") + ")."
                return False, ihalf_edge, error_msg

            # Mark all visited half edges so that they are not processed again.
            # Reduces time spent checking half edge around edge.
            half_edgeX = half_edge
            for k in (0,numh):
                ihalf_edgeX = half_edge.Index()
                is_visited[ihalf_edgeX] = True
                half_edgeX = half_edgeX.NextHalfEdgeAroundEdge()

        return True, 0, None


    ## Check data structure cells.
    # - Return true if no problems found.
    # - Return index of problem cell and error message.
    def CheckCells(self):

        # Initialize
        error_msg = None
        icell = 0

        icmax = max(self.CellIndices(),default=-1)
        if (self.MaxHalfEdgeIndex() < icmax):
            error_msg = "Incorrect value (" + str(self.MaxCellIndex()) +\
                ") of _max_cell_index.  Max cell is " + str(icmax) + "."
            return False, icmax, error_msg

        for j in self._cell_dict:
            icell = j

            cell = self.Cell(icell)

            if (cell.Index() != icell):
                error_msg = "Incorrect cell index for cell " +\
                    str(icell) + "."
                return False, icell, error_msg

            half_edge0 = cell.HalfEdge()
            if (half_edge0.Cell() is not cell):
                error_msg = "Incorrect half edge stored in cell " +\
                    str(icell) + "."
                return False, icell, error_msg

            half_edge = half_edge0
            cell_numv = cell.NumVertices()

            for k in range(1,cell_numv):
                half_edge = half_edge.NextHalfEdgeInCell()

                if (half_edge == half_edge0) :
                    error_msg = "Incorrect number of vertices (" +\
                        str(cell_numv) + ") stored in cell " +\
                        str(icell) + ". Counted " + str(k) + " vertices."
                    return False, icell, error_msg


            if (half_edge.NextHalfEdgeInCell() is not half_edge0):
                error_msg = "Incorrect number of vertices (" +\
                    str(cell_numv) + ") stored in cell " +\
                    str(icell) + ".  Cell has more than " + str(cell_numv) +\
                    " vertices."
                return False, icell, error_msg

        return True, 0, None


    ## Check vertices, half edges and cells.
    # - Return False if problem detected.
    def CheckAll(self):

        flag, iv, error_msg = self.CheckVertices()
        if not flag:
            if (error_msg is None):
                error_msg = "Error related to vertex " + str(iv) + "."

            return False, error_msg

        flag, ihalf_edge, error_msg = self.CheckHalfEdges()
        if not flag:
            if (error_msg is None):
                half_edge = self.HalfEdge(ihalf_edge)
                error_msg = "Error related to half edge " +\
                    half_edge.IndexAndEndpointsStr(",") + "."

            return False, error_msg

        flag, icell, error_msg = self.CheckCells()
        if not flag:
            if (error_msg is None):
                error_msg = "Error related to cell " + str(icell) + "."

            return False, error_msg

        # Passed all checks.
        return True, None


    ## Check if mesh cells are consistently oriented.
    # - Returns true if all adjacent cells are consistently oriented.
    # - Returns index of half edge where half_edge.Cell() and
    #   half_edge.NextHalfEdgeAroundEdge.Cell() have opposite orientations.
    def CheckOrientation(self):

        for ihalf_edge in self.HalfEdgeIndices():
            half_edge = self.HalfEdge(ihalf_edge)

            half_edgeB = half_edge.NextHalfEdgeAroundEdge()

            if (half_edge is half_edgeB):
                # Skip boundary half edge.
                continue

            if (half_edge.FromVertexIndex() == half_edgeB.FromVertexIndex()):
                # Cells half_edge.Cell() and half_edgeB.Cell() have
                #   inconsistent orientations.
                return False, ihalf_edge

        return True, 0


    ## Check manifold edge property.
    # - Return true if all edges have 2 or fewer incident cells.
    # - Returns index of non-manifold edge.
    def CheckManifoldEdges(self):

        for ihalf_edge in self.HalfEdgeIndices():
            half_edge = self.HalfEdge(ihalf_edge)

            nume = half_edge.CountNumHalfEdgesAroundEdge()
            if (nume >= 3):
                return False, ihalf_edge

        return True, 0


    ## Check manifold vertex property.
    # - Return true if cells on each vertex form a fan.
    # - Returns index o non-manifold vertex.
    def CheckManifoldVertices(self):

        for iv in self.VertexIndices():
            v = self.Vertex(iv)

            numh = v.NumHalfEdgesFrom()

            if (numh == 0):
                continue

            half_edge0 = v.KthHalfEdgeFrom(0)

            half_edge = half_edge0.PrevHalfEdgeAroundVertex(iv)

            num_cells = 1
            while (half_edge is not half_edge0 and\
                    not(half_edge.IsBoundary()) and num_cells <= numh):
                num_cells += 1
                half_edge = half_edge.PrevHalfEdgeAroundVertex(iv)

            if (num_cells != numh):
                return False, iv

        return True, 0


    ## Check manifold vertex and edge properties.
    #  - Returns ( flagv, flag_halfe, iv, ihalf_edge )
    #    where flagv is true if mesh has a non-manifold half vertex,
    #    flag_halfe is true if mesh has a non-manifold half edge
    #    iv is the index of a non-manifold vertex, if one exists,
    #    ihalf_edge is the index of a non-manifold edge, if one exist.
    def CheckManifold(self):

        flag_halfe, ihalf_edge = self.CheckManifoldEdges()
        flag_halfe = not flag_halfe

        flagv, iv = self.CheckManifoldVertices()
        flagv = not flagv

        return flagv, flag_halfe, iv, ihalf_edge
