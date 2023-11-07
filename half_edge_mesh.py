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
# - Version 0.0.2

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

import half_edge_mesh_elements


## Half edge mesh data structure.
#  - Initiazed with vertex, half edge and cell classes.
#  - These classes should be derived from HMESH_VERTEX_BASE,
#      HMESH_HALF_EDGE_BASE, and HMESH_CELL_BASE.
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
    # - Returns vertex, even if vertex already exists.
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


    ## Add new half edge with unassigned index.
    def _AddNewHalfEdge(self, cell, vfrom):
        ihalf_edge = self.MaxHalfEdgeIndex() + 1
        half_edge = self._AddHalfEdge(ihalf_edge, cell, vfrom)
        return half_edge


    def _AddAndLinkHalfEdge(self, ihalf_edge, cell, vfrom, vto, hprev):

        half_edge = self._AddHalfEdge(ihalf_edge, cell, vfrom)

        half_edge.prev_half_edge_in_cell = hprev
        if (hprev is not None):
            hprev.next_half_edge_in_cell = half_edge

        half_edgeB = vto.FindHalfEdgeTo(vfrom.Index())
        if (half_edgeB is None):
            # Check whether there is a half edge around the edges
            #   with same orientation as half_edge.
            half_edgeC = vfrom.FindHalfEdgeTo(vto.Index())
            if (half_edgeC is None):
                half_edge.next_half_edge_around_edge = half_edge
            else:
                self._LinkHalfEdgesAroundEdge(half_edgeC, half_edge)
        else:
            self._LinkHalfEdgesAroundEdge(half_edgeB, half_edge)

        vfrom.half_edge_from.append(half_edge)

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


    ## Move boundary half edge to half_edge_from[0] for each vertex
    #    in array cell_vertex[].
    #  - Added: 12-02-2024 - RW
    #  @param cell_vertex[] List of cell vertex indices.
    #  @pre   Each vertex should already have been created.
    def _MoveBoundaryHalfEdgeToHalfEdgeFrom0(self, cell_vertex):

        for i in range(0,len(cell_vertex)):
            iv = cell_vertex[i]
            v = self.Vertex(iv)

            if (v is None):
                raise Exception\
                    ("Programming error." +\
                    " Attempt to access non-existent vertex with identifier " +\
                    str(iv) +\
                    "\n  in HALF_EDGE_MESH_BASE::_MoveBoundaryHalfEdgeFrom0().")

            v.MoveBoundaryHalfEdgeToHalfEdgeFrom0()


    ## Remove half edge from the half_edge_from list of its from_vertex.
    #  - Calls _MoveBoundaryHalfEdgeToHalfEdgeFrom0(), if necessary,
    #    to ensure that first half edge in vertex list is a boundary edge.
    def _RemoveHalfEdgeFromVertexList(self, half_edge0):
        v0 = half_edge0.FromVertex()
        list_length = v0.NumHalfEdgesFrom()
        ilast = list_length-1

        for k in range(0,list_length):
            half_edge = v0.KthHalfEdgeFrom(k)

            if (half_edge0 is half_edge):
                v1 = half_edge.ToVertex();
                if (k != ilast):
                    # Replace half_edge0 with last entry.
                    v0.half_edge_from[k] = v0.half_edge_from[ilast]

                v0.half_edge_from.pop()

                # Deleting a half_edge can create a boundary half edge
                #   at v0 or v1.
                v0.MoveBoundaryHalfEdgeToHalfEdgeFrom0()
                v1.MoveBoundaryHalfEdgeToHalfEdgeFrom0()

                return


    ## Unlink half edge from linked list of half edges around edge.
    def _UnlinkHalfEdgeFromHalfEdgesAroundEdge(self, half_edge):
        prev_half_edge_around_edge = half_edge.PrevHalfEdgeAroundEdge()
        prev_half_edge_around_edge.next_half_edge_around_edge =\
            half_edge.NextHalfEdgeAroundEdge()

        # Link half_edge to self.
        half_edge.next_half_edge_around_edge = half_edge


    ## Relink half edges in cell.
    #  - Overwrites previous links.
    def _RelinkHalfEdgesInCell(self, hprev, hnext):
        hprev.next_half_edge_in_cell = hnext
        hnext.prev_half_edge_in_cell = hprev


    ## Add cell with index icell.
    #  - Raises an exception if some cell already has index icell.
    #  - Use MaxCellIndex() to find an unused cell index.
    def _AddCell(self, icell):
        cell = self.CELL_TYPE()
        cell.index = icell
        self._cell_dict[icell] = cell
        self._max_cell_index = max(self._max_cell_index, icell)

        return cell

    ## Delete vertex.
    def _DeleteVertex(self, v):
        if (v is None):
            # Can't delete None.
            return

        iv = v.Index()
        self._vertex_dict.pop(iv,0)


    ## Delete half edge.
    def _DeleteHalfEdge(self, half_edge):

        self._UnlinkHalfEdgeFromHalfEdgesAroundEdge(half_edge)
        self._RemoveHalfEdgeFromVertexList(half_edge)
        half_edge.next_half_edge_in_cell = None
        half_edge.prev_half_edge_in_cell = None
        ihalf_edge = half_edge.Index()
        self._half_edge_dict.pop(ihalf_edge,0)


    ## Delete half edges in cell.
    def _DeleteCellHalfEdges(self, cell):
        numh = cell.NumVertices()
        half_edge = cell.HalfEdge()
        for k in range(0,numh):
            next_half_edge_in_cell = half_edge.NextHalfEdgeInCell()
            self._DeleteHalfEdge(half_edge)
            half_edge = next_half_edge_in_cell


    ## Delete cell
    def _DeleteCell(self, cell):
        self._DeleteCellHalfEdges(cell)
        icell = cell.Index()
        self._cell_dict.pop(icell,0)


    ## Split half edge.
    def _SplitHalfEdge(self, half_edge, vsplit):
        cell = half_edge.Cell()
        next_half_edge_in_cell = half_edge.NextHalfEdgeInCell()
        new_half_edge = self._AddNewHalfEdge(cell, vsplit)
        self._RelinkHalfEdgesInCell(half_edge, new_half_edge)
        self._RelinkHalfEdgesInCell(new_half_edge, next_half_edge_in_cell)
        vsplit.half_edge_from.append(new_half_edge)


    ## Split edge.
    #  - Return split vertex.
    def _SplitEdge(self, half_edge0):
        numh = half_edge0.CountNumHalfEdgesAroundEdge()

        # Create split vertex.
        vsplit = self.AddNewVertex()
        self.SetCoord(vsplit.Index(), half_edge0.ComputeMidpointCoord())

        # Split half_edge0.
        next_half_edge_around_edge = half_edge0.NextHalfEdgeAroundEdge()
        self._SplitHalfEdge(half_edge0, vsplit)
        new_half_edge0 = half_edge0.NextHalfEdgeInCell()
        self._UnlinkHalfEdgeFromHalfEdgesAroundEdge(half_edge0)

        half_edge = next_half_edge_around_edge
        for k in range(1,numh):
            next_half_edge_around_edge = \
                half_edge.NextHalfEdgeAroundEdge()
            self._SplitHalfEdge(half_edge, vsplit)
            new_half_edge = half_edge.NextHalfEdgeInCell()
            self._UnlinkHalfEdgeFromHalfEdgesAroundEdge(half_edge)
            if (half_edge.FromVertex() is half_edge0.FromVertex()):
                self._LinkHalfEdgesAroundEdge(half_edge0, half_edge)
                self._LinkHalfEdgesAroundEdge(new_half_edge0, new_half_edge)
            else:
                self._LinkHalfEdgesAroundEdge(half_edge0, new_half_edge)
                self._LinkHalfEdgesAroundEdge(new_half_edge0, half_edge)

            half_edge = next_half_edge_around_edge

        return vsplit


    ## Split the cell containing half_edge into quadrilaterals.
    #  @pre Cell has an even number of vertices.
    #  - Add vertex vsplit to cell "interior" and edge from vsplit
    #    to every other cell vertex starting at half_edge.FromVertex().
    #  - Returns vertex vsplit.
    def _SplitCellIntoQuads(self, half_edge0):
        numv = half_edge0.Cell().NumVertices()
        if (numv%2 == 1):
            raise Exception\
                ("Programming error. Cannot split a cell with odd number of vertices into quadrilaterals.")

        # Store cell vertices in a list.
        ivlist = []
        half_edge = half_edge0
        for i in range(0,numv):
            ivlist = ivlist + [ half_edge.FromVertexIndex() ]
            half_edge = half_edge.NextHalfEdgeInCell()

        # Store cell centroid.
        cell_centroid = half_edge0.Cell().ComputeCentroid()

        # Delete cell.
        icell = half_edge.CellIndex()
        self.DeleteCell(icell)

        # Create new vertex.
        vsplit = self.AddNewVertex()
        ivsplit = vsplit.Index()
        self.SetCoord(ivsplit, cell_centroid)

        # Create quadrilaterals.
        for i in range(0,numv):
            if (i%2 == 0):
                i1 = (i+1)%numv
                i2 = (i+2)%numv
                quad_vlist = [ ivsplit, ivlist[i], ivlist[i1], ivlist[i2] ]
                self.AddNewCell(quad_vlist)

        return vsplit


    # *** Public methods. ***

    # *** Get methods. ***

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


    ## Return list of edges.
    #  - Edges are represented by half_edges, one half_edge per edge.
    def GetEdgeList(self):
        elist = []
        for ihalf_edge in self.HalfEdgeIndices():
            half_edge = self.HalfEdge(ihalf_edge)
            if (half_edge is half_edge.MinIndexHalfEdgeAroundEdge()):
                elist.append(half_edge)

        return elist

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

    ## Return half edge (v0,v1) or (v1,v0) if it exists
    #  - Return None if no edge found.
    def FindEdge(self, v0, v1):
        half_edge = v0.FindHalfEdgeTo(v1.Index())

        if not(half_edge is None):
            return half_edge

        half_edge = v1.FindHalfEdgeTo(v0.Index())

        return half_edge


    # *** Count methods. ***

    ## Count number of isolated vertices.
    #  - Isolated vertices are not in any mesh cell.
    def CountNumIsolatedVertices(self):
        num_isolated_vertices = 0
        for iv in self.VertexIndices():
            v = self.Vertex(iv);
            if (v is None):
                # Shouldn't happen, but just in case.
                continue;

            if (v.NumHalfEdgesFrom() == 0):
                num_isolated_vertices = num_isolated_vertices+1

        return num_isolated_vertices


    ## Count number of edges.
    def CountNumEdges(self):
        num_edges = 0;
        for ihalf_edge in self.HalfEdgeIndices():
            half_edge = self.HalfEdge(ihalf_edge)
            if (half_edge is None):
                # Shouldn't happen, but just in case.
                continue

            min_index_half_edge = half_edge.MinIndexHalfEdgeAroundEdge()
            if (half_edge is min_index_half_edge):
                num_edges = num_edges+1

        return num_edges


    ## Count number of boundary edges.
    def CountNumBoundaryEdges(self):
        num_boundary_edges = 0
        for ihalf_edge in self.HalfEdgeIndices():
            half_edge = self.HalfEdge(ihalf_edge)
            if (half_edge is None):
                # Shouldn't happen, but just in case.
                continue
            if (half_edge.IsBoundary()):
                num_boundary_edges = num_boundary_edges+1

        return num_boundary_edges


    ## Count number of cells with a given number of vertices.
    def CountNumCellsOfSize(self, numv):
        num_cells = 0
        for icell in self.CellIndices():
            cell = self.Cell(icell)
            if (cell is None):
                # Shouldn't happen, but just in case.
                continue

            if (cell.NumVertices() == numv):
                num_cells = num_cells+1

        return num_cells


    ## Count number of cells with number of vertices
    #    greater than or equal to.
    def CountNumCellsOfSizeGE(self, numv):
        num_cells = 0
        for icell in self.CellIndices():
            cell = self.Cell(icell)
            if (cell is None):
                # Shouldn't happen, but just in case.
                continue

            if (cell.NumVertices() >= numv):
                num_cells = num_cells+1

        return num_cells


    ## Count number of triangles.
    def CountNumTriangles(self):
        return self.CountNumCellsOfSize(3)

    ## Count number of quadrilaterals.
    def CountNumQuads(self):
        return self.CountNumCellsOfSize(4)

    ## Count number of pentagons.
    def CountNumPentagons(self):
        return self.CountNumCellsOfSize(5)


    # *** Methods to add vertices ***

    ## Add vertex with index iv.
    #  @param iv Vertex index. Should not already be in VertexIndices().
    #  - Returns vertex.
    def AddVertex(self, iv):
        if (iv in self._vertex_dict):
            raise Exception("Illegal argument to AddVertex(). Vertex " + str(iv) + " already exists.")

        return self._CreateVertex(iv)


    ## Add new vertex with unassigned index.
    #  - Returns vertex.
    def AddNewVertex(self):
        iv = self.MaxVertexIndex()+1;
        return self._CreateVertex(iv)


    ## Add numv vertices ith indices from 0 to numv-1.
    #  @pre Half edge mesh has no vertices.
    def AddVertices(self, numv):
        if (len(self._vertex_dict) > 0):
            raise Exception("Error. Cannot call AddVertices() if mesh already has some vertices.")

        for iv in range(0,numv):
            self._CreateVertex(iv)


    # *** Methods to add or delete cells. ***

    ## Add cell with index icell.
    #  - Returns cell.
    #  @param icell Cell index. Should not already be in CellIndices().
    #  @param cell_vertex[] List of the cell vertex indices.
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

        # Link last half edge (hprev) and first half edge (half_edge0)
        self._LinkHalfEdgesInCell(hprev, half_edge0)

        # - This call must be AFTER _LinkHalfEdgesInCell(hprev, half_edge0)
        self._MoveBoundaryHalfEdgeToHalfEdgeFrom0(cell_vertex)

        if (len(cell_vertex) != cell.NumVertices()):
            raise Exception("Error in AddCell().  Incorrect number of vertices in cell.")

        return cell


    ## Add new cell with unassigned index.
    #  - Returns cell.
    #  @param cell_vertex[] List of the cell vertex indices.
    def AddNewCell(self, cell_vertex):
        icell = self.MaxCellIndex()+1;
        cell = self.AddCell(icell, cell_vertex)

        return cell


    ## Delete cell with index icell.
    #  @param icell Cell index. Should be in CellIndices().
    #  @param cell_vertex[] List of the cell vertex indices.
    def DeleteCell(self, icell):
        cell = self.Cell(icell)
        self._DeleteCell(cell)


    # *** Set coordinates. ***

    ## Set coordinates of iv to coord[].
    def SetCoord(self, iv, coord):
        v = self._vertex_dict[iv]

        if (v == None):
            raise Exception("Error in SetCoord(). Vertex " + str(iv) + " does not exist.")

        if (len(coord) < v.Dimension()):
            raise Exception("Error in SetCoord(). Argument coord[] has too few coordinates.")

        for ic in range(0,v.Dimension()):
            v.coord[ic] = coord[ic]


    # *** Split edge/cell. ***

    ## Split edge at midpoint.
    #  - Returns new vertex.
    #  - Splits all half edges around edge containing ihalf_edgeA.
    def SplitEdge(self, ihalf_edgeA):
        half_edgeA = self.HalfEdge(ihalf_edgeA)
        if (half_edgeA is None):
            raise Exception\
                ("Programming error. Argument to SplitEdge is not a half edge index.")

        vsplit = self._SplitEdge(half_edgeA)

        return vsplit


    ## Split the cell containing ihalf_edgeA into quadrilaterals.
    #  @pre Cell has an even number of vertices.
    #  - Add vsplit in cell "interior" and edge from vsplit
    #    to every other cell vertex starting at half_edge.FromVertex().
    def SplitCellIntoQuads(self, ihalf_edgeA):
        half_edgeA = self.HalfEdge(ihalf_edgeA)
        if (ihalf_edgeA is None):
            raise Exception\
                ("Programming error. Argument to SplitCellIntoQuads is not a half edge index.")

        if (half_edgeA.Cell().NumVertices()%2 == 1):
            raise Exception\
                ("Programming error. Cell passed to SplitEdgeIntoQuads() contains an odd number of vertices.")

        vsplit = self._SplitCellIntoQuads(half_edgeA)

        return vsplit


    # *** Check routines. ***

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
                        + ". List incorrectly includes half edge " + str(half_edge.Index()) + "."
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
    #    where flagv is true if each vertex passes manifold check, and
    #    flag_halfe is true if each edge passes manifold check.
    #    iv is the index of a non-manifold vertex, if one exists,
    #    ihalf_edge is the index of a non-manifold edge, if one exist.
    def CheckManifold(self):
        flagv, iv = self.CheckManifoldVertices()
        flag_halfe, ihalf_edge = self.CheckManifoldEdges()
        return flagv, flag_halfe, iv, ihalf_edge


    ## Check vertex index.
    #  - Return true if iv is an index of some mesh vertex.
    def CheckVertexIndex(self, iv):

        # Initialize.
        error_msg = None

        if (not bool(self._vertex_dict)):
            error_msg = "Illegal negative vertex index: " + str(iv) + "\n"
            error_msg = error_msg + "  Vertex indices cannot be negative.\n"
            return False, error_msg

        if (self.Vertex(iv) is None):
            error_msg = "No vertex has index: " + str(iv) + "\n"
            return False, error_msg

        return True, None


    # *** Print error message routines. ***

    ## Print mesh data structure error and error message.
    def PrintErrorMessage(self, out, error_msg):
        out.write("Error detected in mesh data structure.\n")
        if not(error_msg is None):
            out.write(error_msg + "\n")


    ## Print non-manifold vertex message.
    def PrintNonManifoldVertex(self, out, error_prefix, iv):
        out.write(error_prefix + " Non-manifold vertex " + str(iv) + ".\n")

    ## Print non-manifold edge message.
    def PrintNonManifoldEdge(self, out, error_prefix, ihalf_edge):
        out.write(error_prefix + " Non-manifold edge (" +\
                    self.HalfEdge(ihalf_edge).EndpointsStr(",") + ").\n")

    ## Print not oriented.
    def PrintNotOriented(self, out, error_prefix, ihalf_edge):
        out.write(error_prefix + " Inconsistent orientation of cells incident on edge (" +\
                    self.HalfEdge(ihalf_edge).EndpointsStr(",") + ").\n")
