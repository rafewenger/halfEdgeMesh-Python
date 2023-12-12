## \file half_edge_mesh_edit.py
# Extension of half edge mesh supporting mesh edit operations.
# Supports SplitCell, JoinTwoCells, and CollapseEdge.

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
# - Version 0.4.0

import half_edge_mesh


## Half edge mesh editor class.
#  - Inherits from HALF_EDGE_MESH_BASE.
#  - Initialized with vertex, half edge and cell classes.
#  - These classes should be derived from HMESH_VERTEX_BASE,
#      HMESH_HALF_EDGE_BASE and HMESH_CELL_BASE.
class HALF_EDGE_MESH_EDIT_BASE(half_edge_mesh.HALF_EDGE_MESH_BASE):

    # Private/protected member functions.

    def __init__(self, classV, classHE, classC):
        super(HALF_EDGE_MESH_EDIT_BASE,self).__init__(classV, classHE, classC)

        ## Create a dictionary storing boolean flags indicating
        #   visited vertices.
        self._vertex_visited_dict = dict()

    def _RemoveHalfEdgeFromCell(self, half_edge):
        cell = half_edge.Cell()
        prev_half_edge_in_cell = half_edge.PrevHalfEdgeInCell()
        next_half_edge_in_cell = half_edge.NextHalfEdgeInCell()
        prev_half_edge_in_cell.next_half_edge_in_cell =\
            next_half_edge_in_cell
        next_half_edge_in_cell.prev_half_edge_in_cell =\
            prev_half_edge_in_cell
        cell.num_vertices = cell.num_vertices - 1
        if (cell.HalfEdge() == half_edge):
            cell.half_edge = next_half_edge_in_cell

        # Unset prev_half_edge_in_cell and next_half_edge_in_cell.
        half_edge.prev_half_edge_in_cell = None
        half_edge.next_half_edge_in_cell = None


    ## Move half edges in vA.half_edge_from[] to vB.half_edge_from[].
    #  - Clear vA.half_edge_from[]
    def _MoveVertexHalfEdgeFromList(self, vA, vB):
        for k in range(0,vA.NumHalfEdgesFrom()):
            half_edge = vA.KthHalfEdgeFrom(k)
            half_edge.from_vertex = vB

        # Add vA.half_edge_from[] to vB.half_edge_from[].
        vB.half_edge_from.extend(vA.half_edge_from)

        vA.half_edge_from.clear()


    ## Merge two linked lists of half edges around an edge to form
    #   a single linked list.
    def _MergeHalfEdgesAroundEdge(self, half_edgeA, half_edgeB):
        next_half_edge_around_edgeA = half_edgeA.NextHalfEdgeAroundEdge()
        next_half_edge_around_edgeB = half_edgeB.NextHalfEdgeAroundEdge()
        half_edgeA.next_half_edge_around_edge =\
            next_half_edge_around_edgeB
        half_edgeB.next_half_edge_around_edge = \
            next_half_edge_around_edgeA


    ## Merge linked lists of half edges incident on (vA,vC) and (vB,vC)
    #    for some vertex vC.
    def _MergeHalfEdgesIncidentOnVertices(self, vA, vB):

        # Get list of vertices that are neighbors of both vA and vB.
        common_vneighbors_list = \
            self.GetListOfCommonVertexNeighbors(vA, vB)

        for vC in common_vneighbors_list:
            half_edgeA = self.FindEdge(vA, vC)
            half_edgeB = self.FindEdge(vB, vC)
            if (half_edgeA == None) or (half_edgeB == None):
                raise Exception\
                    ("Programming error. Unable to find half edge with given endpoints.")

            self._MergeHalfEdgesAroundEdge(half_edgeA, half_edgeB)


    # Public member functions.

    # *** vertex_visited flag functions ***

    ## Return vertex_visited flag for vertex iv.
    #  - Returns False if iv is not in vertex_visited dictionary.
    def IsVertexVisited(self, iv):
        return self._vertex_visited_dict.get(iv,False)

    ## Set vertex_visited to flag for vertex iv.
    def SetVertexVisitedFlag(self, iv, flag):
        self._vertex_visited_dict[iv] = flag

    ## Set vertex_visited to false for vertex iv to false.
    def ClearVertexVisitedFlag(self, iv):
        self.SetVertexVisitedFlag(iv, False)


    ## Set vertex_visited to flag in all neighbors of vertex iv.
    def SetVisitedFlagsInAdjacentVertices(self, iv, flag):
        v = self.Vertex(iv)
        for k in range(0,v.NumHalfEdgesFrom()):
            half_edgeA = v.KthHalfEdgeFrom(k)
            ivA = half_edgeA.ToVertexIndex()
            self.SetVertexVisitedFlag(ivA, flag)

            half_edgeB = half_edgeA.PrevHalfEdgeInCell()
            ivB = half_edgeB.FromVertexIndex()
            if (ivB != iv):
                # Set visited flags for half_edgeB.FromVertex() in case
                #   of boundary edges or cells with arbitrary orientations.
                self.SetVertexVisitedFlag(ivB, flag)


    ## Set vertex_visited to False in all neighbors of vertex iv.
    def ClearVisitedFlagsInAdjacentVertices(self, iv):
        self.SetVisitedFlagsInAdjacentVertices(iv, False)


    ## Set _visited_flag to flag in all vertices of cell icell.
    def SetVisitedFlagsInAllCellVertices(self, icell, flag):
        cell = self.Cell(icell)
        half_edge = cell.HalfEdge()
        for i in range(cell.NumVertices()):
            self.SetVertexVisitedFlag(half_edge.FromVertexIndex(), flag)
            half_edge = half_edge.NextHalfEdgeInCell()

    ## Set _visited_flag to False in all vertices of cell icell.
    def ClearVisitedFlagsInAllCellVertices(self, icell):
        self.SetVisitedFlagsInAllCellVertices(icell, False)


    ## Set _visited_flag to flag in all vertices of all cells
    #    incident on iv.
    def SetVisitedFlagsInAllVerticesOfAllIncidentCells(self, iv, flag):
        v = self.Vertex(iv)
        for k in range(v.NumHalfEdgesFrom()):
            icell = v.KthHalfEdgeFrom(k).CellIndex()
            self.SetVisitedFlagsInAllCellVertices(icell, flag)


    ## Set _visited_flag to false in all vertices of all cells
    #    incident on ivA.
    def ClearVisitedFlagsInAllVerticesOfAllIncidentCells(self, iv):
        self.SetVisitedFlagsInAllVerticesOfAllIncidentCells(iv, False)


    # *** Get vertex list functions ***

    ## Get list of vertices that are common neighbors of vA and vB.
    #  - Runs in time proportional to number of edges incident on vA
    #    plus the number of edges incident on vB.
    def GetListOfCommonVertexNeighbors(self, vA, vB):
        ivA = vA.Index()
        ivB = vB.Index()
        common_vneighbors_list = []
        self.ClearVisitedFlagsInAdjacentVertices(ivB)
        self.SetVisitedFlagsInAdjacentVertices(ivA, True)

        for k in range(0, vB.NumHalfEdgesFrom()):
            half_edge = vB.KthHalfEdgeFrom(k)
            vC = half_edge.ToVertex()
            ivC = vC.Index()
            if (self.IsVertexVisited(ivC)):
                # vC is a neighbor of vA and vB.
                common_vneighbors_list.append(vC)

                # Clear flag so vC is only appended once
                self.ClearVertexVisitedFlag(ivC)

            prev_half_edge_in_cell = half_edge.PrevHalfEdgeInCell()
            vC = prev_half_edge_in_cell.FromVertex()
            ivC = vC.Index()
            if (self.IsVertexVisited(ivC)):
                # vC is a neighbor of vA and vB
                common_vneighbors_list.append(vC)

                # Clear flag so vC is only appended once
                self.ClearVertexVisitedFlag(ivC)

        # DEBUG:
        # print("common neighbors of " + str(ivA) + " and " + str(ivB) + ": ", end="")
        # for v in common_vneighbors_list:
        #     print(str(v.Index()) + " ", end="")
        # print()

        return common_vneighbors_list


    ## Get list of consecutive cell vertices, starting at half_edgeA.FromVertex()
    #    until half_edgeB.FromVertex().
    #  - Vertices are represented by their indices.
    #  @pre half_edgeA and half_edgeB are in the same cell.
    def GetListOfConsecutiveCellVertices(self, half_edgeA, half_edgeB):
        if (half_edgeA.Cell() != half_edgeB.Cell()):
            raise Exception\
                ("Programming error detected in GetListOfConsecutiveVertices.\n" +\
                    "  Half edges are in different cells.")

        cell = half_edgeA.Cell()
        listA = []
        half_edge = half_edgeA
        # Check that kount is less than cell.NumVertices(), just in case
        #   data structure is corrupted.
        kount = 0
        while (half_edge != half_edgeB) and (kount <= cell.NumVertices()):
            listA.append(half_edge.FromVertexIndex())
            half_edge = half_edge.NextHalfEdgeInCell()
            kount = kount + 1

        if (kount > cell.NumVertices()):
            # Something went wrong.
            icell = cell.Index()
            raise Exception\
                ("Programming error detected in GetListOfConsecutiveVertices.\n" +\
                    "  Error getting vertices in cell " + str(icell) + ".")

        listA.append(half_edgeB.FromVertexIndex())

        return listA


    ## Get list of all cell vertices, starting at half_edge.FromVertex().
    #  - Vertices are represented by their indices.
    #  - Vertices are in order around cell, starting at half_edge.FromVertex().
    def GetListOfAllCellVertices(self, half_edge):
        cell_vlist = self.GetListOfConsecutiveCellVertices\
                        (half_edge, half_edge.PrevHalfEdgeInCell())
        return cell_vlist


    # *** Collapse/join/split functions ***

    ## Collapse edge, mapping two vertices to a single vertex.
    #  - Returns vertex formed by merging two vertices.
    #  - Does not delete any vertices.
    def CollapseEdge(self, ihalf_edge0):
        half_edge0 = self.HalfEdge(ihalf_edge0)
        if (half_edge0 is None):
            # Can't collapse a half edge that doesn't exist.
            return None;

        # Don't collapse half_edge0 if its two endoints (vA,vB) are
        #   in some cell, but edge (vA,vB) is not in the cell.
        if (self.IsIllegalEdgeCollapseH(ihalf_edge0)):
            return None

        vA = half_edge0.FromVertex()
        vB = half_edge0.ToVertex()
        midpoint_coord = half_edge0.ComputeMidpointCoord()

        max_num_half_edges_around_edge =\
            half_edge0.FromVertex().NumHalfEdgesFrom() +\
            half_edge0.ToVertex().NumHalfEdgesFrom()

        half_edges_incident_on_triangles = []
        half_edges_not_incident_on_triangles = []

        # Get list of half edges around edge incident on triangles
        #  and another list of half edges around edge not incident on triangles
        half_edges_incident_on_triangles = []
        half_edges_not_incident_on_triangles = []
        half_edge = half_edge0
        for i in range(max_num_half_edges_around_edge):
            cell = half_edge.Cell()
            if (cell.NumVertices() <= 3):
                half_edges_incident_on_triangles.append(half_edge)
            else:
                half_edges_not_incident_on_triangles.append(half_edge)

            half_edge = half_edge.NextHalfEdgeAroundEdge()
            if (half_edge == half_edge0):
                break

        # Delete triangles incident on edge.
        for half_edge in half_edges_incident_on_triangles:
            cell = half_edge.Cell()
            self.DeleteCell(cell.Index())

        # Merge linked lists of half edges incident on (vA,vC) and (vB,vC)
        #    for some vertex vC.
        self._MergeHalfEdgesIncidentOnVertices(vA, vB)

        for half_edge in half_edges_not_incident_on_triangles:
            self._RemoveHalfEdgeFromVertexList(half_edge)
            self._RemoveHalfEdgeFromCell(half_edge)
            self._half_edge_dict.pop(half_edge.Index(),0)

        # Move all half edges from vA to be from vB.
        # - Moves vA.half_edge_from to vB.half_edge_from
        self._MoveVertexHalfEdgeFromList(vA, vB)

        # Ensure that vB.half_edge_from[0] is boundary_edge.
        vB.MoveBoundaryHalfEdgeToHalfEdgeFrom0()

        # Set vB to midpoint of (vA,vB).
        vB.coord = midpoint_coord

        return vB


    ## Split cell into two cells.
    #  - Split cell at from vertices of ihalf_edgeA and ihalf_edgeB.
    #  - Return half edge representing split edge.
    #  @pre ihalf_edgeA and ihalf_edgeB are in the same cell.
    def SplitCell(self, ihalf_edgeA, ihalf_edgeB):
        half_edgeA = self.HalfEdge(ihalf_edgeA)
        half_edgeB = self.HalfEdge(ihalf_edgeB)
        if (half_edgeA is None) or (half_edgeB is None):
            # Can't split if some half_edge is none.
            return None;

        # Don't split if edges share vertices.
        if (self.IsIllegalSplitCell(ihalf_edgeA, ihalf_edgeB)):
            return None

        listA = self.GetListOfConsecutiveCellVertices(half_edgeA, half_edgeB)
        listB = self.GetListOfConsecutiveCellVertices(half_edgeB, half_edgeA)
        vA = half_edgeA.FromVertex()
        vB = half_edgeB.FromVertex()

        # Delete cell.
        cell = half_edgeA.Cell()
        self.DeleteCell(cell.Index())

        self.AddNewCell(listA)
        self.AddNewCell(listB)

        split_half_edge = self.FindEdge(vA, vB)

        if (split_half_edge == None):
            raise Exception\
                ("Programming error detected in SplitCell.\n" +\
                    "  Unable to find half edge representing split edge.")

        return split_half_edge


    ## Joint two cells sharing an edge.
    #  - Returns new joined cell.
    def JoinTwoCells(self, ihalf_edgeA):
        half_edgeA = self.HalfEdge(ihalf_edgeA)
        if (half_edgeA is None):
            raise Exception("Programming error. Argument to JoinTwoCells is not a half edge index.")

        if (self.IsIllegalJoinTwoCells(ihalf_edgeA)):
            return None

        half_edgeB = half_edgeA.NextHalfEdgeAroundEdge()
        if (half_edgeB.NextHalfEdgeAroundEdge() != half_edgeA):
            raise Exception("Programming error. Half edge passed to JoinTwoCells is in an edge shared by three or more cells.")

        listA = self.GetListOfAllCellVertices(half_edgeA.NextHalfEdgeInCell())
        listB = self.GetListOfAllCellVertices(half_edgeB.NextHalfEdgeInCell())

        # Remove last element of each list, since those elements
        #   are in the other list.
        listA.pop()
        listB.pop()
        listC = listA + listB

        icellA = half_edgeA.CellIndex()
        icellB = half_edgeB.CellIndex()
        self.DeleteCell(icellA)
        self.DeleteCell(icellB)

        cellC = self.AddNewCell(listC)
        return cellC


    # *** Functions to triangulate cells. ***

    ## Triangulate cell from vertex half_edge0.FromVertexIndex().
    #  - Split cell into triangles, each triangle incident
    #    on half_edge0.FromVertexIndex().
    def TriangulateCellFromVertex(self, ihalf_edge0):
        half_edge0 = self.HalfEdge(ihalf_edge0)
        if (half_edge0 == None):
            raise Exception\
                ("Programming error. Argument to TriangulateCellFromVertex() is not a half edge index.")

        cell0 = half_edge0.Cell()
        if (cell0.IsTriangle()):
            # cell0 is already a triangle.
            return

        cell0_vlist = self.GetListOfAllCellVertices(half_edge0)

        self.DeleteCell(cell0.Index())

        numv = len(cell0_vlist)
        for i in range(1,numv-1):
            triangle_vlist =\
                [ cell0_vlist[0], cell0_vlist[i], cell0_vlist[i+1] ]
            self.AddNewCell(triangle_vlist)


    # *** Functions to check potential collapse/join/split operations. ***

    ## Return True if edge collapse is illegal.
    #  - Edge collapse (vA,vB) is illegal if some cell contains
    #    both vA and vB but not edge (vA,vB).
    #  - Version that takes two vertices.
    #  - NOTE: Function suffix is 'V' (for vertex index arguments).
    #  - Takes time proportional to the sum of the number of vertices
    #    in all cells incident on vA or vB.
    def IsIllegalEdgeCollapseV(self, ivA, ivB):
        vA = self.Vertex(ivA)
        vB = self.Vertex(ivB)
        if (vA.NumHalfEdgesFrom() > vB.NumHalfEdgesFrom()):
            # Swap vA and vB to reduce number of cells processed.
            return self.IsIllegalEdgeCollapseV(ivB, ivA)
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
                    half_edge = half_edge.NextHalfEdgeInCell()

            return False


    ## Return True if edge collapse is illegal.
    # - NOTE: Function suffix is 'H' (for half_edge index argument).
    def IsIllegalEdgeCollapseH(self, ihalf_edge):
        half_edge = self.HalfEdge(ihalf_edge)
        return self.IsIllegalEdgeCollapseV\
                    (half_edge.FromVertexIndex(), half_edge.ToVertexIndex())


    ## Return True if edge collapse changes mesh topology.
    #  - NOTE: Function suffix is 'H' (for half_edge index argument)
    def DoesEdgeCollapseChangeMeshTopologyH(self, ihalf_edge):
        flag_hole, iv = self.FindTriangleHole(ihalf_edge)
        if (flag_hole):
            return True

        if (self.IsInteriorEdgeWithBoundaryVertices(ihalf_edge)):
            return True

        half_edge = self.HalfEdge(ihalf_edge)
        icell = half_edge.CellIndex()
        if (self.IsIsolatedTriangle(icell)):
            return True

        if (self.IsInTetrahedron(icell)):
            return True

        return False


    ## Return True if split cell is illegal.
    # - Split cell is illegal
    #     if half_edgeA and half_edgeB are in different cells or
    #     if half_edgeA.FromVertex() and half_edgeB.FromVertex()
    #       are adjacent vertices.
    def IsIllegalSplitCell(self, ihalf_edgeA, ihalf_edgeB):
        half_edgeA = self.HalfEdge(ihalf_edgeA)
        half_edgeB = self.HalfEdge(ihalf_edgeB)
        if not(half_edgeA.Cell() is half_edgeB.Cell()):
            return True

        if (half_edgeA is half_edgeB):
            return True

        if (half_edgeA.FromVertex() is half_edgeB.ToVertex()):
            return True

        if (half_edgeA.ToVertex() is half_edgeB.FromVertex()):
            return True

        vA = half_edgeA.FromVertex()
        vB = half_edgeB.FromVertex()
        half_edge = self.FindEdge(vA,vB)
        if not(half_edge is None):
            # Some cell in mesh intersects half_edgeA.Cell() in (vA,vB)
            #   but (vA,vB) is not an edge of half_edgeA.Cell().
            # Split would create edge (vA,vB) in 3 or more cells.
            return True

        return False


    ## Return True if triangulate cell from vertex changes mesh topology
    # - Triangulate cell from vertex changes mesh topology
    #     if some triangulation diagonal is already a mesh edge.
    def DoesTriangulateCellFromVertexChangeTopology(self, ihalf_edgeA):
        NUM_VERT_PER_TRIANGLE = 3;
        half_edgeA = self.HalfEdge(ihalf_edgeA)
        cellA = half_edgeA.Cell()
        vA = half_edgeA.FromVertex()

        if (cellA.NumVertices() <= NUM_VERT_PER_TRIANGLE):
            # Triangulate cell from vertex does nothing.
            return False

        half_edgeB = half_edgeA.NextHalfEdgeInCell()
        for i in range(2,cellA.NumVertices()-1):
            half_edgeB = half_edgeB.NextHalfEdgeInCell()
            vB = half_edgeB.FromVertex()
            half_edge = self.FindEdge(vA,vB)
            if not(half_edge is None):
                # Some cell in mesh intersects half_edgeA.Cell() in (vA,vB)
                #   but (vA,vB) is not an edge of half_edgeA.Cell().
                # Triangulation would create edge (vA,vB) in 3 or more cells.
                return True

        return False


    ## Return True if join two cells is illegal.
    #  - Join cells is illegal if half_edge is a boundary half edge
    #    or more than two cells are incident on the edge
    #    or some endpoint of half edge has degree 2.
    def IsIllegalJoinTwoCells(self, ihalf_edge):
        TWO = 2
        half_edge = self.HalfEdge(ihalf_edge)
        if (half_edge.IsBoundary()):
            return True

        ivfrom = half_edge.FromVertexIndex()
        ivto = half_edge.ToVertexIndex()
        if not(self.IsVertexIncidentOnMoreThanTwoEdges(ivfrom)):
            return True

        if not(self.IsVertexIncidentOnMoreThanTwoEdges(ivto)):
            return True

        half_edgeX = half_edge.NextHalfEdgeAroundEdge()
        if not(half_edge is half_edgeX.NextHalfEdgeAroundEdge()):
            # More than two cells are incident on edge
            #  (half_edge.FromVertex(), half_edge.ToVertex()).
            return True

        num_shared_vertices = self.CountNumVerticesSharedByTwoCells\
            (half_edge.CellIndex(), half_edgeX.CellIndex())
        if (num_shared_vertices > TWO):
                # Cells share more than two vertices.
                return True

        # Join is LEGAL
        return False


    ## Return True if half edge endpoints and v are in a mesh triangle.
    def IsInTriangle(self, ihalf_edge0, iv):
        half_edge0 = self.HalfEdge(ihalf_edge0)
        v = self.Vertex(iv)
        # Cannot have more than max_numh half edges around an edge.
        max_numh = half_edge0.FromVertex().NumHalfEdgesFrom() +\
                    half_edge0.ToVertex().NumHalfEdgesFrom()

        half_edge = half_edge0
        for k in range(max_numh):
            if (half_edge.Cell().IsTriangle()):
                prev_half_edge = half_edge.PrevHalfEdgeInCell()

                if (prev_half_edge.FromVertex() is v):
                    return True

            half_edge = half_edge.NextHalfEdgeAroundEdge()
            if (half_edge is half_edge0):
                break

        return False


    ## Return True if both endpoints (vfrom,vto) of half_edge
    #    are neighbors of some vertex vC, but (vfrom, vto, vC)
    #    is not a mesh triangle.
    #  - Returns also ivC, the index of the third vertex vC.
    def FindTriangleHole(self, ihalf_edge0):
        half_edge0 = self.HalfEdge(ihalf_edge0)
        vA = half_edge0.FromVertex()
        vB = half_edge0.ToVertex()

        # Get list of vertices that are neighbors of both vA and vB.
        common_vneighbors_list = \
            self.GetListOfCommonVertexNeighbors(vA, vB)

        # Cannot have more than max_numh half edges around an edge.
        max_numh = half_edge0.FromVertex().NumHalfEdgesFrom() +\
            half_edge0.ToVertex().NumHalfEdgesFrom()

        for vC in common_vneighbors_list:
            ivC = vC.Index()
            half_edge = half_edge0
            for k in range(max_numh):
                if not(self.IsInTriangle(half_edge.Index(), ivC)):
                    return True, ivC

                half_edge = half_edge.NextHalfEdgeAroundEdge()
                if (half_edge == half_edge0):
                    break

        return False, 0


    ## Return true if half edge is in interior but both vertices are on boundary.
    def IsInteriorEdgeWithBoundaryVertices(self, ihalf_edge):
        half_edge = self.HalfEdge(ihalf_edge)
        if half_edge.IsBoundary():
            # Half edge is on boundary.
            return False

        if (half_edge.FromVertex().IsBoundary() and\
            half_edge.ToVertex().IsBoundary()):
            # Half edge is in mesh interior, but connects boundary vertices.
            return True

        return False


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
            half_edge = half_edge.NextHalfEdgeInCell()

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
        iv2 = half_edge0.PrevHalfEdgeInCell().FromVertexIndex()

        # Cannot have more than max_numh half edges around an edge.
        max_numh = half_edge0.FromVertex().NumHalfEdgesFrom() +\
                    half_edge0.ToVertex().NumHalfEdgesFrom()

        half_edge = half_edge0
        for k in range(max_numh):
            half_edge = half_edge.NextHalfEdgeAroundEdge()
            if (half_edge == half_edge0):
                break

            cell = half_edge.Cell()

            # Check if cell and vertex iv2 form a tetrahedron.
            if cell.IsTriangle():
                prev_half_edge = half_edge.PrevHalfEdgeInCell()
                next_half_edge = half_edge.NextHalfEdgeInCell()
                iprev = prev_half_edge.Index()
                inext = next_half_edge.Index()

                if self.IsInTriangle(iprev, iv2) and\
                    self.IsInTriangle(inext, iv2):
                    # cell0, cell, and two triangles form a tetrahedron.
                    return True

        return False


    ## Return true if vertex is incident on more than two edges.
    def IsVertexIncidentOnMoreThanTwoEdges(self, iv):
        v = self.Vertex(iv)
        TWO = 2
        if (v.NumHalfEdgesFrom() > TWO):
            return True

        if not(v.IsBoundary()):
            return False

        if (v.NumHalfEdgesFrom() == TWO):
            # Boundary vertex in two cells must have at least
            #   three incident edges.
            return True
        else:
            # Boundary vertex is in just one cell, and has exactly
            #   two incident edges.
            return False


    ## Count number of vertices shared by two cells.
    def CountNumVerticesSharedByTwoCells(self, icellA, icellB):
        num_shared_vertices = 0
        self.ClearVisitedFlagsInAllCellVertices(icellB)
        self.SetVisitedFlagsInAllCellVertices(icellA, True)

        cellB = self.Cell(icellB)
        half_edgeB = cellB.HalfEdge()
        for k in range(0, cellB.NumVertices()):
            iv = half_edgeB.FromVertexIndex()
            if (self.IsVertexVisited(iv)):
                num_shared_vertices = num_shared_vertices+1

            half_edgeB = half_edgeB.NextHalfEdgeInCell()

        return num_shared_vertices
