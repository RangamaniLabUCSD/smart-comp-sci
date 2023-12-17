"""
Functions to create meshes for mechanotransduction example
"""

from typing import Tuple
import pathlib
import numpy as np
import sympy as sym
import dolfin as d
from mpi4py import MPI
from sympy.parsing.sympy_parser import parse_expr
from sympy.solvers.solveset import solveset_real
from smart.mesh_tools import implicit_curve, compute_curvature, gmsh_to_dolfin

def create_3dcell(
    outerExpr: str = "",
    innerExpr: str = "",
    hEdge: float = 0,
    hInnerEdge: float = 0,
    interface_marker: int = 12,
    outer_marker: int = 10,
    inner_vol_tag: int = 2,
    outer_vol_tag: int = 1,
    comm: MPI.Comm = d.MPI.comm_world,
    verbose: bool = False,
    return_curvature: bool = False,
    nanopillars: Tuple[float, float, float] = "",
    thetaExpr: str = "",
) -> Tuple[d.Mesh, d.MeshFunction, d.MeshFunction]:
    """
    Creates an axisymmetric mesh, with the bounding curve defined in
    terms of r and z. (e.g. unit circle defined by "r**2 + (z-1)**2 == 1")
    It is assumed that substrate is present at z = 0, so if the curve extends
    below z = 0 , there is a sharp cutoff.
    Can include one compartment inside another compartment

    Args:
        outerExpr: String implicitly defining an r-z curve for the outer surface
        innerExpr: String implicitly defining an r-z curve for the inner surface
        hEdge: maximum mesh size at the outer edge
        hInnerEdge: maximum mesh size at the edge
            of the inner compartment
        interface_marker: The value to mark facets on the interface with
        outer_marker: The value to mark facets on the outer ellipsoid with
        inner_vol_tag: The value to mark the inner ellipsoidal volume with
        outer_vol_tag: The value to mark the outer ellipsoidal volume with
        comm: MPI communicator to create the mesh with
        verbose: If true print gmsh output, else skip
        nanopillars: tuple with nanopillar radius, height, spacing
    Returns:
        Tuple (mesh, facet_marker, cell_marker)
    """
    import gmsh

    if outerExpr == "":
        ValueError("Outer surface is not defined")

    rValsOuter, zValsOuter = implicit_curve(outerExpr)

    if not innerExpr == "":
        rValsInner, zValsInner = implicit_curve(innerExpr)
        zMid = np.mean(zValsInner)
        ROuterVec = np.sqrt(rValsOuter**2 + (zValsOuter - zMid) ** 2)
        RInnerVec = np.sqrt(rValsInner**2 + (zValsInner - zMid) ** 2)
        maxOuterDim = max(ROuterVec)
        maxInnerDim = max(RInnerVec)
    else:
        zMid = np.mean(zValsOuter)
        ROuterVec = np.sqrt(rValsOuter**2 + (zValsOuter - zMid) ** 2)
        maxOuterDim = max(ROuterVec)
    if np.isclose(hEdge, 0):
        hEdge = 0.1 * maxOuterDim
    if np.isclose(hInnerEdge, 0):
        hInnerEdge = 0.2 * maxOuterDim if innerExpr == "" else 0.2 * maxInnerDim
    # Create the two axisymmetric body mesh using gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", int(verbose))
    gmsh.model.add("3dcell")
    # first add outer body and revolve
    if thetaExpr != "":
        num_theta = 81
        if "rect" in thetaExpr:
            try:
                AR = float(thetaExpr[4:])
                a_rect = np.sqrt(np.pi * AR)
                b_rect = np.sqrt(np.pi / AR)
            except ValueError:
                raise ValueError("For rectangular pattern, " 
                    "thetaExpr must be rectAR, where AR is a number specifying " 
                    "the aspect ratio of the rectangle")
            theta_crit = np.arctan2(b_rect, a_rect)
            theta_incr = 2*np.pi / (num_theta-1)
            num_per_eighth = int(np.floor((num_theta-1) / 8))
            theta_range1 = np.linspace(0.0, theta_crit-theta_incr/4, num_per_eighth)
            theta_range2 = np.linspace(theta_crit+theta_incr/4, 
                                       np.pi-theta_crit-theta_incr/4, 2*num_per_eighth)
            theta_range3 = np.linspace(np.pi-theta_crit+theta_incr/4, 
                                       np.pi+theta_crit-theta_incr/4, 2*num_per_eighth)
            theta_range4 = np.linspace(np.pi+theta_crit+theta_incr/4, 
                                       2*np.pi-theta_crit-theta_incr/4, 2*num_per_eighth)
            theta_range5 = np.linspace(2*np.pi-theta_crit+theta_incr/4, 
                                       2*np.pi, num_per_eighth)
            thetaVec = np.concatenate([theta_range1, theta_range2, theta_range3, theta_range4, theta_range5])
            thetaVec[-1] = 0.0 # for exactness, replace 2*pi with 0.0
        else:
            thetaVec = np.linspace(0.0, 2*np.pi, num_theta)
            thetaVec[-1] = 0.0 # for exactness, replace 2*pi with 0.0
            thetaExprSym = parse_expr(thetaExpr)
    else:
        thetaVec = [0]
    cell_plane_tag = []
    bottom_point_list = []
    outer_spline_list = []
    edge_surf_list = []
    edge_segments = []
    outer_shape = []
    top_point = gmsh.model.occ.add_point(0.0, 0.0, zValsOuter[0])
    all_points_list = [top_point]
    if thetaExpr != "":
        for j in range(len(thetaVec)):
            if j == (len(thetaVec)-1):
                outer_spline_list.append(outer_spline_list[0])
                bottom_point_list.append(bottom_point_list[0])
            else:
                if "rect" in thetaExpr:
                    if np.abs(np.tan(thetaVec[j])) <= b_rect / a_rect:
                        curScale = np.abs(a_rect / (2*np.cos(thetaVec[j])))
                    else:
                        curScale = np.abs(b_rect / (2*np.sin(thetaVec[j])))
                else:
                    curScale = thetaExprSym.subs({"theta": thetaVec[j]})
                    curScale = float(curScale)
                outer_tag_list = []
                for i in range(len(rValsOuter)):
                    if i == 0:
                        outer_tag_list.append(top_point)
                    else:
                        cur_tag = gmsh.model.occ.add_point(
                            rValsOuter[i]*curScale*np.cos(thetaVec[j]), 
                            rValsOuter[i]*curScale*np.sin(thetaVec[j]), 
                            zValsOuter[i])
                        outer_tag_list.append(cur_tag)
                        all_points_list.append(cur_tag)
                outer_spline_list.append(gmsh.model.occ.add_spline(outer_tag_list))
                bottom_point_list.append(outer_tag_list[-1])
            if np.isclose(zValsOuter[-1], 0):  # then include substrate at z=0
                # origin_tag = gmsh.model.occ.add_point(0, 0, 0)
                # symm_axis_tag = gmsh.model.occ.add_line(origin_tag, outer_tag_list[0])
                # if j > 0:
                #     bottom_tag_prev = bottom_tag
                # bottom_tag = gmsh.model.occ.add_line(origin_tag, outer_tag_list[-1])
                # outer_loop_tag = gmsh.model.occ.add_curve_loop(
                #     [outer_spline_tag, symm_axis_tag, bottom_tag]
                # )
                if j > 0:
                    edge_tag = gmsh.model.occ.add_line(bottom_point_list[j-1], bottom_point_list[j])
                    edge_loop_tag = gmsh.model.occ.add_curve_loop(
                        [outer_spline_list[j], outer_spline_list[j-1], edge_tag])
                    edge_surf_list.append(gmsh.model.occ.add_bspline_filling(edge_loop_tag, type="Curved"))
                    edge_segments.append(edge_tag)
            else:
                symm_axis_tag = gmsh.model.occ.add_line(outer_tag_list[0], outer_tag_list[-1])
                outer_loop_tag = gmsh.model.occ.add_curve_loop([outer_spline_list[j], symm_axis_tag])
        # half_bottom_len = np.floor(len(bottom_point_list)/2)
        # bottom_spline1 = gmsh.model.occ.add_spline(bottom_point_list[0:int(half_bottom_len+1)])
        # bottom_spline2 = gmsh.model.occ.add_spline(bottom_point_list[int(half_bottom_len):])
        bottom_loop = gmsh.model.occ.add_curve_loop(edge_segments)
        bottom_surf = gmsh.model.occ.add_plane_surface([bottom_loop])
        # cell_surf = gmsh.model.occ.add_bezier_surface(all_points_list, len(rValsOuter))
        # cell_outer_surf = gmsh.model.occ.fuse([edge_surf_list[0]], edge_surf_list[1:])
        cur_surf_loop = gmsh.model.occ.add_surface_loop([bottom_surf, *edge_surf_list])
        outer_shape = gmsh.model.occ.add_volume([cur_surf_loop])
        outer_shape = [[(3, outer_shape)]]

        # prevScale = float(thetaExprSym.subs({"theta": thetaVec[j-1]}))
        # pt1 = gmsh.model.occ.add_point(
        #     rValsOuter[-1]*prevScale*np.cos(thetaVec[j-1]), 
        #     rValsOuter[-1]*prevScale*np.sin(thetaVec[j-1]), 0.0)
        # pt2 = gmsh.model.occ.add_point(
        #     rValsOuter[-1]*curScale*np.cos(thetaVec[j]), 
        #     rValsOuter[-1]*curScale*np.sin(thetaVec[j]), 0.0)
        # edge_tag = gmsh.model.occ.add_line(pt1, pt2)
        # edge_loop_tag = gmsh.model.occ.add_curve_loop(
        #     [outer_spline_tag, outer_spline_tag_prev, edge_tag])
        # edge_surf_tag = gmsh.model.occ.add_surface_filling(edge_loop_tag)
        # bottom_loop = gmsh.model.occ.add_curve_loop([bottom_tag, bottom_tag_prev, edge_tag])
        # bottom_plane = gmsh.model.occ.add_plane_surface([bottom_loop])
        # cur_surf_loop = gmsh.model.occ.add_surface_loop(
        #     [edge_surf_tag, cell_plane_tag[j], cell_plane_tag[j-1], bottom_plane])
        # outer_shape.append(gmsh.model.occ.add_volume([cur_surf_loop]))
    else:
        outer_tag_list = []
        for i in range(len(rValsOuter)):
            if i == 0:
                outer_tag_list.append(top_point)
            else:
                cur_tag = gmsh.model.occ.add_point(rValsOuter[i], 0.0, zValsOuter[i])
                outer_tag_list.append(cur_tag)
        outer_spline = gmsh.model.occ.add_spline(outer_tag_list)
        origin_tag = gmsh.model.occ.add_point(0, 0, 0)
        symm_axis_tag = gmsh.model.occ.add_line(origin_tag, outer_tag_list[0])
        bottom_tag = gmsh.model.occ.add_line(origin_tag, outer_tag_list[-1])
        outer_loop_tag = gmsh.model.occ.add_curve_loop(
            [outer_spline, symm_axis_tag, bottom_tag]
        )
        cell_plane_tag = gmsh.model.occ.add_plane_surface([outer_loop_tag])
        outer_shape.append(gmsh.model.occ.revolve(
            [(2, cell_plane_tag)], 0, 0, 0, 0, 0, 1, 2 * np.pi))
    
    # if len(outer_shape) > 1:
    #     outer_shape_fused = [[(3, outer_shape[0])]]
    #     for i in [1]:#range(1, len(outer_shape)):
    #         # outer_shape_list.append((3, outer_shape[i]))
    #         outer_shape_fused = gmsh.model.occ.fuse(outer_shape_fused[0], [(3, outer_shape[i])])
    #     # outer_shape = gmsh.model.occ.fuse([outer_shape_list[0]], outer_shape_list[1:])
    #     outer_shape = outer_shape_fused
    outer_shape = outer_shape[0]
    if nanopillars != "":
        nanopillar_rad = nanopillars[0]
        nanopillar_height = nanopillars[1]
        nanopillar_spacing = nanopillars[2]
        zero_idx = np.nonzero(zValsOuter <= 0.0)
        num_pillars = 2*np.floor(rValsOuter[zero_idx]/nanopillar_spacing) + 1
        rMax = nanopillar_spacing * np.floor(rValsOuter[zero_idx]/nanopillar_spacing)
        test_coords = np.linspace(-rMax[0], rMax[0], int(num_pillars[0]))
        xTest, yTest = np.meshgrid(test_coords, test_coords)
        rTest = np.sqrt(xTest**2 + yTest**2)
        keep_logic = rTest <= rValsOuter[zero_idx]-nanopillar_rad-0.1
        xTest, yTest = xTest[keep_logic], yTest[keep_logic]
        for i in range(len(xTest)):
            cyl_tag = gmsh.model.occ.add_cylinder(
                xTest[i], yTest[i], 0.0, 0, 0, nanopillar_height-nanopillar_rad, nanopillar_rad)
            cap_tag = gmsh.model.occ.add_sphere(
                xTest[i], yTest[i], nanopillar_height-nanopillar_rad, nanopillar_rad)
            cur_pillar = gmsh.model.occ.fuse([(3,cyl_tag)], [(3,cap_tag)])
            (outer_shape, outer_shape_map) = gmsh.model.occ.cut(outer_shape, cur_pillar[0])
            outer_shape_list = []
            for j in range(len(outer_shape_map)):
                if outer_shape_map[j]!=[]:
                    outer_shape_list.append(outer_shape_map[j][0])
            outer_shape = outer_shape_list

    outer_shape_tags = []
    for i in range(len(outer_shape)):
        if outer_shape[i][0] == 3:  # pull out tags associated with 3d objects
            outer_shape_tags.append(outer_shape[i][1])
    assert len(outer_shape_tags) == 1  # should be just one 3D body from the full revolution

    if innerExpr == "":
        # No inner shape in this case
        gmsh.model.occ.synchronize()
        gmsh.model.add_physical_group(3, outer_shape_tags, tag=outer_vol_tag)
        facets = gmsh.model.getBoundary([(3, outer_shape_tags[0])])
        # facets = gmsh.model.getBoundary(outer_shape)
        # assert (
        #     len(facets) == 2
        # )  # 2 boundaries because of bottom surface at z = 0, both belong to PM
        facet_tags = []
        for i in range(len(facets)):
            facet_tags.append(facets[i][1])
        gmsh.model.add_physical_group(2, facet_tags, tag=outer_marker)
    else:
        # Add inner shape
        inner_tag_list = []
        for i in range(len(rValsInner)):
            cur_tag = gmsh.model.occ.add_point(rValsInner[i], 0, zValsInner[i])
            inner_tag_list.append(cur_tag)
        inner_spline_tag = gmsh.model.occ.add_spline(inner_tag_list)
        symm_inner_tag = gmsh.model.occ.add_line(inner_tag_list[0], inner_tag_list[-1])
        inner_loop_tag = gmsh.model.occ.add_curve_loop([inner_spline_tag, symm_inner_tag])
        inner_plane_tag = gmsh.model.occ.add_plane_surface([inner_loop_tag])
        inner_shape = gmsh.model.occ.revolve([(2, inner_plane_tag)], 0, 0, 0, 0, 0, 1, 2 * np.pi)
        inner_shape_tags = []
        for i in range(len(inner_shape)):
            if inner_shape[i][0] == 3:  # pull out tags associated with 3d objects
                inner_shape_tags.append(inner_shape[i][1])
        assert len(inner_shape_tags) == 1  # should be just one 3D body from the full revolution

        # Create interface between 2 objects
        two_shapes, (outer_shape_map, inner_shape_map) = gmsh.model.occ.fragment(
            [(3, outer_shape_tags[0])], [(3, inner_shape_tags[0])]
        )
        gmsh.model.occ.synchronize()

        # Get the outer boundary
        outer_shell = gmsh.model.getBoundary(two_shapes, oriented=False)
        # assert (
        #     len(outer_shell) == 2
        # )  # 2 boundaries because of bottom surface at z = 0, both belong to PM
        # Get the inner boundary
        inner_shell = gmsh.model.getBoundary(inner_shape_map, oriented=False)
        # assert len(inner_shell) == 1
        # Add physical markers for facets
        gmsh.model.add_physical_group(
            outer_shell[0][0], [outer_shell[0][1], outer_shell[1][1]], tag=outer_marker
        )
        gmsh.model.add_physical_group(inner_shell[0][0], [inner_shell[0][1]], tag=interface_marker)

        # Physical markers for
        all_volumes = [tag[1] for tag in outer_shape_map]
        inner_volume = [tag[1] for tag in inner_shape_map]
        outer_volume = []
        for vol in all_volumes:
            if vol not in inner_volume:
                outer_volume.append(vol)
        gmsh.model.add_physical_group(3, outer_volume, tag=outer_vol_tag)
        gmsh.model.add_physical_group(3, inner_volume, tag=inner_vol_tag)

    def meshSizeCallback(dim, tag, x, y, z, lc):
        # mesh length is hEdge at the PM and hInnerEdge at the inner membrane
        # between these, the value is interpolated based on the relative distance
        # between the two membranes.
        # Inside the inner shape, the value is interpolated between hInnerEdge
        # and lc3, where lc3 = max(hInnerEdge, 0.2*maxInnerDim)
        # if innerRad=0, then the mesh length is interpolated between
        # hEdge at the PM and 0.2*maxOuterDim in the center
        rCur = np.sqrt(x**2 + y**2)
        RCur = np.sqrt(rCur**2 + (z - zMid) ** 2)
        outer_dist = np.sqrt((rCur - rValsOuter) ** 2 + (z - zValsOuter) ** 2)
        np.append(outer_dist, z)  # include the distance from the substrate
        dist_to_outer = min(outer_dist)
        if innerExpr == "":
            lc3 = 0.2 * maxOuterDim
            dist_to_inner = RCur
            in_outer = True
        else:
            inner_dist = np.sqrt((rCur - rValsInner) ** 2 + (z - zValsInner) ** 2)
            dist_to_inner = min(inner_dist)
            inner_idx = np.argmin(inner_dist)
            inner_rad = RInnerVec[inner_idx]
            R_rel_inner = RCur / inner_rad
            lc3 = max(hInnerEdge, 0.2 * maxInnerDim)
            in_outer = R_rel_inner > 1
        lc1 = hEdge
        lc2 = hInnerEdge
        if in_outer:
            lcTest = lc1 + (lc2 - lc1) * (dist_to_outer) / (dist_to_inner + dist_to_outer)
        else:
            lcTest = lc2 + (lc3 - lc2) * (1 - R_rel_inner)
        return lcTest

    gmsh.model.mesh.setSizeCallback(meshSizeCallback)
    # set off the other options for mesh size determination
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    # this changes the algorithm from Frontal-Delaunay to Delaunay,
    # which may provide better results when there are larger gradients in mesh size
    gmsh.option.setNumber("Mesh.Algorithm", 5)

    gmsh.model.mesh.generate(3)
    rank = MPI.COMM_WORLD.rank
    tmp_folder = pathlib.Path(f"tmp_3dcell_{rank}")
    tmp_folder.mkdir(exist_ok=True)
    gmsh_file = tmp_folder / "3dcell.msh"
    gmsh.write(str(gmsh_file))
    gmsh.finalize()

    # return dolfin mesh of max dimension (parent mesh) and marker functions mf2 and mf3
    dmesh, mf2, mf3 = gmsh_to_dolfin(str(gmsh_file), tmp_folder, 3, comm)
    # remove tmp mesh and tmp folder
    gmsh_file.unlink(missing_ok=False)
    tmp_folder.rmdir()
    # return dolfin mesh, mf2 (2d tags) and mf3 (3d tags)
    if return_curvature:
        if innerExpr == "":
            facet_list = [outer_marker]
            cell_list = [outer_vol_tag]
        else:
            facet_list = [outer_marker, interface_marker]
            cell_list = [outer_vol_tag, inner_vol_tag]
        kappa_mf = compute_curvature(dmesh, mf2, mf3, facet_list, cell_list)
        return (dmesh, mf2, mf3, kappa_mf)
    else:
        return (dmesh, mf2, mf3)


def create_2Dcell(
    outerExpr: str = "",
    innerExpr: str = "",
    hEdge: float = 0,
    hInnerEdge: float = 0,
    interface_marker: int = 12,
    outer_marker: int = 10,
    inner_tag: int = 2,
    outer_tag: int = 1,
    comm: MPI.Comm = d.MPI.comm_world,
    verbose: bool = False,
    half_cell: bool = True,
    return_curvature: bool = False,
    axisymm: bool = True,
) -> Tuple[d.Mesh, d.MeshFunction, d.MeshFunction]:
    """
    Creates a 2D mesh of a cell profile, with the bounding curve defined in
    terms of r and z (e.g. unit circle would be "r**2 + (z-1)**2 == 1)
    It is assumed that substrate is present at z = 0, so if the curve extends
    below z = 0 , there is a sharp cutoff.
    If half_cell = True, only have of the contour is constructed, with a
    left zero-flux boundary at r = 0.
    Can include one compartment inside another compartment.
    Recommended for use with the axisymmetric feature of SMART.

    Args:
        outerExpr: String implicitly defining an r-z curve for the outer surface
        innerExpr: String implicitly defining an r-z curve for the inner surface
        hEdge: maximum mesh size at the outer edge
        hInnerEdge: maximum mesh size at the edge
            of the inner compartment
        interface_marker: The value to mark facets on the interface with
        outer_marker: The value to mark facets on edge of the outer ellipse with
        inner_tag: The value to mark the inner ellipse surface with
        outer_tag: The value to mark the outer ellipse surface with
        comm: MPI communicator to create the mesh with
        verbose: If true print gmsh output, else skip
        half_cell: If true, consider r=0 the symmetry axis for an axisymm shape
    Returns:
        Tuple (mesh, facet_marker, cell_marker)
    """
    import gmsh

    if outerExpr == "":
        ValueError("Outer surface is not defined")

    if return_curvature:
        # create full mesh for curvature analysis and then map onto half mesh
        # if half_cell_with_curvature is True
        half_cell_with_curvature = half_cell
        half_cell = False

    rValsOuter, zValsOuter = implicit_curve(outerExpr)

    if not innerExpr == "":
        rValsInner, zValsInner = implicit_curve(innerExpr)
        zMid = np.mean(zValsInner)
        ROuterVec = np.sqrt(rValsOuter**2 + (zValsOuter - zMid) ** 2)
        RInnerVec = np.sqrt(rValsInner**2 + (zValsInner - zMid) ** 2)
        maxOuterDim = max(ROuterVec)
        maxInnerDim = max(RInnerVec)
    else:
        zMid = np.mean(zValsOuter)
        ROuterVec = np.sqrt(rValsOuter**2 + (zValsOuter - zMid) ** 2)
        maxOuterDim = max(ROuterVec)
    if np.isclose(hEdge, 0):
        hEdge = 0.1 * maxOuterDim
    if np.isclose(hInnerEdge, 0):
        hInnerEdge = 0.2 * maxOuterDim if innerExpr == "" else 0.2 * maxInnerDim
    # Create the 2D mesh using gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", int(verbose))
    gmsh.model.add("2DCell")
    # first add outer body and revolve
    outer_tag_list = []
    for i in range(len(rValsOuter)):
        cur_tag = gmsh.model.occ.add_point(rValsOuter[i], 0, zValsOuter[i])
        outer_tag_list.append(cur_tag)
    outer_spline_tag = gmsh.model.occ.add_spline(outer_tag_list)
    if not half_cell:
        outer_tag_list2 = []
        for i in range(len(rValsOuter)):
            cur_tag = gmsh.model.occ.add_point(-rValsOuter[i], 0, zValsOuter[i])
            outer_tag_list2.append(cur_tag)
        outer_spline_tag2 = gmsh.model.occ.add_spline(outer_tag_list2)
    if np.isclose(zValsOuter[-1], 0):  # then include substrate at z=0
        if half_cell:
            origin_tag = gmsh.model.occ.add_point(0, 0, 0)
            symm_axis_tag = gmsh.model.occ.add_line(origin_tag, outer_tag_list[0])
            bottom_tag = gmsh.model.occ.add_line(origin_tag, outer_tag_list[-1])
            outer_loop_tag = gmsh.model.occ.add_curve_loop(
                [outer_spline_tag, bottom_tag, symm_axis_tag]
            )
        else:
            bottom_tag = gmsh.model.occ.add_line(outer_tag_list[-1], outer_tag_list2[-1])
            outer_loop_tag = gmsh.model.occ.add_curve_loop(
                [outer_spline_tag, outer_spline_tag2, bottom_tag]
            )
    else:
        if half_cell:
            symm_axis_tag = gmsh.model.occ.add_line(outer_tag_list[0], outer_tag_list[-1])
            outer_loop_tag = gmsh.model.occ.add_curve_loop([outer_spline_tag, symm_axis_tag])
        else:
            outer_loop_tag = gmsh.model.occ.add_curve_loop([outer_spline_tag, outer_spline_tag2])
    cell_plane_tag = gmsh.model.occ.add_plane_surface([outer_loop_tag])

    if innerExpr == "":
        # No inner shape in this case
        gmsh.model.occ.synchronize()
        gmsh.model.add_physical_group(2, [cell_plane_tag], tag=outer_tag)
        facets = gmsh.model.getBoundary([(2, cell_plane_tag)])
        facet_tag_list = []
        for i in range(len(facets)):
            facet_tag_list.append(facets[i][1])
        if half_cell:  # if half, set symmetry axis to 0 (no flux)
            xmin, ymin, zmin = (-hInnerEdge / 10, -hInnerEdge / 10, -1)
            xmax, ymax, zmax = (hInnerEdge / 10, hInnerEdge / 10, max(zValsOuter) + 1)
            all_symm_bound = gmsh.model.occ.get_entities_in_bounding_box(
                xmin, ymin, zmin, xmax, ymax, zmax, dim=1
            )
            symm_bound_markers = []
            for i in range(len(all_symm_bound)):
                symm_bound_markers.append(all_symm_bound[i][1])
            gmsh.model.add_physical_group(1, symm_bound_markers, tag=0)
        gmsh.model.add_physical_group(1, facet_tag_list, tag=outer_marker)
    else:
        # Add inner shape
        inner_tag_list = []
        for i in range(len(rValsInner)):
            cur_tag = gmsh.model.occ.add_point(rValsInner[i], 0, zValsInner[i])
            inner_tag_list.append(cur_tag)
        inner_spline_tag = gmsh.model.occ.add_spline(inner_tag_list)
        if half_cell:
            symm_inner_tag = gmsh.model.occ.add_line(inner_tag_list[0], inner_tag_list[-1])
            inner_loop_tag = gmsh.model.occ.add_curve_loop([inner_spline_tag, symm_inner_tag])
        else:
            inner_tag_list2 = []
            for i in range(len(rValsInner)):
                cur_tag = gmsh.model.occ.add_point(-rValsInner[i], 0, zValsInner[i])
                inner_tag_list2.append(cur_tag)
            inner_spline_tag2 = gmsh.model.occ.add_spline(inner_tag_list2)
            inner_loop_tag = gmsh.model.occ.add_curve_loop([inner_spline_tag, inner_spline_tag2])
        inner_plane_tag = gmsh.model.occ.add_plane_surface([inner_loop_tag])
        cell_plane_list = [cell_plane_tag]
        inner_plane_list = [inner_plane_tag]

        outer_volume = []
        inner_volume = []
        all_volumes = []
        inner_marker_list = []
        outer_marker_list = []
        for i in range(len(cell_plane_list)):
            cell_plane_tag = cell_plane_list[i]
            inner_plane_tag = inner_plane_list[i]
            # Create interface between 2 objects
            two_shapes, (outer_shape_map, inner_shape_map) = gmsh.model.occ.fragment(
                [(2, cell_plane_tag)], [(2, inner_plane_tag)]
            )
            gmsh.model.occ.synchronize()

            # Get the outer boundary
            outer_shell = gmsh.model.getBoundary(two_shapes, oriented=False)
            for i in range(len(outer_shell)):
                outer_marker_list.append(outer_shell[i][1])
            # Get the inner boundary
            inner_shell = gmsh.model.getBoundary(inner_shape_map, oriented=False)
            for i in range(len(inner_shell)):
                inner_marker_list.append(inner_shell[i][1])
            for tag in outer_shape_map:
                all_volumes.append(tag[1])
            for tag in inner_shape_map:
                inner_volume.append(tag[1])

            for vol in all_volumes:
                if vol not in inner_volume:
                    outer_volume.append(vol)

        # Add physical markers for facets
        if half_cell:  # if half, set symmetry axis to 0 (no flux)
            xmin, ymin, zmin = (-hInnerEdge / 10, -hInnerEdge / 10, -1)
            xmax, ymax, zmax = (hInnerEdge / 10, hInnerEdge / 10, max(zValsOuter) + 1)
            all_symm_bound = gmsh.model.occ.get_entities_in_bounding_box(
                xmin, ymin, zmin, xmax, ymax, zmax, dim=1
            )
            symm_bound_markers = []
            for i in range(len(all_symm_bound)):
                symm_bound_markers.append(all_symm_bound[i][1])
            # note that this first call sets the symmetry axis to tag 0 and
            # this is not overwritten by the next calls to add_physical_group
            gmsh.model.add_physical_group(1, symm_bound_markers, tag=0)
        gmsh.model.add_physical_group(1, outer_marker_list, tag=outer_marker)
        gmsh.model.add_physical_group(1, inner_marker_list, tag=interface_marker)

        # Physical markers for "volumes"
        gmsh.model.add_physical_group(2, outer_volume, tag=outer_tag)
        gmsh.model.add_physical_group(2, inner_volume, tag=inner_tag)

    def meshSizeCallback(dim, tag, x, y, z, lc):
        # mesh length is hEdge at the PM and hInnerEdge at the inner membrane
        # between these, the value is interpolated based on the relative distance
        # between the two membranes.
        # Inside the inner shape, the value is interpolated between hInnerEdge
        # and lc3, where lc3 = max(hInnerEdge, 0.2*maxInnerDim)
        # if innerRad=0, then the mesh length is interpolated between
        # hEdge at the PM and 0.2*maxOuterDim in the center
        rCur = np.sqrt(x**2 + y**2)
        RCur = np.sqrt(rCur**2 + (z - zMid) ** 2)
        outer_dist = np.sqrt((rCur - rValsOuter) ** 2 + (z - zValsOuter) ** 2)
        np.append(outer_dist, z)  # include the distance from the substrate
        dist_to_outer = min(outer_dist)
        if innerExpr == "":
            lc3 = 0.2 * maxOuterDim
            dist_to_inner = RCur
            in_outer = True
        else:
            inner_dist = np.sqrt((rCur - rValsInner) ** 2 + (z - zValsInner) ** 2)
            dist_to_inner = min(inner_dist)
            inner_idx = np.argmin(inner_dist)
            inner_rad = RInnerVec[inner_idx]
            R_rel_inner = RCur / inner_rad
            lc3 = max(hInnerEdge, 0.2 * maxInnerDim)
            in_outer = R_rel_inner > 1
        lc1 = hEdge
        lc2 = hInnerEdge
        if in_outer:
            lcTest = lc1 + (lc2 - lc1) * (dist_to_outer) / (dist_to_inner + dist_to_outer)
        else:
            lcTest = lc2 + (lc3 - lc2) * (1 - R_rel_inner)
        return lcTest

    gmsh.model.mesh.setSizeCallback(meshSizeCallback)
    # set off the other options for mesh size determination
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    # this changes the algorithm from Frontal-Delaunay to Delaunay,
    # which may provide better results when there are larger gradients in mesh size
    gmsh.option.setNumber("Mesh.Algorithm", 5)

    gmsh.model.mesh.generate(2)
    rank = MPI.COMM_WORLD.rank
    tmp_folder = pathlib.Path(f"tmp_2DCell_{rank}")
    tmp_folder.mkdir(exist_ok=True)
    gmsh_file = tmp_folder / "2DCell.msh"
    gmsh.write(str(gmsh_file))
    gmsh.finalize()

    # return dolfin mesh of max dimension (parent mesh) and marker functions mf2 and mf3
    dmesh, mf2, mf3 = gmsh_to_dolfin(str(gmsh_file), tmp_folder, 2, comm)
    # remove tmp mesh and tmp folder
    gmsh_file.unlink(missing_ok=False)
    tmp_folder.rmdir()
    # return dolfin mesh, mf2 (2d tags) and mf3 (3d tags)
    if return_curvature:
        if innerExpr == "":
            facet_list = [outer_marker]
            cell_list = [outer_tag]
        else:
            facet_list = [outer_marker, interface_marker]
            cell_list = [outer_tag, inner_tag]
        if half_cell_with_curvature:  # will likely not work in parallel...
            dmesh_half, mf2_half, mf3_half = create_2Dcell(
                outerExpr,
                innerExpr,
                hEdge,
                hInnerEdge,
                interface_marker,
                outer_marker,
                inner_tag,
                outer_tag,
                comm,
                verbose,
                half_cell=True,
                return_curvature=False,
            )
            kappa_mf = compute_curvature(
                dmesh, mf2, mf3, facet_list, cell_list, half_mesh_data=(dmesh_half, mf2_half),
                axisymm=axisymm,
            )
            (dmesh, mf2, mf3) = (dmesh_half, mf2_half, mf3_half)
        else:
            kappa_mf = compute_curvature(dmesh, mf2, mf3, facet_list, cell_list, axisymm=axisymm)
        return (dmesh, mf2, mf3, kappa_mf)
    else:
        return (dmesh, mf2, mf3)
