"""
Functions to create meshes for mechanotransduction example
"""

from typing import Tuple
import pathlib
import numpy as np
import dolfin as d
from mpi4py import MPI
from sympy.parsing.sympy_parser import parse_expr
from smart.mesh_tools import implicit_curve, compute_curvature, gmsh_to_dolfin, facet_topology

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
    thetaExpr: str = "",
    use_tmp: bool = False,
    dcurv: float = 0.1,
    sym_fraction: float = 1.0,
) -> Tuple[d.Mesh, d.MeshFunction, d.MeshFunction]:
    """
    Creates a 3d cell mesh.
    The inner contour (e.g. nucleus or other organelle) is defined
    implicitly through innerExpr, which is rotated about the z axis
    to form an axisymmetric shape (e.g. unit sphere centered at (0, 2) 
    defined by innerExpr = "r**2 + (z-2)**2 - 1")
    The outer cell contour is defined in terms of 
    cylindrical coordinates r, z, and theta.
    It is assumed that r can be expressed as a function of z and theta.
    If r = r1(z)T(theta), r1 is defined implicitly by outerExpr or innerExpr 
    (e.g. circle with radius 5 defined by outerExpr = "r**2 + z**2 - 25")
    and T is defined by thetaExpr (for an axisymmetric geometry, thetaExpr = "1")
    It is assumed that substrate is present at z = 0, so if the curve extends
    below z = 0 , there is a sharp cutoff.
    
    Some relevant examples that include theta dependence:
    * thetaExpr = "T0 + T1*cos(5*theta)", where T0 and T1 are user-defined numbers, 
        describes a five-pointed star geometry
    * thetaExpr = "a*b/sqrt((b*cos(theta))**2 + (a*sin(theta))**2)", where a and b
        are user-defined numbers, describes an ellipse contact region
    
    Special cases:
    * thetaExpr = "rectAR", where AR is a number representing an aspect ratio,
        defines a rectangular geometry at the cell contact region.
    * thetaExpr not given or thetaExpr = "" -> define axisymmetric 3d geometry by
        rotating the gmsh object

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
        return_curvature: If true, return curvatures as a vertex mesh function
        thetaExpr: String defining the theta dependence of the outer shape
        dcurv: curving parameter for rectangular shape
    Returns:
        Tuple (mesh, facet_marker, cell_marker)
    Or, if return_curvature = True, Returns:
        Tuple (mesh, facet_marker, cell_marker, curvature_marker)
    """
    import gmsh

    if outerExpr == "":
        ValueError("Outer surface is not defined")

    if sym_fraction == 0.0:
        ValueError("Use create_2dcell for axisymmetric mesh")
    elif not np.isclose(1.0/sym_fraction, round(1.0/sym_fraction)):
        ValueError("sym_fraction must 1, 1/2, 1/3, 1/4, etc.")
    elif sym_fraction != 1.0 and return_curvature:
        print("Warning: Curvature calculations not supported for meshes"
                      "with sym_fraction not 1.0, not returning curvature")
        return_curvature = False

    rValsOuter, zValsOuter = implicit_curve(outerExpr)
    zMax = max(zValsOuter)

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
            theta_range1 = np.linspace(0.0, theta_crit-theta_incr, num_per_eighth)
            theta_range2 = np.linspace(theta_crit, 
                                       np.pi-theta_crit-theta_incr, 2*num_per_eighth)
            theta_range3 = np.linspace(np.pi-theta_crit, 
                                       np.pi+theta_crit-theta_incr, 2*num_per_eighth)
            theta_range4 = np.linspace(np.pi+theta_crit, 
                                       2*np.pi-theta_crit-theta_incr, 2*num_per_eighth)
            theta_range5 = np.linspace(2*np.pi-theta_crit, 
                                       2*np.pi, num_per_eighth+1)
            thetaVec = np.concatenate([theta_range1, theta_range2, theta_range3, 
                                       theta_range4, theta_range5])
            if sym_fraction < 1:
                thetaVec = thetaVec[thetaVec <= 2*np.pi*sym_fraction-theta_incr]
                thetaVec = np.concatenate([thetaVec, [2*np.pi*sym_fraction]])
            else:
                thetaVec[-1] = 0.0 # for exactness, replace 2*pi with 0.0
            scaleVec = []
            for j in range(len(thetaVec)):
                if np.abs(np.tan(thetaVec[j])) <= b_rect / a_rect:
                    scaleVec.append(np.abs(a_rect / (2*np.cos(thetaVec[j]))))
                else:
                    scaleVec.append(np.abs(b_rect / (2*np.sin(thetaVec[j]))))
                if dcurv > 0:
                    # smooth corners
                    xCur = np.cos(thetaVec[j])*scaleVec[j]
                    yCur = np.sin(thetaVec[j])*scaleVec[j]
                    if np.abs(xCur) > (a_rect/2 - dcurv) and np.abs(yCur) > (b_rect/2 - dcurv):
                        x1 = np.sign(xCur)*(a_rect/2 - dcurv)
                        y1 = np.sign(yCur)*(b_rect/2 - dcurv)
                        aPoly = dcurv**2 * (1 + np.tan(thetaVec[j])**2)
                        bPoly = 2*dcurv*np.tan(thetaVec[j])*(x1*np.tan(thetaVec[j])-y1)
                        cPoly = np.tan(thetaVec[j])*(-2*x1*y1+x1**2*np.tan(thetaVec[j]))+y1**2-dcurv**2
                        root1 = (-bPoly + np.sqrt(bPoly**2 - 4*aPoly*cPoly))/(2*aPoly)
                        root2 = (-bPoly - np.sqrt(bPoly**2 - 4*aPoly*cPoly))/(2*aPoly)
                        if x1 > 0:
                            rootVal = root1
                        else:
                            rootVal = root2
                        if np.abs(rootVal) <= 1:
                            cosCur = rootVal
                            thetaCircle = np.sign(yCur)*np.arccos(cosCur)
                            xCircle = x1 + dcurv*np.cos(thetaCircle)
                            yCircle = y1 + dcurv*np.sin(thetaCircle)
                            scaleVec[j] *= np.sqrt(xCircle**2+yCircle**2)/np.sqrt(xCur**2 + yCur**2)
                        else:
                            raise ValueError(f"Unable to smooth corners for {thetaExpr}")
            # define test values for smoothing below
            xVals = rValsOuter[-1]*np.multiply(np.array(scaleVec), np.cos(np.array(thetaVec)))
            yVals = rValsOuter[-1]*np.multiply(np.array(scaleVec), np.sin(np.array(thetaVec)))
            
        else:
            thetaVec = np.linspace(0.0, 2*np.pi*sym_fraction, round((num_theta-1)*sym_fraction)+1)
            if sym_fraction == 1:
                thetaVec[-1] = 0.0 # for exactness, replace 2*pi with 0.0
            thetaExprSym = parse_expr(thetaExpr)
            scaleVec = []
            for j in range(len(thetaVec)):
                scaleVec.append(float(thetaExprSym.subs({"theta": thetaVec[j]})))
            # define test values for smoothing below
            thetaTest = np.linspace(0.0, 2*np.pi, num_theta)
            scaleVecTest = []
            for j in range(len(thetaTest)):
                scaleVecTest.append(float(thetaExprSym.subs({"theta": thetaTest[j]})))
            xVals = rValsOuter[-1]*np.multiply(np.array(scaleVecTest), np.cos(np.array(thetaTest)))
            yVals = rValsOuter[-1]*np.multiply(np.array(scaleVecTest), np.sin(np.array(thetaTest)))
        dtheta = min(2*np.pi/(num_theta-1), hEdge/max(xVals), hEdge/max(yVals))
        # define a smooth shape (ellipsoid) to average out sharp edges
        AR = (np.max(xVals) - np.min(xVals)) / (np.max(yVals) - np.min(yVals))
        a = np.sqrt(AR)
        b = 1/np.sqrt(AR)
        thetaExprSmooth = f"{a}*{b}/sqrt(({b}*cos(theta))**2 + ({a}*sin(theta))**2)"
        thetaExprSmooth = parse_expr(thetaExprSmooth)
        scaleVecSmooth = []
        for j in range(len(thetaVec)):
            scaleVecSmooth.append(float(thetaExprSmooth.subs({"theta": thetaVec[j]})))
    else:
        thetaVec = [0]
        dtheta = hEdge/max(rValsOuter)
    cell_plane_tag = []
    bottom_point_list = []
    outer_spline_list = []
    edge_surf_list = []
    edge_segments = []
    top_point = gmsh.model.occ.add_point(0.0, 0.0, zValsOuter[0])
    all_points_list = [top_point]
    if thetaExpr != "":
        # define theta dependence
        for j in range(len(thetaVec)):
            if j == (len(thetaVec)-1) and sym_fraction == 1:
                outer_spline_list.append(outer_spline_list[0])
                bottom_point_list.append(bottom_point_list[0])
            else:
                outer_tag_list = []
                for i in range(len(rValsOuter)):
                    if i == 0:
                        outer_tag_list.append(top_point)
                    else:
                        # average out sharp edges (cell becomes more and more ellipsoidal
                        # away from the substrate)
                        zScale1 = (zMax - zValsOuter[i])/zMax
                        zScale2 = zValsOuter[i]/zMax
                        scaleCur = np.sqrt(zScale1*scaleVec[j]**2 + zScale2*scaleVecSmooth[j]**2)
                        # scaleCur = zScale1*scaleVec[j] + zScale2*scaleVecSmooth[j]
                        xCur = scaleCur*rValsOuter[i]*np.cos(thetaVec[j])
                        yCur = scaleCur*rValsOuter[i]*np.sin(thetaVec[j])
                        zCur = zValsOuter[i]
                        cur_tag = gmsh.model.occ.add_point(xCur, yCur, zCur)
                        outer_tag_list.append(cur_tag)
                        all_points_list.append(cur_tag)
                outer_spline_list.append(gmsh.model.occ.add_spline(outer_tag_list))
                bottom_point_list.append(outer_tag_list[-1])
            if j > 0:
                edge_tag = gmsh.model.occ.add_line(bottom_point_list[j-1], bottom_point_list[j])
                edge_loop_tag = gmsh.model.occ.add_curve_loop(
                    [outer_spline_list[j], outer_spline_list[j-1], edge_tag])
                edge_surf_list.append(gmsh.model.occ.add_bspline_filling(edge_loop_tag, type="Curved"))
                edge_segments.append(edge_tag)
        if sym_fraction < 1:
            origin_tag = gmsh.model.occ.add_point(0, 0, 0)
            bottom_line1 = gmsh.model.occ.add_line(origin_tag, bottom_point_list[0])
            bottom_line2 = gmsh.model.occ.add_line(origin_tag, bottom_point_list[-1])
            bottom_loop = gmsh.model.occ.add_curve_loop([bottom_line1, bottom_line2, *edge_segments])
            bottom_surf = gmsh.model.occ.add_plane_surface([bottom_loop])
            origin_line = gmsh.model.occ.add_line(origin_tag, top_point)
            front_loop = gmsh.model.occ.add_curve_loop([origin_line, outer_spline_list[0], bottom_line1])
            front_surf = gmsh.model.occ.add_plane_surface([front_loop])
            back_loop = gmsh.model.occ.add_curve_loop([origin_line, outer_spline_list[-1], bottom_line2])
            back_surf = gmsh.model.occ.add_plane_surface([back_loop])
            cur_surf_loop = gmsh.model.occ.add_surface_loop([bottom_surf, front_surf, back_surf, *edge_surf_list])
            outer_shape = gmsh.model.occ.add_volume([cur_surf_loop])
            outer_shape = [(3, outer_shape)]
        else:
            # now define total outer shape from edge_segments and edge_surf_list    
            bottom_loop = gmsh.model.occ.add_curve_loop(edge_segments)
            bottom_surf = gmsh.model.occ.add_plane_surface([bottom_loop])
            cur_surf_loop = gmsh.model.occ.add_surface_loop([bottom_surf, *edge_surf_list])
            outer_shape = gmsh.model.occ.add_volume([cur_surf_loop])
            outer_shape = [(3, outer_shape)]
    else:
        # rotate shape 2*pi in the case of no theta dependence
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
        outer_shape = gmsh.model.occ.revolve([(2, cell_plane_tag)], 0, 0, 0, 0, 0, 1, 2 * np.pi*sym_fraction)

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
        inner_shape = gmsh.model.occ.revolve([(2, inner_plane_tag)], 
                                             0, 0, 0, 0, 0, 1, 2 * np.pi * sym_fraction)
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
    if use_tmp:
        tmp_folder = pathlib.Path(f"/root/tmp/tmp_3dcell_{rank}")
    else:
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
    # fix mf2 labels in certain cases
    for f in d.facets(dmesh):
        topology, cellIndices = facet_topology(f, mf3)
        if topology == "boundary":
            mf2.set_value(f.index(), outer_marker)
        elif topology == "interface":
            mf2.set_value(f.index(), interface_marker)
        else:
            mf2.set_value(f.index(), 0)

    # if applicable, exclude no flux surfaces from boundary
    if sym_fraction < 1:
        for f in d.facets(dmesh):
            # calculate current angle theta
            theta_cur = np.arctan2(f.midpoint().y(), f.midpoint().x())
            if (np.abs(theta_cur-2*np.pi*sym_fraction) < dtheta/10 or 
                np.abs(theta_cur) < dtheta/10):
                mf2.set_value(f.index(), 0)

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
    use_tmp: bool = False,
    nanopillar: bool = False,
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
        use_tmp: argument to use tmp directory for gmsh file creation
        nanopillar: if true, add a central nanopillar 
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
            if nanopillar:
                contact_tag = gmsh.model.occ.add_point(0.3, 0, 0)
                nanopillar_tag_list = [contact_tag]
                xCur = 0.3
                zCur = 0.0
                thetaCur = 0.0
                while xCur > 0:
                    if zCur < 3.0:
                        zCur += 0.01
                    else:
                        thetaCur += np.pi/20
                        zCur = 3.0 + 0.3*np.sin(thetaCur)
                        xCur = 0.3*np.cos(thetaCur)
                    cur_tag = gmsh.model.occ.add_point(xCur, 0, zCur)
                    nanopillar_tag_list.append(cur_tag)
                left_tag = gmsh.model.occ.add_point(0, 0, 3.3)
                nanopillar_tag_list.append(left_tag)
                nanopillar_spline_tag = gmsh.model.occ.add_spline(nanopillar_tag_list)
                symm_axis_tag = gmsh.model.occ.add_line(left_tag, outer_tag_list[0])
                bottom_tag = gmsh.model.occ.add_line(contact_tag, outer_tag_list[-1])
                outer_loop_tag = gmsh.model.occ.add_curve_loop(
                    [outer_spline_tag, bottom_tag, nanopillar_spline_tag, symm_axis_tag]
                )
            else:
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
        if nanopillar:
            logical_idx = np.logical_or(rValsInner > 0.4, zValsInner > 3.4)
            rValsInner = rValsInner[logical_idx]
            zValsInner = zValsInner[logical_idx]
            rValsInner = np.append(rValsInner, 0.4)
            zValsInner = np.append(zValsInner, zValsInner[-1])
        for i in range(len(rValsInner)):
            cur_tag = gmsh.model.occ.add_point(rValsInner[i], 0, zValsInner[i])
            inner_tag_list.append(cur_tag)
        inner_spline_tag = gmsh.model.occ.add_spline(inner_tag_list)
        if half_cell:
            if nanopillar:
                contact_tag = gmsh.model.occ.add_point(0.4, 0, zValsInner[-1])
                nanopillar_tag_list = [contact_tag]
                xCur = 0.4
                zCur = zValsInner[-1]
                thetaCur = 0.0
                while xCur > 0:
                    if zCur < 3.0:
                        zCur += 0.01
                    else:
                        thetaCur += np.pi/20
                        zCur = 3.0 + 0.4*np.sin(thetaCur)
                        xCur = 0.4*np.cos(thetaCur)
                    cur_tag = gmsh.model.occ.add_point(xCur, 0, zCur)
                    nanopillar_tag_list.append(cur_tag)
                left_tag = gmsh.model.occ.add_point(0, 0, 3.4)
                nanopillar_tag_list.append(left_tag)
                nanopillar_spline_tag_inner = gmsh.model.occ.add_spline(nanopillar_tag_list)
                symm_inner_tag = gmsh.model.occ.add_line(left_tag, inner_tag_list[0])
                inner_loop_tag = gmsh.model.occ.add_curve_loop(
                    [inner_spline_tag, nanopillar_spline_tag_inner, symm_inner_tag]
                )
            else:
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
    if use_tmp:
        tmp_folder = pathlib.Path(f"/root/tmp/tmp_2DCell_{rank}")
    else:
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
