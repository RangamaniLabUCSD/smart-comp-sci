import sys
import argparse
from pathlib import Path
import dolfin as d
import numpy as np
from smart import mesh_tools

from mech_parser_args import (
    add_preprocess_mech_mesh_arguments,
    Shape,
    shape2symfraction,
)


# Define expressions to be used for geometry definitions.


outerExpr15 = "(1 - z**4/(1000+z**4)) * (r**2 + z**2) + 0.4*(r**2 + (z+15)**2)*z**4 / (10 + z**4) - 225"
innerExpr15 = "(r/5.3)**2 + ((z-7.8/2)/2.4)**2 - 1"
outerExpr13 = "(1 - z**4/(2000+z**4)) * (r**2 + z**2) + 0.4*(r**2 + (z+9.72)**2)*z**4 / (15 + z**4) - 169"
innerExpr13 = "(r/5.3)**2 + ((z-9.6/2)/2.4)**2 - 1"
outerExpr10 = "(1 - z**4/(35+z**4)) * (r**2 + z**2) + 0.94*(r**2 + (z+.01)**2)*z**4 / (30 + z**4) - 100"
innerExpr10 = "(r/5.3)**2 + ((z-5)/2.4)**2 - 1"
sphCapOuter = "r**2 + z**2 - 100"
sphCapInner = "r**2 + (z-4)**2 - 9"


def shape2theta(shape: Shape) -> str:
    # Based on thetaStr = ["", "", "", star_str1, star_str1, star_str1, "rect0.6", "rect0.6", "rect0.6"]
    # Other values are
    # star_str1 = "0.98 + 0.2814*cos(5*theta)"
    # star_str2 = "0.95 + 0.4416*cos(5*theta)"
    # star_str2 = "0.9 + 0.6164*cos(5*theta)"
    if Shape[shape].value == "circle":
        return ""
    elif Shape[shape].value == "start":
        return "0.98 + 0.2814*cos(5*theta)"
    elif Shape[shape].value == "rect":
        return "rect0.6"
    else:
        raise ValueError(f"Invalid shape {shape}")


def main(mesh_folder: str, shape: Shape, hEdge=0.6, hInnerEdge=0.6):
    here = Path(__file__).parent.absolute()
    sys.path.append((here / ".." / "utils").as_posix())

    import spread_cell_mesh_generation as mesh_gen

    sym_fraction = shape2symfraction(shape)

    # initialize current mesh
    if sym_fraction == 0:
        cell_mesh, facet_markers, cell_markers, curv_markers = mesh_gen.create_2Dcell(
            outerExpr=outerExpr13,
            innerExpr=innerExpr13,
            hEdge=hEdge,
            hInnerEdge=hInnerEdge,
            half_cell=True,
            return_curvature=True,
        )
    else:
        cell_mesh, facet_markers, cell_markers = mesh_gen.create_3dcell(
            outerExpr=outerExpr13,
            innerExpr=innerExpr13,
            hEdge=hEdge,
            hInnerEdge=hInnerEdge,
            thetaExpr=shape2theta(shape),
        )
        # nanopillars=[0.5, 2.0, 2.0])
        for f in d.facets(cell_mesh):
            topology, cellIndices = mesh_tools.facet_topology(f, cell_markers)
            if topology == "boundary":
                facet_markers.set_value(f.index(), 10)

        curv_markers = mesh_tools.compute_curvature(
            cell_mesh, facet_markers, cell_markers, [10, 12], [1, 2]
        )

        # if applicable, define symmetries of current model
        if sym_fraction < 1:
            for c in d.cells(cell_mesh):
                # calculate current angle theta
                theta_cur = np.arctan2(c.midpoint().y(), c.midpoint().x())
                if theta_cur > 2 * np.pi * sym_fraction or theta_cur < 0.0:
                    cell_markers.set_value(c.index(), 0)
                    for f in d.facets(c):
                        facet_markers.set_value(f.index(), 0)

    mesh_folder = Path(mesh_folder)
    mesh_folder.mkdir(exist_ok=True)
    mesh_file = mesh_folder / "spreadCell_mesh.h5"
    mesh_tools.write_mesh(cell_mesh, facet_markers, cell_markers, mesh_file)
    # save curvatures for reference
    curv_file_name = mesh_folder / "curvatures.xdmf"
    with d.XDMFFile(str(curv_file_name)) as curv_file:
        curv_file.write(curv_markers)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    add_preprocess_mech_mesh_arguments(parser)
    args = vars(parser.parse_args())
    raise SystemExit(main(**args))
