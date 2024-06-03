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
    elif Shape[shape].value == "star":
        return "0.98 + 0.2814*cos(5*theta)"
    elif Shape[shape].value == "rect":
        return "rect0.6"
    else:
        raise ValueError(f"Invalid shape {shape}")


def refine(mesh, cell_markers, facet_markers, num_refinements, curv_markers=None):
    if num_refinements > 0:
        print(
            f"Original mesh has {mesh.num_cells()} cells, "
            f"{mesh.num_facets()} facets and "
            f"{mesh.num_vertices()} vertices"
        )
        d.parameters["refinement_algorithm"] = "plaza_with_parent_facets"
        if curv_markers is not None:
            f_orig = d.Function(d.FunctionSpace(mesh, "CG", 1))
            f_orig.set_allow_extrapolation(True)
            f_orig.vector()[:] = curv_markers.array()
        for _ in range(num_refinements):
            mesh = d.adapt(mesh)
            cell_markers = d.adapt(cell_markers, mesh)
            facet_markers = d.adapt(facet_markers, mesh)

        print(
            f"Refined mesh has {mesh.num_cells()} cells, "
            f"{mesh.num_facets()} facets and "
            f"{mesh.num_vertices()} vertices"
        )
        if curv_markers is not None:
            adapted_curv_markers = d.MeshFunction("double", mesh, 0)
            f = d.interpolate(f_orig, d.FunctionSpace(mesh, "CG", 1))
            adapted_curv_markers.array()[:] = f.vector()[:]
            return mesh, cell_markers, facet_markers, adapted_curv_markers

    return mesh, cell_markers, facet_markers, curv_markers


def main(
    mesh_folder: str,
    shape: Shape,
    hEdge: float = 0.6,
    hInnerEdge: float = 0.6,
    num_refinements: int = 0,
    full_3d=False,
):
    here = Path(__file__).parent.absolute()
    sys.path.append((here / ".." / "utils").as_posix())

    import spread_cell_mesh_generation as mesh_gen

    if full_3d:
        sym_fraction = 1
    else:
        sym_fraction = shape2symfraction(shape)

    # initialize current mesh
    if sym_fraction == 0:
        cell_mesh, facet_markers, cell_markers = mesh_gen.create_2Dcell(
            outerExpr=outerExpr13,
            innerExpr=innerExpr13,
            hEdge=hEdge,
            hInnerEdge=hInnerEdge,
            half_cell=True,
        )
        cell_mesh, cell_markers, facet_markers, _ = refine(
            cell_mesh,
            cell_markers,
            facet_markers,
            num_refinements,
        )
    else:
        print(f"h={hEdge},{hInnerEdge} theta={shape2theta(shape)} sym_fraction={sym_fraction}")
        cell_mesh, facet_markers, cell_markers = mesh_gen.create_3dcell(
            outerExpr=outerExpr13,
            innerExpr=innerExpr13,
            hEdge=hEdge,
            hInnerEdge=hInnerEdge,
            thetaExpr=shape2theta(shape),
            sym_fraction=sym_fraction,
        )

        cell_mesh, cell_markers, facet_markers, _ = refine(
            cell_mesh, cell_markers, facet_markers, num_refinements
        )

    mesh_folder = Path(mesh_folder)
    mesh_folder.mkdir(exist_ok=True, parents=True)
    mesh_file = mesh_folder / "spreadCell_mesh.h5"
    mesh_tools.write_mesh(cell_mesh, facet_markers, cell_markers, mesh_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    add_preprocess_mech_mesh_arguments(parser)
    args = vars(parser.parse_args())
    raise SystemExit(main(**args))
