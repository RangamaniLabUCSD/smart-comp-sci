import sys
import argparse
from pathlib import Path
import dolfin as d
import numpy as np
from smart import mesh_tools

from phosphorylation_parser_args import add_phosphorylation_preprocess_arguments


def refinement2hEdge(refinement: int, hEdge: float = 0.2) -> float:
    return hEdge / 2**refinement


def main(
    mesh_folder: str,
    curRadius: float = 1.0,
    hEdge: float = 0.2,
    num_refinements: int = 0,
    axisymmetric: bool = False,
    rect: bool = False,
):
    print(f"Generating mesh with hEdge={hEdge}")
    # initialize current mesh
    if rect:
        xSize, ySize, zSize = curRadius * 10, curRadius * 10, curRadius
        cell_mesh = d.BoxMesh(
            d.Point(0.0, 0.0, 0.0), d.Point(xSize, ySize, zSize), 10, 10, 2
        )
        for _ in range(num_refinements):
            cell_mesh = d.refine(cell_mesh, redistribute=False)

        cell_markers = d.MeshFunction("size_t", cell_mesh, 3, 1)
        facet_markers = d.MeshFunction("size_t", cell_mesh, 2, 0)
        z0 = d.CompiledSubDomain("near(x[2], 0.0)")
        z1 = d.CompiledSubDomain("near(x[2], zSize)", zSize=zSize)
        facet_markers.set_all(0)
        z0.mark(facet_markers, 10)
        z1.mark(facet_markers, 10)

        # for f in d.facets(cell_mesh):
        #     x, y, z = f.midpoint()[:]
        #     if np.isclose(z, 0.) or np.isclose(z, zSize):
        #             # or np.isclose(y, 0.) or np.isclose(y, ySize)\
        #             # or np.isclose(z, 0.) or np.isclose(z, zSize):
        #         facet_markers[f] = 10
    elif axisymmetric:
        cell_mesh, facet_markers, cell_markers = mesh_tools.create_2Dcell(
            outerExpr=f"r**2 + (z-({curRadius}+1))**2 - {curRadius}**2",
            innerExpr="",
            hEdge=refinement2hEdge(num_refinements, hEdge),
            hInnerEdge=refinement2hEdge(num_refinements, hEdge),
        )
    else:
        cell_mesh, facet_markers, cell_markers = mesh_tools.create_spheres(
            curRadius,
            0,
            hEdge=refinement2hEdge(num_refinements, hEdge),
            hInnerEdge=refinement2hEdge(num_refinements, hEdge),
        )

    mesh_folder = Path(mesh_folder)
    mesh_folder.mkdir(exist_ok=True, parents=True)
    mesh_file = mesh_folder / "DemoSphere.h5"
    with d.XDMFFile((mesh_folder / "mesh.xdmf").as_posix()) as f:
        f.write(cell_mesh)
    with d.XDMFFile((mesh_folder / "facet_markers.xdmf").as_posix()) as f:
        f.write(facet_markers)
    with d.XDMFFile((mesh_folder / "cell_markers.xdmf").as_posix()) as f:
        f.write(cell_markers)
    mesh_tools.write_mesh(cell_mesh, facet_markers, cell_markers, mesh_file)
    print("Saved mesh to", mesh_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    add_phosphorylation_preprocess_arguments(parser)
    args = vars(parser.parse_args())
    raise SystemExit(main(**args))
