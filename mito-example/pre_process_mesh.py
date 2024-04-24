from pathlib import Path
import argparse
import dolfin as d

from smart import mesh_tools
import mito_parser_args
import numpy as np


def main(
    input_mesh_file,
    output_mesh_file,
    input_curv_file,
    output_curv_file,
    num_refinements,
    single_compartment_im=False,
):
    if not Path(input_mesh_file).is_file():
        raise FileNotFoundError(f"File {input_mesh_file} does not exist")

    mito_mesh = d.Mesh(Path(input_mesh_file).as_posix())
    cell_markers = d.MeshFunction("size_t", mito_mesh, 3, mito_mesh.domains())
    facet_markers = d.MeshFunction("size_t", mito_mesh, 2, mito_mesh.domains())

    if single_compartment_im:
        facet_markers.array()[np.where(facet_markers.array() == 11)[0]] = 12
        print(
            f"Remove cristae markers, single compartment im is {single_compartment_im}"
        )

    if num_refinements > 0:
        print(
            f"Original mesh has {mito_mesh.num_cells()} cells, "
            f"{mito_mesh.num_facets()} facets and "
            f"{mito_mesh.num_vertices()} vertices"
        )
        d.parameters["refinement_algorithm"] = "plaza_with_parent_facets"
        for _ in range(num_refinements):
            mito_mesh = d.adapt(mito_mesh)
            cell_markers = d.adapt(cell_markers, mito_mesh)
            facet_markers = d.adapt(facet_markers, mito_mesh)
        print(
            f"Refined mesh has {mito_mesh.num_cells()} cells, "
            f"{mito_mesh.num_facets()} facets and "
            f"{mito_mesh.num_vertices()} vertices"
        )

    Path(output_mesh_file).parent.mkdir(exist_ok=True, parents=True)
    mesh_tools.write_mesh(mito_mesh, facet_markers, cell_markers, output_mesh_file)
    curvature = d.MeshFunction("double", mito_mesh, str(input_curv_file))
    curvature.array()[np.where(curvature.array() > 1e9)[0]] = 0
    with d.XDMFFile(str(output_curv_file)) as curv_file:
        curv_file.write(curvature)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    mito_parser_args.add_preprocess_mito_mesh_arguments(parser)
    args = vars(parser.parse_args())
    raise SystemExit(main(**args))
