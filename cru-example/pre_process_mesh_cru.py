from pathlib import Path
import argparse
import dolfin as d

from smart import mesh_tools
import cru_args

def main(input_mesh_file, output_mesh_file, num_refinements):
    if not Path(input_mesh_file).is_file():
        raise FileNotFoundError(f"File {input_mesh_file} does not exists")

    cru_mesh = d.Mesh(Path(input_mesh_file).as_posix())
    cell_markers = d.MeshFunction("size_t", cru_mesh, 3, cru_mesh.domains())
    facet_markers = d.MeshFunction("size_t", cru_mesh, 2, cru_mesh.domains())

    if num_refinements > 0:
        print(
            f"Original mesh has {cru_mesh.num_cells()} cells, "
            f"{cru_mesh.num_facets()} facets and "
            f"{cru_mesh.num_vertices()} vertices"
        )
        d.parameters["refinement_algorithm"] = "plaza_with_parent_facets"
        for _ in range(num_refinements):
            cru_mesh = d.adapt(cru_mesh)
            cell_markers = d.adapt(cell_markers, cru_mesh)
            facet_markers = d.adapt(facet_markers, cru_mesh)
            psd_domain = d.adapt(psd_domain, cru_mesh)
            int_domain = d.adapt(int_domain, cru_mesh)
        print(
            f"Refined mesh has {cru_mesh.num_cells()} cells, "
            f"{cru_mesh.num_facets()} facets and "
            f"{cru_mesh.num_vertices()} vertices"
        )

    Path(output_mesh_file).parent.mkdir(exist_ok=True, parents=True)
    mesh_tools.write_mesh(cru_mesh, facet_markers, cell_markers, output_mesh_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    cru_args.add_preprocess_cru_mesh_arguments(parser)
    args = vars(parser.parse_args())
    raise SystemExit(main(**args))