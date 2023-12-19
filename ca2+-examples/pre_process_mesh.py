from pathlib import Path
import argparse
import dolfin as d

from smart import mesh_tools
import ca2_parser_args

def main(input_mesh_file, output_mesh_file):
    if not Path(input_mesh_file).is_file():
        raise FileNotFoundError(f"File {input_mesh_file} does not exists")
    

    spine_mesh = d.Mesh(Path(input_mesh_file).as_posix())
    cell_markers = d.MeshFunction("size_t", spine_mesh, 3, spine_mesh.domains())
    facet_markers = d.MeshFunction("size_t", spine_mesh, 2, spine_mesh.domains())

    # for i in range(len(facet_array)):
    #     if (
    #         facet_array[i] == 11
    #     ):  # this indicates PSD; in this case, set to 10 to indicate it is a part of the PM
    #         facet_array[i] = 10

    Path(output_mesh_file).parent.mkdir(exist_ok=True, parents=True)
    mesh_tools.write_mesh(spine_mesh, facet_markers, cell_markers, output_mesh_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    ca2_parser_args.add_preprocess_mesh_arguments(parser)
    args = vars(parser.parse_args())
    raise SystemExit(main(**args))