from pathlib import Path

import dolfin as d

from smart import mesh_tools


mesh_file = Path.cwd().parent / "meshes" / "2spine_PM10_PSD11_ERM12_cyto1_ER2.xml"
print(f"Load mesh {mesh_file}")
spine_mesh = d.Mesh(mesh_file.as_posix())
cell_markers = d.MeshFunction("size_t", spine_mesh, 3, spine_mesh.domains())
facet_markers = d.MeshFunction("size_t", spine_mesh, 2, spine_mesh.domains())


mesh_folder = Path("ellipseSpine_mesh")
mesh_folder.mkdir(exist_ok=True)
new_mesh_file = mesh_folder / "ellipseSpine_mesh.h5"

mesh_tools.write_mesh(spine_mesh, facet_markers, cell_markers, new_mesh_file)
print(f"Mesh saved to {new_mesh_file}")
