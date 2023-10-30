from pathlib import Path
import smart
from IPython import embed
import argparse
import dolfin as d
d.parameters["refinement_algorithm"] = "plaza_with_parent_facets"


def convert_from_h5(in_path: Path | str, out_path: Path | str):
    in_path = Path(in_path).absolute()
    out_path = Path(out_path).absolute()
    assert in_path.suffix == ".h5"
    assert out_path.suffix == ".xdmf"
    mesh = d.Mesh(d.MPI.comm_world)
    with d.HDF5File(d.MPI.comm_world, str(in_path), "r") as hdf5, d.XDMFFile(d.MPI.comm_world, str(out_path)) as xdmf:
        hdf5.read(mesh, "/mesh", False)
        cf = d.MeshFunction("size_t", mesh, mesh.topology().dim(), value=0)
        print(mesh.topology().dim(), mesh.topology().dim()-1)
        ff = d.MeshFunction("size_t", mesh, mesh.topology().dim()-1, value=0)
        hdf5.read(cf, "/mf2")
        from IPython import embed
        embed()
        hdf5.read(ff, "/mf1")

        xdmf.write(mesh)
        xdmf.write(cf)
        xdmf.write(ff)


outerExpr = "(1 - z**4/(1000+z**4)) * (r**2 + z**2) + 0.4*(r**2 + (z+15)**2)*z**4 / (10 + z**4) - 225"
innerExpr = "(r/5.3)**2 + ((z-7.8/2)/2.4)**2 - 1"

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--hEdge",
                    dest="hEdge", default=0.3, type=float,
                    help="Maximum mesh size for coarsest mesh at outer edge")
parser.add_argument("--hInnerEdge",
                    dest="hInnerEdge", default=0.5, type=float,
                    help="Maximum mesh size for coarsest mesh at inner edge")
parser.add_argument("--num-refinements",
                    dest="num_ref", default=0, type=int,
                    help="Number of mesh refinements", choices=[0, 1, 2, 3])
parser.add_argument("--outerExpr", default=outerExpr, type=str,
                    help="Expression for outer boundary as an `r-z` curve")
parser.add_argument("--innerExpr", default=innerExpr, type=str,
                    help="Expression for inner boundary as an `r-z` curve")
parser.add_argument("--mesh-dir", dest="mesh_dir", default="spreadCell_mesh",
                    type=str, help="Directory for outputting mesh files")
args = parser.parse_args()

mesh_folder = Path(args.mesh_dir)
mesh_folder.mkdir(exist_ok=True)


if d.MPI.comm_world.rank == 0:
    cell_mesh, facet_markers, cell_markers = smart.mesh_tools.create_2Dcell(
        outerExpr=outerExpr,
        innerExpr=innerExpr,
        hEdge=args.hEdge,
        hInnerEdge=args.hInnerEdge,
        comm=d.MPI.comm_self,
        verbose=True)

    mesh_file = mesh_folder / "spreadCell_mesh_0.h5"
    smart.mesh_tools.write_mesh(
        cell_mesh, facet_markers, cell_markers, mesh_file)

d.MPI.comm_world.Barrier()


mesh = d.Mesh(d.MPI.comm_world)
with d.HDF5File(mesh.mpi_comm(), str((mesh_folder / "spreadCell_mesh_0.h5").absolute()), "r") as hdf5:
    hdf5.read(mesh, "/mesh", False)
    cf = d.MeshFunction("size_t", mesh, mesh.topology().dim(), value=0)
    cf.rename("cell_marker", "")
    ff = d.MeshFunction("size_t", mesh, mesh.topology().dim()-1, value=0)
    ff.rename("facet_marker", "")
    hdf5.read(cf, "/mf2")
    hdf5.read(ff, "/mf1")


for i in range(args.num_ref):
    mesh = d.adapt(mesh)
    cf = d.adapt(cf, mesh)
    ff = d.adapt(cf, mesh)
    with d.HDF5File(mesh.mpi_comm(), str((mesh_folder / f"spreadCell_mesh_{i}.h5").absolute()), "w") as hdf5:
        hdf5.write(mesh, "/mesh")
        hdf5.write(cf, "/mf2")
        hdf5.write(ff, "/mf1")
    convert_from_h5(
        mesh_folder / f"spreadCell_mesh_{i}.h5", f"spreadCell_mesh_{i}.xdmf")
