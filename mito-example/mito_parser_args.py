import argparse
from pathlib import Path

def add_run_mito_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--mesh-file",
        type=Path,
        default=Path("mito_mesh.h5")
    )
    parser.add_argument(
        "--curv-file",
        type=Path,
        default=Path("mito_curv.xdmf")
    )
    parser.add_argument(
        "-o",
        "--outdir",
        type=Path,
        default=Path("results")
    )
    parser.add_argument(
        "-dt",
        "--time-step",
        type=float,
        default=2e-5,
    )
    parser.add_argument(
        "--curv-dep",
        type=float,
        default=0,
    )
    parser.add_argument(
        "--D",
        type=float,
        default=15.0
    )
    parser.add_argument(
        "--single-compartment-im",
        action="store_true",
    )


def add_preprocess_mito_mesh_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--input-mesh-file",
        type=Path,
        default=Path.cwd().parent / "meshes" / "mito1_mesh.xml"
    )
    parser.add_argument(
        "--input-curv-file",
        type=Path,
        default=Path.cwd().parent / "meshes" / "mito1_curvature.xml"
    )
    parser.add_argument(
        "--output-mesh-file",
        type=Path,
        default=Path("mito_mesh.h5")
    )
    parser.add_argument(
        "--output-curv-file",
        type=Path,
        default=Path.cwd().parent / "meshes" / "mito_curv.xdmf"
    )
    parser.add_argument(
        "-n",
        "--num-refinements",
        type=int,
        default=0,
    )
    parser.add_argument(
        "--single-compartment-im",
        action="store_true",
    )