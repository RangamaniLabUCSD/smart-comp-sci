import argparse
from pathlib import Path

def add_run_cru_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--mesh-file",
        type=Path,
        default=Path.cwd() / "cru_mesh" /  "cru_mesh.h5"
    )
    parser.add_argument(
        "-o",
        "--outdir",
        type=Path,
        default=Path("results_cru")
    )
    parser.add_argument(
        "-dt",
        "--time-step",
        type=float,
        default=1e-4,
    )
    parser.add_argument(
        "--enforce-mass-conservation",
        action='store_true',
        default=True,
    )
    parser.add_argument(
        "--no-serca",
        action='store_true',
        default=False,
    )


def add_preprocess_cru_mesh_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--input-mesh-file",
        type=Path,
        default=Path.cwd().parent / "meshes_local" / "CRU_mesh.xml"
    )
    parser.add_argument(
        "--output-mesh-file",
        type=Path,
        default=Path("cru_mesh.h5")
    )
    parser.add_argument(
        "-n",
        "--num-refinements",
        type=int,
        default=0,
    )