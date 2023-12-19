import argparse
from pathlib import Path

def add_run_dendritic_spine_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--mesh-file",
        type=Path,
        default=Path.cwd() / "ellipseSpine_mesh" /  "ellipseSpine_mesh.h5"
    )
    parser.add_argument(
        "-o",
        "--outdir",
        type=Path,
        default=Path("results_2spines")
    )
    parser.add_argument(
        "-n",
        "--num-refinements",
        type=int,
        default=0,
    )
    parser.add_argument(
        "-dt",
        "--time-step",
        type=float,
        default=0.0002,
    )


def add_preprocess_mesh_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--input-mesh-file",
        type=Path,
        default=Path.cwd().parent / "meshes_local" / "1spine_PM10_PSD11_ERM12_cyto1_ER2.xml"
    )
    parser.add_argument(
        "--output-mesh-file",
        type=Path,
        default=Path("ellipseSpine_mesh.h5")
    )