import argparse
from pathlib import Path

def add_run_dendritic_spine_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--mesh-file",
        type=Path,
        default=Path.cwd() / "spine_mesh" /  "spine_mesh.h5"
    )
    parser.add_argument(
        "-o",
        "--outdir",
        type=Path,
        default=Path("results_dendritic_spine")
    )
    parser.add_argument(
        "-dt",
        "--time-step",
        type=float,
        default=0.0002,
    )


def add_preprocess_spine_mesh_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--input-mesh-file",
        type=Path,
        default=Path.cwd().parent / "meshes_local" / "1spine_PM10_PSD11_ERM12_cyto1_ER2_coarser.xml"
    )
    parser.add_argument(
        "--output-mesh-file",
        type=Path,
        default=Path("spine_mesh.h5")
    )
    parser.add_argument(
        "-n",
        "--num-refinements",
        type=int,
        default=0,
    )


def add_dendritic_spine_postprocess_arguments(
    parser: argparse.ArgumentParser,
) -> None:
    parser.add_argument("-i", "--results-folder", type=Path, default="./results")
    parser.add_argument("-o", "--output-folder", type=Path, default="./")
    parser.add_argument(
        "-s",
        "--skip-if-processed",
        action="store_true",
        default=False,
        help=(
            "Skip loading results from results folder "
            "if processed results are found in the output folder"
        )
    )
    parser.add_argument(
        "--use-tex",
        action="store_true",
        default=False,
        help="Use LaTex rendering for figures",
    )
    parser.add_argument("-f", "--format", type=str, default="png", help="Format of images")