import argparse
from pathlib import Path


def add_phosphorylation_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--mesh-folder",
        type=Path,
        default=Path.cwd().parent / "meshes-phosphorylation",
    )
    parser.add_argument(
        "-o", "--outdir", type=Path, default=Path("results_phosphorylation")
    )
    parser.add_argument(
        "-dt",
        "--time-step",
        type=float,
        default=0.01,
    )
    parser.add_argument(
        "--curRadius",
        type=float,
        default=1.0,
    )
    parser.add_argument(
        "--diffusion",
        type=float,
        default=10.0,
        help="Diffusion coefficient"
    )
    parser.add_argument(
        "--axisymmetric",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--rect",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--write-checkpoint",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--comparison-results-folder",
        type=Path,
        default="",
    )
    parser.add_argument(
        "--comparison-mesh-folder",
        type=Path,
        default="",
    )


def add_phosphorylation_preprocess_arguments(
    parser: argparse.ArgumentParser,
) -> None:
    parser.add_argument(
        "--mesh-folder",
        type=Path,
        default=Path.cwd().parent / "meshes-phosphorylation",
    )
    parser.add_argument("--curRadius", type=float, default=1.0)
    parser.add_argument("--hEdge", type=float, default=0.2)
    parser.add_argument("--num-refinements", type=int, default=0)
    parser.add_argument(
        "--axisymmetric",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--rect",
        action="store_true",
        default=False,
    )

def add_phosphorylation_postprocess(
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