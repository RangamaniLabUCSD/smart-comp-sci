import argparse
from pathlib import Path
from enum import Enum


class Shape(str, Enum):
    circle = "circle"
    star = "star"
    rect = "rect"


def shape2symfraction(shape: Shape):
    if Shape[shape].value == "circle":
        return 0
    elif Shape[shape].value == "star":
        return 1 / 10
    elif Shape[shape].value == "rect":
        return 1 / 4
    else:
        raise ValueError(f"Invald shape {shape}")


# EVals = [0.1, 5.7, 7e7, 0.1, 5.7, 7e7, 0.1, 5.7, 7e7]
# shapeStr = ["circle","circle","circle", "star", "star", "star", "rect", "rect", "rect"]
# z_cutoff = 1e-4*np.ones(9)
# thetaStr = ["", "", "", star_str1, star_str1, star_str1, "rect0.6", "rect0.6", "rect0.6"]
# symmList specifications: if zero, consider axisymmetric.
# if 1, no symmetries.
# If a fraction, specifies the symmetry of the shape
# (e.g. 1/2 means the shape is symmetric about the x=0 axis, 1/4 means symmetric about x=0 and y=0)
# symmList = [0, 0, 0, 1/2, 1/2, 1/2, 1/4, 1/4, 1/4]


def add_mechanotransduction_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--mesh-folder",
        type=Path,
        default=Path.cwd().parent / "meshes_local" / "spreadCell_mesh_circle_R13smooth",
    )
    parser.add_argument(
        "-o", "--outdir", type=Path, default=Path("results_mechanotransduction")
    )
    parser.add_argument(
        "-dt",
        "--time-step",
        type=float,
        default=0.01,
    )
    parser.add_argument(
        "--e-val",
        type=float,
        default=7e7,
    ) # stiffness for glass coverslip
    parser.add_argument(
        "--z-cutoff",
        type=float,
        default=1e-4,
    ) # stimulation only below z-cutoff
    parser.add_argument(
        "--axisymmetric",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--well-mixed",
        action="store_true",
        default=False,
    )


def add_preprocess_mech_mesh_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--mesh-folder",
        type=Path,
        default=Path.cwd().parent / "meshes_local" / "spreadCell_mesh_circle_R13smooth",
    )
    parser.add_argument(
        "--shape", type=str, default="circle", choices=Shape._member_names_
    )
    parser.add_argument("--hEdge", type=float, default=0.6)
    parser.add_argument("--hInnerEdge", type=float, default=0.6)
    parser.add_argument("--num-refinements", type=int, default=0)
    parser.add_argument("--full-3d", action="store_true", default=False)


def add_mechanotransduction_postprocess_arguments(
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