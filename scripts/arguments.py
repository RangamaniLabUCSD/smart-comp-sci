import sys
from pathlib import Path
import argparse


here = Path(__file__).parent.absolute()


def dendritic_spine_example(parser: argparse.ArgumentParser):
    sys.path.insert(0, (here / ".." / "dendritic-spine-example").as_posix())
    import dendritic_spine_args

    dendritic_spine_args.add_run_dendritic_spine_arguments(parser)


def preprocess_spine_mesh(parser: argparse.ArgumentParser):
    sys.path.insert(0, (here / ".." / "dendritic-spine-example").as_posix())
    import dendritic_spine_args

    dendritic_spine_args.add_preprocess_spine_mesh_arguments(parser)


def postprocess_spine_mesh(parser: argparse.ArgumentParser):
    sys.path.insert(0, (here / ".." / "dendritic-spine-example").as_posix())
    import dendritic_spine_args

    dendritic_spine_args.add_dendritic_spine_postprocess_arguments(parser)


def cru_example(parser: argparse.ArgumentParser):
    sys.path.insert(0, (here / ".." / "cru-example").as_posix())
    import cru_args

    cru_args.add_run_cru_arguments(parser)


def preprocess_cru_mesh(parser: argparse.ArgumentParser):
    sys.path.insert(0, (here / ".." / "cru-example").as_posix())
    import cru_args

    cru_args.add_preprocess_cru_mesh_arguments(parser)


def mechanotransduction_example(parser: argparse.ArgumentParser):
    sys.path.insert(0, (here / ".." / "mechanotransduction-example").as_posix())
    import mech_parser_args

    mech_parser_args.add_mechanotransduction_arguments(parser)


def mechanotransduction_postprocess(parser: argparse.ArgumentParser):
    sys.path.insert(0, (here / ".." / "mechanotransduction-example").as_posix())
    import mech_parser_args

    mech_parser_args.add_mechanotransduction_postprocess_arguments(parser)


def preprocess_mech_mesh(parser: argparse.ArgumentParser):
    sys.path.insert(0, (here / ".." / "mechanotransduction-example").as_posix())
    import mech_parser_args

    mech_parser_args.add_preprocess_mech_mesh_arguments(parser)


def mito_example(parser: argparse.ArgumentParser):
    sys.path.insert(0, (here / ".." / "mito-example").as_posix())
    import mito_parser_args

    mito_parser_args.add_run_mito_arguments(parser)


def preprocess_mito_mesh(parser: argparse.ArgumentParser):
    sys.path.insert(0, (here / ".." / "mito-example").as_posix())
    import mito_parser_args

    mito_parser_args.add_preprocess_mito_mesh_arguments(parser)


def phosphorylation_example(parser: argparse.ArgumentParser):
    sys.path.insert(0, (here / ".." / "phosphorylation-example").as_posix())
    import phosphorylation_parser_args

    phosphorylation_parser_args.add_phosphorylation_arguments(parser)


def phosphorylation_preprocess(parser: argparse.ArgumentParser):
    sys.path.insert(0, (here / ".." / "phosphorylation-example").as_posix())
    import phosphorylation_parser_args

    phosphorylation_parser_args.add_phosphorylation_preprocess_arguments(parser)

def postprocess_phosphorylation(parser: argparse.ArgumentParser):
    sys.path.insert(0, (here / ".." / "phosphorylation-example").as_posix())
    import phosphorylation_parser_args

    phosphorylation_parser_args.add_phosphorylation_postprocess(parser)



def setup_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # Root parser
    parser.add_argument(
        "-n",
        "--ntasks",
        default=1,
        type=int,
        help="Number of cores to use for each program",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Just print the command and do not run it",
    )


    cg = parser.add_argument_group("cluster", "Options relevant for clusters")
    cg.add_argument(
        "-p",
        "--partition",
        default="defq",
        type=str,
        help="Which partition to use on the cluster",
    )
    cluster = cg.add_mutually_exclusive_group()
    cluster.add_argument(
        "--submit-ex3",
        action="store_true",
        help="Add this flag if you want to submit the job on the ex3 cluster",
    )
    cluster.add_argument(
        "--submit-tscc",
        action="store_true",
        help="Add this flag if you want to submit the job on the hopper cluster at TSCC",
    )

    subparsers = parser.add_subparsers(dest="command")

    # Convert notebooks
    subparsers.add_parser("convert-notebooks", help="Convert notebooks to python files")

    # Dendritic spine
    preprocess_spine_mesh_parser = subparsers.add_parser(
        "dendritic-spine-preprocess", help="Preprocess mesh for dendritic spine example"
    )
    preprocess_spine_mesh(preprocess_spine_mesh_parser)

    dendritic_spine_parser = subparsers.add_parser(
        "dendritic-spine", help="Run dendritic spine example"
    )
    dendritic_spine_example(dendritic_spine_parser)

    postprocess_spine_mesh_parser = subparsers.add_parser(
        "dendritic-spine-postprocess", help="postprocess mesh for dendritic spine example"
    )
    postprocess_spine_mesh(postprocess_spine_mesh_parser)

    # CRU (calcium release unit)
    preprocess_cru_mesh_parser = subparsers.add_parser(
        "cru-preprocess", help="Preprocess mesh for CRU example"
    )
    preprocess_cru_mesh(preprocess_cru_mesh_parser)

    cru_parser = subparsers.add_parser(
        "cru", help="Run CRU example"
    )
    cru_example(cru_parser)

    # Mechanotransduction example
    preprocess_spine_mesh_parser = subparsers.add_parser(
        "mechanotransduction-preprocess", help="Preprocess mesh for mechanotransduction example"
    )
    preprocess_mech_mesh(preprocess_spine_mesh_parser)

    mechanotransduction_parser = subparsers.add_parser(
        "mechanotransduction", help="Run mechanotransduction example"
    )
    mechanotransduction_example(mechanotransduction_parser)

    mechanotransduction_postprocess_parser = subparsers.add_parser(
        "mechanotransduction-postprocess", help="Postprocess mechanotransduction example"
    )
    mechanotransduction_postprocess(mechanotransduction_postprocess_parser)


    # Mito example
    preprocess_mito_mesh_parser = subparsers.add_parser(
        "mito-preprocess", help="Preprocess mesh for mito example"
    )
    preprocess_mito_mesh(preprocess_mito_mesh_parser)

    mito_parser = subparsers.add_parser("mito", help="Run mito example")
    mito_example(mito_parser)

    # Phosphorylation example
    phosphorylation_preprocess_parser = subparsers.add_parser(
        "phosphorylation-preprocess",
        help="Preprocess mesh for phosphorylation example",
    )
    phosphorylation_preprocess(phosphorylation_preprocess_parser)

    phosphorylation_parser = subparsers.add_parser(
        "phosphorylation", help="Run phosphorylation example"
    )
    phosphorylation_example(phosphorylation_parser)

    postprocess_phosphorylation_parser = subparsers.add_parser(
        "phosphorylation-postprocess", help="Postprocess phosphorylation example"
    )
    postprocess_phosphorylation(postprocess_phosphorylation_parser)
    return parser