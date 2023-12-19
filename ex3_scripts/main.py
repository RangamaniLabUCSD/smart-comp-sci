import sys
from textwrap import dedent
from pathlib import Path
import argparse
import shlex
import tempfile
import time
import subprocess as sp

ex3_template = dedent(
    """#!/bin/bash
#SBATCH --job-name="{job_name}"
#SBATCH --partition=defq
#SBATCH --time=3-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --output=%j-%x-stdout.txt
#SBATCH --error=%j-%x-stderr.txt


module use /cm/shared/ex3-modules/202309a/defq/modulefiles
module load python-fenics-dolfin-2019.2.0.dev0

SCRATCH_DIRECTORY=/global/D1/homes/${{USER}}/smart-comp-sci/{job_name}/${{SLURM_JOBID}}
mkdir -p ${{SCRATCH_DIRECTORY}}
echo "Scratch directory: ${{SCRATCH_DIRECTORY}}"

echo 'Run command: {python} {script} --outdir "${{SCRATCH_DIRECTORY}}" {args}'
{python} {script} --outdir "${{SCRATCH_DIRECTORY}}" {args}
# Move log file to results folder
mv ${{SLURM_JOBID}}-* ${{SCRATCH_DIRECTORY}}
"""
)


here = Path(__file__).parent.absolute()


def add_argument_dendritic_spine_example(parser: argparse.ArgumentParser):
    sys.path.insert(0, (here / ".." / "ca2+-examples").as_posix())
    import ca2_parser_args

    ca2_parser_args.add_run_dendritic_spine_arguments(parser)


def add_argument_preprocess_spine_mesh(parser: argparse.ArgumentParser):
    sys.path.insert(0, (here / ".." / "ca2+-examples").as_posix())
    import ca2_parser_args

    ca2_parser_args.add_preprocess_mesh_arguments(parser)


def run_preprocess_spine_mesh(
    input_mesh_file: Path,
    output_mesh_file: Path,
    dry_run: bool,
    num_refinements: int,
    **kwargs,
):
    args = [
        "--input-mesh-file",
        Path(input_mesh_file).as_posix(),
        "--output-mesh-file",
        Path(output_mesh_file).as_posix(),
        "--num-refinements",
        num_refinements,
    ]
    script = (
        (here / ".." / "ca2+-examples" / "pre_process_mesh.py")
        .absolute()
        .resolve()
        .as_posix()
    )
    args = list(map(str, args))
    if dry_run:
        print(f"Run command: {sys.executable} {script} {' '.join(args)}")
        return

    sp.run([sys.executable, script, *args])


def run_dendritic_spine_example(
    mesh_file: Path,
    outdir: Path,
    time_step: float,
    dry_run: bool = False,
    submit_ex3: bool = False,
    **kwargs,
):
    args = [
        "--mesh-file",
        Path(mesh_file).as_posix(),
        "--time-step",
        time_step,
    ]
    if not submit_ex3:
        args.extend(["--outdir", Path(outdir).as_posix()])

    script = (
        (here / ".." / "ca2+-examples" / "dendritic_spine.py")
        .absolute()
        .resolve()
        .as_posix()
    )
    # Turn all arguments into strings
    args = list(map(str, args))
    args_str = shlex.join(args)
    if dry_run:
        print(f"Run command: {sys.executable} {script} {args_str}")
        return

    if submit_ex3:
        job_file = Path("tmp_job.sbatch")
        job_file.write_text(
            ex3_template.format(
                job_name="dendritic-spine",
                python=sys.executable,
                script=script,
                args=args_str,
            )
        )

        sp.run(["sbatch", job_file])
        job_file.unlink()
    else:
        sp.run([sys.executable, script, *args])


def convert_notebooks(dry_run: bool = False, **kwargs):
    import jupytext

    for example_folder in (here / "..").iterdir():
        # Only look in folders that contains the word example
        if "example" not in example_folder.stem:
            continue

        # Loop over all files in that folder
        for f in example_folder.iterdir():
            if not f.suffix == ".ipynb":
                continue

            print(f"Convert {f} to {f.with_suffix('.py')}")
            if dry_run:
                continue

            text = jupytext.reads(f.read_text())
            jupytext.write(text, f.with_suffix(".py"))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Just print the command and do not run it",
    )
    parser.add_argument(
        "--submit-ex3",
        action="store_true",
        help="Add this flag if you want to submit the job on the ex3 cluster",
    )

    subparsers = parser.add_subparsers(dest="command")
    subparsers.add_parser("convert-notebooks", help="Convert notebooks to python files")

    preprocess_spine_mesh = subparsers.add_parser(
        "preprocess-spine-mesh", help="Preprocess mesh for dendritic spine example"
    )
    add_argument_preprocess_spine_mesh(preprocess_spine_mesh)

    dendritic_spine = subparsers.add_parser(
        "dendritic-spine", help="Run dendritic spine example"
    )
    add_argument_dendritic_spine_example(dendritic_spine)

    args = vars(parser.parse_args())
    # if args[""]

    if args["command"] == "convert-notebooks":
        print("Convert notebook")
        convert_notebooks(**args)

    elif args["command"] == "preprocess-spine-mesh":
        print("Run preprocess dendritic spine mesh")
        run_preprocess_spine_mesh(**args)

    elif args["command"] == "dendritic-spine":
        print("Run dendritic spine example")
        run_dendritic_spine_example(**args)


if __name__ == "__main__":
    raise SystemExit(main())
