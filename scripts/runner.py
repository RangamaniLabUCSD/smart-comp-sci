import sys
import json
from textwrap import dedent
from pathlib import Path
import subprocess as sp

ex3_template = dedent(
    """#!/bin/bash
#SBATCH --job-name="{job_name}"
#SBATCH --partition={partition}
#SBATCH --time=6-00:00:00
#SBATCH --ntasks={ntasks}
#SBATCH --output=%j-%x-stdout.txt
#SBATCH --error=%j-%x-stderr.txt

module purge
module use /cm/shared/ex3-modules/latest/modulefiles
module load  gcc-10.1.0
module load libgfortran-5.0.0
. /home/henriknf/local/src/spack/share/spack/setup-env.sh
spack env activate fenics-dev-{partition}

SCRATCH_DIRECTORY=/global/D1/homes/${{USER}}/smart-comp-sci/{job_name}/${{SLURM_JOBID}}
mkdir -p ${{SCRATCH_DIRECTORY}}
echo "Scratch directory: ${{SCRATCH_DIRECTORY}}"

echo 'Run command: python {script} --outdir "${{SCRATCH_DIRECTORY}}" {args}'
echo 'Partition: {partition}'
echo 'ntasks: {ntasks}'
mpirun -n {ntasks} python {script} --outdir "${{SCRATCH_DIRECTORY}}" {args}
# Move log file to results folder
mv ${{SLURM_JOBID}}-* ${{SCRATCH_DIRECTORY}}
"""
)

tscc_template = dedent(
    """#!/bin/bash
#SBATCH --job-name="{job_name}"
#SBATCH --partition=condo
#SBATCH --time=100:00:00
#SBATCH --ntasks={ntasks}
#SBATCH --output=%j-%x-stdout.txt
#SBATCH --error=%j-%x-stderr.txt
#SBATCH --account=csd786
#SBATCH --qos=condo

module load singularitypro/3.11
module load mpich/ge/gcc/64/3.4.2
cd /tscc/nfs/home/eafrancis/gitrepos/smart-comp-sci/scripts

SCRATCH_DIRECTORY=/tscc/lustre/ddn/scratch/eafrancis/smart-comp-sci-data/{job_name}
mkdir -p ${{SCRATCH_DIRECTORY}}
echo "Scratch directory: ${{SCRATCH_DIRECTORY}}"

echo 'Run command in container: python {script} {args}'
singularity exec --bind $HOME:/root/shared,\
$TMPDIR:/root/tmp,$SCRATCH_DIRECTORY:/root/scratch \
/tscc/nfs/home/eafrancis/smart-newmeshview.sif \
python3 {script} {args}

# Move log file to results folder
mv ${{SLURM_JOBID}}-* ${{SCRATCH_DIRECTORY}}
"""
)

here = Path(__file__).parent.absolute()


def run(
    args,
    dry_run: bool,
    script: str,
    submit_ex3: bool,  # FIXME: Change this and next param to enum
    submit_tscc: bool,
    job_name: str = "",
    ntasks: int = 1,
    partition: str = "defq",
):
    in_args = list(map(str, args))
    args_str = " ".join(in_args)
    if dry_run:
        print(f"Run command: {sys.executable} {script} {args_str}")
        return

    if submit_ex3:
        template = ex3_template
    elif submit_tscc:
        template = tscc_template
    else:
        sp.run(["mpirun", "-n", str(ntasks), sys.executable, script, *in_args])
        return

    job_file = Path("tmp_job.sbatch")
    job_file.write_text(
        template.format(
            job_name=job_name,
            script=script,
            args=args_str,
            ntasks=ntasks,
            partition=partition,
        )
    )
    sp.run(["sbatch", job_file.as_posix()])
    job_file.unlink()


def preprocess_phosphorylation_mesh(
    mesh_folder: Path,
    curRadius: float,
    hEdge: float,
    num_refinements: int,
    axisymmetric: bool,
    rect: bool,
    dry_run: bool,
    **kwargs,
):
    args = [
        "--mesh-folder",
        Path(mesh_folder).as_posix(),
        "--curRadius",
        curRadius,
        "--hEdge",
        hEdge,
        "--num-refinements",
        num_refinements,
    ]
    if axisymmetric:
        args.append("--axisymmetric")

    if rect:
        args.append("--rect")

    script = (
        (here / ".." / "phosphorylation-example" / "pre_process_mesh.py")
        .absolute()
        .resolve()
        .as_posix()
    )
    run(
        args=args,
        dry_run=dry_run,
        script=script,
        submit_ex3=False,
        submit_tscc=False,
    )


def phosphorylation_example(
    mesh_folder: Path,
    outdir: Path,
    time_step: float,
    curRadius: float,
    axisymmetric: bool,
    rect: bool,
    diffusion: float = 10.0,
    dry_run: bool = False,
    submit_ex3: bool = False,
    submit_tscc: bool = False,
    ntasks: int = 1,
    partition: str = "defq",
    write_checkpoint: bool = False,
    comparison_results_folder: Path = "",
    comparison_mesh_folder: Path = "",
    **kwargs,
):
    args = [
        "--mesh-folder",
        Path(mesh_folder).as_posix(),
        "--time-step",
        time_step,
        "--curRadius",
        curRadius,
        "--diffusion",
        diffusion,
        "--comparison-results-folder",
        comparison_results_folder,
        "--comparison-mesh-folder",
        comparison_mesh_folder,
    ]
    if axisymmetric:
        args.append("--axisymmetric")

    if rect:
        args.append("--rect")

    if write_checkpoint:
        args.append("--write-checkpoint")

    if submit_ex3 is False:
        args.extend(["--outdir", Path(outdir).as_posix()])

    script = (
        (here / ".." / "phosphorylation-example" / "phosphorylation.py")
        .absolute()
        .resolve()
        .as_posix()
    )
    run(
        job_name="phosphorylation",
        args=args,
        dry_run=dry_run,
        script=script,
        submit_ex3=submit_ex3,
        submit_tscc=submit_tscc,
        ntasks=ntasks,
        partition=partition,
    )


def postprocess_phosphorylation(
    results_folder: Path,
    output_folder: Path,
    skip_if_processed: bool = False,
    use_tex: bool = False,
    dry_run: bool = False,
    format: str = "png",
    **kwargs,
):
    args = [
        "--results-folder",
        Path(results_folder).as_posix(),
        "--output-folder",
        Path(output_folder).as_posix(),
        "--format",
        format,
    ]
    if skip_if_processed:
        args.append("--skip-if-processed")
    if use_tex:
        args.append("--use-tex")
    script = (
        (here / ".." / "phosphorylation-example" / "postprocess.py")
        .absolute()
        .resolve()
        .as_posix()
    )
    run(
        args=args,
        dry_run=dry_run,
        script=script,
        submit_ex3=False,
        submit_tscc=False,
    )


def preprocess_mito_mesh(
    input_mesh_file: Path,
    output_mesh_file: Path,
    input_curv_file: Path,
    output_curv_file: Path,
    dry_run: bool,
    num_refinements: int,
    single_compartment_im: bool,
    **kwargs,
):
    args = [
        "--input-mesh-file",
        Path(input_mesh_file).as_posix(),
        "--output-mesh-file",
        Path(output_mesh_file).as_posix(),
        "--input-curv-file",
        Path(input_curv_file).as_posix(),
        "--output-curv-file",
        Path(output_curv_file).as_posix(),
        "--num-refinements",
        num_refinements,
    ]

    if single_compartment_im:
        args.append("--single-compartment-im")

    script = (
        (here / ".." / "mito-example" / "pre_process_mesh.py")
        .absolute()
        .resolve()
        .as_posix()
    )
    run(
        args=args,
        dry_run=dry_run,
        script=script,
        submit_ex3=False,
        submit_tscc=False,
    )


def mito_example(
    mesh_file: Path,
    curv_file: Path,
    outdir: Path,
    time_step: float,
    curv_dep: float,
    D: float,
    dry_run: bool = False,
    submit_ex3: bool = False,
    submit_tscc: bool = False,
    ntasks: int = 1,
    partition: str = "defq",
    **kwargs,
):
    args = [
        "--mesh-file",
        Path(mesh_file).as_posix(),
        "--curv-file",
        Path(curv_file).as_posix(),
        "--time-step",
        time_step,
        "--curv-dep",
        curv_dep,
        "--D",
        D,
    ]

    if not submit_ex3:
        args.extend(["--outdir", Path(outdir).as_posix()])

    script = (
        (here / ".." / "mito-example" / "mito_atp_dynamics.py")
        .absolute()
        .resolve()
        .as_posix()
    )
    run(
        job_name="mito",
        args=args,
        dry_run=dry_run,
        script=script,
        submit_ex3=submit_ex3,
        submit_tscc=submit_tscc,
        ntasks=ntasks,
        partition=partition,
    )


def pre_preprocess_mech_mesh(
    mesh_folder: Path,
    shape: str,
    hEdge: float,
    hInnerEdge: float,
    num_refinements: int,
    dry_run: bool,
    full_3d: bool,
    **kwargs,
):
    args = [
        "--mesh-folder",
        Path(mesh_folder).as_posix(),
        "--shape",
        shape,
        "--hEdge",
        hEdge,
        "--hInnerEdge",
        hInnerEdge,
        "--num-refinements",
        num_refinements,
    ]

    if full_3d:
        args.append("--full-3d")

    script = (
        (here / ".." / "mechanotransduction-example" / "pre_process_mesh.py")
        .absolute()
        .resolve()
        .as_posix()
    )
    run(
        args=args,
        dry_run=dry_run,
        script=script,
        submit_ex3=False,
        submit_tscc=False,
    )


def mechanotransduction_example(
    mesh_folder: Path,
    outdir: Path,
    time_step: float,
    e_val: float,
    z_cutoff: float,
    axisymmetric: bool,
    well_mixed: bool,
    dry_run: bool = False,
    submit_ex3: bool = False,
    submit_tscc: bool = False,
    ntasks: int = 1,
    partition: str = "defq",
    **kwargs,
):
    args = [
        "--mesh-folder",
        Path(mesh_folder).as_posix(),
        "--time-step",
        time_step,
        "--e-val",
        e_val,
        "--z-cutoff",
        z_cutoff,
    ]
    if axisymmetric:
        args.append("--axisymmetric")
    if well_mixed:
        args.append("--well-mixed")

    if submit_ex3 is False:
        args.extend(["--outdir", Path(outdir).as_posix()])

    script = (
        (here / ".." / "mechanotransduction-example" / "mechanotransduction.py")
        .absolute()
        .resolve()
        .as_posix()
    )
    run(
        job_name="mechanotransduction",
        args=args,
        dry_run=dry_run,
        script=script,
        submit_ex3=submit_ex3,
        submit_tscc=submit_tscc,
        ntasks=ntasks,
        partition=partition,
    )


def postprocess_mechanotransduction(
    results_folder: Path,
    output_folder: Path,
    skip_if_processed: bool = False,
    use_tex: bool = False,
    dry_run: bool = False,
    format: str = "png",
    **kwargs,
):
    args = [
        "--results-folder",
        Path(results_folder).as_posix(),
        "--output-folder",
        Path(output_folder).as_posix(),
        "--format",
        format,
    ]
    if skip_if_processed:
        args.append("--skip-if-processed")
    if use_tex:
        args.append("--use-tex")

    script = (
        (here / ".." / "mechanotransduction-example" / "postprocess.py")
        .absolute()
        .resolve()
        .as_posix()
    )
    run(
        args=args,
        dry_run=dry_run,
        script=script,
        submit_ex3=False,
        submit_tscc=False,
    )


def preprocess_cru_mesh(
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
        (here / ".." / "cru-example" / "pre_process_mesh_cru.py")
        .absolute()
        .resolve()
        .as_posix()
    )
    run(
        args=args,
        dry_run=dry_run,
        script=script,
        submit_ex3=False,
        submit_tscc=False,
    )


def cru_example(
    mesh_file: Path,
    outdir: Path,
    time_step: float,
    no_serca: bool,
    dry_run: bool = False,
    submit_ex3: bool = False,
    submit_tscc: bool = False,
    **kwargs,
):
    args = [
        "--mesh-file",
        Path(mesh_file).as_posix(),
        "--time-step",
        time_step,
    ]

    if no_serca:
        args.append("--no-serca")

    if submit_ex3 is False:
        args.extend(["--outdir", Path(outdir).as_posix()])

    script = (here / ".." / "cru-example" / "cru.py").absolute().resolve().as_posix()
    run(
        job_name="cru",
        args=args,
        dry_run=dry_run,
        script=script,
        submit_ex3=submit_ex3,
        submit_tscc=submit_tscc,
    )


def preprocess_spine_mesh(
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
        (here / ".." / "dendritic-spine-example" / "pre_process_mesh.py")
        .absolute()
        .resolve()
        .as_posix()
    )
    run(
        args=args,
        dry_run=dry_run,
        script=script,
        submit_ex3=False,
        submit_tscc=False,
    )


def dendritic_spine_postprocess(
    results_folder: Path,
    output_folder: Path,
    skip_if_processed: bool = False,
    use_tex: bool = False,
    dry_run: bool = False,
    format: str = "png",
    **kwargs,
):
    args = [
        "--results-folder",
        Path(results_folder).as_posix(),
        "--output-folder",
        Path(output_folder).as_posix(),
        "--format",
        format,
    ]
    if skip_if_processed:
        args.append("--skip-if-processed")
    if use_tex:
        args.append("--use-tex")

    script = (
        (here / ".." / "dendritic-spine-example" / "postprocess.py")
        .absolute()
        .resolve()
        .as_posix()
    )
    run(
        args=args,
        dry_run=dry_run,
        script=script,
        submit_ex3=False,
        submit_tscc=False,
    )


def dendritic_spine_example(
    mesh_file: Path,
    outdir: Path,
    time_step: float,
    dry_run: bool = False,
    submit_ex3: bool = False,
    submit_tscc: bool = False,
    ntasks: int = 1,
    partition: str = "defq",
    **kwargs,
):
    args = [
        "--mesh-file",
        Path(mesh_file).as_posix(),
        "--time-step",
        time_step,
    ]

    if submit_ex3 is False:
        args.extend(["--outdir", Path(outdir).as_posix()])

    script = (
        (here / ".." / "dendritic-spine-example" / "dendritic_spine.py")
        .absolute()
        .resolve()
        .as_posix()
    )
    run(
        job_name="dendritic_spine",
        args=args,
        dry_run=dry_run,
        script=script,
        submit_ex3=submit_ex3,
        submit_tscc=submit_tscc,
        ntasks=ntasks,
        partition=partition,
    )


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
