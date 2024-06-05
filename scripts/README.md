# Scripts

In this folder we provide an entry point for running all the notebooks. To see all the different options you can do

```bash
$ python3 main.py --help
usage: main.py [-h] [--dry-run] [--submit-ex3] [--submit-saga]
               {convert-notebooks,dendritic-spine-preprocess,dendritic-spine,mechanotransduction-preprocess,mechanotransduction,mechanotransduction-postprocess,mito-preprocess,mito,phosphorylation-preprocess,phosphorylation,phosphorylation-postprocess}
               ...

positional arguments:
  {convert-notebooks,dendritic-spine-preprocess,dendritic-spine,mechanotransduction-preprocess,mechanotransduction,mechanotransduction-postprocess,mito-preprocess,mito,phosphorylation-preprocess,phosphorylation,phosphorylation-postprocess}
    convert-notebooks   Convert notebooks to python files
    dendritic-spine-preprocess
                        Preprocess mesh for dendritic spine example
    dendritic-spine     Run dendritic spine example
    mechanotransduction-preprocess
                        Preprocess mesh for mechanotransduction example
    mechanotransduction
                        Run mechanotransduction example
    mechanotransduction-postprocess
                        Postprocess mechanotransduction example
    mito-preprocess     Preprocess mesh for mito example
    mito                Run mito example
    phosphorylation-preprocess
                        Preprocess mesh for phosphorylation example
    phosphorylation     Run phosphorylation example
    phosphorylation-postprocess
                        Postprocess phosphorylation example

optional arguments:
  -h, --help            show this help message and exit
  --dry-run             Just print the command and do not run it
  --submit-ex3          Add this flag if you want to submit the job on the ex3
                        cluster
  --submit-saga         Add this flag if you want to submit the job on the
                        saga cluster
```

There are currently 6 ways to execute the scripts. Note that options 3, 4, and 5 can be readily adapted to run on other HPC clusters as needed.

1. By running the script without any additional flags, e.g
   ```bash
   python3 main.py phosphorylation
   ```
   will run the script directly as a normal script
2. You can submit a job to the [`ex3` cluster](https://www.ex3.simula.no) by passing the `--submit-ex3` flag, e.g
   ```bash
   python3 main.py --submit-ex3 phosphorylation
   ```
3. You can submit a job to the [`SAGA` supercomputer](https://documentation.sigma2.no/hpc_machines/saga.html) by passing the `--submit-saga` flag, e.g
   ```bash
   python3 main.py --submit-saga phosphorylation
   ```
4. You can submit a job to the [`TSCC` cluster at UCSD](https://www.sdsc.edu/services/hpc/tscc/) by passing the `--submit-tscc` flag, e.g
   ```bash
   python3 main.py --submit-tscc phosphorylation
   ```
5. You can display the command that will be run without actually running the script by passing the flag `--dry-run`, e.g
   ```bash
   python3 main.py --dry-run phosphorylation
   ```
6. You can navigate to the example folders and run the notebooks directly using `jupyter`

## Setup environment

All the code in this repository depends on [`smart`](https://rangamanilabucsd.github.io/smart), which in turn depends on the [development version of legacy FEniCs](https://bitbucket.org/fenics-project/dolfin/src/master/). While `smart` is a pure Python package and can be easily installed with `pip` (i.e `python3 -m pip install fenics-smart`), the development version of FEniCs can be tricky to install, and we recommend to use [Docker](https://www.docker.com) for running the code locally, and [Spack](https://spack.io) for running on HPC systems.

### Set up environment using Docker

We provide a pre-built docker image containing both the development version of FEniCS and smart which you can pull using

```bash
docker pull ghcr.io/rangamanilabucsd/smart:v2.2.2
```

If you prefer to run the code in jupyter notebooks we also provide a docker image for that

```bash
docker pull ghcr.io/rangamanilabucsd/smart-lab:v2.2.2
```

You can read more about how to initialize a container and running the code in the [smart documentation](https://rangamanilabucsd.github.io/smart/README.html#installation).

### Set up environment using Spack

In order to setup an environment using Spack you need to first make sure to have working `C` and `fortran` compiler installed on your system. On the `ex3` computing cluster which we have used for running the experiments in the paper, you can setup a these compilers using the following commands

```bash
module use /cm/shared/ex3-modules/latest/modulefiles
module load  gcc-10.1.0
module load libgfortran-5.0.0
```

Next you need to clone Spack and activate the `spack` environment

```bash
git clone https://github.com/spack/spack.git
. ./spack/share/spack/setup-env.sh
```

The next step is to create a new Spack environment and install the dependencies.
If you know in advance that you want to submit the jobs on a specific node of the cluster, it might be beneficial to run an interactive job and install the dependencies while logged in to the node.
To start an interactive job on ex3 on the `defq` partition you can e.g run the command

```bash
srun --partition=defq --nodes=1 --ntasks-per-node=1 --time=05:00:00 --pty bash -i
```

Now we will create a new environment

```bash
spack env create fenics-dev-defq
spack env activate fenics-dev-defq
```

and we will add FEniCS development version

```bash
spack add fenics@=master%gcc@=10.1.0 + python ^python@3.10 py-pip
```

Note that we also specify the version of `gcc` that we loaded initially. We also add `pip` so that we can add additional packages (such as `smart`). Finally, you need to do

```bash
spack concretize
spack install
```

which will take a few hours to complete.

Once this is complete you can install `smart` and the additional packages needed to run the code. We have compiled the exact dependencies used when generating the figures in the paper in the file [`requirements-ex3.txt`](requirements-ex3.txt), which can be installed with

```bash
python3 -m pip install -r requirements-ex3.txt`
```

Now you would need to modify the template in [`runner.py`](runner.py) to match the specifications of your system.

## Running the scripts

All the scripts are available as Jupyter notebooks. If you want to run the examples using the `main.py` script in the following folder, you need to convert the notebooks to python files first.

### Convert notebooks to python files

In order to run the scripts on the cluster we need to first convert the notebooks to Python files.
To do this we will use [jupytext](https://jupytext.readthedocs.io/en/latest/) which is part of the requirements.
To convert all the notebooks into python files you can do

```bash
python3 main.py convert-notebooks
```

inside this folder (called `ex3_scripts`).

### Examples

Examples have a pre-processing, running and postprocessing step which needs to be run in this order.
The pre-processing step typically involves setting up the mesh and markers, while the post-processing step will generate figures.
For easy generation of figures, Jupyter notebooks can be run in each folder after either running examples locally or downloading the[``SMART Analysis data'' Zenodo dataset](https://zenodo.org/doi/10.5281/zenodo.11252054).

#### Mechanotransduction example

**Pre-process**

```bash
python3 main.py mechanotransduction-preprocess --mesh-folder meshes-mechanotransduction --shape circle --hEdge 0.3 --hInnerEdge 0.3 --num-refinements 0
```

**Running**

```bash
python3 main.py mechanotransduction --mesh-folder meshes-mechanotransduction --time-step 0.01 --e-val 70000000.0 --z-cutoff 70000000.0 --outdir results-mechanotransduction
```

**Post-process**

```bash
python3 main.py mechanotransduction-postprocess --results-folder results-mechanotransduction --output-folder figures-mechanotransduction
```

#### Phosphorylation example

**Pre-process**

```bash
python3 main.py phosphorylation-preprocess --mesh-folder meshes-phosphorylation -curRadius 1.0 --hEdge 0.2 --num-refinements 0
```

**Running**

```bash
python3 main.py phosphorylation --mesh-folder meshes-phosphorylation/DemoSphere.h5 --time-step 0.01 --curRadius 1.0 --outdir results_phosphorylation
```

**Post-process**

```bash
python3 main.py phosphorylation-postprocess --results-folder results-phosphorylation --output-folder figures-phosphorylation --format png
```

#### Dendritic spine example

To run this example, the mesh must be first downloaded from the [``SMART Demo Meshes" Zenodo dataset](https://zenodo.org/records/10480304).

**Pre-process**

```bash
python3 main.py dendritic-spine-preprocess --input-mesh-file meshes_local/1spine_PM10_PSD11_ERM12_cyto1_ER2.xml  --output-mesh-file meshes-dendritic-spine/spine_mesh.h5 --num-refinements 0
```

**Running**

```bash
python3 main.py dendritic-spine --mesh-folder meshes-dendritic-spine/spine_mesh.h5 --time-step 0.0002 --outdir results_dendritic-spine
```

**Post-process**

```bash
python3 main.py dendritic-spine-postprocess --results-folder results_dendritic-spine --output-folder figures-dendritic-spine --format png
```

#### CRU example

To run this example, the mesh must be first downloaded from the [``SMART Demo Meshes" Zenodo dataset](https://zenodo.org/records/10480304).

**Pre-process**

```bash
python3 main.py cru-preprocess --input-mesh-file meshes_local/CRU_mesh.xml --output-mesh-file meshes-cru/cru_mesh.h5
```

**Running**

```bash
python3 main.py cru --mesh-file meshes-cru/cru_mesh.h5 --outdir cru-results
```

#### Mito example

To run this example, the mesh must be first downloaded from the [``SMART Demo Meshes" Zenodo dataset](https://zenodo.org/records/10480304).
Currently, the curvature-sensitive distribution is not used, but we test with vs. without restricting inner membrane species to the cristae.

**Pre-process**

```bash
python3 main.py mito-preprocess --input-mesh-file meshes_local/mito1_coarser2_mesh.xml --output-mesh-file meshes-mito/mito_mesh.h5 --input-curv-file meshes_local/mito1_coarser2_curvature.xml --output-curv-file meshes-mito/mito_curv.xdmf --single-compartment-im
```

**Running**

```bash
python3 mito --mesh-file meshes-mito/mito_mesh.h5  --curv-file meshes-mito/mito_curv.xdmf --outdir mito-results --single-compartment-im
```
