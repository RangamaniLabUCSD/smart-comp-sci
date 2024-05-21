# smart-comp-sci

Examples and numerical tests for SMART scientific computing paper.

To run the code in this repository, it is necessary to install [SMART (Spatial Modeling Algorithms for Reaction and Transport)](https://github.com/RangamaniLabUCSD/smart.git).
See more info about running the code and reproducing the results in the [scripts](scripts) folder.

## Installation

To run the scripts, we advice usage of docker, and the following base image
`ghcr.io/scientificcomputing/fenics-gmsh:2024-02-19`, which after installation of docker, can be started with

```bash
docker run -ti -v $(pwd):/root/shared -w /root/shared -p 8888:8888 --name smart-comp-sci ghcr.io/scientificcomputing/fenics-gmsh:2024-05-13@sha256:e067aeefedc074230a4155e1aa0b0d0c8e3049ee982e43e580fefc2f39fe8175
```

This should preferably be started from the root of this git repo, as `-v` shared the current directory on your computer with the docker container.

> [!NOTE]
> This uses a Docker image build for AMD (not ARM), which means that users of ARM processors. We are aiming to fix this.

This will launch a terminal with [FEniCS](https://bitbucket.org/fenics-project/dolfin/src/master/) installed.
To install the compatible version of SMART, call

```bash
python3 -m pip install git+https://github.com/RangamaniLabUCSD/smart.git@development
```

If the current repository is in your root, call

```bash
python3 -m pip install -e .[lab] -U
```

To run notebooks in your browser, call

```bash
jupyter lab --ip 0.0.0.0 --no-browser --allow-root
```

All meshes can be downloaded from the [``SMART Demo Meshes" Zenodo dataset](https://zenodo.org/records/10480304).
To run any of the Jupyter notebook versions of the examples, these meshes should be present in the main folder of the local repository.
Alternatively, the paths can be provided when running in Python scripts as described in the README in the `scripts` folder.

The output from all analyses are freely available from the [``SMART Analysis data'' Zenodo dataset]().

Shield: [![CC BY-SA 4.0][cc-by-sa-shield]][cc-by-sa]

This work is licensed under a
[Creative Commons Attribution-ShareAlike 4.0 International License][cc-by-sa].

[![CC BY-SA 4.0][cc-by-sa-image]][cc-by-sa]

[cc-by-sa]: http://creativecommons.org/licenses/by-sa/4.0/
[cc-by-sa-image]: https://licensebuttons.net/l/by-sa/4.0/88x31.png
[cc-by-sa-shield]: https://img.shields.io/badge/License-CC%20BY--SA%204.0-lightgrey.svg