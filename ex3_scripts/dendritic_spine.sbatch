#!/bin/bash
#SBATCH --job-name="dendritic_spine"
#SBATCH --partition=defq
#SBATCH --time=3-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --output=%j-%x-stdout.txt
#SBATCH --error=%j-%x-stderr.txt


module use /cm/shared/ex3-modules/202309a/defq/modulefiles
module load python-fenics-dolfin-2019.2.0.dev0

LIBRARY_ROOT=$(realpath $(pwd)/..)
EXAMPLEDIR=$LIBRARY_ROOT/ca2+-examples
export OMPI_MCA_btl='^ofi'  # James told me to set this to remove a warning

SCRATCH_DIRECTORY=/global/D1/homes/${USER}/smart-comp-sci/dendritic_spine/${SLURM_JOBID}
mkdir -p ${SCRATCH_DIRECTORY}
echo "Scratch directory: ${SCRATCH_DIRECTORY}"
export RESULTSDIR=$SCRATCH_DIRECTORY

echo "Run command: ${LIBRARY_ROOT}/venv/bin/python3 ${LIBRARY_ROOT}/ca2+-examples/dendritic_spine.py"
# Assuming that you installed the software in a virutal environment in ${LIBRARY_ROOT}/venv
RESULTSDIR=$SCRATCH_DIRECTORY EXAMPLEDIR=$EXAMPLEDIR ${LIBRARY_ROOT}/venv/bin/python3 ${LIBRARY_ROOT}/ca2+-examples/dendritic_spine.py
cp slurm-output/${SLURM_JOBID}-* ${SCRATCH_DIRECTORY}
