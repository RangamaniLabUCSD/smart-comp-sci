# ex3


## Setup environment
Before submitting any scripts we need to set up the environment on the cluster. First we load the installed modules
```
module use /cm/shared/ex3-modules/202309a/defq/modulefiles
module load python-fenics-dolfin-2019.2.0.dev0
```
Next we create a python virtual environment in the root of the repo. 
```
python3 -m venv venv
```
We activate the virtual environment
```
. venv/bin/activate
```
and install the rest of the dependencies that are not already installed
```
python3 -m pip install -r requirements.txt
```

## Convert notebooks to python files
In order to run the scripts on the cluster we need to first convert the notebooks to python files. To do this we will use `jupytext` which is part of the requirements. For example you can use
```
jupytext ca2+-examples/dendritic_spine.ipynb --to py
```

## Submit job
```
sbatch dendritic_spine.sbatch 
```