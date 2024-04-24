# Create meshes
python3 main.py mito-preprocess -n 0 --input-mesh-file mito1_mesh.xml --input-curv-file mito1_curvature.xml --output-mesh-file --output-mesh-file meshes-mito/mito-mesh.h5 --output-curv-file meshes-mito/mito-curv.xdmf
python3 main.py mito-preprocess -n 1 --input-mesh-file mito1_mesh.xml --input-curv-file mito1_curvature.xml --output-mesh-file --output-mesh-file meshes-mito/mito-mesh-refined.h5 --output-curv-file meshes-mito/mito-curv-refined.xdmf

# Run simulations
for ntasks in 1 2 4 6 8 10 12 14 16
    mpirun -n ${ntasks} python3 main.py mito  --mesh-file meshes-mito/mito-mesh.h5 --curv-file meshes-mito/mito-curv.xdmf -o results-mito/base-${ntasks}
    sleep 5
    mpirun -n ${ntasks} python3 main.py mito  --mesh-file meshes-mito/mito-mesh-refined.h5 --curv-file meshes-mito/mito-curv-refined.xdmf -o results-mito/refined-${ntasks}