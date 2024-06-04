# # Dendritic spine model
echo "Run dendritic spine example"
python3 main.py --submit-tscc --ntasks 1 dendritic-spine --mesh-file /root/scratch/meshes/spine_mesh.h5 --outdir /root/scratch/dendritic-spine-results

# # Run cru example with vs. without serca
python3 main.py --submit-tscc --ntasks 1 cru --mesh-file /root/scratch/meshes/cru_mesh.h5 --outdir /root/scratch/cru-results
sleep 10
python3 main.py --submit-tscc --ntasks 1 cru --mesh-file /root/scratch/meshes/cru_mesh.h5 --outdir /root/scratch/cru-results-noserca --no-serca