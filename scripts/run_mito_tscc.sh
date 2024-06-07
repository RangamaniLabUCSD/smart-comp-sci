echo "Running case 1, uniform distribution in inner membrane"
python3 main.py --ntasks 1 --submit-tscc mito --mesh-file /root/scratch/meshes/mito_mesh.h5  --curv-file /root/scratch/meshes/mito_curv.xdmf --outdir /root/scratch/single-compartment-im --single-compartment-im
sleep 10 # avoid opening mesh file at the same time in other cases

echo "Running case 2, cristae-localized species"
python3 main.py --ntasks 1 --submit-tscc mito --mesh-file /root/scratch/meshes/mito_mesh.h5  --curv-file /root/scratch/meshes/mito_curv.xdmf --outdir /root/scratch/cristae-loc
sleep 10

echo "Running fast diffusion, D = 150.0"
python3 main.py --ntasks 1 --submit-tscc mito --mesh-file /root/scratch/meshes/mito_mesh.h5  --curv-file /root/scratch/meshes/mito_curv.xdmf --outdir /root/scratch/D150 --single-compartment-im --D 150
