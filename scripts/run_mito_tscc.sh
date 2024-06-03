echo "Running case 1, curv dep = -20"
python3 main.py --ntasks 1 --submit-tscc mito --mesh-file /root/scratch/meshes/mito_mesh.h5  --curv-file /root/scratch/meshes/mito_curv.xdmf --outdir /root/scratch/curvdepneg20 --curv-dep=-20
sleep 10 # avoid opening mesh file at the same time in other cases

echo "Running case 2, curv dep = -10"
python3 main.py --ntasks 1 --submit-tscc mito --mesh-file /root/scratch/meshes/mito_mesh.h5  --curv-file /root/scratch/meshes/mito_curv.xdmf --outdir /root/scratch/curvdepneg10 --curv-dep=-10
sleep 10

echo "Running case 3, curv dep = 0"
python3 main.py --ntasks 1 --submit-tscc mito --mesh-file /root/scratch/meshes/mito_mesh.h5  --curv-file /root/scratch/meshes/mito_curv.xdmf --outdir /root/scratch/curvdep0 --curv-dep=0
sleep 10

echo "Running case 4, curv dep = 10"
python3 main.py --ntasks 1 --submit-tscc mito --mesh-file /root/scratch/meshes/mito_mesh.h5  --curv-file /root/scratch/meshes/mito_curv.xdmf --outdir /root/scratch/curvdep10 --curv-dep=10
sleep 10

echo "Running case 5, curv dep = 20"
python3 main.py --ntasks 1 --submit-tscc mito --mesh-file /root/scratch/meshes/mito_mesh.h5  --curv-file /root/scratch/meshes/mito_curv.xdmf --outdir /root/scratch/curvdep20 --curv-dep=20

echo "Running fast diffusion, D = 150.0"
python3 main.py --ntasks 1 --submit-tscc mito --mesh-file /root/scratch/meshes/mito_mesh.h5  --curv-file /root/scratch/meshes/mito_curv.xdmf --outdir /root/scratch/D150 --curv-dep=0 --D 150
