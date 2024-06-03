# # circle sims
python3 main.py --ntasks 1 --submit-tscc mechanotransduction --mesh-folder /root/scratch/meshes/circle_hEdge_0.1_full3dfalse --time-step 0.01 --e-val 0.1 --outdir /root/scratch/circle_E0.1_full3dfalse --axisymmetric
sleep 10 # sleep between calls to avoid opening the same h5 file at the same time
python3 main.py --ntasks 1 --submit-tscc mechanotransduction --mesh-folder /root/scratch/meshes/circle_hEdge_0.1_full3dfalse --time-step 0.01 --e-val 5.7 --outdir /root/scratch/circle_E5.7_full3dfalse --axisymmetric
sleep 10
python3 main.py --ntasks 1 --submit-tscc mechanotransduction --mesh-folder /root/scratch/meshes/circle_hEdge_0.1_full3dfalse --time-step 0.01 --e-val 70000000 --outdir /root/scratch/circle_E70000000_full3dfalse --axisymmetric
python3 main.py --ntasks 1 --submit-tscc mechanotransduction --mesh-folder /root/scratch/meshes/circle_hEdge_0.3_full3dtrue --time-step 0.01 --e-val 0.1 --outdir /root/scratch/circle_E0.1_full3dtrue
sleep 10
python3 main.py --ntasks 1 --submit-tscc mechanotransduction --mesh-folder /root/scratch/meshes/circle_hEdge_0.3_full3dtrue --time-step 0.01 --e-val 5.7 --outdir /root/scratch/circle_E5.7_full3dtrue
sleep 10
python3 main.py --ntasks 1 --submit-tscc mechanotransduction --mesh-folder /root/scratch/meshes/circle_hEdge_0.3_full3dtrue --time-step 0.01 --e-val 70000000 --outdir /root/scratch/circle_E70000000_full3dtrue

# rect sims
python3 main.py --ntasks 1 --submit-tscc mechanotransduction --mesh-folder /root/scratch/meshes/rect_hEdge_0.3_full3dfalse --time-step 0.01 --e-val 0.1 --outdir /root/scratch/rect_E0.1_full3dfalse
sleep 10
python3 main.py --ntasks 1 --submit-tscc mechanotransduction --mesh-folder /root/scratch/meshes/rect_hEdge_0.3_full3dfalse --time-step 0.01 --e-val 5.7 --outdir /root/scratch/rect_E5.7_full3dfalse
sleep 10
python3 main.py --ntasks 1 --submit-tscc mechanotransduction --mesh-folder /root/scratch/meshes/rect_hEdge_0.3_full3dfalse --time-step 0.01 --e-val 70000000 --outdir /root/scratch/rect_E70000000_full3dfalse
python3 main.py --ntasks 1 --submit-tscc mechanotransduction --mesh-folder /root/scratch/meshes/rect_hEdge_0.3_full3dtrue --time-step 0.01 --e-val 0.1 --outdir /root/scratch/rect_E0.1_full3dtrue
sleep 10
python3 main.py --ntasks 1 --submit-tscc mechanotransduction --mesh-folder /root/scratch/meshes/rect_hEdge_0.3_full3dtrue --time-step 0.01 --e-val 5.7 --outdir /root/scratch/rect_E5.7_full3dtrue
sleep 10
python3 main.py --ntasks 1 --submit-tscc mechanotransduction --mesh-folder /root/scratch/meshes/rect_hEdge_0.3_full3dtrue --time-step 0.01 --e-val 70000000 --outdir /root/scratch/rect_E70000000_full3dtrue

# star sims
python3 main.py --ntasks 1 --submit-tscc mechanotransduction --mesh-folder /root/scratch/meshes/star_hEdge_0.3_full3dfalse --time-step 0.01 --e-val 0.1 --outdir /root/scratch/star_E0.1_full3dfalse
sleep 10
python3 main.py --ntasks 1 --submit-tscc mechanotransduction --mesh-folder /root/scratch/meshes/star_hEdge_0.3_full3dfalse --time-step 0.01 --e-val 5.7 --outdir /root/scratch/star_E5.7_full3dfalse
sleep 10
python3 main.py --ntasks 1 --submit-tscc mechanotransduction --mesh-folder /root/scratch/meshes/star_hEdge_0.3_full3dfalse --time-step 0.01 --e-val 70000000 --outdir /root/scratch/star_E70000000_full3dfalse
python3 main.py --ntasks 1 --submit-tscc mechanotransduction --mesh-folder /root/scratch/meshes/star_hEdge_0.3_full3dtrue --time-step 0.01 --e-val 0.1 --outdir /root/scratch/star_E0.1_full3dtrue
sleep 10
python3 main.py --ntasks 1 --submit-tscc mechanotransduction --mesh-folder /root/scratch/meshes/star_hEdge_0.3_full3dtrue --time-step 0.01 --e-val 5.7 --outdir /root/scratch/star_E5.7_full3dtrue
sleep 10
python3 main.py --ntasks 1 --submit-tscc mechanotransduction --mesh-folder /root/scratch/meshes/star_hEdge_0.3_full3dtrue --time-step 0.01 --e-val 70000000 --outdir /root/scratch/star_E70000000_full3dtrue
