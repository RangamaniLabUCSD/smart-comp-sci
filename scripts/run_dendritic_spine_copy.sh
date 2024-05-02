# Download meshes at https://zenodo.org/records/10480304
python3 main.py dendritic-spine-preprocess --input-mesh-file /root/shared/gitrepos/smart-comp-sci/meshes/1spine_PM10_PSD11_ERM12_cyto1_ER2_coarser.xml --output-mesh-file meshes-dendritic-spine/1spine_mesh_coarser.h5 --num-refinements 0
python3 main.py dendritic-spine-preprocess --input-mesh-file /root/shared/gitrepos/smart-comp-sci/meshes/1spine_PM10_PSD11_ERM12_cyto1_ER2_coarser.xml --output-mesh-file meshes-dendritic-spine/1spine_mesh_coarser_refined_1.h5 --num-refinements 1
# python3 main.py dendritic-spine-preprocess --input-mesh-file meshes-dendritic-spine/1spine_PM10_PSD11_ERM12_cyto1_ER2_coarser.xml --output-mesh-file meshes-dendritic-spine/1spine_mesh_coarser_refined_2.h5 --num-refinements 2

partition=defq
# # # Base case
for timestep in 0.002 0.001 0.0005 0.00025
do
    python3 main.py dendritic-spine --mesh-file meshes-dendritic-spine/1spine_mesh_coarser.h5 --time-step $timestep --outdir results_dendritic_spine/coarse_dt$timestep
    sleep 5
    python3 main.py dendritic-spine --mesh-file meshes-dendritic-spine/1spine_mesh_coarser_refined_1.h5 --time-step $timestep --outdir results_dendritic_spine/refined1_dt$timestep
    sleep 5
    # python3 main.py dendritic-spine --mesh-file meshes-dendritic-spine/1spine_mesh_coarser_refined_2.h5 --time-step $timestep
    # sleep 5
done