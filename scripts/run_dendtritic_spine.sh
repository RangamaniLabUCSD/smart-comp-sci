# Download meshes at https://zenodo.org/records/10480304
# python3 main.py dendritic-spine-preprocess --input-mesh-file meshes-dendritic-spine/1spine_PM10_PSD11_ERM12_cyto1_ER2_coarser.xml --output-mesh-file meshes-dendritic-spine/1spine_mesh_coarser.h5 --num-refinements 0
# python3 main.py dendritic-spine-preprocess --input-mesh-file meshes-dendritic-spine/1spine_PM10_PSD11_ERM12_cyto1_ER2_coarser.xml --output-mesh-file meshes-dendritic-spine/1spine_mesh_coarser_refined_1.h5 --num-refinements 1
# python3 main.py dendritic-spine-preprocess --input-mesh-file meshes-dendritic-spine/1spine_PM10_PSD11_ERM12_cyto1_ER2_coarser.xml --output-mesh-file meshes-dendritic-spine/1spine_mesh_coarser_refined_2.h5 --num-refinements 2
# python3 main.py dendritic-spine-preprocess --input-mesh-file meshes-dendritic-spine/1spine_PM10_PSD11_ERM12_cyto1_ER2.xml --output-mesh-file meshes-dendritic-spine/1spine_mesh.h5 --num-refinements 0
# python3 main.py dendritic-spine-preprocess --input-mesh-file meshes-dendritic-spine/2spine_PM10_PSD11_ERM12_cyto1_ER2.xml --output-mesh-file meshes-dendritic-spine/2spine_mesh.h5 --num-refinements 0
partition=defq
# # # Base case
for ntasks in 1
# 2 4 6 8 10 12 14 16
do
    # for timestep in  0.002 0.001 0.0005 0.00025 0.000125
    for timestep in  0.000125
    do
        python3 main.py --submit-ex3 --ntasks $ntasks --partition=$partition dendritic-spine --mesh-file meshes-dendritic-spine/1spine_mesh_coarser.h5 --time-step $timestep 
        sleep 5
        python3 main.py --submit-ex3 --ntasks $ntasks --partition=$partition dendritic-spine --mesh-file meshes-dendritic-spine/1spine_mesh_coarser_refined_1.h5 --time-step $timestep
        # sleep 5
        # python3 main.py --submit-ex3 --ntasks $ntasks --partition=$partition dendritic-spine --mesh-file meshes-dendritic-spine/1spine_mesh_coarser_refined_2.h5 --time-step $timestep
        # sleep 5
        # python3 main.py --submit-ex3 --ntasks $ntasks --partition=$partition dendritic-spine --mesh-file meshes-dendritic-spine/1spine_mesh.h5 --time-step $timestep # -o results-dendritic-spine/1spine_timestep_${timestep}
        # sleep 5
        # python3 main.py --submit-ex3 --partition=$partition dendritic-spine --mesh-file meshes-dendritic-spine/2spine_mesh.h5 --time-step $timestep -o results-dendritic-spine/2spine_timestep_${timestep}
        # sleep 5
    done
done