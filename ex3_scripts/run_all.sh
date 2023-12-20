python3 main.py preprocess-spine-mesh --input-mesh-file /home/henriknf/local/src/smart-comp-sci/meshes_local/1spine_PM10_PSD11_ERM12_cyto1_ER2.xml --output-mesh-file meshes/ellipse1Spine_mesh_0.h5 --num-refinements 0
python3 main.py preprocess-spine-mesh --input-mesh-file /home/henriknf/local/src/smart-comp-sci/meshes_local/1spine_PM10_PSD11_ERM12_cyto1_ER2.xml --output-mesh-file meshes/ellipse1Spine_mesh_1.h5 --num-refinements 1
python3 main.py preprocess-spine-mesh --input-mesh-file /home/henriknf/local/src/smart-comp-sci/meshes_local/1spine_PM10_PSD11_ERM12_cyto1_ER2.xml --output-mesh-file meshes/ellipse1Spine_mesh_2.h5 --num-refinements 2
# # Base case
python3 main.py --submit-ex3 dendritic-spine --mesh-file meshes/ellipse1Spine_mesh_0.h5  --time-step 0.0002
# # Run spatial convergence
python3 main.py --submit-ex3 dendritic-spine --mesh-file meshes/ellipse1Spine_mesh_1.h5 --time-step 0.0002
python3 main.py --submit-ex3 dendritic-spine --mesh-file meshes/ellipse1Spine_mesh_2.h5 --time-step 0.0002
# Run temporal convergence
python3 main.py --submit-ex3 dendritic-spine --mesh-file meshes/ellipse1Spine_mesh_0.h5 --time-step 0.0004
python3 main.py --submit-ex3 dendritic-spine --mesh-file meshes/ellipse1Spine_mesh_0.h5 --time-step 0.0001
python3 main.py --submit-ex3 dendritic-spine --mesh-file meshes/ellipse1Spine_mesh_0.h5 --time-step 0.00005
python3 main.py --submit-ex3 dendritic-spine --mesh-file meshes/ellipse1Spine_mesh_0.h5 --time-step 0.00001
