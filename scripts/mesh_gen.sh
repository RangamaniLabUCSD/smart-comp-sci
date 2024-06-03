# # echo "Saving mesh for mito example..."
# python3 main.py mito-preprocess --input-mesh-file /root/shared/gitrepos/smart-comp-sci/meshes/mito1_coarser2_mesh.xml \
# --output-mesh-file /root/scratch/write-meshes/mito/meshes/mito_mesh.h5 --input-curv-file /root/shared/gitrepos/smart-comp-sci/meshes/mito1_coarser2_curvature.xml \
# --output-curv-file /root/scratch/write-meshes/mito/meshes/mito_curv.xdmf --single-compartment-im

# # echo "Saving mesh for dendritic spine example..."
# python3 main.py dendritic-spine-preprocess --input-mesh-file /root/shared/gitrepos/smart-comp-sci/meshes/1spine_PM10_PSD11_ERM12_cyto1_ER2_coarser.xml \
# --output-mesh-file /root/scratch/write-meshes/dendritic_spine/meshes/spine_mesh.h5

# # echo "Saving mesh for cru example..."
# python3 main.py cru-preprocess --input-mesh-file /root/shared/gitrepos/smart-comp-sci/meshes/CRU_mesh.xml \
# --output-mesh-file /root/scratch/write-meshes/cru/meshes/cru_mesh.h5

echo "Saving meshes for mechanotransduction example"
# Create meshes with symm considerations
# echo "Writing axisymmetric circle mesh..."
# python3 main.py mechanotransduction-preprocess --shape circle --hEdge 0.1 --hInnerEdge 0.1 --mesh-folder /root/scratch/write-meshes/mechanotransduction/meshes/circle_hEdge_0.1_full3dfalse
# echo "Writing quarter rect mesh..."
# python3 main.py mechanotransduction-preprocess --shape rect --hEdge 0.3 --hInnerEdge 0.3 --mesh-folder /root/scratch/write-meshes/mechanotransduction/meshes/rect_hEdge_0.3_full3dfalse
echo "Writing partial star mesh..."
python3 main.py mechanotransduction-preprocess --shape star --hEdge 0.3 --hInnerEdge 0.3 --mesh-folder /root/shared/write-meshes/mechanotransduction/meshes/star_hEdge_0.3_full3dfalse
# Create full 3d meshes
# echo "Writing full circle mesh..."
# python3 main.py mechanotransduction-preprocess --shape circle --hEdge 0.3 --hInnerEdge 0.3 --full-3d --mesh-folder /root/scratch/write-meshes/mechanotransduction/meshes/circle_hEdge_0.3_full3dtrue
# echo "Writing full rect mesh..."
# python3 main.py mechanotransduction-preprocess --shape rect --hEdge 0.3 --hInnerEdge 0.3 --full-3d --mesh-folder /root/scratch/write-meshes/mechanotransduction/meshes/rect_hEdge_0.3_full3dtrue
echo "Writing full star mesh..."
python3 main.py mechanotransduction-preprocess --shape star --hEdge 0.3 --hInnerEdge 0.3 --full-3d --mesh-folder /root/shared/write-meshes/mechanotransduction/meshes/star_hEdge_0.3_full3dtrue
echo "Done."