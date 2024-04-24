# Create meshes
for refinement in 0 1 2 3 4 5 6 7 8
do
    python3 main.py phosphorylation-preprocess --num-refinements $refinement --rect --curRadius 2.0 --mesh-folder "meshes-phosphorylation/rect_refined_${refinement}"
done

# # Diffusion, Refinements, timestep
# for timestep in .64 .32 .16 .08 .04 .02 .01 #.02 .04 .08 .16 .32 .64
# do
#     for refinement in 0 1 2 3 4 5 6 7 8
#     do
#         for diffusion in 0.1 10.0 100.0
#         do
#             python3 main.py phosphorylation --time-step $timestep --rect --diffusion $diffusion --curRadius 2.0 --mesh-folder "meshes-phosphorylation/rect_refined_${refinement}" --outdir "results-phosphorylation2/rect_refined_${refinement}_D_${diffusion}_dt_${timestep}"
#         done
#     done
# done