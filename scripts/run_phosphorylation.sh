# Create meshes
# for refinement in 0 1 2 3 4 5 6 7 8
# do
#     python3 main.py phosphorylation-preprocess --num-refinements $refinement --rect --curRadius 2.0 --mesh-folder "meshes-phosphorylation/rect_refined_${refinement}"
# done

partition=defq

# Diffusion, Refinements, timestep
for timestep in .64 .32 .16 .08 .04 .02 .01
do
    for refinement in 0 1 2 3
    do
        for diffusion in 10.0 100.0 1000.0
        do
            python3 main.py --partition=$partition --submit-ex3 phosphorylation --time-step $timestep --rect --diffusion $diffusion --curRadius 2.0 --mesh-folder "meshes-phosphorylation/rect_refined_${refinement}"
            sleep 5
        done
    done
done