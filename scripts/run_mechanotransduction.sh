partition=defq

shape=circle
timestep=0.01
eval=70000000

for refinement in 0 1 2 3
do
    python3 main.py --submit-ex3 --partition=$partition mechanotransduction --mesh-folder meshes-mechanotransduction/${shape}_refined_${refinement} --time-step $timestep --e-val $eval --axisymmetric --well-mixed
    sleep 5
done


