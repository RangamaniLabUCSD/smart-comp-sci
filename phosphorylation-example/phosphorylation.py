# # Example 3: Protein phosphorylation and diffusion in 3D cell geometry
#
# Here, we implement the model of protein phosphorylation at the cell membrane and diffusion in the cytosol from [Meyers, Craig, and Odde 2006, Current Biology](https://doi.org/10.1016/j.cub.2006.07.056).
#
# This model geometry consists of 2 domains - one surface and one volume:
# - plasma membrane (PM) - cell surface
# - cytosol - intracellular volume
#
# In this case, we only model the case of a spherical cell, where the cytosol corresponds to the interior of the sphere and the PM corresponds to the surface of the sphere.
#
# This model includes a single species, A, which is phosphorylated at the cell membrane. The unphosphorylated form of A ($A_{dephos}$) can be computed from mass conservation; everywhere $c_{A_{phos}} + c_{A_{dephos}} = c_{Tot}$, which is a constant in both time and space if the phosphorylated vs. unphosphorylated forms have the same diffusion coefficient.
#
# There are two reactions - one in the PM and other in the cytosol. At the membrane, $A_{dephos}$ is phosphorylated by a first-order reaction with rate $k_{kin}$, and in the cytosolic volume, $A_{phos}$ is dephosphorylated by a first order reaction with rate $k_p$. The resulting equations are:
#
# $$
# \frac{\partial{c_{A_{phos}}}}{\partial{t}} = D_{A_{phos}} \nabla ^2 c_{A_{phos}} - k_p c_{A_{phos}} \quad \text{in} \; \Omega_{Cyto}\\
# \text{B.C.:} \quad D_{A_{phos}}  (\textbf{n} \cdot \nabla c_{A_{phos}})  = k_{kin} c_{A_{dephos}} \quad \text{on} \; \Gamma_{PM}
# $$
#
# where we note that $c_{A_{dephos}} = c_{Tot} - c_{A_{phos}}$ in the boundary condition due to mass conservation.
#
# In this file, we test this model over multiple cell sizes and compare the results to analytical predictions. Please note that because we are testing several different geometries, this file may take an hour or more to complete execution.

# +
import dolfin as d
import numpy as np
import pathlib
import argparse
import logging

from smart import config, mesh, model, mesh_tools
from smart.units import unit

from smart import config, mesh, model, mesh_tools, visualization
from smart.model_assembly import (
    Compartment,
    Parameter,
    Reaction,
    Species,
    SpeciesContainer,
    ParameterContainer,
    CompartmentContainer,
    ReactionContainer,
)

from matplotlib import pyplot as plt
from phosphorylation_parser_args import add_phosphorylation_arguments


# -

# We will set the logging level to `INFO`. This will display some output during the simulation. If you want to get even more output you could set the logging level to `DEBUG`.

smart_logger = logging.getLogger("smart")
smart_logger.setLevel(logging.DEBUG)
logger = logging.getLogger("example3")
logger.setLevel(logging.INFO)
logger.info("Starting phosphorylation example")
d.parameters["form_compiler"]["quadrature_degree"] = 4 

# -

parser = argparse.ArgumentParser()
add_phosphorylation_arguments(parser)
args = vars(parser.parse_args())
timer = d.Timer("phosphorylation-example")

if not pathlib.Path.exists(args["mesh_folder"]):
    args["rect"] = True
    here = pathlib.Path(__file__).parent
    import sys
    sys.path.insert(0, (here / ".." / "scripts").as_posix())
    import runner as main_run
    args["mesh_folder"] = here / "mesh"
    main_run.preprocess_phosphorylation_mesh(
        mesh_folder = args["mesh_folder"],
        curRadius = 2.0,
        hEdge = 0.2,
        num_refinements = 0,
        axisymmetric = False,
        rect = args["rect"],
        dry_run = False)
    args["outdir"] = here / "results"

# Futhermore, you could also save the logs to a file by attaching a file handler to the logger as follows.
#
# ```
# file_handler = logging.FileHandler("filename.log")
# file_handler.setFormatter(logging.Formatter(smart.config.base_format))
# logger.addHandler(file_handler)
# ```

# Now, we define various units used in this problem.

uM = unit.uM
um = unit.um
molecule = unit.molecule
sec = unit.sec
dimensionless = unit.dimensionless
D_unit = um**2 / sec
flux_unit = molecule / (um**2 * sec)
vol_unit = uM

# ## Generate model
# Next we generate the model, which consists of four containers - compartments, species, reactions, and parameters.
#
# ### Compartments
# As described above, the two compartments are the cytosol ("Cyto") and the plasma membrane ("PM"). These are initialized by calling:
# ```
# compartment_var = Compartment(name, dimensionality, compartment_units, cell_marker)
# ```
# where
# - name: string naming the compartment
# - dimensionality: topological dimensionality (i.e. 3 for Cyto, 2 for PM)
# - compartment_units: length units for the compartment (um for both here)
# - cell_marker: integer marker value identifying each compartment in the parent mesh

if args["axisymmetric"]:
    parent_dim = 2
else:
    parent_dim = 3

Cyto = Compartment("Cyto", parent_dim, um, 1)
PM = Compartment("PM", parent_dim-1, um, 10)

# Create a compartment container.

cc = CompartmentContainer()
cc.add([PM, Cyto])

# ### Species
# In this case, we have a single species, "A", which exists in the cytosol. A single species is initialized by calling:
# ```
# species_var = Species(
#             name, initial_condition, concentration_units,
#             D, diffusion_units, compartment_name, group (opt)
#         )
# ```
# where
# - name: string naming the species
# - initial_condition: initial concentration for this species (can be an expression given by a string to be parsed by sympy - the only unknowns in the expression should be x, y, and z)
# - concentration_units: concentration units for this species (μM here)
# - D: diffusion coefficient
# - diffusion_units: units for diffusion coefficient (μm<sup>2</sup>/sec here)
# - compartment_name: each species should be assigned to a single compartment ("Cyto", here)
# - group (opt): for larger models, specifies a group of species this belongs to;
#             for organizational purposes when there are multiple reaction modules

Aphos = Species("Aphos", 0.1, vol_unit, args["diffusion"], D_unit, "Cyto")

# Create a species container.

sc = SpeciesContainer()
sc.add([Aphos])

# ### Parameters and Reactions
# Parameters and reactions are generally defined together, although the order does not strictly matter. Parameters are specified as:
# ```
# param_var = Parameter(name, value, unit, group (opt), notes (opt), use_preintegration (opt))
# ```
# where
# - name: string naming the parameter
# - value: value of the given parameter
# - unit: units associated with given value
# - group (optional): optional string placing this reaction in a reaction group; for organizational purposes when there are multiple reaction modules
# - notes (optional): string related to this parameter
# - use_preintegration (optional): in the case of a time-dependent parameter, uses preintegration in the solution process
#
# Reactions are specified by a variable number of arguments (arguments are indicated by (opt) are either never
# required or only required in some cases, for more details see notes below and API documentation):
# ```
# reaction_var = Reaction(
#                 name, lhs, rhs, param_map,
#                 eqn_f_str (opt), eqn_r_str (opt), reaction_type (opt), species_map,
#                 explicit_restriction_to_domain (opt), group (opt), flux_scaling (opt)
#             )
# ```
# - name: string naming the reaction
# - lhs: list of strings specifying the reactants for this reaction
# - rhs: list of strings specifying the products for this reaction
#     ***NOTE: the lists "reactants" and "products" determine the stoichiometry of the reaction;
#        for instance, if two A's react to give one B, the reactants list would be ["A","A"],
#        and the products list would be ["B"]
# - param_map: relationship between the parameters specified in the reaction string and those given
#               in the parameter container. By default, the reaction parameters are "kon" and "koff" when
#               a system obeys simple mass action. If the forward rate is given by a parameter "k1" and the
#               reverse rate is given by "k2", then param_map = {"on":"k1", "off":"k2"}
# - eqn_f_str: For systems not obeying simple mass action, this string specifies the forward reaction rate
#              By default, this string is "on*{all reactants multiplied together}"
# - eqn_r_str: For systems not obeying simple mass action, this string specifies the reverse reaction rate
#              By default, this string is "off*{all products multiplied together}"
# - reaction_type (opt): either "custom" or "mass_action" (default is "mass_action") [never a required argument]
# - species_map: same format as param_map; required if other species not listed in reactants or products appear in the
#             reaction string
# - explicit_restriction_to_domain: string specifying where the reaction occurs; required if the reaction is not
#                                   constrained by the reaction string (e.g., if production occurs only at the boundary,
#                                   as it does here, but the species being produced exists through the entire volume)
# - group (opt): string placing this reaction in a reaction group; for organizational purposes when there are multiple reaction modules
# - flux_scaling (opt): in certain cases, a given reactant or product may experience a scaled flux (for instance, if we assume that
#                 some of the molecules are immediately sequestered after the reaction); in this case, to signify that this flux
#                 should be rescaled, we specify ''flux_scaling = {scaled_species: scale_factor}'', where scaled_species is a
#                 string specifying the species to be scaled and scale_factor is a number specifying the rescaling factor

Atot = Parameter("Atot", 1.0, vol_unit)
# Phosphorylation of Adephos at the PM
kkin = Parameter("kkin", 50.0, 1 / sec)
curRadius = args["curRadius"]  # first radius value to test
# vol to surface area ratio of the cell (overwritten for each cell size)
if args["rect"]:
    VolSA = Parameter("VolSA", curRadius / 20, um)
else:
    VolSA = Parameter("VolSA", curRadius / 3, um)
r1 = Reaction(
    "r1",
    [],
    ["Aphos"],
    param_map={"kon": "kkin", "Atot": "Atot", "VolSA": "VolSA"},
    eqn_f_str="kon*VolSA*(Atot - Aphos)",
    species_map={"Aphos": "Aphos"},
    explicit_restriction_to_domain="PM",
)
# Dephosphorylation of Aphos in the cytosol
kp = Parameter("kp", 10.0, 1 / sec)
r2 = Reaction(
    "r2",
    ["Aphos"],
    [],
    param_map={"kon": "kp"},
    eqn_f_str="kp*Aphos",
    species_map={"Aphos": "Aphos"},
)

# Create parameter and reaction containers.

# +
pc = ParameterContainer()
pc.add([Atot, kkin, VolSA, kp])

rc = ReactionContainer()
rc.add([r1, r2])
# -

# ## Create/load in mesh
#
# In SMART we have different levels of meshes. Here, for our first mesh, we specify a sphere of radius 1.
#
# $$
# \Omega: r \in [0, 1] \subset \mathbb{R}^3\\
# \text{where} \qquad r = \sqrt{x^2 + y^2 + z^2}
# $$
#
# which will serve as our parent mesh, giving the overall cell geometry.
#
# Different domains can be specified within this parent mesh by assigning marker values to cells (3D) or facets (2D) within the mesh. A subdomain within the parent mesh, defined by a region which shares the same marker value, is referred to as a child mesh.
#
# Here, we have two child meshes corresponding to the 2 compartments specified in the compartment container. As defined above, "PM" is a 2D compartment defined by facets with marker value 10 and "Cyto" is a 3D compartment defined by cells with marker value 1. These subdomains are defined by:
# - $\Omega_{Cyto}: r \in [0, 1) \subset \mathbb{R}^3$
# - $\Gamma_{PM}: r=1 \subset \mathbb{R}^3$
#
# We generate the parent mesh with appropriate markers using gmsh in the function `mesh_tools.create_spheres`

# Load mesh
mesh_file = args["mesh_folder"] / "DemoSphere.h5"
parent_mesh = mesh.ParentMesh(
    str(mesh_file),
    mesh_filetype="hdf5",
    name="parent_mesh",
)

if args["comparison_results_folder"] != pathlib.Path(""):
    comm = parent_mesh.mpi_comm
    comparison_mesh = d.Mesh(comm)
    mesh_filename = f'{str(args["comparison_mesh_folder"])}/DemoSphere.h5'
    hdf5 = d.HDF5File(comparison_mesh.mpi_comm(), mesh_filename, "r")
    hdf5.read(comparison_mesh, "/mesh", False)
    cell_markers = d.MeshFunction("size_t", comparison_mesh, 2, value=0)
    hdf5.read(cell_markers, f"/mf2")
    hdf5.close()
    cyto_comparison_mesh = d.create_meshview(cell_markers, 1)
    Aphos_sol = d.Function(d.FunctionSpace(cyto_comparison_mesh, "P", 1))
    comparison_file = d.XDMFFile(f"{str(args['comparison_results_folder'])}/Aphos.xdmf")
    time_ref = np.loadtxt(f"{str(args['comparison_results_folder'])}/tvec.txt")

# ## Initialize model and solver
#
# Now we modify the solver configuration for this problem. In the solver config, we set the final t as 1 s, the initial dt at .01 s (without any additional specifications, this will be the time step for the whole simulation), and the time precision (number of digits after the decimal point to round to) as 6.

config_cur = config.Config()
model_cur = model.Model(pc, sc, cc, rc, config_cur, parent_mesh)
config_cur.flags.update(
    {
        "axisymmetric_model": args["axisymmetric"],
    }
)
config_cur.solver.update(
    {
        "final_t": 1.0,
        "initial_dt": args["time_step"],
        "time_precision": 6,
    }
)


# Now we initialize the model and solver.

model_cur.initialize()

# ## Solve model and store output
#
# We create XDMF files where we will store the output and store model information in a .pkl file.

# +
# Write initial condition(s) to file
results = dict()
result_folder = args["outdir"]
result_folder.mkdir(exist_ok=True, parents=True)
for species_name, species in model_cur.sc.items:
    results[species_name] = d.XDMFFile(
        model_cur.mpi_comm_world, str(result_folder / f"{species_name}.xdmf")
    )
    if args["write_checkpoint"]:
        results[species_name].write_checkpoint(model_cur.sc[species_name].u["u"], 
                                                   "u", model_cur.t, d.XDMFFile.Encoding.HDF5, False)
    else:
        results[species_name].write(model_cur.sc[species_name].u["u"], model_cur.t)


# model_cur.to_pickle(result_folder / "model_cur.pkl")

if model_cur.mpi_comm_world.rank == 0:
    import json

    # Dump config to results folder
    (result_folder / "config.json").write_text(
        json.dumps(
            {
                "solver": config_cur.solver.__dict__,
                "flags": config_cur.flags.__dict__,
                "reaction_database": config_cur.reaction_database,
                "mesh_file": str(args["mesh_folder"]),
                "outdir": str(args["outdir"]),
                "time_step": args["time_step"],
                "curRadius": args["curRadius"],
                "diffusion": args["diffusion"],
                "hmin": parent_mesh.dolfin_mesh.hmin(),
                "hmax": parent_mesh.dolfin_mesh.hmax()
            }
        )
    )

# -

# We now run the solver until t reaches final_t, recording the average Aphos concentration at each time point. Then plot the final concentration profile using pyvista.

# +
# Set loglevel to warning in order not to pollute notebook output
smart_logger.setLevel(logging.WARNING)

# save integration measure and volume for computing average Aphos at each time step
dx = d.Measure("dx", domain=model_cur.cc["Cyto"].dolfin_mesh)
if args["axisymmetric"]:
    x = d.SpatialCoordinate(model_cur.cc['Cyto'].dolfin_mesh)
    volume = d.assemble_mixed(1.0*x[0]*dx)
else:
    volume = d.assemble_mixed(1.0 * dx)
# Solve
avg_Aphos = [Aphos.initial_condition]
L2vec = [0]

while True:
    logger.info(f"Solve for time step {model_cur.t}")
    # Solve the system
    model_cur.monolithic_solve()
    # Save results for post processing
    for species_name, species in model_cur.sc.items:
        if args["write_checkpoint"]:
            results[species_name].write_checkpoint(model_cur.sc[species_name].u["u"], 
                                                   "u", model_cur.t, d.XDMFFile.Encoding.HDF5, True)
        else:
            results[species_name].write(model_cur.sc[species_name].u["u"], model_cur.t)
    # compute average Aphos concentration at each time step
    if args["axisymmetric"]:
        int_val = d.assemble_mixed(x[0]*model_cur.sc['Aphos'].u['u']*dx)
    else:
        int_val = d.assemble_mixed(model_cur.sc["Aphos"].u["u"] * dx)
    avg_Aphos.append(int_val / volume)

    if args["comparison_results_folder"] != pathlib.Path(""):
        test_idx = np.argmin(np.abs(time_ref - float(model_cur.t)))
        comparison_file.read_checkpoint(Aphos_sol, "u", test_idx)
        Aphos_sol.set_allow_extrapolation(True)
        Aphos_proj = d.Function(sc["Aphos"].V)
        vals = Aphos_proj.vector().get_local()
        coords = sc["Aphos"].V.tabulate_dof_coordinates()
        for i in range(len(coords)):
            vals[i] = Aphos_sol(coords[i])
        Aphos_proj.vector().set_local(vals)
        Aphos_proj.vector().apply("insert")
        L2vec.append(np.sqrt(d.assemble_mixed((Aphos_proj-model_cur.sc["Aphos"].u["u"])**2 *dx)))

    # End if we've passed the final time
    if model_cur.t >= model_cur.final_t:
        break
# visualization.plot(model_cur.sc["Aphos"].u["u"])
# -

# L2 error
xvec = d.SpatialCoordinate(cc["Cyto"].dolfin_mesh)
if args["axisymmetric"]:
    r = d.sqrt(xvec[0]**2 + xvec[1]**2 + (xvec[2]-(curRadius+1))**2)
else:
    r = d.sqrt(xvec[0]**2 + xvec[1]**2 + (xvec[2])**2)
k_kin = kkin.value
k_p = kp.value
cT = Atot.value
D = Aphos.D
if args["rect"]:
    k_kin = kkin.value * VolSA.value
    zSize = curRadius
    zFactor = zSize / np.sqrt(D/k_p)
    Az = (k_kin*zSize/(2*D))*np.exp(zFactor) / ((zFactor/2)*(np.exp(zFactor)-1) + (k_kin*zSize/(2*D))*(1+np.exp(zFactor)))
    Bz = (k_kin*zSize/(2*D)) / ((zFactor/2)*(np.exp(zFactor)-1) + (k_kin*zSize/(2*D))*(1+np.exp(zFactor)))
    expFactor = np.sqrt(k_p/(D))
    sol = cT * ((Az*d.exp(-expFactor*xvec[2]) + Bz*d.exp(expFactor*xvec[2])))
else:
    thieleMod = curRadius / np.sqrt(D/k_p)
    C1 = k_kin*cT*curRadius**2/((3*D*(np.sqrt(k_p/D)-(1/curRadius)) + k_kin*curRadius)*np.exp(thieleMod) +
                                (3*D*(np.sqrt(k_p/D)+(1/curRadius))-k_kin*curRadius)*np.exp(-thieleMod))
    sol = C1*(d.exp(r/np.sqrt(D/k_p))-d.exp(-r/np.sqrt(D/k_p)))/r
    
L2norm = np.sqrt(d.assemble_mixed((sol-model_cur.sc["Aphos"].u["u"])**2 *dx))


logger.info("Done with solve loop")
timer.stop()
timings = d.timings(
    d.TimingClear.keep,
    [d.TimingType.wall, d.TimingType.user, d.TimingType.system],
).str(True)
(result_folder / f"timings_rank{model_cur.mpi_comm_world.rank}.txt").write_text(timings)
print(timings)

if model_cur.mpi_comm_world.size > 1:
    d.MPI.comm_world.Barrier()

if model_cur.mpi_comm_world.rank == 0: 
    # We plot the average Aphos over time.
    fig, ax = plt.subplots()
    ax.plot(model_cur.tvec, avg_Aphos)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Aphos concentration (μM)")
    fig.savefig(result_folder / "avg_Aphos.png")
    np.savetxt(result_folder / "avg_Aphos.txt", avg_Aphos)
    np.savetxt(result_folder / "tvec.txt", model_cur.tvec)
    np.savetxt(result_folder / "L2norm.txt", np.array([L2norm]))
    np.savetxt(result_folder / "L2vec.txt", np.array(L2vec))
    (result_folder / "timings.txt").write_text(timings)