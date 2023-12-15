#!/usr/bin/env python
# coding: utf-8

# # ATP dynamics in realistic mitochondridal geometries
# 
# Here, we implement the model presented in [Garcia et al, Scientific Reports](https://www.biorxiv.org/content/10.1101/2022.08.16.500715v2.full#ref-39), which describes production of ATP in mitochondria.
# 
# The geometry in this model is divided into 4 domains - two volumes and two surfaces:
# - outer mitochondrial membrane (OM)
# - intermembrane space (IMS)
# - Inner membrane (IM)
# - Matrix (MAT) (volume inside the mitochondrial matrix)

# In[1]:


import dolfin as d
import sympy as sym
import numpy as np
import pathlib
import logging
import gmsh  # must be imported before pyvista if dolfin is imported first

from smart import config, mesh, model, mesh_tools, visualization
from smart.units import unit
from smart.model_assembly import (
    Compartment,
    Parameter,
    Reaction,
    Species,
    SpeciesContainer,
    ParameterContainer,
    CompartmentContainer,
    ReactionContainer,
    sbmodel_from_locals
)

from matplotlib import pyplot as plt
import matplotlib.image as mpimg
from matplotlib import rcParams

logger = logging.getLogger("smart")
logger.setLevel(logging.DEBUG)


# First, we define the various units for the inputs and define universal constants.

# In[2]:


# Aliases - base units
mM = unit.mM
um = unit.um
molecule = unit.molecule
sec = unit.sec
dimensionless = unit.dimensionless
# Aliases - units used in model
D_unit = um**2 / sec
flux_unit = mM * um / sec
vol_unit = mM
surf_unit = molecule / um**2

F = 9.649e4 # Faraday's constant (C/mol)
R = 8.315 # gas constant (J/mol-K)


# ## Model generation
# 
# We first define compartments and the add them to a compartment container. Note that we can specify nonadjacency for surfaces in the model, which is not required, but can speed up the solution process.

# In[3]:


IMS = Compartment("IMS", 3, um, 1)
OM = Compartment("OM", 2, um, 10)
Mat = Compartment("Mat", 3, um, 2)
IM = Compartment("IM", 2, um, 12)
OM.specify_nonadjacency(["IM", "Mat"])
IM.specify_nonadjacency(["OM"])

cc = CompartmentContainer()
cc.add([IMS, OM, Mat, IM])


# Define all model species and and place in species container.

# In[4]:


# ATP synthase states
E_tot = Parameter("E_tot", 267/1.545, surf_unit)
E_Mat = Species("E_Mat", 0.0, surf_unit, 0.0, D_unit, "IM")
E_IMS = Species("E_IMS", 1.0*E_tot.value, surf_unit, 0.0, D_unit, "IM")
E_Mat_H3Star = Species("E_Mat_H3Star", 0.0, surf_unit, 0.0, D_unit, "IM")
E_Mat_H3S = Species("E_Mat_H3S", 0.0, surf_unit, 0.0, D_unit, "IM")
E_Mat_H3 = Species("E_Mat_H3", 0.0, surf_unit, 0.0, D_unit, "IM")
# ANT states
L_tot = Parameter("L_tot", 16471/1.545, surf_unit)
L = Species("L", 1.0*L_tot.value, surf_unit, 0.0, D_unit, "IM")
TL = Species("TL", 0.0, surf_unit, 0.0, D_unit, "IM")
LTsp = Species("LTsp", 0.0, surf_unit, 0.0, D_unit, "IM")
DL = Species("DL", 0.0, surf_unit, 0.0, D_unit, "IM")
LD = Species("LD", 0.0, surf_unit, 0.0, D_unit, "IM")
TLD = Species("TLD", 0.0, surf_unit, 0.0, D_unit, "IM")
DLT = Species("DLT", 0.0, surf_unit, 0.0, D_unit, "IM")
DLD = Species("DLD", 0.0, surf_unit, 0.0, D_unit, "IM") # DLD + DLD' in original model
# DLDPrime = Species("DLDPrime", 0.0, surf_unit, 0.0, D_unit, "IM")
# TLT = Species("TLT", 0.0, surf_unit, 0.0, D_unit, "IM") # TLT + TLT' in original model
# TLTPrime = Species("TLTPrime", 0.0, surf_unit, 0.0, D_unit, "IM")
# ATP/ADP in matrix and IMS
D_Mat = Species("D_Mat", 0.8*2.0, vol_unit, 15.0, D_unit, "Mat")
T_Mat = Species("T_Mat", 13.0, vol_unit, 15.0, D_unit, "Mat")
T_IMS = Species("T_IMS", 6.5, vol_unit, 15.0, D_unit, "IMS")
D_IMS = Species("D_IMS", 0.1, vol_unit, 15.0, D_unit, "IMS")


# Define mitochondrial conditions (temperature, voltage, pH)

# In[5]:


# pH_Mat = 7.6 # matrix pH
# pH_c = 7.2 # cristae pH
T = 310 # body temperature (K)
dPsi = 180 # membrane voltage (mV)
dPsiB_IMS = -50 # phase boundary potential bulk IMS -> IM (mV) (I think this is reversed)
dPsiB_Mat = 0.0 # phase boundary potential from IM -> Mat (mV)
dPsi_m = dPsi + dPsiB_IMS - dPsiB_Mat
F_RT = F / (1000*R*T) # F/(RT) (1/1000)*C/K (1/mV)


# Now, we define parameters and reactions for each portion of the model. First, we define the ATP synthase dynamics:

# In[6]:


# E1: Movement of proton binding site in ATP synthase
k_16 = Parameter("k_16", 100.0*np.exp(3*F_RT*dPsi_m/2), 1/sec) # movement of proton binding site Mat->IMS side
k_61 = Parameter("k_61", 4.98e7*np.exp(-3*F_RT*dPsi_m/2), 1/sec) # movement of proton binding site IMS->Mat side
E1 = Reaction("E1", ["E_Mat"], ["E_IMS"], {"on": "k_16", "off": "k_61"})

# E2: bind/release of 3 protons in IMS
k_65 = Parameter("k_65", 3969, 1/sec) # binding rate of protons in IMS
k_56 = Parameter("k_56", 2.75e5*np.exp(3*F_RT*dPsiB_IMS), 1/sec) # release of protons in IMS
E2 = Reaction("E2", ["E_IMS"], [],
              param_map={"E_tot":"E_tot", "k_56":"k_56", "k_65":"k_65"},
              species_map={"E_Mat":"E_Mat", "E_IMS":"E_IMS", "E_Mat_H3Star":"E_Mat_H3Star", "E_Mat_H3S":"E_Mat_H3S", "E_Mat_H3":"E_Mat_H3"},
              eqn_f_str="k_65*E_IMS - k_56*(E_tot-E_Mat-E_IMS-E_Mat_H3Star-E_Mat_H3S-E_Mat_H3)")

# E3: movement of 3 protons from IMS to Matrix
k_54 = Parameter("k_54", 100.0, 1/sec)
k_45 = Parameter("k_45", 100.0, 1/sec)
E3 = Reaction("E3", [], ["E_Mat_H3Star"],
              param_map={"E_tot":"E_tot", "k_54":"k_54", "k_45":"k_45"},
              species_map={"E_Mat":"E_Mat", "E_IMS":"E_IMS", "E_Mat_H3Star":"E_Mat_H3Star", "E_Mat_H3S":"E_Mat_H3S", "E_Mat_H3":"E_Mat_H3"},
              eqn_f_str="k_54*(E_tot-E_Mat-E_IMS-E_Mat_H3Star-E_Mat_H3S-E_Mat_H3) - k_45*E_Mat_H3Star")

# E4: movement of 3 protons from IMS to matrix without producing ATP
k_52 = Parameter("k_52", 1e-20, 1/sec)
k_25 = Parameter("k_25", 5.85e-30, 1/sec)
E4 = Reaction("E4", [], ["E_Mat_H3"],
              param_map={"E_tot":"E_tot", "k_52":"k_52", "k_25":"k_25"},
              species_map={"E_Mat":"E_Mat", "E_IMS":"E_IMS", "E_Mat_H3Star":"E_Mat_H3Star", "E_Mat_H3S":"E_Mat_H3S", "E_Mat_H3":"E_Mat_H3"},
              eqn_f_str="k_52*(E_tot-E_Mat-E_IMS-E_Mat_H3Star-E_Mat_H3S-E_Mat_H3) - k_25*E_Mat_H3")

# E5: binding/unbinding of ADP-P to ATP synthase
k_43 = Parameter("k_43", 2e3, 1/(mM*sec))
k_34 = Parameter("k_34", 100.0, 1/sec)
E5 = Reaction("E5", ["E_Mat_H3Star", "D_Mat"], ["E_Mat_H3S"],
              {"on": "k_43", "off": "k_34"})

# E6: ATP production
k_32 = Parameter("k_32", 5e3, 1/sec)
k_23 = Parameter("k_23", 5e3, 1/(mM*sec))
E6 = Reaction("E6", ["E_Mat_H3S"], ["E_Mat_H3", "T_Mat"],
              {"on": "k_32", "off": "k_23"})

# E7: detachment of protons in the matrix
k_21 = Parameter("k_21", 40.0*np.exp(3*F_RT*dPsiB_Mat), 1/sec) # unbinding rate of protons in matrix
k_12 = Parameter("k_12", 25.0, 1/sec) # binding rate of protons in matrix
E7 = Reaction("E7", ["E_Mat_H3"], ["E_Mat"], {"on": "k_21", "off": "k_12"})


# Now, we define the reactions associated with adenine nucleotide transporters (ANTs):

# In[7]:


# L1: Define all kinetic constants for ANTs
koff_Tm = Parameter("koff_Tm", 4e4, 1/sec) # unbinding of matrix ATP
kon_Tm = Parameter("kon_Tm", 6.4e6/1000, 1/(mM*sec)) # binding of matrix ATP
koff_Ti = Parameter("koff_Ti", 200.0, 1/sec) # unbinding of IMS ATP
kon_Ti = Parameter("kon_Ti", 4e5/1000, 1/(mM*sec)) # binding of IMS ATP
koff_Dm = Parameter("koff_Dm", 4e4, 1/sec) # unbinding of matrix ADP
kon_Dm = Parameter("kon_Dm", 4e6/1000, 1/(mM*sec)) # binding of matrix ADP
koff_Di = Parameter("koff_Di", 100.0, 1/sec) # unbinding of IMS ADP
kon_Di = Parameter("kon_Di", 4e6/1000, 1/(mM*sec)) # binding of IMS ADP
k_p = Parameter("k_p", 9.2, 1/sec) # ATP transport Mat->IMS (productive)
k_cp = Parameter("k_cp", 0.35, 1/sec) # ATP transport IMS -> Mat (nonproductive)

# L1: Binding of matrix ATP to free ANT (L)
L1 = Reaction("L1", ["L", "T_Mat"], ["LTsp"], {"on": "kon_Tm", "off": "koff_Tm"})

# L2: Binding of matrix ADP to free ANT (L)
L2 = Reaction("L2", ["L", "D_Mat"], ["LD"], {"on": "kon_Dm", "off": "koff_Dm"})

# L3: Binding of IMS ATP to free ANT (L)
L3 = Reaction("L3", ["L", "T_IMS"], ["TL"], {"on": "kon_Ti", "off": "koff_Ti"})

# L4: Binding of IMS ADP to free ANT (L)
L4 = Reaction("L4", ["L", "D_IMS"], ["DL"], {"on": "kon_Di", "off": "koff_Di"})

# L5: Binding of matrix ATP to TL
L5 = Reaction("L5", ["TL", "T_Mat"], [],
              param_map={"L_tot":"L_tot", "kon_Tm":"kon_Tm", "koff_Tm":"koff_Tm"},
              species_map={"L":"L", "LTsp":"LTsp", "LD":"LD", "TL":"TL", "DL":"DL", "TLD":"TLD", "DLT":"DLT", "DLD":"DLD", "T_Mat":"T_Mat"},
              eqn_f_str="kon_Tm*TL*T_Mat - koff_Tm*(L_tot-L-LTsp-LD-TL-DL-TLD-DLT-DLD)")

# L6: Binding of matrix ADP to TL
L6 = Reaction("L6", ["TL", "D_Mat"], ["TLD"], {"on": "kon_Dm", "off": "koff_Dm"})

# L7: Binding of IMS ATP to LT
L7 = Reaction("L7", ["LTsp", "T_IMS"], [],
              param_map={"L_tot":"L_tot", "kon_Ti":"kon_Ti", "koff_Tm":"koff_Ti"},
              species_map={"L":"L", "LTsp":"LTsp", "LD":"LD", "TL":"TL", "DL":"DL", "TLD":"TLD", "DLT":"DLT", "DLD":"DLD", "T_IMS":"T_IMS"},
              eqn_f_str="kon_Ti*LTsp*T_IMS - koff_Ti*(L_tot-L-LTsp-LD-TL-DL-TLD-DLT-DLD)")

# L8: Binding of IMS ADP to LT
L8 = Reaction("L8", ["LTsp", "D_IMS"], ["DLT"], {"on": "kon_Di", "off": "koff_Di"})

# L9: Binding of matrix ATP to DL
L9 = Reaction("L9", ["DL", "T_Mat"], ["DLT"], {"on": "kon_Tm", "off": "koff_Tm"})

# L10: Binding of matrix ADP to DL
L10 = Reaction("L10", ["DL", "D_Mat"], ["DLD"], {"on": "kon_Dm", "off": "koff_Dm"})

# L11: Binding of IMS ATP to LD
L11 = Reaction("L11", ["LD", "T_IMS"], ["TLD"], {"on": "kon_Ti", "off": "koff_Ti"})

# L12: Binding of IMS ADP to LD
L12 = Reaction("L12", ["LD", "D_IMS"], ["DLD"], {"on": "kon_Di", "off": "koff_Di"})

# L13: Exchange of ADP for ATP (can go either way)
L13 = Reaction("L13", ["DLT"], ["TLD"], {"on": "k_p", "off": "k_cp"})


# Finally, consider the export of ATP into the cytosol:

# In[8]:


k_vdac = Parameter("k_vdac", 1e6/1000, 1/(mM*sec))
VDAC = Parameter("VDAC", 1e4, surf_unit)
T_cyto = Parameter("T_cyto", 6.5, mM)
V1 = Reaction("V1", ["T_IMS"],[],
              explicit_restriction_to_domain="OM",
              param_map={"k_vdac":"k_vdac", "VDAC":"VDAC", "T_cyto":"T_cyto"},
              species_map={"T_IMS":"T_IMS"},
              eqn_f_str="k_vdac * VDAC * (T_cyto - T_IMS)")


# ## Create and load in mesh
# 
# Here, we consider an "ellipsoid-in-an-ellipsoid" geometry. The inner ellipsoid represents the ER and the volume between the ER boundary and the boundary of the outer ellipsoid represents the cytosol.

# In[9]:


cur_dir = pathlib.Path.cwd()
parent_dir = cur_dir.parent
mito_mesh = d.Mesh(f"{str(parent_dir)}/meshes/mito1_mesh.xml")
cell_markers = d.MeshFunction("size_t", mito_mesh, 3, mito_mesh.domains())
facet_markers_orig = d.MeshFunction("size_t", mito_mesh, 2, mito_mesh.domains())
facet_markers = d.MeshFunction("size_t", mito_mesh, 2, mito_mesh.domains())
facet_array = facet_markers.array()[:]
for i in range(len(facet_array)):
    if facet_array[i] == 11: # this indicates cristae
        facet_array[i] = 12 # consider cristae the same as inner membrane for now
    elif facet_array[i] > 1e9: # unassigned
        facet_array[i] = 0
# facet_markers.array()[np.where(facet_markers.array() > 1e9)[0]] = 0 # set unassigned facets to 0

# Write mesh and meshfunctions to file
mesh_folder = pathlib.Path("mesh")
mesh_folder.mkdir(exist_ok=True)
mesh_path = mesh_folder / "mito1.h5"
mesh_tools.write_mesh(
    mito_mesh, facet_markers, cell_markers, filename=mesh_path
)
parent_mesh = mesh.ParentMesh(
    mesh_filename=str(mesh_path),
    mesh_filetype="hdf5",
    name="parent_mesh",
)
visualization.plot_dolfin_mesh(mito_mesh, cell_markers, facet_markers, clip_logic=False)


# Initialize model and solver.

# In[10]:


pc, sc, cc, rc = sbmodel_from_locals(locals().values())
config_cur = config.Config()
config_cur.flags.update({"allow_unused_components": True})
model_cur = model.Model(pc, sc, cc, rc, config_cur, parent_mesh)
dt = .001
config_cur.solver.update(
    {
        "final_t": 0.1,
        "initial_dt": dt,
        "time_precision": 6,
    }
)
model_cur.initialize()


# Initialize XDMF files for saving results, save model information to .pkl file, then solve the system until `model_cur.t > model_cur.final_t`

# In[11]:


# Write initial condition(s) to file
results = dict()
result_folder = pathlib.Path(f"results")
result_folder.mkdir(exist_ok=True)
for species_name, species in model_cur.sc.items:
    results[species_name] = d.XDMFFile(
        model_cur.mpi_comm_world, str(result_folder / f"{species_name}.xdmf")
    )
    results[species_name].parameters["flush_output"] = True
    results[species_name].write(model_cur.sc[species_name].u["u"], model_cur.t)
model_cur.to_pickle("model_cur.pkl")

# Set loglevel to warning in order not to pollute notebook output
# logger.setLevel(logging.WARNING)

concVec = np.array([sc["T_IMS"].initial_condition])
IMSMesh = model_cur.cc["IMS"].dolfin_mesh
dx_IMS = d.Measure("dx", domain=IMSMesh)
OMMesh = model_cur.cc["OM"].dolfin_mesh
dx_OM = d.Measure("dx", domain=OMMesh)
volume_IMS = d.assemble(1.0*dx_IMS)
sa_OM = d.assemble(1.0*dx_OM)
volume_cyto = 0.306
cyto_convert = 1/(6.0221408e5 * volume_cyto)
# Solve
displayed = False
while True:
    T_cyto_prev = pc["T_cyto"].value
    T_cyto_flux = d.assemble(k_vdac.value * VDAC.value * (T_cyto_prev - sc["T_IMS"].u["u"]) * dx_OM)
    pc["T_cyto"].value = T_cyto_prev + T_cyto_flux * dt * cyto_convert
    # Solve the system
    model_cur.monolithic_solve()
    # update estimate for T_cyto
    T_cyto_flux = d.assemble(k_vdac.value * VDAC.value * (pc["T_cyto"].value - sc["T_IMS"].u["u"]) * dx_OM)
    pc["T_cyto"].value = T_cyto_prev + T_cyto_flux * dt * cyto_convert
    pc["T_cyto"].value_vector = np.vstack((pc["T_cyto"].value_vector, [float(model_cur.t), pc["T_cyto"].value]))
    pc["T_cyto"].dolfin_constant.assign(pc["T_cyto"].value)
    # Save results for post processing
    for species_name, species in model_cur.sc.items:
        results[species_name].write(model_cur.sc[species_name].u["u"], model_cur.t)
    int_val = d.assemble(model_cur.sc["T_IMS"].u["u"]*dx_IMS)
    curConc = np.array([int_val / volume_IMS])
    concVec = np.concatenate((concVec, curConc))
    np.savetxt(result_folder / f"tvec.txt", np.array(model_cur.tvec).astype(np.float32))
    np.savetxt(result_folder / f"T_cyto.txt", np.array(pc["T_cyto"].value_vector).astype(np.float32))
    # if model_cur.t > .025 and not displayed:  # display first time after .025 s
    #     visualization.plot(model_cur.sc['Ca'].u['u'])
    #     displayed = True
    # End if we've passed the final time
    if model_cur.t >= model_cur.final_t:
        break


# Plot concentration over time.

# In[ ]:


concVec[0] = sc["T_Mat"].initial_condition
plt.plot(model_cur.tvec, concVec)
plt.xlabel("Time (s)")
plt.ylabel("ATP concentration (mM)")
plt.title("SMART simulation")

