# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.15.2
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Mechanotransduction model in a mammalian cell
#
# Here, we implement the model originally presented in [Sun et al. 2016, Biophysical Journal], then expanded by [Scott et al. 2021, PNAS] and [Francis and Rangamani 2023, BioRxiv]. The goal is to describe the nuclear translocation of the mechanosensitive transcriptional regulators YAP/TAZ and MRTF, given some mechanical input from the substrate.
#
# The geometry in this model is divided into 4 domains - two volumes and two surfaces:
# - plasma membrane (PM)
# - Cytosol (Cyto)
# - Nuclear membrane (NM)
# - Interior of the nucleus (Nuc)
#
# This model has 19 species after accounting for mass conservation where possible.
# ```

# +
import argparse
import dolfin as d
import sympy as sym
import numpy as np
import pathlib
import logging
import gmsh # must be imported before pyvista if dolfin is imported first

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
# -

# Aliases - base units
uM = unit.uM
um = unit.um
molecule = unit.molecule
sec = unit.sec
dimensionless = unit.dimensionless
# Aliases - units used in model
D_unit = um**2 / sec
flux_unit = uM * um / sec
vol_unit = uM
surf_unit = molecule / um**2
# stiffness units
kPa = unit.kilopascal

# +
# cell_rad = 5.0
# nuc_rad = 0.5*cell_rad
outerExpr = "(1 - z**4/(1000+z**4)) * (r**2 + z**2) + 0.4*(r**2 + (z+15)**2)*z**4 / (10 + z**4) - 225"
innerExpr = "(r/5.3)**2 + ((z-7.8/2)/2.4)**2 - 1"

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--hEdge",
                    dest="hEdge", default=0.3, type=float,
                    help="Maximum mesh size for coarsest mesh at outer edge")
parser.add_argument("--hInnerEdge",
                    dest="hInnerEdge", default=0.5, type=float,
                    help="Maximum mesh size for coarsest mesh at inner edge")
parser.add_argument("--num_refinements",
                    dest="num-ref", default=0, type=int,
                    help="Number of mesh refinements",choices=[0,1,2,3])

args = parser.parse_args()
# cell_mesh, facet_markers, cell_markers = mesh_tools.create_axisymm(outerExpr = outerExpr, innerExpr = innerExpr,
#                                                                    hEdge=0.6, hInnerEdge=0.6)
cell_mesh, facet_markers, cell_markers = mesh_tools.create_2Dcell(outerExpr = outerExpr, innerExpr = innerExpr,
                                                                   hEdge=args.hEdge, hInnerEdge=args.hInnerEdge)
from IPython import embed;embed()

exit()

# cell_mesh, facet_markers, cell_markers = mesh_tools.create_spheres(outerRad = cell_rad, innerRad = nuc_rad,
#                                                                    hEdge=0.6, hInnerEdge=0.6)
mesh_folder = pathlib.Path("spreadCell_mesh")
mesh_folder.mkdir(exist_ok=True)
mesh_file = mesh_folder / "spreadCell_mesh.h5"
mesh_tools.write_mesh(cell_mesh, facet_markers, cell_markers, mesh_file)

parent_mesh = mesh.ParentMesh(
    mesh_filename=str(mesh_file),
    mesh_filetype="hdf5",
    name="parent_mesh",
)

dx = d.Measure("dx", domain=cell_mesh, subdomain_data=cell_markers)
ds = d.Measure("ds", domain=cell_mesh, subdomain_data=facet_markers)
x = d.SpatialCoordinate(cell_mesh)
# vol_cyto_orig = (4/3) * np.pi * (cell_rad**3 - nuc_rad**3)
# vol_nuc_orig = (4/3) * np.pi * nuc_rad**3
# sa_pm_orig = 4 * np.pi * cell_rad**2
vol_cyto_cur = d.assemble(1.0*x[0]*dx(1)) # actual cell volume (not used currently)
vol_cyto = 2300 # used to compute pFAK and RhoA fluxes
# vol_nuc = d.assemble(1.0*x[0]*dx(2))
sa_pm_cur = d.assemble(1.0*x[0]*ds(10)) # actual pm surface area (not used currently)
sa_pm = 1260 # used to compute pFAK and RhoA fluxes
# -

Cyto = Compartment("Cyto", 2, um, 1)
PM = Compartment("PM", 1, um, 10)
Nuc = Compartment("Nuc", 2, um, 2)
NM = Compartment("NM", 1, um, 12)
PM.specify_nonadjacency(['NM', 'Nuc'])
NM.specify_nonadjacency(['PM'])

E = Species("E", 10.0, dimensionless, 0.0, D_unit, "PM") # variable for substrate stiffness stimulation at the boundary
pFAK = Species("pFAK", 0.3, vol_unit, 10.0, D_unit, "Cyto")
RhoA_GDP = Species("RhoA_GDP", 1.0, vol_unit, 1.0, D_unit, "Cyto")
RhoA_GTP = Species("RhoA_GTP", 33.6, surf_unit, 0.3, D_unit, "PM")
ROCK_A = Species("ROCK_A", 0.0, vol_unit, 75.0, D_unit, "Cyto")
mDia_A = Species("mDia_A", 0.0, vol_unit, 1.0, D_unit, "Cyto")
Myo_A = Species("Myo_A", 1.5, vol_unit, 0.8, D_unit, "Cyto")
LIMK_A = Species("LIMK_A", 0.1, vol_unit, 10.0, D_unit, "Cyto")
Cofilin_NP = Species("Cofilin_NP", 1.8, vol_unit, 10.0, D_unit, "Cyto")
FActin = Species("FActin", 17.9, vol_unit, 0.6, D_unit, "Cyto")
GActin = Species("GActin", 482.4, vol_unit, 13.37, D_unit, "Cyto")
LaminA = Species("LaminA", 0.0, surf_unit, 0.001, D_unit, "NM")
NPC_A = Species("NPC_A", 0.0, surf_unit, 0.001, D_unit, "NM")
YAPTAZ = Species("YAPTAZ", 0.7, vol_unit, 19.0, D_unit, "Cyto") #non-phosphorylated in the cytosol
YAPTAZ_phos = Species("YAPTAZ_phos", 0.2, vol_unit, 19.0, D_unit, "Cyto") #phosphorylated in the cytosol
YAPTAZ_nuc = Species("YAPTAZ_nuc", 0.7, vol_unit, 19.0, D_unit, "Nuc")
MRTF = Species("MRTF", 0.0, vol_unit, 5.0, D_unit, "Cyto") #free MRTF in the cytosol, very rough estimate for diffusion
MRTF_bound = Species("MRTF_bound", 0.65, vol_unit, 5.0, D_unit, "Cyto") #free MRTF in the cytosol, very rough estimate for diffusion
MRTF_nuc = Species("MRTF_nuc", 0.301, vol_unit, 5.0, D_unit, "Nuc")

# +
# now define reactions and parameters across 3 modules
# module A: Substrate stiffness -> initial signaling events
# a1: FAK phosphorylation at the membrane
FAK_tot = Parameter("FAK_tot", 1.0, vol_unit)
k_f = Parameter("k_f", 0.015, 1/sec)
k_sf = Parameter("k_sf", 0.379, 1/sec)
C = Parameter("C", 3.25, dimensionless)
cytoConvert = Parameter("cytoConvert", vol_cyto/sa_pm, um)
a1 = Reaction("a1", [], ["pFAK"], 
              param_map={"FAK_tot":"FAK_tot", "k_f":"k_f", "k_sf":"k_sf", "C":"C", "cytoConvert":"cytoConvert"}, 
              species_map={"E":"E", "pFAK":"pFAK"},
              explicit_restriction_to_domain="PM", 
              eqn_f_str="cytoConvert*(FAK_tot-pFAK)*(k_f + k_sf*E/(E+C))")

# a2: FAK dephosphorylation throughout the cytosol
k_df = Parameter("k_df", 0.035, 1/sec)
a2 = Reaction("a2", ["pFAK"], [], 
              param_map={"k_df":"k_df"}, 
              species_map={"pFAK":"pFAK"}, 
              eqn_f_str="k_df*pFAK")

# a3: RhoA activation and deactivation
k_fkrho = Parameter("k_fkrho", 0.0168, 1/sec)
gammaConst = Parameter("gammaConst", 77.56, uM**(-5))
# n = Parameter("n", 5, dimensionless)
k_drho = Parameter("k_drho", 0.625, 1/sec)
a3 = Reaction("a3", ["RhoA_GDP"], ["RhoA_GTP"], 
              param_map={"k_fkrho":"k_fkrho", "gammaConst":"gammaConst", "k_drho":"k_drho", 
                         "cytoConvert":"cytoConvert"}, 
              species_map={"pFAK":"pFAK", "RhoA_GDP":"RhoA_GDP", "RhoA_GTP":"RhoA_GTP"}, 
              eqn_f_str="RhoA_GDP*cytoConvert*k_fkrho*(gammaConst*pFAK**5 + 1)",
              eqn_r_str="k_drho*RhoA_GTP")

#a4: ROCK activation by RhoA_GTP at PM
k_rrho = Parameter("k_rrho", 0.648, 1/(sec*uM))
ROCK_tot = Parameter("ROCK_tot", 1.0, uM)
a4 = Reaction("a4", [], ["ROCK_A"], 
              param_map={"k_rrho":"k_rrho", "ROCK_tot":"ROCK_tot"}, 
              species_map={"ROCK_A":"ROCK_A", "RhoA_GTP":"RhoA_GTP"}, 
              eqn_f_str="k_rrho*RhoA_GTP*(ROCK_tot-ROCK_A)")

#a5: Deactivation of ROCK in the cytosol
k_drock = Parameter("k_drock", 0.8, 1/sec)
a5 = Reaction("a5", ["ROCK_A"], [], 
              param_map={"k_drock":"k_drock"}, 
              species_map={"ROCK_A":"ROCK_A"}, 
              eqn_f_str="k_drock*ROCK_A")

# +
#Module B: cytoskeletal signaling
#b1: mDia activation by RhoA_GTP at PM
k_mrho = Parameter("k_mrho", 0.002, 1/(sec*uM))
mDia_tot = Parameter("mDia_tot", 0.8, uM)
b1 = Reaction("b1", [], ["mDia_A"], 
              param_map={"k_mrho":"k_mrho", "mDia_tot":"mDia_tot"}, 
              species_map={"mDia_A":"mDia_A", "RhoA_GTP":"RhoA_GTP"}, 
              eqn_f_str="k_mrho*RhoA_GTP*(mDia_tot-mDia_A)")

#b2: Deactivation of mDia in the cytosol
k_dmdia = Parameter("k_dmdia", 0.005, 1/sec)
b2 = Reaction("b2", ["mDia_A"], [], 
              param_map={"k_dmdia":"k_dmdia"}, 
              species_map={"mDia_A":"mDia_A"}, 
              eqn_f_str="k_dmdia*mDia_A")

#b3: activation and deactivation of myosin
Myo_tot = Parameter("Myo_tot", 5.0, vol_unit)
k_mr = Parameter("k_mr", 0.03, 1/sec)
ROCK_B = Parameter("ROCK_B", 0.3, vol_unit)
epsilon = Parameter("epsilon", 36.0, 1/uM)
sc1 = Parameter("sc1", 20.0, 1/uM)
k_dmy = Parameter("k_dmy", 0.067, 1/sec)
b3 = Reaction("b3", [], ["Myo_A"], 
              param_map={"Myo_tot":"Myo_tot", "k_mr":"k_mr", "ROCK_B":"ROCK_B", "epsilon":"epsilon", "sc1":"sc1", "k_dmy":"k_dmy"}, 
              species_map={"Myo_A":"Myo_A", "ROCK_A":"ROCK_A"}, 
              eqn_f_str="(Myo_tot-Myo_A)*k_mr*(1 + epsilon*(ROCK_A/2)*(tanh(sc1*(ROCK_A-ROCK_B)) + 1))\
                         - k_dmy*Myo_A")

#b4: activation and deactivation of LIMK
LIMK_tot = Parameter("LIMK_tot", 2.0, vol_unit)
k_lr = Parameter("k_lr", 0.07, 1/sec)
tau = Parameter("tau", 55.49, 1/uM)
k_dl = Parameter("k_dl", 2.0, 1/sec)
b4 = Reaction("b4", [], ["LIMK_A"], 
              param_map={"LIMK_tot":"LIMK_tot", "k_lr":"k_lr", "ROCK_B":"ROCK_B", "tau":"tau", "sc1":"sc1", "k_dl":"k_dl"}, 
              species_map={"LIMK_A":"LIMK_A", "ROCK_A":"ROCK_A"}, 
              eqn_f_str="(LIMK_tot-LIMK_A)*k_lr*(1 + tau*(ROCK_A/2)*(tanh(sc1*(ROCK_A-ROCK_B)) + 1))\
                         - k_dl*LIMK_A")

#b5: dephos. and phos. of Cofilin
Cofilin_tot = Parameter("Cofilin_tot", 2.0, vol_unit)
k_turnover = Parameter("k_turnover", 0.04, 1/sec)
k_catCof = Parameter("k_catCof", 0.34, 1/sec)
k_mCof = Parameter("k_mCof", 4.0, uM)
b5 = Reaction("b5", [], ["Cofilin_NP"], 
              param_map={"Cofilin_tot":"Cofilin_tot", "k_turnover":"k_turnover", "k_catCof":"k_catCof", "k_mCof":"k_mCof"}, 
              species_map={"LIMK_A":"LIMK_A", "Cofilin_NP":"Cofilin_NP"}, 
              eqn_f_str="(Cofilin_tot-Cofilin_NP)*k_turnover - k_catCof*LIMK_A*Cofilin_NP/(k_mCof + Cofilin_NP)")

#b6: actin polymerization and depolymerization
# Actin_tot = Parameter("Actin_tot", 500.0, uM)
k_ra = Parameter("k_ra", 0.4, 1/sec)
alpha = Parameter("alpha", 50.0, 1/uM)
mDia_B = Parameter("mDia_B", 0.165, uM)
k_dep = Parameter("k_dep", 3.5, 1/sec)
k_fc1 = Parameter("k_fc1", 4.0, 1/(uM*sec))
b6 = Reaction("b6", ["GActin"], ["FActin"], 
              param_map={"k_ra":"k_ra", "alpha":"alpha", "sc1":"sc1", "mDia_B":"mDia_B", "k_dep":"k_dep", "k_fc1":"k_fc1"}, 
              species_map={"FActin":"FActin", "GActin":"GActin", "mDia_A":"mDia_A", "Cofilin_NP":"Cofilin_NP"}, 
              eqn_f_str="GActin*k_ra*(1 + alpha*(mDia_A/2)*(tanh(sc1*(mDia_A-mDia_B)) + 1))\
                         - (k_dep + k_fc1*Cofilin_NP)*FActin")

# +
#Module C: nucleo-cytoplasmic transport
#c1: YAP/TAZ dephos. and phos. in the cytosol
# YAPTAZ_tot = Parameter("YAPTAZ_tot", 1385000, molecule)
k_CN = Parameter("k_CN", 0.56, 1/sec)
k_CY = Parameter("k_CY", 0.00076, 1/(sec*uM**2))
k_NC = Parameter("k_NC", 0.14, 1/sec)
# concToNum_cyto = Parameter("concToNum_cyto", 602.214*vol_cyto, molecule/uM)
# concToNum_nuc = Parameter("concToNum_nuc", 602.214*vol_nuc, molecule/uM)
c1 = Reaction("c1", ["YAPTAZ_phos"], ["YAPTAZ"], 
              param_map={"k_CN":"k_CN", "k_CY":"k_CY", "k_NC":"k_NC"}, 
              species_map={"YAPTAZ":"YAPTAZ", "YAPTAZ_phos":"YAPTAZ_phos", "FActin":"FActin", "Myo_A":"Myo_A"}, 
              eqn_f_str="YAPTAZ_phos*(k_CN + k_CY*FActin*Myo_A) - k_NC*YAPTAZ")

#c2: MRTF binding/unbinding from G-actin in the cytosol
# MRTF_tot = Parameter("YAPTAZ_tot", 1000000, molecule)
k_offMRTF = Parameter("k_offMRTF", 1.0, 1/sec)
k_onMRTF = Parameter("k_onMRTF", 2.3e-5, 1/(sec*uM**2))
c2 = Reaction("c2", ["MRTF_bound"], ["MRTF"], 
              param_map={"k_offMRTF":"k_offMRTF", "k_onMRTF":"k_onMRTF"}, 
              species_map={"MRTF_bound":"MRTF_bound", "GActin":"GActin", "MRTF":"MRTF"}, 
              eqn_f_str="MRTF_bound*k_offMRTF - k_onMRTF*GActin**2 * MRTF")

#c3: LaminA dephos. and phos.
LaminA_tot = Parameter("LaminA_tot", 3500.0, surf_unit)
k_fl = Parameter("k_fl", 0.46, 1/sec)
p = Parameter("p", 9.0e-6, kPa/uM**2.6)
C_LaminA = Parameter("C_LaminA", 100.0, kPa)
k_rl = Parameter("k_rl", 0.001, 1/sec)
c3 = Reaction("c3", [], ["LaminA"], 
              param_map={"LaminA_tot":"LaminA_tot", "k_fl":"k_fl", "p":"p", "C_LaminA":"C_LaminA", "k_rl":"k_rl"}, 
              species_map={"LaminA":"LaminA", "FActin":"FActin"}, 
              eqn_f_str="(LaminA_tot - LaminA)*k_fl*p*FActin**2.6/(C_LaminA + p*FActin**2.6) - k_rl*LaminA")

#c4: opening and closing of NPCs
NPC_tot = Parameter("NPC_tot", 6.5, surf_unit)
k_fNPC = Parameter("k_fNPC", 2.8e-7, 1/(sec*uM**2 *surf_unit))
k_rNPC = Parameter("k_rNPC", 8.7, 1/sec)
c4 = Reaction("c4", [], ["NPC_A"], 
              param_map={"NPC_tot":"NPC_tot", "k_fNPC":"k_fNPC", "k_rNPC":"k_rNPC"}, 
              species_map={"NPC_A":"NPC_A", "LaminA":"LaminA", "FActin":"FActin", "Myo_A":"Myo_A"}, 
              eqn_f_str="(NPC_tot - NPC_A)*k_fNPC*LaminA*FActin*Myo_A - k_rNPC*NPC_A")

#c5: nuclear translocation of YAP/TAZ
k_insolo = Parameter("k_insolo", 1.0, surf_unit/(sec*uM))
k_in2 = Parameter("k_in2", 10.0, 1/(sec*uM))
k_out = Parameter("k_out", 1.0, surf_unit/(sec*uM))
c5 = Reaction("c5", ["YAPTAZ"], ["YAPTAZ_nuc"], 
              param_map={"k_insolo":"k_insolo", "k_in2":"k_in2", "k_out":"k_out"}, 
              species_map={"YAPTAZ":"YAPTAZ", "YAPTAZ_nuc":"YAPTAZ_nuc", "NPC_A":"NPC_A"}, 
              eqn_f_str="YAPTAZ*(k_insolo + k_in2*NPC_A) - k_out*YAPTAZ_nuc")

#c6: nuclear translocation of MRTF
k_insoloMRTF = Parameter("k_insoloMRTF", 2.708, surf_unit/(sec*uM))
k_in2MRTF = Parameter("k_in2MRTF", 8.5, 1/(sec*uM))
k_outMRTF = Parameter("k_outMRTF", 1.0, surf_unit/(sec*uM))
c6 = Reaction("c6", ["MRTF"], ["MRTF_nuc"], 
              param_map={"k_insoloMRTF":"k_insoloMRTF", "k_in2MRTF":"k_in2MRTF", "k_outMRTF":"k_outMRTF"}, 
              species_map={"MRTF":"MRTF", "MRTF_nuc":"MRTF_nuc", "NPC_A":"NPC_A"}, 
              eqn_f_str="MRTF*(k_insoloMRTF + k_in2MRTF*NPC_A) - k_outMRTF*MRTF_nuc")

# +

pc, sc, cc, rc = sbmodel_from_locals(locals().values())
# -

configCur = config.Config()
configCur.flags.update({
    "allow_unused_components": True,
    "axisymmetric_model": True})
model_cur = model.Model(pc, sc, cc, rc, configCur, parent_mesh)
configCur.solver.update(
    {
        "final_t": 10000.0,
        "initial_dt": 1.0,
        "time_precision": 6,
        "use_snes": True
    }
)

EVals = [0.1, 5.7, 7e7]
for i in range(len(EVals)):
    sc['E'].initial_condition = EVals[i]
    model_cur.initialize(initialize_solver=True)
    # Write initial condition(s) to file
    results = dict()
    result_folder = pathlib.Path(f"results_R15CellAxisymm_3DStim_E{EVals[i]}kPa")
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
    # Solve
    displayed = False
    while True:
        # Solve the system
        model_cur.monolithic_solve()
        if model_cur.idx_nl[-1] in [0,1]:
            dt_scale=1.4
        elif model_cur.idx_nl[-1] in [2,3,4]:
            dt_scale=1.1
        elif model_cur.idx_nl[-1] in [5,6,7,8,9,10]:
            dt_scale=1.0
        # decrease time step
        elif model_cur.idx_nl[-1] in [11,12,13,14,15,16,17,18,19,20]:
            dt_scale=0.8
        elif model_cur.idx_nl[-1] >= 20:
            dt_scale=0.5
        if model_cur.idx_l[-1] <= 5 and dt_scale>=1.0:
            dt_scale *= 1.05
        if model_cur.idx_l[-1] >= 10:
            dt_scale = min(dt_scale*0.8, 0.8)
        model_cur.set_dt(float(model_cur.dt)*dt_scale)

        # Save results for post processing
        for species_name, species in model_cur.sc.items:
            results[species_name].write(model_cur.sc[species_name].u["u"], model_cur.t)
        # End if we've passed the final time
        if model_cur.t >= model_cur.final_t:
            break
    np.savetxt(result_folder / f"tvec.txt", np.array(model_cur.tvec).astype(np.float32))
