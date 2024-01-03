import dolfin as d
import smart
import h5py
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import tkinter as tk
import os
import re

from tkinter import filedialog
from tkinter import messagebox
root = tk.Tk()
root.withdraw()

def analyze_all(mesh_file="", results_path="", display=True):
    """
    Function for post-processing of XDMF files written from SMART simulations
    Currently relies on loading in the hdf5 file that stores the dolfin mesh
    with associated facet and cell markers for domains.
    If the mesh_file and results_path are not provided, a dialog
    box opens to allow the user to specify the respective paths.
    If display=True, then a plot displays with normalized plots
    of all variables in the results folder over time.
    Currently, this function iterates over all hdf5 files in the results folder
    and returns the spatial average of each variable over time.

    Args:
        mesh_file: path to hdf5 file with dolfin mesh and mesh functions
        results_path: path to folder with SMART results (XDMF files and all associated HDF5s)
    
    Returns:
        tVec: list of time vectors from each results file
        results_stored: list of vectors with each variable spatially averaged at each time point
    """
    if mesh_file == "":
        messagebox.showinfo(title="Load mesh file.",
                            message="Please select the mesh file.")
        mesh_file = filedialog.askopenfilename()
    parent_mesh = smart.mesh.ParentMesh(
        mesh_filename=mesh_file,
        mesh_filetype="hdf5",
        name="parent_mesh",)

    comm = parent_mesh.mpi_comm
    dmesh = parent_mesh.dolfin_mesh
    dim = parent_mesh.dimensionality
    mf_cell = d.MeshFunction("size_t", dmesh, dim, value=0)
    mf_facet = d.MeshFunction("size_t", dmesh, dim-1, value=0)
    hdf5 = d.HDF5File(comm, mesh_file, "r")
    hdf5.read(mf_cell, f"/mf{dim}")
    hdf5.read(mf_facet, f"/mf{dim-1}")
    hdf5.close()

    if results_path == "":
        messagebox.showinfo(title="Load results directory",
                            message="Please select the results directory")
        results_path = filedialog.askdirectory()

    results_file_list = []
    tVec = []
    for file in os.listdir(results_path):
        # check the files which are end with specific extension
        if file.endswith(".h5"):
            results_file_list.append(file)
        if file.endswith(".xdmf") and tVec==[]:
            xdmf_file = open(f"{results_path}/{file}", "r")
            xdmf_string = xdmf_file.read()
            found_pattern = re.findall(r"Time Value=\"?[^\s]+", xdmf_string)
            for i in range(len(found_pattern)):
                tVec.append(float(found_pattern[i][12:-1]))

    cell_vals = np.unique(mf_cell.array())
    cell_vals = cell_vals[np.logical_and(cell_vals !=0, cell_vals < 1e9)]
    facet_vals = np.unique(mf_facet.array())
    facet_vals = facet_vals[np.logical_and(facet_vals !=0, facet_vals < 1e9)]
    child_meshes = []
    child_mesh_len = []
    for i in range(len(cell_vals)+len(facet_vals)):
        if i < len(cell_vals):
            mesh = d.MeshView.create(mf_cell, cell_vals[i])
        else:
            mesh = d.MeshView.create(mf_facet, facet_vals[i-len(cell_vals)])
        child_meshes.append(mesh)
        child_mesh_len.append(len(mesh.coordinates()))

    results_stored = []
    for j in range(len(results_file_list)):
        cur_file = h5py.File(f"{results_path}/{results_file_list[j]}", "r")
        try:
            test_array = cur_file["VisualisationVector"]["0"][:]
        except:
            results_stored.append([])
            continue
        find_mesh = len(test_array) == np.array(child_mesh_len)
        if len(np.nonzero(find_mesh)[0]) > 1:
            ValueError("Could not identify submesh")
        else:
            cur_mesh = child_meshes[np.nonzero(find_mesh)[0][0]]

        Vcur = d.FunctionSpace(cur_mesh, "P", 1)
        dvec = d.Function(Vcur)
        time_pt_tags = cur_file["VisualisationVector"].keys()
        num_time_points = len(time_pt_tags)
        dof_map = d.dof_to_vertex_map(Vcur)[:]
        dx_cur = d.Measure("dx", domain=cur_mesh)
        vol_cur = d.assemble(1.0*dx_cur)
        var_avg = []

        for i in range(num_time_points-1):
            try:
                cur_array = cur_file["VisualisationVector"][str(i)][:]
            except:
                var_avg.append(0.0)
                continue
            # array matches mesh ordering; reorder according to dof mapping for Vcur
            cur_array = cur_array[dof_map]
            dvec.vector().set_local(cur_array)
            dvec.vector().apply("insert")
            var_avg.append(d.assemble(dvec*dx_cur)/vol_cur)
            print(f"Done with time step {i} for file {j}")

        results_stored.append(var_avg)
        if var_avg == []:
            continue
        else:
            var_avg = np.array(var_avg)
            max_idx = min([len(var_avg), len(tVec)])
            if display:
                plt.plot(tVec[0:max_idx], var_avg[0:max_idx]/max(var_avg), label=results_file_list[j][0:-3])
    if display:
        plt.legend()
        plt.show()

    return (tVec, results_stored)