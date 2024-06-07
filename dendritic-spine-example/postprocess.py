from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd
import json
from dataclasses import dataclass
from typing import Any
from itertools import cycle
import numpy as np
import re
import matplotlib.pyplot as plt

import dendritic_spine_args


ntasks_pattern = re.compile("ntasks: (?P<n>\d+)")


@dataclass
class Data:
    timings_: list[dict[str, Any]]
    petsc_timings_: list[dict[str, Any]]
    config: dict[str, Any]
    t: np.ndarray
    concVec: np.ndarray
    gradVec: np.ndarray
    stderr: str
    stdout: str
    ntasks: int
    folder: Path = Path("")

    @property
    def mesh(self) -> str:
        return Path(self.config["mesh_file"]).stem

    @property
    def dt(self) -> float:
        return self.config["solver"]["initial_dt"]

    @property
    def timings(self):
        if not hasattr(self, "_timings"):
            self._timings = pd.DataFrame(self.timings_)
            # Now update all dtypes
            self._timings["reps"] = self._timings["reps"].astype(int)
            self._timings["wall avg"] = self._timings["wall avg"].astype(float)
            self._timings["wall tot"] = self._timings["wall tot"].astype(float)
            self._timings["usr avg"] = self._timings["usr avg"].astype(float)
            self._timings["usr tot"] = self._timings["usr tot"].astype(float)
            self._timings["sys avg"] = self._timings["sys avg"].astype(float)
            self._timings["sys tot"] = self._timings["sys tot"].astype(float)
            self._timings["name"] = self._timings["name"].astype(str)

        return self._timings

    @property
    def petsc_timings(self):
        return pd.DataFrame(self.petsc_timings_)

    @property
    def num_refinements(self) -> int:
        if "refined" in self.mesh:
            return int(self.mesh.split("_")[-1])
        return 0

    @property
    def total_run_time(self) -> float:
        try:
            return self.timings[self.timings["name"] == "dendritic-spine-example"][
                "wall tot"
            ].values[0]
        except IndexError:
            return np.nan

    def to_json(self) -> dict[str, Any]:
        return {
            "t": self.t.tolist(),
            "concVec": self.concVec.tolist(),
            "gradVec": self.gradVec.tolist(),
            "config": self.config,
            "timings_": self.timings_,
            "stderr": self.stderr,
            "stdout": self.stdout,
            "ntasks": self.ntasks,
            "petsc_timings_": self.petsc_timings_,
        }


def parse_ntasks(stdout: str) -> int:
    for line in stdout.splitlines():
        if m := re.match(ntasks_pattern, line):
            return int(m.group("n"))

    return 1


def parse_timings(timings: str) -> dict[str, Any]:
    f = lambda x: len(x) > 0 and "|" not in x

    header = list(map(str.strip, filter(f, timings.splitlines()[0].split("  "))))
    header[0] = "name"

    data = []
    for item_str in timings.splitlines()[2:]:
        item = list(map(str.strip, filter(f, item_str.split("  "))))
        data.append(dict(zip(header, item)))

    data.append(dict(zip(header, item)))
    return data


def load_all_data(main_path: Path):
    all_data = []
    for folder in (f for f in Path(main_path).iterdir() if f.is_dir()):
        try:
            data = load_data(folder=folder)
            print(data.mesh)
            print(f"Load data from folder {folder}")
        except FileNotFoundError as e:
            print(f"Skipping folder {folder}, due to {e}")
            continue

        all_data.append(data)
    return all_data


def plot_data(all_data: list[Data], output_folder, format: str = "png"):
    data = [d for d in all_data if d.ntasks == 1]
    fig, ax = plt.subplots(2, 4, sharex=True, sharey="row", figsize=(15, 8))

    mesh2index = {
        "1spine_mesh_coarser": 0,
        "1spine_mesh_coarser_refined_1": 1,
        "1spine_mesh_coarser_refined_2": 2,
        "1spine_mesh": 3,
    }

    dts = list(sorted({d.dt for d in data}))
    dts2color = {d: c for d, c in zip(dts, cycle(plt.cm.tab10.colors))}

    timings = [np.zeros_like(dts) for _ in mesh2index]
    lines = []
    labels = []
    for d in sorted(data, key=lambda x: x.dt):
        index = mesh2index[d.mesh]
        timings[index][dts.index(d.dt)] = d.total_run_time

        (l,) = ax[0, index].plot(d.t, d.concVec, label=d.dt, color=dts2color[d.dt])
        if index == 0:
            lines.append(l)
            labels.append(f"{d.dt:.2e}")
        ax[1, index].plot(d.t, d.gradVec, label=d.dt, color=dts2color[d.dt])

    for k, v in mesh2index.items():
        ax[0, v].set_title(" ".join(k.split("_")))

    ax[0, 0].set_ylabel("Average Cytosolic calcium (μM)")
    ax[1, 0].set_ylabel("Average Gradient of Cytosolic calcium (μM)")

    lgd = fig.legend(
        lines, labels, title="Time step", loc="center right", bbox_to_anchor=(1.1, 0.5)
    )
    fig.subplots_adjust(right=0.99)
    fig.savefig(
        (output_folder / "results.png").with_suffix(f".{format}"),
        bbox_extra_artists=(lgd,),
        bbox_inches="tight",
    )

    # Plot timings
    fig_t, ax_t = plt.subplots()

    x = np.arange(0, len(dts))
    ax_t.bar(
        x - 0.25,
        timings[0],
        width=0.25,
        label=" ".join(list(mesh2index.keys())[0].split("_")),
    )
    ax_t.bar(
        x, timings[1], width=0.25, label=" ".join(list(mesh2index.keys())[1].split("_"))
    )
    ax_t.bar(
        x + 0.25,
        timings[2],
        width=0.25,
        label=" ".join(list(mesh2index.keys())[2].split("_")),
    )
    ax_t.set_yscale("log")
    ax_t.set_xticks(x)
    ax_t.set_xticklabels(dts)
    ax_t.legend()
    ax_t.grid()
    fig_t.savefig((output_folder / "timings.png").with_suffix(f".{format}"))


def plot_refinement_study(all_data: list[Data], output_folder, format: str = "png"):
    data = sorted(
        [d for d in all_data if "coarser" in d.mesh], key=lambda x: x.num_refinements
    )
    dts = list(sorted({d.dt for d in data}))
    num_refinements = list(sorted({d.num_refinements for d in data}))
    # Find index where we have all refinements
    ref_data = [None] * len(num_refinements)
    for d in data:
        print(d.dt, d.num_refinements)#, d.folder.stem)
        # index = num_refinements.index(d.num_refinements)
        # if np.isclose(d.dt, 0.00025):
        #     ref_data[index] = [d.t, d.concVec]

    plotted = set()
    fig, ax = plt.subplots(1, len(num_refinements), sharex=True, sharey="row", figsize=(8, 3))
    # fig_err, ax_err = plt.subplots(1, len(num_refinements), sharex=True, sharey="row", figsize=(8, 2))
    lines = []
    labels = []
    for d in data:
        # index = dts.index(d.dt)
        index = num_refinements.index(d.num_refinements)
        # key = (d.num_refinements, index)
        key = (d.dt, index)
        if key in plotted:
            continue
        plotted.add(key)
        # (l,) = ax[0, index].plot(d.t, d.concVec, label=d.num_refinements)
        (l,) = ax[index].plot(d.t, d.concVec, label=d.dt)
        # interpCur = np.interp(ref_data[index][0], d.t, d.concVec)
        # ax_err[index].plot(ref_data[index][0], 100*(interpCur-ref_data[index][1])/ref_data[index][1])
        # ax[1, index].plot(d.t, d.gradVec, label=d.num_refinements)
        # ax[0, index].set_title(f"dt = {d.dt}")
        # if index == 3:
        #     lines.append(l)
        #     labels.append(d.num_refinements)
        ax[index].set_title(f"{d.num_refinements} refinements")
        ax[index].set_xlabel("Time (s)")
        if index == 2:
            lines.append(l)
            labels.append(d.dt)
    ax[0].set_ylabel("Average Cytosolic calcium (μM)")
    # ax[1, 0].set_ylabel("Average Gradient of Cytosolic calcium (μM)")
    lgd = fig.legend(lines, labels, title="Time step (s)")
    # lgd = fig.legend(
    #     lines,
    #     labels,
    #     title="Number of refinements",
    #     loc="center right",
    #     bbox_to_anchor=(1.1, 0.5),
    # )
    # fig.subplots_adjust(right=0.9)
    fig.savefig(
        (output_folder / "refinement_study.png").with_suffix(f".{format}"),
        bbox_extra_artists=(lgd,),
        bbox_inches="tight",
    )


def load_timings(folder: Path) -> list[dict[str, Any]]:
    timings = (folder / "timings.txt").read_text()

    # # Read total run time from the start and end timestamp from the logs
    # logs = next(folder.glob(f"{folder.name}*stderr.txt")).read_text()
    # start_time = datetime.datetime.fromisoformat(logs.splitlines()[0].split(",")[0].strip())
    # end_time = datetime.datetime.fromisoformat(logs.splitlines()[-1].split(",")[0].strip())
    # total_run_time = (end_time - start_time).total_seconds()

    f = lambda x: len(x) > 0 and "|" not in x

    header = list(map(str.strip, filter(f, timings.splitlines()[0].split("  "))))
    header[0] = "name"

    data = []
    for item_str in timings.splitlines()[2:]:
        item = list(map(str.strip, filter(f, item_str.split("  "))))
        data.append(dict(zip(header, item)))

    # item = ["Total run time", 1] + [total_run_time] * (len(header) - 2)
    data.append(dict(zip(header, item)))
    return data


def load_data(folder: Path = Path("82094")) -> Data:
    config_file = folder / "config.json"
    if not config_file.is_file():
        raise FileNotFoundError(config_file)

    t_file = folder / "tvec.txt"
    if not t_file.is_file():
        raise FileNotFoundError(t_file)

    t = np.loadtxt(t_file)
    concVec = np.load(folder / "concVec.npy")
    gradVec = np.load(folder / "gradVec.npy")

    if (petsc_timings_file := folder / "petsc_timings.py").is_file():
        petsc_timings_ = {}
        exec(petsc_timings_file.read_text(), petsc_timings_)
        stages = petsc_timings_.pop("Stages")["Main Stage"]
        petsc_timings = []
        for k, v in stages.items():
            petsc_timings.append({"name": k, **v[0]})

    else:
        petsc_timings = []

    config = json.loads(config_file.read_text())
    timings = load_timings(folder=folder)
    if (stdout_file := folder / f"{folder.name}-dendritic_spine-stdout.txt").is_file():
        stdout = stdout_file.read_text()
        ntasks = parse_ntasks(stdout=stdout)
    else:
        stdout = ""
        ntasks = 1

    if (stderr_file := folder / f"{folder.name}-dendritic_spine-stderr.txt").is_file():
        stderr = stderr_file.read_text()
    else:
        stderr = ""

    return Data(
        timings_=timings,
        config=config,
        t=t,
        concVec=concVec,
        gradVec=gradVec,
        stderr=stderr,
        stdout=stdout,
        ntasks=ntasks,
        folder=folder,
        petsc_timings_=petsc_timings,
    )


def plot_timings_stacked(all_data: list[Data], output_folder, format: str = "png"):
    data_refined0 = [
        d for d in all_data if "coarser" in d.mesh and d.num_refinements == 0
    ]
    data_refined1 = [
        d for d in all_data if "coarser" in d.mesh and d.num_refinements == 1
    ]
    data_refined2 = [
        d for d in all_data if "coarser" in d.mesh and d.num_refinements == 2
    ]

    timings_refined0 = [d.timings for d in data_refined0]
    timings_refined0_concat = pd.concat(tuple(timings_refined0))
    timings_refined0_group = timings_refined0_concat.groupby("name")
    timings_refined0_mean = timings_refined0_group.mean()
    timings_refined0_std = timings_refined0_group.std()

    petsc_timings_refined0 = [d.petsc_timings for d in data_refined0]
    petsc_timings_refined0_concat = pd.concat(tuple(petsc_timings_refined0))
    petsc_timings_refined0_group = petsc_timings_refined0_concat.groupby("name")
    petsc_timings_refined0_mean = petsc_timings_refined0_group.mean()
    petsc_timings_refined0_std = petsc_timings_refined0_group.std()

    timings_refined1 = [d.timings for d in data_refined1]
    timings_refined1_concat = pd.concat(tuple(timings_refined1))
    timings_refined1_group = timings_refined1_concat.groupby("name")
    timings_refined1_mean = timings_refined1_group.mean()
    timings_refined1_std = timings_refined1_group.std()

    petsc_timings_refined1 = [d.petsc_timings for d in data_refined1]
    petsc_timings_refined1_concat = pd.concat(tuple(petsc_timings_refined1))
    petsc_timings_refined1_group = petsc_timings_refined1_concat.groupby("name")
    petsc_timings_refined1_mean = petsc_timings_refined1_group.mean()
    petsc_timings_refined1_std = petsc_timings_refined1_group.std()

    timings_refined2 = [d.timings for d in data_refined2]
    timings_refined2_concat = pd.concat(tuple(timings_refined2))
    timings_refined2_group = timings_refined2_concat.groupby("name")
    timings_refined2_mean = timings_refined2_group.mean()
    timings_refined2_std = timings_refined2_group.std()

    petsc_timings_refined2 = [d.petsc_timings for d in data_refined2]
    petsc_timings_refined2_concat = pd.concat(tuple(petsc_timings_refined2))
    petsc_timings_refined2_group = petsc_timings_refined2_concat.groupby("name")
    petsc_timings_refined2_mean = petsc_timings_refined2_group.mean()
    petsc_timings_refined2_std = petsc_timings_refined2_group.std()

    petsc_names = list(petsc_timings_refined0_mean.index)
    # names = list(timings_refined0_mean.index)

    # Create bar plot
    fig, ax = plt.subplots()

    keys = [
        "KSPSolve",
        "SNESJacobianEval",
        "SNESFunctionEval",
    ]

    indices = [list(petsc_names).index(key) for key in keys]
    x = np.arange(0, 3)
    width = 0.7
    bottom = np.zeros_like(x)

    colors = plt.cm.tab10.colors

    total_index = timings_refined0_mean.index.tolist().index(
        "dendritic-spine-example [main]"
    )

    total_time_mean = [
        timings_refined0_mean["wall tot"].values[total_index],
        timings_refined1_mean["wall tot"].values[total_index],
        timings_refined2_mean["wall tot"].values[total_index],
    ]
    tot_name = "dendritic-spine-example [main]"
    all_total_times = [np.unique(timings_refined0_group.get_group(tot_name)["wall tot"].values), 
                       np.unique(timings_refined1_group.get_group(tot_name)["wall tot"].values), 
                       np.unique(timings_refined2_group.get_group(tot_name)["wall tot"].values)]
    total_time_std = [
        timings_refined0_std["wall tot"].values[total_index],
        timings_refined1_std["wall tot"].values[total_index],
        timings_refined2_std["wall tot"].values[total_index],
    ]
    total_time_mean = [
        petsc_timings_refined0_mean["time"]["summary"],
        petsc_timings_refined1_mean["time"]["summary"],
        petsc_timings_refined2_mean["time"]["summary"]
    ]
    all_total_times = [
        petsc_timings_refined0_group.get_group("summary").time.values, 
        petsc_timings_refined1_group.get_group("summary").time.values, 
        petsc_timings_refined2_group.get_group("summary").time.values
    ]
    total_time_std = [
        petsc_timings_refined0_std["time"]["summary"],
        petsc_timings_refined1_std["time"]["summary"],
        petsc_timings_refined2_std["time"]["summary"],
    ]
    

    assemble_time = [
        petsc_timings_refined0_mean["time"].values[indices[1]]
        + petsc_timings_refined0_mean["time"].values[indices[2]],
        petsc_timings_refined1_mean["time"].values[indices[1]]
        + petsc_timings_refined1_mean["time"].values[indices[2]],
        petsc_timings_refined2_mean["time"].values[indices[1]]
        + petsc_timings_refined2_mean["time"].values[indices[2]],
    ]
    assemble_time_std = [
        np.sqrt(petsc_timings_refined0_std["time"].values[indices[1]]**2
        + petsc_timings_refined0_std["time"].values[indices[2]]**2),
        np.sqrt(petsc_timings_refined1_std["time"].values[indices[1]]**2
        + petsc_timings_refined1_std["time"].values[indices[2]]**2),
        np.sqrt(petsc_timings_refined2_std["time"].values[indices[1]]**2
        + petsc_timings_refined2_std["time"].values[indices[2]]**2),
    ]
    all_assemble_times =[
        petsc_timings_refined0_group.get_group(petsc_names[indices[1]]).time.values +
        petsc_timings_refined0_group.get_group(petsc_names[indices[2]]).time.values,
        petsc_timings_refined1_group.get_group(petsc_names[indices[1]]).time.values +
        petsc_timings_refined1_group.get_group(petsc_names[indices[2]]).time.values,        
        petsc_timings_refined2_group.get_group(petsc_names[indices[1]]).time.values +
        petsc_timings_refined2_group.get_group(petsc_names[indices[2]]).time.values,
    ]

    ksp_time = [
        petsc_timings_refined0_mean["time"].values[indices[0]],
        petsc_timings_refined1_mean["time"].values[indices[0]],
        petsc_timings_refined2_mean["time"].values[indices[0]],
    ]
    ksp_time_std = [
        petsc_timings_refined0_std["time"].values[indices[0]],
        petsc_timings_refined1_std["time"].values[indices[0]],
        petsc_timings_refined2_std["time"].values[indices[0]],
    ]
    all_ksp_times =[
        petsc_timings_refined0_group.get_group(petsc_names[indices[0]]).time.values,
        petsc_timings_refined1_group.get_group(petsc_names[indices[0]]).time.values,        
        petsc_timings_refined2_group.get_group(petsc_names[indices[0]]).time.values,
    ]

    rest_time = [
        total_time_mean[0] - assemble_time[0] - ksp_time[0],
        total_time_mean[1] - assemble_time[1] - ksp_time[1],
        total_time_mean[2] - assemble_time[2] - ksp_time[2],
    ]
    rest_time_std = [
        np.sqrt(total_time_std[0]**2 - assemble_time_std[0]**2 - ksp_time_std[0]**2),
        np.sqrt(total_time_std[1]**2 - assemble_time_std[1]**2 - ksp_time_std[1]**2),
        np.sqrt(total_time_std[2]**2 - assemble_time_std[2]**2 - ksp_time_std[2]**2),
    ]
    all_rest_times =[
        all_total_times[0] - all_assemble_times[0] - all_ksp_times[0],
        all_total_times[1] - all_assemble_times[1] - all_ksp_times[1],
        all_total_times[2] - all_assemble_times[2] - all_ksp_times[2],
    ]
    all_timing_points = [all_assemble_times, all_ksp_times, all_rest_times]

    lines = []
    labels = []
    bottom = np.zeros(len(x))
    # y = np.zeros_like(x)
    for i, (label, yi, yerr) in enumerate(
        [
            ("Assemble", assemble_time, assemble_time_std),
            ("KSP solve", ksp_time, ksp_time_std),
            ("Other", rest_time, rest_time_std),
        ]
    ):
        y = np.divide(yi, total_time_mean)
        # yerr = np.divide(
        #     np.sqrt(
        #         np.square(np.divide(yi, total_time_mean)) * np.square(total_time_std)
        #         + np.square(y) * np.square(np.divide(total_time_std, total_time_mean))
        #     ),
        #     total_time_mean,
        # )
        yerr = np.divide(yerr, total_time_mean)
        l = ax.bar(
            x,
            y,
            width=width,
            yerr=yerr,
            capsize=5,
            color=colors[i],
            bottom=bottom,
        )
        for j in range(len(x)):
            x_cur = np.linspace(x[j]-0.08, x[j]+0.08, len(all_timing_points[i][j]))
            sum_timings = all_timing_points[0][j]
            for k in range(1,i+1):
                sum_timings = sum_timings + all_timing_points[k][j]
            ax.plot(x_cur, np.divide(sum_timings, all_total_times[j]),
                    "o", markersize=4, color=colors[i], markeredgecolor="k")
        lines.append(l)
        labels.append(label)
        bottom += y

    ax.set_xticks(x)
    ax.set_xticklabels(["standard", "fine", "extra fine"])
    ax.set_yticks(np.arange(0, 1.1, 0.2))
    ax.set_ylabel("Percentage of total time")
    ax.set_yticklabels([f"{i:.0%}" for i in np.arange(0, 1.1, 0.2)])
    ax.spines.top.set_visible(False)
    ax.spines.right.set_visible(False)
    # ax.spines.bottom.set_visible(False)
    # ax.spines.left.set_visible(False)

    # # Put legend outside of the plot to the right
    fig.subplots_adjust(right=0.9)
    lgd = fig.legend(
        lines,
        labels,
        loc="center left",
        bbox_to_anchor=(0.9, 0.5),
        # title="Timings",
        title_fontsize="large",
    )

    fig.savefig(
        (output_folder / "timings_profile_stacked.png").with_suffix(f".{format}"),
        bbox_extra_artists=(lgd,),
        bbox_inches="tight",
    )


def plot_timings(all_data: list[Data], output_folder, format: str = "png"):
    data_refined0 = [
        d for d in all_data if "coarser" in d.mesh and d.num_refinements == 0
    ]
    data_refined1 = [
        d for d in all_data if "coarser" in d.mesh and d.num_refinements == 1
    ]
    data_refined2 = [
        d for d in all_data if "coarser" in d.mesh and d.num_refinements == 2
    ]

    timings_refined0 = [d.timings for d in data_refined0]
    timings_refined0_concat = pd.concat(tuple(timings_refined0))
    timings_refined0_mean = timings_refined0_concat.groupby("name").mean()
    timings_refined0_std = timings_refined0_concat.groupby("name").std()

    petsc_timings_refined0 = [d.petsc_timings for d in data_refined0]
    petsc_timings_refined0_concat = pd.concat(tuple(petsc_timings_refined0))
    petsc_timings_refined0_mean = petsc_timings_refined0_concat.groupby("name").mean()
    petsc_timings_refined0_std = petsc_timings_refined0_concat.groupby("name").std()

    timings_refined1 = [d.timings for d in data_refined1]
    timings_refined1_concat = pd.concat(tuple(timings_refined1))
    timings_refined1_mean = timings_refined1_concat.groupby("name").mean()
    timings_refined1_std = timings_refined1_concat.groupby("name").std()

    petsc_timings_refined1 = [d.petsc_timings for d in data_refined1]
    petsc_timings_refined1_concat = pd.concat(tuple(petsc_timings_refined1))
    petsc_timings_refined1_mean = petsc_timings_refined1_concat.groupby("name").mean()
    petsc_timings_refined1_std = petsc_timings_refined1_concat.groupby("name").std()

    timings_refined2 = [d.timings for d in data_refined2]
    timings_refined2_concat = pd.concat(tuple(timings_refined2))
    timings_refined2_mean = timings_refined2_concat.groupby("name").mean()
    timings_refined2_std = timings_refined2_concat.groupby("name").std()

    petsc_timings_refined2 = [d.petsc_timings for d in data_refined2]
    petsc_timings_refined2_concat = pd.concat(tuple(petsc_timings_refined2))
    petsc_timings_refined2_mean = petsc_timings_refined2_concat.groupby("name").mean()
    petsc_timings_refined2_std = petsc_timings_refined2_concat.groupby("name").std()

    petsc_names = list(petsc_timings_refined0_mean.index)
    # names = list(timings_refined0_mean.index)

    # Create bar plot
    fig, ax = plt.subplots(1, 2, figsize=(12, 6))

    keys = [
        "KSPSolve",
        "SNESJacobianEval",
        "SNESFunctionEval",
    ]

    indices = [list(petsc_names).index(key) for key in keys]
    x = np.arange(0, 3)
    width = 0.25
    bottom = np.zeros_like(x)

    colors = plt.cm.tab10.colors

    total_index = timings_refined0_mean.index.tolist().index(
        "dendritic-spine-example [main]"
    )

    total_time_mean = [
        timings_refined0_mean["wall tot"].values[total_index],
        timings_refined1_mean["wall tot"].values[total_index],
        timings_refined2_mean["wall tot"].values[total_index],
    ]
    total_time_std = [
        timings_refined0_std["wall tot"].values[total_index],
        timings_refined1_std["wall tot"].values[total_index],
        timings_refined2_std["wall tot"].values[total_index],
    ]

    lines = []
    labels = []
    for j, (index, key) in enumerate(zip(indices, keys)):
        y = np.zeros_like(x)
        err = np.zeros_like(x)

        counts = np.zeros_like(x)
        counts_err = np.zeros_like(x)
        for k, (timings_mean, timings_std) in enumerate(
            [
                (petsc_timings_refined0_mean, petsc_timings_refined0_std),
                (petsc_timings_refined1_mean, petsc_timings_refined1_std),
                (petsc_timings_refined2_mean, petsc_timings_refined2_std),
            ]
        ):
            y[k] = float(timings_mean["time"].values[index])
            err[k] = float(timings_std["time"].values[index])

            counts[k] = float(timings_mean["count"].values[index])
            counts_err[k] = float(timings_std["count"].values[index])

        print(key, y, bottom)
        l = ax[0].bar(
            x + (j - 1) * width,
            y,
            width=width,
            yerr=err,
            capsize=5,
            color=colors[j],
        )
        lines.append(l)
        labels.append(key)

        ax[1].bar(
            x + (j - 1) * width,
            counts,
            width=width,
            yerr=counts_err,
            capsize=5,
            color=colors[j],
        )

    for axi in ax:
        axi.set_xticks(x)
        # axi.set_xlabel("Number of refinements")
        axi.set_xticklabels(["standard", "fine", "extra fine"])

    (l,) = ax[0].plot(
        [-width, width],
        [total_time_mean[0], total_time_mean[0]],
        "-.",
        color="k",
    )
    lines.append(l)
    labels.append("Total run time (coarse)")
    (l,) = ax[0].plot(
        [1 - width, 1 + width],
        [total_time_mean[1], total_time_mean[1]],
        "--",
        color="k",
    )
    lines.append(l)
    labels.append("Total run time (fine)")
    (l,) = ax[0].plot(
        [2 - width, 2 + width], [total_time_mean[2], total_time_mean[2]], ":", color="k"
    )
    lines.append(l)
    labels.append("Total run time (finest)")

    # Put legend outside of the plot to the right
    fig.subplots_adjust(right=0.9)
    lgd = fig.legend(
        lines,
        labels,
        loc="center left",
        bbox_to_anchor=(0.9, 0.5),
        title="Timings",
        title_fontsize="large",
    )
    ax[0].set_yscale("log")
    ax[0].set_ylabel("Time (s)")
    ax[1].set_ylabel("Number of calls")

    fig.savefig(
        (output_folder / "timings_profile.png").with_suffix(f".{format}"),
        bbox_extra_artists=(lgd,),
        bbox_inches="tight",
    )

    num_dofs = [49194, 323328, 2299090]

    x = np.array(num_dofs)
    fig, ax = plt.subplots()
    ax.errorbar(
        num_dofs,
        total_time_mean,
        yerr=total_time_std,
        linewidth=2,
        fmt="-",
        label="Total run time",
    )
    ax.plot(num_dofs[0], timings_refined0_mean["wall tot"].values[total_index],
            "o", color="b", markersize=4)
    ax.plot(num_dofs[1], timings_refined1_mean["wall tot"].values[total_index],
            "o", color="b", markersize=4)
    ax.plot(num_dofs[2], timings_refined2_mean["wall tot"].values[total_index],
            "o", color="b", markersize=4)
    ax.plot(
        num_dofs,
        5e-4 * x * np.log(x),
        "--",
        color="k",
        linewidth=2,
        label="$\mathcal{O}(N\mathrm{log}N)$",
    )
    ax.plot(num_dofs, 5e-3 * x, ":", color="k", linewidth=2, label="$\mathcal{O}(N)$")
    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.spines.right.set_visible(False)
    ax.spines.top.set_visible(False)
    ax.set_xlabel("Number of degrees of freedom (N)")
    ax.set_ylabel("Total time (s)")
    ax.legend()
    fig.savefig(
        (output_folder / "total_time_vs_dofs.png").with_suffix(f".{format}"),
        bbox_inches="tight",
    )
    print(np.diff(total_time_mean) / np.diff(num_dofs))

    print("Percentage of total time (assembly)")
    print(
        (
            "standard: ",
            (petsc_timings_refined0_mean["time"].values[indices[1]]
            + petsc_timings_refined0_mean["time"].values[indices[2]])
            / timings_refined0_mean["wall tot"].values[total_index],
        )
    )

    print(
        (
            "fine",
            (petsc_timings_refined1_mean["time"].values[indices[1]]
            + petsc_timings_refined1_mean["time"].values[indices[2]])
            / timings_refined1_mean["wall tot"].values[total_index],
        )
    )

    print(
        (
            "extra fine",
            (petsc_timings_refined2_mean["time"].values[indices[1]]
            + petsc_timings_refined2_mean["time"].values[indices[2]])
            / timings_refined2_mean["wall tot"].values[total_index],
        )
    )

    print("Percentage of total time (KSP)")
    print(
        "standard",
        petsc_timings_refined0_mean["time"].values[indices[0]]
        / timings_refined0_mean["wall tot"].values[total_index],
    )
    print(
        "fine: ",
        petsc_timings_refined1_mean["time"].values[indices[0]]
        / timings_refined1_mean["wall tot"].values[total_index],
    )
    print(
        "extra fine: ",
        petsc_timings_refined2_mean["time"].values[indices[0]]
        / timings_refined2_mean["wall tot"].values[total_index],
    )

    print("Setup time")
    index = timings_refined0_mean.index.tolist().index("Setup [main]")
    print(
        "standard: ",
        timings_refined0_mean["wall tot"].values[index]
        / timings_refined0_mean["wall tot"].values[total_index],
    )
    print(
        "fine: ",
        timings_refined1_mean["wall tot"].values[index]
        / timings_refined1_mean["wall tot"].values[total_index],
    )
    print(
        "extra fine: ",
        timings_refined2_mean["wall tot"].values[index]
        / timings_refined2_mean["wall tot"].values[total_index],
    )


def plot_linf_error(
    all_data: list[Data], output_folder, format: str = "png", subdomains=[]
):
    here = Path(__file__).parent
    import dolfin
    import sys

    sys.path.append((here / ".." / "utils").as_posix())

    import smart_analysis

    try:
        import ufl_legacy as ufl
    except ImportError:
        import ufl

    # find data indices associated with coarsest and finest meshes (with min dt)
    num_refinements = []
    dt = []
    for i in range(len(all_data)):
        num_refinements.append(all_data[i].num_refinements)
        dt.append(all_data[i].dt)
    finest_logic = np.logical_and(
        np.array(dt) == min(dt), np.array(num_refinements) == max(num_refinements)
    )
    finest_idx = np.nonzero(finest_logic)[0]
    if len(finest_idx) != 1:
        raise ValueError("Could not find the finest mesh case")
    else:
        finest_idx = finest_idx[0]
    coarsest_logic = np.logical_and(
        np.array(dt) == min(dt), np.array(num_refinements) == min(num_refinements)
    )
    coarsest_idx = np.nonzero(coarsest_logic)[0]
    if len(coarsest_idx) != 1:
        raise ValueError("Could not find the coarsest mesh case")
    else:
        coarsest_idx = coarsest_idx[0]
    coarsest_data = all_data[coarsest_idx]
    finest_data = all_data[finest_idx]
    max_errs_file = output_folder / "max_errs.txt"
    l2_errs_file = output_folder / "l2_errs.txt"
    l1_errs_file = output_folder / "l1_errs.txt"
    avg_coarsest_file = output_folder / "avg_coarsest.txt"
    avg_finest_file = output_folder / "avg_finest.txt"

    if (
        not max_errs_file.is_file()
        or not l2_errs_file.is_file()
        or not l1_errs_file.is_file()
        or not avg_coarsest_file.is_file()
        or not avg_finest_file.is_file()
    ):
        # pull out mesh files from data structure
        mesh_file_coarsest = str(
            here
            / ".."
            / "scripts"
            / "meshes-dendritic-spine"
            / f"{coarsest_data.mesh}.h5"
        )
        mesh_file_finest = str(
            here
            / ".."
            / "scripts"
            / "meshes-dendritic-spine"
            / f"{finest_data.mesh}.h5"
        )
        # Load solutions
        coarsest_solutions = smart_analysis.load_solution(
            mesh_file_coarsest, coarsest_data.folder / "Ca.h5", 0
        )

        finest_solutions = smart_analysis.load_solution(
            mesh_file_finest, finest_data.folder / "Ca.h5", 0
        )

        u_coarsest = next(iter(coarsest_solutions))
        u_finest = next(iter(finest_solutions))

        V_coarsest = u_coarsest.function_space()
        V_finest = u_finest.function_space()
        u_coarsest_interp = dolfin.Function(V_finest)
        u_err = dolfin.Function(V_finest)

        if len(subdomains) > 0:
            coarsest_mf = dolfin.MeshFunction(
                "size_t", V_coarsest.mesh(), V_coarsest.mesh().topology().dim(), 0
            )
            finest_mf = dolfin.MeshFunction(
                "size_t", V_finest.mesh(), V_finest.mesh().topology().dim(), 0
            )
            for k, subdomain in enumerate(subdomains):
                if len(subdomain) == 6:
                    meshes = [V_coarsest.mesh(), V_finest.mesh()]
                    mfs = [coarsest_mf, finest_mf]
                    for j in range(len(meshes)):
                        mesh = meshes[j]
                        mf = mfs[j]
                        # then defines a box to specify region of integration [x0, y0, z0, x1, y1, z1]
                        for c in dolfin.cells(mesh):
                            xCur = c.midpoint().x()
                            yCur = c.midpoint().y()
                            zCur = c.midpoint().z()
                            if (
                                xCur > subdomain[0]
                                and xCur < subdomain[3]
                                and yCur > subdomain[1]
                                and yCur < subdomain[4]
                                and zCur > subdomain[2]
                                and zCur < subdomain[5]
                            ):
                                mf[c] = k + 1
            dx_coarsest = dolfin.Measure(
                "dx", V_coarsest.mesh(), subdomain_data=coarsest_mf
            )
            dx_finest = dolfin.Measure("dx", V_finest.mesh(), subdomain_data=finest_mf)
            vol_coarsest = [dolfin.assemble(1.0 * dx_coarsest)]
            vol_finest = [dolfin.assemble(1.0 * dx_finest)]
            for k in range(len(subdomains)):
                vol_coarsest.append(dolfin.assemble(1.0 * dx_coarsest(k + 1)))
                vol_finest.append(dolfin.assemble(1.0 * dx_finest(k + 1)))
        else:
            dx_coarsest = dolfin.Measure("dx", V_coarsest.mesh())
            vol_coarsest = [dolfin.assemble(1.0 * dx_coarsest)]
            dx_finest = dolfin.Measure("dx", V_finest.mesh())
            vol_finest = [dolfin.assemble(1.0 * dx_finest)]

        u_err_fname = output_folder / "u_err.xdmf"
        u_err_fname.unlink(missing_ok=True)
        u_err_fname.with_suffix(".h5").unlink(missing_ok=True)
        err_file = dolfin.XDMFFile(dolfin.MPI.comm_world, str(u_err_fname))
        err_file.parameters["flush_output"] = True

        avg_coarsest = np.zeros([len(coarsest_data.t), len(subdomains) + 1])
        avg_finest = np.zeros([len(finest_data.t), len(subdomains) + 1])
        avg_coarsest[0][0] = dolfin.assemble(u_coarsest * dx_coarsest) / vol_coarsest[0]
        avg_finest[0][0] = dolfin.assemble(u_finest * dx_finest) / vol_finest[0]
        for k in range(len(subdomains)):
            avg_coarsest[0][k + 1] = (
                dolfin.assemble(u_coarsest * dx_coarsest(k + 1)) / vol_coarsest[k + 1]
            )
            avg_finest[0][k + 1] = (
                dolfin.assemble(u_finest * dx_finest(k + 1)) / vol_finest[k + 1]
            )

        print("Computing errors")
        max_errs = [0.0]
        l2_errs = [0.0]
        l1_errs = [0.0]
        i = 1
        for u_coarsest, u_finest, t in zip(
            coarsest_solutions, finest_solutions, coarsest_data.t
        ):
            avg_coarsest[i][0] = (
                dolfin.assemble(u_coarsest * dx_coarsest) / vol_coarsest[0]
            )
            avg_finest[i][0] = dolfin.assemble(u_finest * dx_finest) / vol_finest[0]
            for k in range(len(subdomains)):
                avg_coarsest[i][k + 1] = (
                    dolfin.assemble(u_coarsest * dx_coarsest(k + 1))
                    / vol_coarsest[k + 1]
                )
                avg_finest[i][k + 1] = (
                    dolfin.assemble(u_finest * dx_finest(k + 1)) / vol_finest[k + 1]
                )

            for j, point in enumerate(V_finest.tabulate_dof_coordinates()):
                u_coarsest_interp.vector()[j] = u_coarsest(point)
                u_err.vector()[j] = u_coarsest_interp.vector()[j] - u_finest.vector()[j]

            err_file.write(u_err, t)

            max_errs.append(max(np.abs(u_err.vector()[:])))
            l2_errs.append(
                np.sqrt(
                    dolfin.assemble((u_coarsest_interp - u_finest) ** 2 * dx_finest)
                )
            )
            l1_errs.append(
                dolfin.assemble(
                    ufl.algebra.Abs(u_coarsest_interp - u_finest) * dx_finest
                )
            )
            print(f"Processed error data {i+1} of {len(coarsest_data.t)}")
            i += 1

        # Save as text files
        np.savetxt(max_errs_file, max_errs)
        np.savetxt(l2_errs_file, l2_errs)
        np.savetxt(l1_errs_file, l1_errs)
        np.savetxt(avg_coarsest_file, avg_coarsest)
        np.savetxt(avg_finest_file, avg_finest)

    max_errs = np.loadtxt(max_errs_file)
    l2_errs = np.loadtxt(l2_errs_file)
    l1_errs = np.loadtxt(l1_errs_file)
    avg_coarsest = np.loadtxt(avg_coarsest_file)
    avg_finest = np.loadtxt(avg_finest_file)

    # Plot errors in three subplots
    fig, ax = plt.subplots(3, 1)  # , figsize=(8, 12))
    ax[0].plot(max_errs)
    ax[0].set_title("Max error")
    ax[1].plot(l2_errs)
    ax[1].set_title("L2 error")
    ax[2].plot(l1_errs)
    ax[2].set_title("L1 error")
    fig.savefig((output_folder / "errors.png").with_suffix(f".{format}"))

    # Plot calcium conc avg differences
    if len(subdomains) > 0:
        cases = []
        coarse_maxes = []
        fine_maxes = []
        coarse_avgs = []
        fine_avgs = []
        for i in range(len(subdomains)):
            cases.append(f"Region{i}")
            coarse_maxes.append(max(avg_coarsest[:, i + 1]))
            fine_maxes.append(max(avg_finest[:, i + 1]))
            coarse_avgs.append(
                np.trapz(avg_coarsest[:, i + 1], coarsest_data.t)
                / np.ptp(coarsest_data.t)
            )
            fine_avgs.append(
                np.trapz(avg_finest[:, i + 1], finest_data.t) / np.ptp(finest_data.t)
            )
        cases = tuple(cases)
        maxes = {"coarse": tuple(coarse_maxes), "fine": tuple(fine_maxes)}
        avgs = {"coarse": tuple(coarse_avgs), "fine": tuple(fine_avgs)}

        x = np.arange(len(cases))  # the label locations
        width = 0.25  # the width of the bars

        fig, ax = plt.subplots(2, 1, figsize=(4, 6))
        multiplier = 0
        for attribute, measurement in maxes.items():
            offset = width * multiplier
            ax[0].bar(x + offset, measurement, width, label=attribute)
            multiplier += 1
        # Add some text for labels, title and custom x-axis tick labels, etc.
        ax[0].set_ylabel("Max calcium (μM)")
        ax[0].set_xticks(x + width / 2, cases)
        ax[0].legend()

        multiplier = 0
        for attribute, measurement in avgs.items():
            offset = width * multiplier
            ax[1].bar(x + offset, measurement, width, label=attribute)
            multiplier += 1
        # Add some text for labels, title and custom x-axis tick labels, etc.
        ax[1].set_ylabel("Avg calcium (μM)")
        ax[1].set_xticks(x + width / 2, cases)
        fig.savefig((output_folder / "summary_ca_error.png").with_suffix(f".{format}"))


def print_dofs(all_data):
    data = sorted(
        [
            d
            for d in all_data
            if "coarser" in d.mesh and np.isclose(d.dt, 0.001) and d.ntasks == 1
        ],
        key=lambda x: x.num_refinements,
    )
    mesh_files = {
        d.num_refinements: (Path(d.config["mesh_file"]), d.folder) for d in data
    }

    import sys

    sys.path.append((Path(__file__).parent / ".." / "utils").as_posix())
    import smart_analysis

    for num_refinement, (mesh_file, results_folder) in mesh_files.items():
        print(f"Load {num_refinement=}")
        solution = next(
            smart_analysis.load_solution(
                mesh_file.as_posix(), results_folder / "Ca.h5", 0
            )
        )
        hmin = solution.function_space().mesh().hmin()
        hmax = solution.function_space().mesh().hmax()
        num_dofs = solution.function_space().dim()
        print(
            f"Mesh: {mesh_file.stem}, {num_refinement=}, {num_dofs=}, {hmin=}, {hmax=}"
        )

    # print(f"Mesh: {d.mesh}")
    # print(f"Number of dofs: {d.concVec.size}")
    # print(f"Number of gradients: {d.gradVec.size}")


def main(
    results_folder: Path,
    output_folder: Path,
    format: str = "png",
    skip_if_processed: bool = False,
    use_tex: bool = False,
) -> int:
    plt.rcParams["text.usetex"] = use_tex

    output_folder.mkdir(exist_ok=True, parents=True)
    results_file = output_folder / "results_dendritic_spine.json"

    if skip_if_processed and results_file.is_file():
        print(f"Load results from {results_file}")
        all_results = [Data(**r) for r in json.loads(results_file.read_text())]
    else:
        print(f"Gather results from {results_folder}")
        all_results = load_all_data(results_folder)
        print(f"Save results to {results_file.absolute()}")
        results_file.write_text(
            json.dumps([r.to_json() for r in all_results], indent=4)
        )
    # plot_data(all_results, output_folder, format=format)
    # plot_refinement_study(all_results, output_folder, format=format)
    # plot_timings(all_results, output_folder, format=format)
    plot_timings_stacked(all_results, output_folder, format=format)
    # plot_linf_error(all_results, output_folder, format=format)
    # print_dofs(all_results)
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    dendritic_spine_args.add_dendritic_spine_postprocess_arguments(parser)
    args = vars(parser.parse_args())
    raise SystemExit(main(**args))
