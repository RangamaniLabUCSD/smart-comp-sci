from __future__ import annotations
# module use /cm/shared/ex3-modules/202309a/defq/modulefiles
# module load python-fenics-dolfin-2019.2.0.dev0
# . ~/local/src/smart-comp-sci/venv/bin/activate
from pathlib import Path
import pandas as pd
import json
from typing import NamedTuple, Any
from itertools import cycle
import numpy as np
import datetime
import matplotlib.pyplot as plt


class Data(NamedTuple):
    timings: pd.DataFrame
    config: dict[str, Any]
    t: np.ndarray
    c: np.ndarray

    @property
    def enforce_mass_conservation(self) -> bool:
        return self.config.get("flags", {}).get("enforce_mass_conservation", False)
    
    @property
    def num_refinements(self) -> int:
        return int(Path(self.config["mesh_file"]).stem.split("_")[-1])
    
    @property
    def dt(self) -> float:
        return self.config["solver"]["initial_dt"]
    
    @property
    def total_run_time(self) -> float:
        return self.timings[self.timings["name"] == "Total run time"]["wall tot"].values[0]


def load_all_data(main_path: str):
    all_data = []
    for folder in (f for f in Path(main_path).iterdir() if f.is_dir()):
        
        try:
            data = load_data(folder=folder)
        except FileNotFoundError as e:
            print(f"Skipping folder {folder}, due to {e}")
            continue
        all_data.append(data)
    return all_data

def is_default(d: Data) -> bool:
    return np.isclose(d.dt, 0.0002) and not d.enforce_mass_conservation and d.num_refinements == 0

def plot_data(data: list[Data]):
    fig, ax = plt.subplots(3, 1, sharex=True, sharey=True, figsize=(8, 8))
    fig_t, ax_t = plt.subplots(1, 3, sharey=True, figsize=(12, 4))    
    linestyles = [cycle(["-", "--", ":", "-."]) for _ in range(3)]

    # Get the default data
    default_data = next(d for d in data if is_default(d))

    # Plot the temporal convergence
    temporal_convergence_data = sorted([default_data] + [di for di in data if not is_default(di) and not np.isclose(di.dt, 0.0002)], key=lambda x: x.dt)
    for d in temporal_convergence_data:
        ax[0].plot(d.t, d.c, linestyle=next(linestyles[0]), label=f'dt={d.dt}')

    x = np.arange(len(temporal_convergence_data))
    ax_t[0].bar(x, [d.total_run_time for d in temporal_convergence_data])
    ax_t[0].set_xticks(x)
    ax_t[0].set_xticklabels([d.dt for d in temporal_convergence_data])
    ax_t[0].set_xlabel("dt")
    ax_t[0].set_title("Temporal convergence")
    ax_t[0].set_ylabel("Time [s]")
        
    # Plot the spatial convergence
    spatial_convergence_data = sorted([default_data] + [di for di in data if not is_default(di) and di.num_refinements > 0], key=lambda x: x.num_refinements)
    for d in spatial_convergence_data:
        ax[1].plot(d.t, d.c, linestyle=next(linestyles[1]), label=f'# refinements = {d.num_refinements}')

    x = np.arange(len(spatial_convergence_data))
    ax_t[1].bar(x, [d.total_run_time for d in spatial_convergence_data])
    ax_t[1].set_xticks(x)
    ax_t[1].set_xticklabels([d.num_refinements for d in spatial_convergence_data])
    ax_t[1].set_xlabel("# refinements")
    ax_t[1].set_title("Spatial convergence")
    

    # Plot mass conservation
    mass_conservation_data = next(d for d in data if d.enforce_mass_conservation)
    ax[2].plot(default_data.t, default_data.c, linestyle=next(linestyles[2]), label='no mass conservation')
    ax[2].plot(mass_conservation_data.t, mass_conservation_data.c, linestyle=next(linestyles[2]), label='mass conservation')
    x = np.arange(2)
    ax_t[2].bar(x, [default_data.total_run_time, mass_conservation_data.total_run_time])
    ax_t[2].set_xticks(x)
    ax_t[2].set_xticklabels(["No", "Yes"])
    ax_t[2].set_xlabel("Mass conservation")
    ax_t[2].set_title("Mass conservation")

    for axi in ax:
        axi.legend()
        axi.set_ylabel("[Ca$^{2+}$]")

    for axi in ax_t:
        axi.grid()
        axi.set_yscale("log")
        axi.set_ylim(1e3, 3e5)

    ax[0].set_title("Temporal convergence")
    ax[1].set_title("Spatial convergence (dt = 0.0002)")
    ax[2].set_title("Mass conservation")
    ax[2].set_xlabel("Time [s]")

    fig.savefig("ca2+results.png")
    fig_t.savefig("timings.png")


        
def load_timings(folder: Path):
    timings = (folder / "timings.txt").read_text()

    # Read total run time from the start and end timestamp from the logs
    logs = next(folder.glob(f"{folder.name}*stderr.txt")).read_text()
    start_time = datetime.datetime.fromisoformat(logs.splitlines()[0].split(",")[0].strip())
    end_time = datetime.datetime.fromisoformat(logs.splitlines()[-1].split(",")[0].strip())
    total_run_time = (end_time - start_time).total_seconds()

    f = lambda x: len(x) > 0 and "|" not in x

    header = list(map(str.strip, filter(f, timings.splitlines()[0].split("  "))))
    header[0] = "name"
    
    data = []
    for item_str in timings.splitlines()[2:]:  
        item = list(map(str.strip, filter(f, item_str.split("  "))))
        data.append(dict(zip(header, item)))

    item = ["Total run time", 1] + [total_run_time] * (len(header) - 2)
    data.append(dict(zip(header, item)))
    return pd.DataFrame(data)



def load_data(folder: Path = Path("82094")) -> Data:

    config_file = folder / "config.json"
    if not config_file.is_file():
        raise FileNotFoundError(config_file)
    
    t_file = folder / "tvec.npy"
    if not t_file.is_file():
        raise FileNotFoundError(t_file)

    t = np.load(t_file)
    c = np.load(folder / "concVec.npy")

    config = json.loads(config_file.read_text())
    timings = load_timings(folder=folder)

    return Data(timings=timings, config=config, t=t, c=c)
    
def main():
    data = load_all_data("/global/D1/homes/henriknf/smart-comp-sci/dendritic-spine")
    plot_data(data)
    

if __name__ == "__main__":
    main()