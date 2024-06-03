from typing import NamedTuple, Dict, Any
from pathlib import Path
import pandas as pd
import argparse
import matplotlib.pyplot as plt
import numpy as np
import re
import json
import phosphorylation_parser_args

ntasks_pattern = re.compile("ntasks: (?P<n>\d+)")
cmap = plt.get_cmap("tab10")
linestyles = ["-", "--", "-.", ":"]
markers = ["", ".", "o"]


class Data(NamedTuple):
    t: np.ndarray
    ss: np.ndarray
    l2: np.ndarray
    config: Dict[str, Any]
    timings_: Dict[int, Dict[str, Any]]
    stderr: str
    stdout: str
    ntasks: int


    @property
    def radius(self):
        return self.config["curRadius"]

    @property
    def dt(self):
        return self.config["time_step"]

    @property
    def diffusion(self):
        return self.config.get("diffusion", 10.0)

    @property
    def axisymmetric(self):
        return "axisymmetric" in self.config["mesh_file"]

    @property
    def timings(self):
        return {k: pd.DataFrame(t) for k, t in self.timings_.items()}

    @property
    def min_run_time(self) -> float:
        """Return the sum of the run time of all ranks."""
        return min(self.run_time_dict.values())
    
    @property
    def mean_run_time(self) -> float:
        """Return the sum of the run time of all ranks."""
        return np.mean(list(self.run_time_dict.values()))

    @property
    def max_run_time(self) -> float:
        """Return the maximum run time of all ranks."""
        return max(self.run_time_dict.values())

    @property
    def run_time_dict(self) -> Dict[int, float]:
        return {
            k: float(v[v["name"] == "phosphorylation-example"].iloc[0]["wall tot"])
            for k, v in self.timings.items()
        }

    @property
    def hmin(self) -> float:
        return self.config.get("hmin", -1.0)

    @property
    def hmax(self) -> float:
        return self.config.get("hmax", -1.0)

    @property
    def refinement(self):
        return int(self.config["mesh_file"].split("refined_")[-1][0])

    def __str__(self):
        return f"radius={self.radius}, dt={self.dt}, refinement={self.refinement}"

    def to_json(self) -> Dict[str, Any]:
        return {
            "t": self.t.tolist(),
            "ss": self.ss.tolist(),
            "l2": float(self.l2),
            "config": self.config,
            "timings_": self.timings_,
            "stderr": self.stderr,
            "stdout": self.stdout,
            "ntasks": self.ntasks,
        }


def parse_ntasks(stdout: str) -> int:
    for line in stdout.splitlines():
        if m := re.match(ntasks_pattern, line):
            return int(m.group("n"))

    return 1


def parse_timings(timings: str) -> Dict[str, Any]:
    f = lambda x: len(x) > 0 and "|" not in x

    header = list(map(str.strip, filter(f, timings.splitlines()[0].split("  "))))
    header[0] = "name"

    data = []
    for item_str in timings.splitlines()[2:]:
        item = list(map(str.strip, filter(f, item_str.split("  "))))
        data.append(dict(zip(header, item)))

    data.append(dict(zip(header, item)))
    return data


def load_results(folder):
    data = []
    for result_folder in Path(folder).iterdir():

        if not result_folder.is_dir():
            continue
        if not (result_folder / "avg_Aphos.txt").exists():
            continue
        stdout = (
            result_folder / f"{result_folder.name}-phosphorylation-stdout.txt"
        ).read_text()
        stderr = (
            result_folder / f"{result_folder.name}-phosphorylation-stderr.txt"
        ).read_text()
        ntasks = parse_ntasks(stdout=stdout)
        config = json.loads((result_folder / "config.json").read_text())
        t = np.loadtxt(result_folder / "tvec.txt")
        ss = np.loadtxt(result_folder / "avg_Aphos.txt")
        l2 = np.loadtxt(result_folder / "L2norm.txt")

        timings = {
            int(f.stem.lstrip("timings_rank")): parse_timings(f.read_text())
            for f in result_folder.glob("timings_rank*")
        }
        data.append(
            Data(
                t=t,
                ss=ss,
                l2=l2,
                config=config,
                timings_=timings,
                stderr=stderr,
                stdout=stdout,
                ntasks=ntasks,
            )
        )
        print(data[-1])

    if len(data) == 0:
        raise RuntimeError(f"No results found in folder {folder}")
    return data


def plot_time_step_vs_error(
    all_results, output_folder, format, axisymmetric=False
):
    radii = sorted({r.radius for r in all_results if r.ntasks == 1})
    diffusions = sorted({r.diffusion for r in all_results if r.ntasks == 1})
    refinements = sorted({r.refinement for r in all_results if r.ntasks == 1})
    time_steps = sorted({r.dt for r in all_results if r.ntasks == 1})
    x = np.arange(len(time_steps))
    width = 0.9 / len(refinements)

    rates = []

    fig, ax = plt.subplots(len(diffusions), 1, sharex=True, figsize=(10, 10))
    fig_t, ax_t = plt.subplots(
        len(radii), len(diffusions), sharex=True, figsize=(10, 10)
    )
    lines = []
    labels = []
    for j, diffusion in enumerate(diffusions):
        for axi in [ax, ax_t]:
            axi[j].set_title(f"D = {diffusion}")

        for k, refinement in enumerate(refinements):
            results = list(
                sorted(
                    filter(
                        lambda d: np.isclose(d.radius, 2.0)
                        and d.refinement == refinement
                        and np.isclose(d.diffusion, diffusion)
                        and d.ntasks == 1
                        and (d.axisymmetric is axisymmetric),
                        all_results,
                    ),
                    key=lambda d: d.dt,
                )
            )

            

            if len(results) == 0:
                continue
    
            dts = np.array([d.dt for d in results])
            l2 = np.array([d.l2 for d in results])
            timings = np.array([d.max_run_time for d in results])

            for r in results:
                rates.append(
                    {
                        "radius": 2.0,
                        "diffusion": diffusion,
                        "refinement": refinement,
                        "hmin": r.hmin,
                        "hmax": r.hmax,
                        "dt": r.dt,
                        "l2": float(r.l2),
                    }
                )

            
            l, = ax[j].loglog(
                dts,
                l2,
                marker="o",
                label=f"h {results[0].hmin:.3f}",
            )
            if j == 0:
                lines.append(l)
                labels.append(f"{results[0].hmin:.3f}")
            ax_t[j].bar(
                x[np.isin(time_steps, dts)] + k * width,
                timings,
                width,
                label=f"refinement {refinement}",
            )
        ax_t[j].set_xlabel("Time step [s]")
        ax_t[j].set_ylabel("Total run time [s]")
        ax_t[j].set_xticks(x + width * len(refinements) / 2)
        ax_t[j].set_xticklabels([f"{t:.1e}" for t in time_steps], rotation=45)
        ax_t[j].set_yscale("log")

        ax[j].set_ylabel(r"$ \| u_e - u \|_2$")

    l, = ax[0].loglog(dts[4:], ((dts[4:] ** 2.1) / dts[-1]) * 0.15, "k--")
    lines.append(l)
    labels.append("$O((\Delta t)^{2.1})$")
    l,  = ax[1].loglog(dts[3:], ((dts[3:] ** 2.6) / dts[-1]) * 0.2, "k:")
    lines.append(l)
    labels.append("$O((\Delta t)^{2.6})$")
    ax[2].loglog(dts[2:], ((dts[2:] ** 2.6) / dts[-1]) * 0.2, "k:")

    ax[-1].set_xlabel("Time step [s]")
    lgd = fig.legend(lines, labels, loc="center right", title="$h$")
    fig.subplots_adjust(right=0.85)

    fig.savefig(
        output_folder / f"time_step_vs_error.{format}",
        bbox_extra_artists=(lgd,),
        bbox_inches="tight",
        dpi=300,
    )
    fig_t.savefig(
        output_folder / f"total_time_vs_error.{format}",
        dpi=300,
    )


def plot_refinement_vs_error(
    all_results, output_folder, format, axisymmetric=False
):
    radii = sorted({r.radius for r in all_results if r.ntasks == 1})
    diffusions = sorted({r.diffusion for r in all_results if r.ntasks == 1})
    refinements = sorted({r.refinement for r in all_results if r.ntasks == 1})
    time_steps = sorted({r.dt for r in all_results if r.ntasks == 1})
    x = np.arange(len(refinements))
    width = 0.9 / len(time_steps)
    rates = []

    fig, ax = plt.subplots(len(diffusions),1, sharex=True, figsize=(10, 10))
    # fig_t, ax_t = plt.subplots(
    #     len(radii), len(diffusions), sharex=True, figsize=(10, 10)
    # )
    lines = []
    labels = []
    for j, diffusion in enumerate(diffusions):
        # if j > 0:
        #     continue
        print(f"D = {diffusion}")
        for axi in [ax]:
            axi[j].set_title(f"D = {diffusion}")

        for k, dt in enumerate(time_steps):
            print(f"dt = {dt}")
            print()
            results = list(
                sorted(
                    filter(
                        lambda d: np.isclose(d.radius, 2.0)
                        and np.isclose(d.dt, dt)
                        and np.isclose(d.diffusion, diffusion)
                        and d.ntasks == 1
                        and (d.axisymmetric is axisymmetric),
                        all_results,
                    ),
                    key=lambda d: d.refinement,
                )
            )

            
            if len(results) == 0:
                continue
    
            hmaxs = np.array([d.hmax for d in results])

            l2 = np.array([d.l2 for d in results])
            timings = np.array([d.max_run_time for d in results])

            for r in results:
                print(len(r.t))
                print("Refinement: ", r.refinement, "(hmin: ", r.hmin, ", hmax: ", r.hmax, ")", "l2: ", r.l2)
                rates.append(
                    {
                        "radius": 2.0,
                        "diffusion": diffusion,
                        "refinement": r.refinement,
                        "hmin": r.hmin,
                        "hmax": r.hmax,
                        "dt": r.dt,
                        "l2": float(r.l2),
                    }
                )

    
            l, = ax[j].loglog(
                hmaxs,
                l2,
                marker="o",
                # label=f"dt={dt:.1e}",
            )
            if j == 0:
                lines.append(l)
                labels.append(f"{dt:.1e}")
   

        ax[j].set_ylabel(r"$ \| u_e - u \|_2$")
    ax[-1].set_xlabel("$h$")

    l, = ax[0].loglog(hmaxs, ((hmaxs ** 2.0) / hmaxs[-1]) * 0.015, "k--")
    lines.append(l)
    labels.append("$O((h)^{2.0})$")
    ax[1].loglog(hmaxs, ((hmaxs ** 2.0) / hmaxs[-1]) * 0.0015, "k--")
    ax[2].loglog(hmaxs, ((hmaxs ** 2.0) / hmaxs[-1]) * 0.00015, "k--")

    lgd = fig.legend(lines, labels, loc="center right", title="Time step [s]")
    fig.subplots_adjust(right=0.85)

    fig.savefig(
        output_folder / f"refinement_vs_error.{format}",
        bbox_extra_artists=(lgd,),
        bbox_inches="tight",
        dpi=300,
    )



def get_convergence_rates(all_results, output_folder, axisymmetric=False):
    radii = sorted({r.radius for r in all_results if r.ntasks == 1})
    diffusions = sorted({r.diffusion for r in all_results if r.ntasks == 1})
    refinements = sorted({r.refinement for r in all_results if r.ntasks == 1})
    time_steps = sorted({r.dt for r in all_results if r.ntasks == 1})

    rates = []


    for j, diffusion in enumerate(diffusions):

        for k, refinement in enumerate(refinements):
            results = list(
                sorted(
                    filter(
                        lambda d: np.isclose(d.radius, 2.0)
                        and d.refinement == refinement
                        and np.isclose(d.diffusion, diffusion)
                        and d.ntasks == 1
                        and (d.axisymmetric is axisymmetric),
                        all_results,
                    ),
                    key=lambda d: d.dt,
                )
            )

            for r in results:
                rates.append(
                    {
                        "diffusion": diffusion,
                        "refinement": refinement,
                        "hmin": r.hmin,
                        "hmax": r.hmax,
                        "dt": r.dt,
                        "l2": float(r.l2),
                    }
                )

    rates_df = pd.DataFrame(rates)
    rates_spatial = []
    rates_temporal = []


    for diffusion in diffusions:
        rates_diffusion = rates_df[rates_df["diffusion"] == diffusion]

        # Spatial convergence rates
        for dt in time_steps:
            rates_dt = rates_diffusion[rates_diffusion["dt"] == dt].sort_values(
                by="hmin", ascending=False
            )
            h = rates_dt.hmin[1:].values
            h_ = rates_dt.hmin[:-1].values

            e = rates_dt.l2[1:].values
            e_ = rates_dt.l2[:-1].values

            r = np.log(e / e_) / np.log(h / h_)
            for i in range(1, len(r)):
                rates_spatial.append(
                    {
                        "diffusion": diffusion,
                        "hmin_i": h[i],
                        "hmin_i-1": h_[i],
                        "l2_i": e[i],
                        "l2_i-1": e_[i],
                        "dt": dt,
                        "rate": r[i],
                    }
                )

        # Temporal convergence rates
        for refinement in refinements:
            rates_ref = rates_diffusion[
                rates_diffusion["refinement"] == refinement
            ].sort_values(by="dt", ascending=False)
            dt = rates_ref.dt.values[1:]
            dt_ = rates_ref.dt.values[:-1]

            e = rates_ref.l2.values[1:]
            e_ = rates_ref.l2.values[:-1]

            r = np.log(e / e_) / np.log(dt / dt_)
            for i in range(1, len(r)):
                rates_temporal.append(
                    {
                        "diffusion": diffusion,
                        "dt_i": dt[i],
                        "dt_i-1": dt_[i],
                        "l2_i": e[i],
                        "l2_i-1": e_[i],
                        "refinement": refinement,
                        "hmin": rates_ref.hmin.values[i],
                        "rate": r[i],
                    }
                )

    rates_spatial_df = pd.DataFrame(rates_spatial)
    rates_temporal_df = pd.DataFrame(rates_temporal)

    rates_spatial_df.to_csv(output_folder / "rates_spatial.csv", index=False)
    rates_temporal_df.to_csv(output_folder / "rates_temporal.csv", index=False)


def plot_convergence_finest(all_results, output_folder, format):
    diffusions = sorted({r.diffusion for r in all_results if r.ntasks == 1})
    hmins = np.array(sorted({r.hmin for r in all_results if r.ntasks == 1}))
    time_steps = np.array(sorted({r.dt for r in all_results if r.ntasks == 1}))

    dt_min = time_steps[0]
    hmin_min = hmins[0]

    # dt
    fig, ax = plt.subplots()
    for diffusion in diffusions:
        results = list(
            sorted(
                filter(
                    lambda d: np.isclose(d.diffusion, diffusion)
                    and d.ntasks == 1
                    and np.isclose(d.hmin, hmin_min),
                    all_results,
                ),
                key=lambda d: d.dt
            )
        )
  
        l2 = np.array([d.l2 for d in results])
        ax.loglog(time_steps, l2, marker="o", label=f"D = {diffusion}")
    ax.loglog(time_steps, 0.1 * (time_steps / time_steps[-1]) ** 2.0, "k--", label="$O((\Delta t)^{2.0})$")
    ax.set_title(f"hmin = {hmin_min}")
    ax.set_xlabel("Time step [s]")
    ax.set_ylabel(r"$ \| u_e - u \|_2$")
    ax.legend(title="Diffusion")
    fig.savefig(
        output_folder / f"convergence_finest_dt.{format}",
        dpi=300,
    )

    # hmin
    fig, ax = plt.subplots()
    for diffusion in diffusions:
        results = list(
            sorted(
                filter(
                    lambda d: np.isclose(d.diffusion, diffusion)
                    and d.ntasks == 1
                    and np.isclose(d.dt, dt_min),
                    all_results,
                ),
                key=lambda d: d.hmin
            )
        )
  
        l2 = np.array([d.l2 for d in results])
        ax.loglog(hmins, l2, marker="o", label=f"D = {diffusion}")

    ax.loglog(hmins, 0.1 * (hmins / hmins[-1]) ** 2.0, "k--", label="$O((h)^{2.0})$")
    ax.set_title(f"dt = {dt_min}")
    ax.set_xlabel("h")
    ax.set_ylabel(r"$ \| u_e - u \|_2$")
    ax.legend(title="Diffusion")
    fig.savefig(
        output_folder / f"convergence_finest_hmin.{format}",
        dpi=300,
    )

def main(
    results_folder: Path,
    output_folder: Path,
    format: str = "png",
    skip_if_processed: bool = False,
    use_tex: bool = False,
) -> int:
    plt.rcParams["text.usetex"] = use_tex

    output_folder.mkdir(exist_ok=True, parents=True)
    results_file = output_folder / "results_phosphorylation.json"
    if skip_if_processed and results_file.is_file():
        print(f"Load results from {results_file}")
        all_results = [Data(**r) for r in json.loads(results_file.read_text())]
    else:
        print(f"Gather results from {results_folder}")
        all_results = load_results(Path(results_folder))
        print(f"Save results to {results_file.absolute()}")
        results_file.write_text(
            json.dumps([r.to_json() for r in all_results], indent=4)
        )

    # Only include results with diffusion greater than 1
    all_results = [r for r in all_results if r.diffusion > 1]

    plot_convergence_finest(all_results, output_folder, format)
    get_convergence_rates(all_results, output_folder)
    plot_time_step_vs_error(all_results, output_folder, format)
    plot_refinement_vs_error(all_results, output_folder, format)
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    phosphorylation_parser_args.add_phosphorylation_postprocess(parser)
    args = vars(parser.parse_args())
    raise SystemExit(main(**args))
