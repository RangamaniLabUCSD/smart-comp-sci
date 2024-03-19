# ## Compare model results to analytical solution and previous results
#
# Here, we plot the steady-state concentration as a function of cell radius, according to the analytical solution and the SMART numerical solution. The analytical solution for the average concentration in the cytosol at steady state is given in Meyers and Odde 2006 and is included here for ease of reference:
#
# $$
# \bigl< c_{A_{phos}} \bigr> = \frac{6C_1}{R} \left[ \frac{\cosh{\Phi}}{\Phi} - \frac{\sinh{\Phi}}{\Phi^2} \right]\\
# \text{where} \quad C_1 = \frac{k_{kin} c_{Tot} R^2}{\left[3D(1/L_{gradient} - 1/R) + k_{kin}R \right] e^\Phi + \left[3D(1/L_{gradient} + 1/R) - k_{kin}R \right] e^{-\Phi}}\\
# \text{and} \quad \Phi = \frac{R}{L_{gradient}} \quad \text{and} \quad L_{gradient} = \sqrt{\frac{D_{A_{phos}}}{k_p}}
# $$
from typing import NamedTuple, Dict, Any
from pathlib import Path
import pandas as pd
import argparse
import matplotlib.pyplot as plt
import numpy as np
import json
import phosphorylation_parser_args


class Data(NamedTuple):
    t: np.ndarray
    ss: np.ndarray
    l2: np.ndarray
    config: Dict[str, Any]
    timings_: Dict[str, Any]

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
        return pd.DataFrame(self.timings_)

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
        }


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
        config = json.loads((result_folder / "config.json").read_text())
        t = np.loadtxt(result_folder / "tvec.txt")
        ss = np.loadtxt(result_folder / "avg_Aphos.txt")
        l2 = np.loadtxt(result_folder / "L2norm.txt")
        timings = parse_timings((result_folder / "timings.txt").read_text())
        data.append(Data(t=t, ss=ss, l2=l2, config=config, timings_=timings))
        print(data[-1])

    if len(data) == 0:
        raise RuntimeError(f"No results found in folder {folder}")
    return data


def analytical_solution(radius, D=10.0):
    k_kin = 50
    k_p = 10
    cT = 1
    thieleMod = radius / np.sqrt(D / k_p)
    C1 = (
        k_kin
        * cT
        * radius**2
        / (
            (3 * D * (np.sqrt(k_p / D) - (1 / radius)) + k_kin * radius)
            * np.exp(thieleMod)
            + (3 * D * (np.sqrt(k_p / D) + (1 / radius)) - k_kin * radius)
            * np.exp(-thieleMod)
        )
    )
    return (6 * C1 / radius) * (
        np.cosh(thieleMod) / thieleMod - np.sinh(thieleMod) / thieleMod**2
    )


def plot_error_analytical_solution_different_radius(all_results, output_folder, format):

    diffusions = sorted({r.diffusion for r in all_results})
    refinements = sorted({r.refinement for r in all_results})

    for diffusion in diffusions:
        for axisymmetric in [True, False]:
            fig_steady, ax_steady = plt.subplots()
            fig_percent, ax_percent = plt.subplots()
            fig_l2, ax_l2 = plt.subplots()
            fig_time, ax_time = plt.subplots()
            for refinement in refinements:
                results = list(
                    sorted(
                        filter(
                            lambda d: d.refinement == refinement
                            and np.isclose(d.dt, 0.01)
                            and np.isclose(d.diffusion, diffusion)
                            and (d.axisymmetric is axisymmetric),
                            all_results,
                        ),
                        key=lambda d: d.radius,
                    )
                )
                radiusVec = np.array([d.radius for d in results])
                cA = analytical_solution(radiusVec, D=diffusion)
                ss_vec = np.array([d.ss[-1] for d in results])
                percentError = 100 * np.abs(ss_vec - cA) / cA
                l2 = np.array([d.l2 for d in results])
                total_run_time = [
                    float(
                        result.timings[
                            result.timings["name"] == "phosphorylation-example"
                        ].iloc[0]["wall tot"]
                    )
                    for result in results
                ]

                ax_steady.plot(
                    radiusVec,
                    ss_vec,
                    linestyle=":",
                    marker="o",
                    label=f"SMART simulation (refinement: {refinement})",
                )
                ax_percent.plot(
                    radiusVec, percentError, label=f"refinement {refinement}"
                )
                ax_l2.semilogy(radiusVec, l2, label=f"refinement {refinement}")
                rmse = np.sqrt(np.mean(percentError**2))
                print(
                    f"RMSE with respect to analytical solution = {rmse:.3f}% for refinement {refinement}, D: {diffusion}, axisymetric: {axisymmetric}"
                )
                ax_time.plot(radiusVec, total_run_time, marker="o")

            radiusTest = np.logspace(0, 1, 100)
            cA_smooth = analytical_solution(radiusTest, D=diffusion)
            ax_steady.plot(radiusTest, cA_smooth, label="Analytical solution")

            ax_l2.set_xlabel("Cell radius (μm)")
            ax_l2.set_ylabel("$ \| u_e - u \|^2$")
            ax_l2.set_title("$\ell^2$ error from analytical solution")
            ax_l2.legend()
            fig_l2.savefig(
                (
                    output_folder
                    / f"error_radius_l2_diffusion_{diffusion}_axisymmetric_{axisymmetric}.{format}"
                )
            )
            plt.close(fig_l2)

            ax_time.set_xlabel("Cell radius (μm)")
            ax_time.set_ylabel("Total run time [s]")
            ax_time.legend()
            fig_time.savefig(
                output_folder
                / f"total_time_radius_diffusion_{diffusion}_axisymmetric_{axisymmetric}.{format}"
            )
            plt.close(fig_time)

            ax_steady.set_xlabel("Cell radius (μm)")
            ax_steady.set_ylabel("Steady state concentration (μM)")
            ax_steady.legend()
            fig_steady.savefig(
                output_folder
                / f"steady_state_radius_diffusion_{diffusion}_axisymmetric_{axisymmetric}.{format}"
            )
            plt.close(fig_steady)

            ax_percent.set_xlabel("Cell radius (μm)")
            ax_percent.set_ylabel("Percent error from analytical solution")
            ax_percent.legend()
            fig_percent.savefig(
                output_folder
                / f"percent_error_radius_diffusion_{diffusion}_axisymmetric_{axisymmetric}.{format}"
            )
            plt.close(fig_percent)


def plot_error_different_refinements(all_results, output_folder, format):
    radii = sorted({r.radius for r in all_results})
    diffusions = sorted({r.diffusion for r in all_results})

    rates = [("Diffusion", "axisymmetric", "radius", "Rate 0-1", "Rate 1-2")]
    rates_l2 = [("Diffusion", "axisymmetric", "radius", "Rate 0-1", "Rate 1-2")]
    for diffusion in diffusions:
        for axisymmetric in [True, False]:

            fig, ax = plt.subplots()
            fig_l2, ax_l2 = plt.subplots()
            fig_time, ax_time = plt.subplots()
            for radius in radii:
                print(radius)
                results = list(
                    sorted(
                        filter(
                            lambda d: np.isclose(d.radius, radius)
                            and np.isclose(d.dt, 0.01)
                            and np.isclose(d.diffusion, diffusion)
                            and (d.axisymmetric is axisymmetric),
                            all_results,
                        ),
                        key=lambda d: d.refinement,
                    )
                )
                refinements = np.array([d.refinement for d in results])
                cA = analytical_solution(radius, D=diffusion)
                ss_vec = np.array([d.ss[-1] for d in results])
                percentError = 100 * np.abs(ss_vec - cA) / cA
                rates.append(
                    (
                        diffusion,
                        axisymmetric,
                        radius,
                        *np.log(percentError[:-1] / percentError[1:]) / np.log(2),
                    )
                )

                ax.plot(refinements, percentError, marker="o", label=radius)
                ax.set_xticks(refinements)

                l2 = np.array([d.l2 for d in results])
                rates_l2.append(
                    (
                        diffusion,
                        axisymmetric,
                        radius,
                        *np.log(l2[:-1] / l2[1:]) / np.log(2),
                    )
                )
                ax_l2.semilogy(refinements, l2, marker="o", label=radius)
                ax_l2.set_xticks(refinements)

                total_run_time = [
                    float(
                        result.timings[
                            result.timings["name"] == "phosphorylation-example"
                        ].iloc[0]["wall tot"]
                    )
                    for result in results
                ]
                ax_time.semilogy(refinements, total_run_time, marker="o", label=radius)
                ax_time.set_xticks(refinements)

            ax.legend(title="Radius")
            ax.set_xlabel("Refinement")
            ax.set_ylabel("Percent error from analytical solution")
            fig.savefig(
                output_folder
                / f"percent_error_refinement_diffusion_{diffusion}_axisymmetric_{axisymmetric}.{format}"
            )
            plt.close(fig)

            ax_l2.set_xlabel("Refinement")
            ax_l2.set_ylabel("$ \| u_e - u \|^2$")
            fig_l2.savefig(
                (
                    output_folder
                    / f"error_refinement_l2_diffusion_{diffusion}_axisymmetric_{axisymmetric}.{format}"
                )
            )
            plt.close(fig_l2)

            ax_time.legend(title="Radius")
            ax_time.set_xlabel("Refinement")
            ax_time.set_ylabel("Total run time [s]")
            fig_time.savefig(
                (
                    output_folder
                    / f"total_time_refinement_diffusion_{diffusion}_axisymmetric_{axisymmetric}.{format}"
                )
            )
            plt.close(fig_time)

    print("Rates L2")
    for r in rates_l2:
        print(r)
    print("Rates percent error")
    for r in rates:
        print(r)


def plot_error_different_timesteps(all_results, output_folder, format):
    radii = sorted({r.radius for r in all_results})
    diffusions = sorted({r.diffusion for r in all_results})

    for diffusion in diffusions:
        for axisymmetric in [True, False]:
            fig, ax = plt.subplots()
            fig_l2, ax_l2 = plt.subplots()
            fig_time, ax_time = plt.subplots()
            for radius in radii:
                print(radius)
                results = list(
                    sorted(
                        filter(
                            lambda d: np.isclose(d.radius, radius)
                            and d.refinement == 0
                            and np.isclose(d.diffusion, diffusion)
                            and (d.axisymmetric is axisymmetric),
                            all_results,
                        ),
                        key=lambda d: d.dt,
                    )
                )
                dts = np.array([d.dt for d in results])
                cA = analytical_solution(radius, D=diffusion)
                ss_vec = np.array([d.ss[-1] for d in results])
                percentError = 100 * np.abs(ss_vec - cA) / cA
                ax.plot(dts, percentError, marker="o", label=radius)
                ax.set_xticks(dts)

                l2 = np.array([d.l2 for d in results])
                ax_l2.semilogy(dts, l2, marker="o", label=radius)
                ax_l2.set_xticks(dts)

                total_run_time = [
                    float(
                        result.timings[
                            result.timings["name"] == "phosphorylation-example"
                        ].iloc[0]["wall tot"]
                    )
                    for result in results
                ]
                ax_time.semilogy(dts, total_run_time, marker="o", label=radius)
                ax_time.set_xticks(dts)

            ax.legend(title="Radius")
            ax.set_xlabel("Time step [s]")
            ax.set_ylabel("Percent error from analytical solution")
            fig.savefig(
                output_folder
                / f"percent_error_timestep_diffusion_{diffusion}_axisymmetric_{axisymmetric}.{format}"
            )
            plt.close(fig)

            ax_l2.set_xlabel("Time step [s]")
            ax_l2.set_ylabel("$ \| u_e - u \|^2$")
            fig_l2.savefig(
                (
                    output_folder
                    / f"error_timestep_l2_diffusion_{diffusion}_axisymmetric_{axisymmetric}.{format}"
                )
            )
            plt.close(fig_l2)

            ax_time.legend(title="Radius")
            ax_time.set_xlabel("Time step [s]")
            ax_time.set_ylabel("Total run time [s]")
            fig_time.savefig(
                output_folder
                / f"total_time_timestep_diffusion_{diffusion}_axisymmetric_{axisymmetric}.{format}"
            )
            plt.close(fig_time)


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
        (output_folder / "results_phosphorylation.json").write_text(
            json.dumps([r.to_json() for r in all_results], indent=4)
        )

    plot_error_analytical_solution_different_radius(all_results, output_folder, format)
    plot_error_different_refinements(all_results, output_folder, format)
    plot_error_different_timesteps(all_results, output_folder, format)
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    phosphorylation_parser_args.add_phosphorylation_postprocess(parser)
    args = vars(parser.parse_args())
    raise SystemExit(main(**args))
