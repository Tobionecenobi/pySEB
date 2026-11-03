#!/usr/bin/env python3
"""SEB chain → pyPRISM solution scattering. See pyprism_seb.md for setup."""

from pathlib import Path
import sys

import numpy as np
import pyPRISM

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import pyseb

# Lengths are in interaction-site diameters; density counts sites, not chains.
n_sites = 8
bond_length = 1.5
site_density = 0.1

# Build scattering sites joined by freely rotating, non-scattering rods.
world = pyseb.World("pyPRISM prototype")
for i in range(n_sites):
    if i == 0:
        graph = world.Add("Point", "s0", "site")
    else:
        world.Link("ThinRod", f"b{i}.end1", f"s{i - 1}.point", "bond")
        world.Link("Point", f"s{i}.point", f"b{i}.end2", "site")
world.Add(graph, "chain")
parameters = {"beta_site": 1.0, "beta_bond": 0.0, "L_bond": bond_length}

# Evaluate SEB on pyPRISM's grid. SEB F(0)=1; discrete-site omega(0)=N.
system = pyPRISM.System(["site"], kT=1.0)
system.domain = pyPRISM.Domain(dr=0.025, length=2048)
q = system.domain.k
omega = n_sites * np.asarray(world.EvaluateFormFactor("chain", parameters, q))
system.omega["site", "site"] = pyPRISM.omega.FromArray(omega, k=q)

# Add repulsion between sites on different chains; internal statistics stay fixed.
system.density["site"] = site_density
system.diameter["site"] = 1.0
system.potential["site", "site"] = pyPRISM.potential.HardSphere()
system.closure["site", "site"] = pyPRISM.closure.PercusYevick()
prism = system.createPRISM()
result = prism.solve(options={"disp": False, "fatol": 1e-9, "maxiter": 200})
residual = np.max(np.abs(prism.cost(result.x)))
if not result.success or not np.isfinite(residual) or residual > 1e-7:
    raise RuntimeError(f"pyPRISM did not converge: {result.message}; residual={residual}")

# S includes intra- and intermolecular scattering: I/V = b² * site_density * S.
structure = pyPRISM.calculate.structure_factor(prism, normalize=True)["site", "site"]
pair = pyPRISM.calculate.pair_correlation(prism)["site", "site"]

if __name__ == "__main__":
    import argparse
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("pyprism-output"))
    output = parser.parse_args().output
    output.mkdir(parents=True, exist_ok=True)
    np.savetxt(output / "scattering.csv",
               np.column_stack((q, omega, structure, site_density * structure)),
               delimiter=",", header="q_sigma,omega_seb,S_site,I_per_volume_b1", comments="")
    np.savetxt(output / "pair_correlation.csv",
               np.column_stack((system.domain.r, pair)),
               delimiter=",", header="r_over_sigma,g_inter", comments="")

    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    axes[0].loglog(q, omega, label="Isolated chain (SEB)")
    axes[0].loglog(q, structure, label="Solution (SEB → pyPRISM)")
    axes[0].set(xlim=(q[0], 20), xlabel=r"$q\sigma$", ylabel="Scattering per site")
    axes[0].legend()
    axes[1].plot(system.domain.r, pair)
    axes[1].set(xlim=(0, 8), xlabel=r"$r/\sigma$", ylabel="Intermolecular g(r)")
    fig.savefig(output / "pyprism_seb.png", dpi=160)
    plt.close(fig)
    print(f"Solver residual: {residual:.2g}. Results: {output.resolve()}")
