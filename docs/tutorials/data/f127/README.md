# F127 SANS data for the corona tutorial

`F127_0_merged_1.txt` is the unchanged experimental curve released by Kexin Dai
and Bradley D. Olsen with *Accelerated small angle neutron scattering algorithms
for polymeric materials*, Soft Matter (2025),
[DOI: 10.1039/D4SM01350F](https://doi.org/10.1039/D4SM01350F).

- Sample: Pluronic F127, 1.0 wt% in D2O.
- [Upstream repository](https://github.com/olsenlabmit/Accelerated-SANS-Code).
- Pinned revision: `7f8a4b2548702f6e47e989266365d3efd771df01`.
- [Pinned data file](https://github.com/olsenlabmit/Accelerated-SANS-Code/blob/7f8a4b2548702f6e47e989266365d3efd771df01/PluronicF-127/Data/F127_0_merged_1.txt).
- Retrieved 2026-09-17; the upstream MIT license is included as `LICENSE`.
- SHA256: `d5879a8c3ffcda8efbd69454f63e0636f28b806c8195489a97019c36fcf00e99`.

The file contains 98 whitespace-separated rows without a header:
q [A^-1], intensity [cm^-1], and supplied intensity standard deviation [cm^-1].
No intensities, uncertainties, or rows have been changed in the stored file.

The tutorial preserves the released fitting script's initial selection
(MATLAB rows `22:end-2`), then selects the higher-q contiguous block after the
single q-order reset. This uses original one-based rows 59–96: 38 points,
q = 0.03511192–0.4641589 A^-1. All other measured rows remain in this file.

## Reproduce the tutorial

The authoritative construction, numerical model, fitting, symbolic checks,
and OpenUSD export code live in
[fit-f127-corona.qmd](../../fit-f127-corona.qmd).
From the repository root, with pySEB and the tutorial dependencies installed:

```sh
quarto render docs/tutorials/fit-f127-corona.qmd --to html --execute --no-cache
```

Only the experimental file is loaded as numerical input. The QMD calculates
its main four-cylinder-arm fit and optional 15-model comparison, and generates
the figures, table, interactive views, and OpenUSD downloads during rendering.
Generated scenes and standalone viewers go under
`docs/tutorials/generated/f127-corona/` and are excluded from version control.
No exploratory scripts, saved fit results, or saved scenes are required.

The tutorial compares effective analytical models on this selected q range.
It does not establish a unique molecular corona structure. Scene conformations
illustrate the models and do not enter the scattering calculation.
