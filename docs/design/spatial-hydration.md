# Spatial hydration in atomic SAXS clouds

This records the initial volume-density implementation. The later
[staged reconciliation](hydration-reconciliation.md) resolves the reference
weighting and grid differences with an explicit Pepsi 3 convention and
separate FoXS water factors. The historical comparisons below use the
original volume-density convention.

Implemented 2026-10-15 as an explicit alternative to atom-centered hydration.
The default remains `AtomicSAXSHydrationModel.AtomCentered`; request
`SpatialShell` before constructing the cloud.

## Geometry and normalization

For inner radius r and thickness t, the shell region is the union of spheres
of radius r+t around positive-occupancy atoms, excluding the union of spheres
of radius r. The default r=t=3 A is inspired by the small-molecule shell
construction in [Grudinin et al. (2017), section 2.3](https://inria.hal.science/hal-01516719).
The inner radius is a center-distance envelope, not an additional offset on
top of each atom's van der Waals radius. No adaptive protein-size rule or
empirical density correction is applied.

Deterministic grid quadrature uses cubic volume elements h^3 (default h=3 A).
Overlaps are counted once. A grid node represents its centered voxel;
inclusion is decided by the node position, so the discretized boundary and
volume require convergence checks. Translation preserves the grid relative
to the molecule; rotation invariance is approached by refinement. Internal
pockets satisfying the distance rule are included. Solvent connectivity and
inter-subunit burial are not modeled.

Each point contributes excess density times volume, multiplied by the
existing water form factor normalized to one at q=0. Signed shell contrast
is supported. Atomic SASA and the old 30 A3-per-exposed-group convention are
not involved. No exposed-atom hydration is added on top of the spatial shell.

Shell points use the same Debye pair and reference-distance machinery as
atomic scatterers, supplying atom-shell/shell-shell interference, attachment
amplitudes, and Guinier moments. Their internal moment includes the normalized
water form factor. The reference named `center` remains the atomic centroid.
Compensated summation reduces cancellation and accumulation error when many
quadrature points contribute to low-q intensities.

Geometry is built even if spatial hydration is initially disabled, and
changing contrast or toggling hydration reuses distances. Grid parameters
require a new cloud. Selecting enabled spatial hydration without shell points
is rejected. Generated points are inspectable in `get_scatterers()` via
positive `hydration_volume` and are exported with coordinates in OpenUSD.
The builder rejects invalid dimensions, empty grids, more than two million
candidate grid nodes, and more than ten thousand shell points. Exact pair
caching remains quadratic in the total atomic-plus-shell count.

## Independent tests

- A non-scattering atom surrounded by the grid shell is compared to the
  analytic uniform spherical-shell amplitude and form factor, multiplied by
  the normalized water form factor. For r=3, outer radius=6 A and h=1 A,
  volume error is 1.3%; normalized amplitude error is below 0.014 over the
  tested q=0.1, 0.3, 0.5 inverse Angstroms. Refinement from 2 to 1 A improves
  both volume and amplitude. These are quadrature tolerances, not experimental
  fitting thresholds.
- A central constant scatterer and two spatial hydration points test all
  interference terms, attachment amplitudes, reference phase factors, and
  positive/negative contrasts against independent closed-form expressions.
- Low-q slopes are checked against the cloud's analytic Guinier moment.
- Tests check union geometry, translation, preservation of the atomic center,
  disabled-hydration equivalence to the original atomic curve, invalid grid
  settings, and rejection of enabling an unbuilt shell.
- An OpenUSD test checks spatial-model metadata, point-volume metadata, and
  preservation of the attachment reference.

## Pepsi comparison: distinguish shape from normalization

The existing [Pepsi benchmark](atomic-saxs-pepsi.md) now accepts
`--spatial-shell` and `--shell-spacing`. It preserves the raw fixed-contrast
comparison and adds a separate diagnostic that matches only the
hydration-induced forward-intensity ratio. This is a sample-dependent scalar
normalization, not a fitted curve or a production default.

For 1TIU model 1:

| Model | Shell points | Shell volume A3 | RMS difference from Pepsi |
|---|---:|---:|---:|
| Atom-centered, contrast 0.0167 e/A3 | — | exposure-weighted | 9.52% |
| Spatial, h=3 A, contrast 0.0167 e/A3 | 667 | 18009 | 32.38% |
| Spatial, h=2 A, contrast 0.0167 e/A3 | 2275 | 18200 | 33.06% |
| Spatial, h=3 A, forward-ratio matched | 667 | 18009 | 1.23% |
| Spatial, h=2 A, forward-ratio matched | 2275 | 18200 | 1.20% |

RMS differences refer to individually I(0)-normalized curves at 101 q points
from 0 to 0.5 inverse Angstroms. The matched contrasts are 0.002860 and
0.002830 e/A3, respectively. The unhydrated baseline is unchanged (0.594% RMS).
At equal nominal contrast, the spatial model increases I(0) by about 62%,
whereas Pepsi increases it by about 9.5%. Refinement changes the geometric
volume by about 1.1%, so grid refinement does not explain this normalization
gap. The new volume-weighted shell is therefore **not Pepsi-equivalent at
the same nominal density contrast**. Its spatial scattering shape is much
closer after isolating that normalization difference.

The tiny-cloud checks are less close to Pepsi's coarse discrete shell:
forward-matched RMS differences on one carbon fall from 27.9% at h=3 A to
8.0% at h=1 A; for two carbon centers 2 A apart they fall from 24.3% to
14.5%. The corresponding raw equal-contrast RMS differences at h=1 A are
16.1% and 19.7%. These results are preserved rather than used to tune a
sample-specific shell grid. Agreement with the independent continuous-shell
solution is the implementation check; Pepsi's discrete shell is a different
quadrature and normalization reference. For crambin at h=3 A the raw RMS
difference is 34.9%, or 2.29% after matching the forward-intensity ratio.

The analytical tests validate the defined shell model, not a unique physical
hydration density for a protein. Do not insert the forward-matched contrasts
as general defaults. An independent derivation of the reference program's
hydration normalization is needed for a direct parameter mapping. FoXS's
different displaced-solvent convention also remains a separate issue.

```bash
.venv-sympy/bin/python Examples/compare_atomic_saxs_pepsi.py \
  --pepsi /path/to/Pepsi-SAXS --spatial-shell --shell-spacing 3 \
  --output /tmp/atomic-pepsi-titin-spatial
```

Metric snapshots are in `spatial-hydration/`; the generated CSVs preserve
the raw and forward-matched curves separately. No production constants were
tuned to the protein comparison.

Validation: the Python suite passes 91 tests with one skipped, including
19 atomic SAXS tests and the spatial OpenUSD check. The rebuilt C++
NumericalSubunitTest and PDBParserTest also pass.
