# SEB → pyPRISM prototype

This example builds an eight-site freely jointed chain in SEB, passes its
intramolecular correlations to pyPRISM, and solves for a solution of those
chains with hard-sphere repulsion between sites on different molecules.
It uses existing SEB functionality; there are no changes to the C++ core.

## Run

Use a separate Python 3.11 environment with pySEB installed or a compatible
local pySEB development build. The following versions support pyPRISM's legacy
build on Linux (a C compiler with OpenMP support and Python headers is needed):

```bash
python3.11 -m venv /tmp/pyseb-prism-env
source /tmp/pyseb-prism-env/bin/activate
python -m pip install 'numpy==1.23.5' 'scipy==1.15.2' 'Cython==0.29.37' \
    'matplotlib==3.10.1' pint setuptools wheel
python -m pip install --no-build-isolation --no-deps 'pyPRISM==1.0.4'
# If a local Python 3.11 SEB build is not already available:
python -m pip install pyseb
python Examples/pyprism_seb.py --output /tmp/pyseb-prism-results
```

Run from the repository root, or launch the example by absolute path. Ordinary
pySEB use does not require pyPRISM. Installing pyPRISM 1.0.4 with default build
isolation on Python 3.11 can fail because its shipped C extension was generated
for older Python versions. Preinstalling Cython and NumPy and disabling build
isolation regenerates it. NumPy 1.23.5 also retains legacy types used by that
extension. The prototype itself does not use pyPRISM's trajectory extension.

The bridge and solver also ran with NumPy 2.2.4 at runtime, using a pyPRISM wheel
built separately with NumPy 1.23.5 and Cython 0.29.37. The conservative environment
above keeps the build and runtime versions consistent.

The script runs from top to bottom: build the chain, evaluate SEB, supply the
result to pyPRISM, solve, and plot. Edit `n_sites`, `bond_length`, and
`site_density` near the top to explore other parameters. Outputs are:

- `pyprism_seb.png`: isolated-chain and solution scattering, plus the
  intermolecular site pair correlation.
- `scattering.csv`: SEB input, solution scattering, and scattering per volume
  for unit scattering length.
- `pair_correlation.csv`: intermolecular site pair correlation versus distance.

## What SEB supplies

Each site is an analytical `Point` with unit scattering length. Neighboring sites
are connected by `ThinRod` subunits with zero scattering length. Only the sites
scatter; the rods provide fixed bond lengths with independently averaged
orientations. This gives discrete freely jointed chains without converting a
continuous density into artificial interaction sites.

For N identical unit-weight sites, SEB's normalized form factor is

```text
F(q) = <sum_i sum_j exp(i q · (r_i - r_j))> / N²
omega(q) = N F(q)
```

The script evaluates `World.EvaluateFormFactor` directly on `system.domain.k`
and supplies `pyPRISM.omega.FromArray(omega, k=q)`. In particular, omega(0)=N
and omega tends to 1 at high q, retaining discrete-site self-scattering.

## What pyPRISM adds

The interaction-site diameter sigma is the length unit. Bond length is 1.5 sigma,
site number density is 0.1 / sigma³, chain number density is 0.1/N / sigma³, and
kT=1. The site diameter specifies intermolecular repulsion; it does not give the
point scatterers a finite scattering radius. These are illustrative reduced
units, not parameters fitted to a particular material.

pyPRISM uses the hard-sphere potential and Percus–Yevick closure to calculate
intermolecular correlations. Its site-density-normalized structure factor
includes both intra- and intermolecular scattering:

```text
S_site(q) = omega(q) + rho_site h_inter(q)
I(q)/V = b² rho_site S_site(q)
```

Here b is the common, q-independent site scattering length; the exported
intensity column uses b=1. In the dilute limit S_site approaches omega, not the
SEB form factor F. Do not multiply F by S_site again: that would double-count
intramolecular scattering. Absolute experimental units require assigning
physical units to sigma and b.

## Validation and limits

The example checks solver success and residuals. The separate regression tests
check zero-angle normalization, the high-q self term, agreement with an
independent bond-separation sum and pyPRISM's built-in `FreelyJointedChain`,
agreement of the two solved structure factors, and the dilute limit at
1e-8 / sigma³. This keeps validation machinery out of the example's main flow.

The default calculation gave maximum absolute errors of approximately
2.4e-12 for both omega and the solved structure factor. The dilute result
matched omega to about 1.1e-7 relative error. Exact errors depend on the numerical
environment. Regression checks can be run with:

```bash
python -m unittest discover -s pyseb/tests/python -p test_pyprism_prototype.py
```

This validates the connection and normalization, not the physical accuracy of
PRISM for arbitrary architectures. The model uses isotropic, fixed ideal-chain
statistics and the usual site-averaging approximation. Neither intrachain
excluded volume nor interaction-induced shape changes are calculated. There is
one site type; extracting AA, AB, and BB correlations for a diblock is a later
step. Grid convergence and interaction-model validation would be needed before
quantitative use for a material.

References: [pyPRISM array input](https://pyprism.readthedocs.io/en/latest/api/pyPRISM.omega.FromArray.html),
[freely jointed chain](https://pyprism.readthedocs.io/en/latest/api/pyPRISM.omega.FreelyJointedChain.html),
[structure-factor conventions](https://pyprism.readthedocs.io/en/latest/api/pyPRISM.calculate.structure_factor.html).
