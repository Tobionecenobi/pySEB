# Atomic SAXS comparison with Pepsi-SAXS

The [later staged hydration validation](hydration-reconciliation.md) resolves
the grid, point-weight normalization, and water-factor differences identified
here. This page preserves the earlier baseline comparisons.

See also the [independent FoXS comparison](atomic-saxs-foxs.md), which
isolates a different displaced-solvent Gaussian convention. CRYSOL has not
been tested because no licensed installation is available.

Measured 2026-10-15 using Pepsi-SAXS 3.0 (May 2020), downloaded from
[Inria](https://files.inria.fr/NanoDFiles/Website/Software/Pepsi-SAXS/Linux/3.0/Pepsi-SAXS-Linux.zip).
The external executable is not bundled. Its help describes academic use only.
Executable SHA256: `65fb5640325213a458279d1b73707cb80d25516635b0519b95c2dee26b97a98f`.

## Protocol

`Examples/compare_atomic_saxs_pepsi.py` extracts the first model's ATOM records
and supplies the identical PDB to both calculators. The fixtures are model 1
of the existing 1TIU heavy-atom ensemble (688 atoms) and
[RCSB 1CRN](https://files.rcsb.org/download/1CRN.pdb) (327 atoms).
Hydrogens are implicit; temperature factors are disabled in pySEB. Both use
Fraser excluded volumes, default water solvent, and no instrument smearing.
No experimental curve or fitting is used. Pepsi expansion order is fixed to
30, with 101 q points from 0 to 0.5 inverse Angstroms.

Each intensity is divided by its own I(0), removing concentration and unit
scales. Relative differences are pySEB/Pepsi minus one; RMS is unweighted
over the entire grid. Vacuum comparisons use Pepsi's Iat output column.
Hydration contrasts are fixed at zero and 5% of bulk water; the latter is
0.0167 electrons/A3 in pySEB. pySEB retains 30 A3 per fully exposed group.
The two programs have different hydration geometries, so equal contrast does
not imply equal shell scattering power or equal physical shell volume.

## Baseline results (before bonded-H distance correction)

| Structure | Vacuum RMS / max | No hydration RMS / max | Hydrated RMS / max | Hydration I(0) ratio, pySEB / Pepsi |
|---|---:|---:|---:|---:|
| 1TIU model 1 | 0.43% / 0.94% | 1.05% / 1.90% | 10.72% / 14.57% | 1.0445 / 1.0950 |
| 1CRN | 0.35% / 0.51% | 1.24% / 2.23% | 8.62% / 12.76% | 1.0424 / 1.1017 |

These are cross-model discrepancies, not experimental goodness-of-fit or
formal acceptance thresholds. Vacuum and excluded-solvent shapes agree
approximately, not bit for bit. Atom-group conventions and numerical
approximations have not been isolated as causes of the remaining differences.
Absolute electron-count/intensity equivalence is not established by this
normalized comparison.

For the hydrated 1TIU case, raising Pepsi's expansion order from 30 to 40
changes intensities by at most 0.025%, much smaller than the observed
cross-model discrepancy.

Pepsi reports a spatial hydration layer with 3 A shell width and solvent
radius; pySEB adds water form factors at exposed atom centers. The observed
hydration differences therefore have an identifiable model-level cause.
Changing only the contrast or volume multiplier cannot establish equivalent
spatial scattering. No constants were tuned to make these results agree.

The atomic core can be reviewed with this approximate independent agreement
and the existing analytic tests as evidence. Hydration should remain opt-in
and explicitly approximate. This benchmark does not validate physical
hydration-shell parameters, a complete solution ensemble, or the titin fit.
If Pepsi-like hydration is required, a spatial shell representation and
further validation are needed before claiming that capability.

## Reproduce

```bash
.venv-sympy/bin/python Examples/compare_atomic_saxs_pepsi.py \
  --pepsi /path/to/Pepsi-SAXS --output /tmp/atomic-pepsi
curl -L --fail -o /tmp/1CRN.pdb https://files.rcsb.org/download/1CRN.pdb
.venv-sympy/bin/python Examples/compare_atomic_saxs_pepsi.py \
  --pepsi /path/to/Pepsi-SAXS --pdb /tmp/1CRN.pdb \
  --output /tmp/atomic-pepsi-crambin
```

Each run saves the exact input PDB, input and executable checksums, Pepsi
curves/logs/stdout, comparison CSVs, and JSON metrics. Metric snapshots are
stored alongside this report in `atomic-saxs-pepsi/`.

The baseline Python suite passed 87 tests (one skipped), including all 15
atomic SAXS tests, using the current `build-python` extension. Its CMake build
was up to date. The C++ NumericalSubunitTest and PDBParserTest also passed
after rebuilding the previously stale `build-docs` library.

## Isolating the differences with simple inputs

`Examples/diagnose_atomic_saxs_pepsi.py --pepsi /path/to/Pepsi-SAXS
--output /tmp/atomic-pepsi-groups` generates isolated atomic groups and
C-C/C-O pairs at 2 and 10 A separations. They are mathematical diagnostics,
not chemically complete molecules. Output includes the exact PDB records,
component intensities, and inferred forward amplitudes/volumes.

Before correcting bonded-H geometry, bare C/N/O/S vacuum normalized shapes
disagreed by at most 0.051% over this grid, whereas CH, CH2, and CH3 differed
by 1.30%, 2.30%, and 2.90%. This identifies the coincident implicit-H
approximation as a separate issue from hydration. Pepsi uses five-Gaussian
atomic/group factors; pySEB uses four-Gaussian neutral-atom factors, so exact
agreement is not expected even after correcting distances.

Using Pepsi's single-group component amplitudes, the isolated pair curves
obey `I = f1^2 + f2^2 + 2*f1*f2*sinc(q*d)` to within 0.009% in normalized
vacuum, excluded-solvent, and solution intensities. All these selected
solution amplitudes have the same positive sign. This test removes the
form-factor convention from the check of pair interference.

At 5% hydration contrast, Pepsi's isolated carbon has 24 surrounding shell
points, shell Rg 4.975 A, and effective hydration volume about 216 A3.
pySEB assigns 30 A3 at the atom center (unit exposure). Pepsi's shell-only
I(0.5)/I(0) is about 0.057. Matching its forward amplitude by increasing the
pySEB volume would still leave the spatial distribution different.

The first correction uses the published spherical-amplitude approximation
`f_group(q) = f_heavy(q) + n*f_H(q)*sinc(q*d)` from
[Grudinin et al. (2017), equation (8), Table 1](https://inria.hal.science/hal-01516719).
No distances are fitted to Pepsi output. The corresponding low-q moment
includes `n*f_H(0)*d^2`, and an independent angular-quadrature test checks
the finite-q amplitude and Guinier moment. This is not a reproduction of
Pepsi's charged/resonance group tables. S-H retains the coincident model
because that table supplies no distance. Hydration geometry remains unchanged.

## Results after the bonded-H correction

| Isolated group | Baseline maximum vacuum shape difference | Corrected |
|---|---:|---:|
| CH | 1.30% | 0.132% |
| CH2 | 2.30% | 0.194% |
| CH3 | 2.90% | 0.243% |
| NH | 0.993% | 0.057% |
| OH | 0.788% | 0.065% |

| Structure | Vacuum RMS | No hydration RMS | Hydrated RMS |
|---|---:|---:|---:|
| 1TIU model 1 | 0.078% | 0.594% | 9.52% |
| 1CRN | 0.193% | 1.936% | 7.47% |

The isolated-group correction improves vacuum scattering in both proteins.
The full unhydrated crambin curve becomes less similar to Pepsi despite
that vacuum improvement: the previous errors partly compensated between
vacuum and excluded solvent. This is why protein-level agreement alone is
not sufficient to validate the components. The next isolated discrepancy
is the solvent volume convention: inferred Pepsi group volumes include
C=16.43, CH=21.69, CH2=26.52, CH3=32.03, O=9.203 A3, versus pySEB's
16.44, 21.59, 26.74, 31.89, 9.13 A3. These have not been copied from binary
outputs into the production profile; an explicitly sourced volume convention
would be needed for such a change.

After the correction the Python suite passes 88 tests (one skipped),
including the new bonded-H orientation-average and Guinier check.
The rebuilt C++ NumericalSubunitTest and PDBParserTest also pass.
Snapshots with `-corrected` in their filenames record the new results;
the original protein snapshots and `groups-baseline` preserve the baseline.
