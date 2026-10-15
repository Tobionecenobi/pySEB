# Hydration reconciliation from atoms to proteins

The staged checks resolve the hydration weight and geometry differences
identified in the initial [spatial-shell experiment](spatial-hydration.md).
No molecular intensity scale or contrast is fitted in these comparisons.
The executables are Pepsi-SAXS 3.0 (May 2020) and IMP/FoXS 2.25.0.

## What changed

Pepsi and FoXS do not implement the same hydration model. There is no single
set of points and weights that should reproduce both at identical parameter
labels. pySEB now supports the water form factors and shell conventions
explicitly, preserving its previous defaults.

### Pepsi

Three missing details were isolated before returning to proteins:

1. The shell grid is anchored to the minimum atomic coordinates with a
   half-cell offset. On an isolated atom it gives 24 points, not the centered
   grid's 26.
2. The exclusion boundary uses the atomic surface radius plus the solvent
   radius. Using radii inferred from excluded volume fails on glycine.
   Using the profile's surface radii reproduces the shell point counts and
   shell radii of gyration through both protein controls.
3. The native hydration basis contains a bulk-density prefactor. At the
   reference settings each 3 A cell has effective weight `27 * 0.334` before
   applying the reported absolute contrast `0.0167` for `--dro 5`.

The density factor is an executable-specific normalization convention,
verified using isolated C, O, CH3, pairs, and proteins. It is not a claim that
the nominal contrast is a unique physical hydration density. The new
`AtomicSAXSHydrationShellConvention.Pepsi3` retains that reference convention;
`VolumeDensity` retains the geometric-volume interpretation. Their stored
`hydration_volume` values have different meanings.

The water factor also needs the O-H separation. `BondedWater` implements the
spherical-amplitude approximation with O-H distance 0.967 A and updates the
corresponding Guinier moment. The method and atomic distances follow
[Grudinin et al. (2017), equation (8), Table 1](https://inria.hal.science/hal-01516719).
Pepsi's five-Gaussian factors differ slightly from pySEB's four-Gaussian
factors; that residual is preserved rather than fitted away.

### FoXS

FoXS's water contribution is *solution water*, not neutral O+2H. The
`FoXSReciprocal` water model subtracts the O and two H displaced-solvent
Gaussians individually, with the 16*pi width used by the actual
[IMP implementation](https://github.com/salilab/imp/blob/2.25.0/modules/saxs/src/FormFactorTable.cpp).
`FoXSFast` implements the normalized water amplitude `exp(-0.115*q^2)`,
corresponding to its common intensity modulation `exp(-0.23*q^2)`.

FoXS's coefficient is `c2 * fw(q) * exposure`. pySEB's mapping is
`contrast * volume_per_exposed_atom = c2 * fw(0)` and the selected water
factor is normalized at q=0. This is a fixed per-exposure mapping, not a
scale inferred from the molecule's I(0).

Two checks distinguish water weights from solvent accessibility:

- With native FoXS exposure fractions, water weights and pair sums agree
  to about 1e-8 relative RMS in reciprocal mode through proteins. Fast-mode
  differences reflect its histogram/interpolation approximations.
- With independently calculated FreeSASA fractions, using the same radii and
  a 1.8 A probe, protein hydration intensities differ by about 1.2–1.3% RMS.
  These differences remain even when the weight formulas agree.

The two-carbon control explains why independent sampling need not duplicate
the native output. For identical expanded spheres separated by 2 A, the
continuous exposed fraction is `1/2+d/(4R) = 0.64804`. FreeSASA gives
0.6485 and 0.6460, while FoXS's default dot sampling gives 0.67114 and
0.65101. The resulting hydration-only intensity differs by about 4.1%.
Changing the water normalization to compensate would be incorrect.
`sasa_use_profile_radii` allows an explicit common-radius calculation without
silently changing the existing ProtOr default.

## Progression and results

`Examples/validate_hydration_ladder.py` runs isolated C/O/CH3 groups, C-C and
C-O pairs, a non-collinear three-atom setup, glycine, then 1TIU and optionally
other proteins. Protein comparisons are gated on the small-control results.
It compares hydration-only terms as well as the complete Pepsi intensity,
and records forward ratios separately from normalized shapes.

| Setup | Pepsi hydration RMS | FoXS reciprocal hydration RMS, common exposure |
|---|---:|---:|
| Isolated C/O/CH3 | 0.063% | below 1e-14 |
| C-C / C-O pairs | 0.063–0.066% | below 1e-14 |
| Three atoms | 0.075% | below 1e-14 |
| Glycine | 0.071% | 2.5e-9 |
| 1TIU model 1 | 0.063% | 1.4e-8 |
| 1CRN | 0.064% | 2.3e-8 |

The FoXS column uses fractional relative RMS; the Pepsi column is percent.
Pepsi hydration forward-intensity ratios are approximately 1.00008 for every
setup, so the close shapes do not depend on a fitted forward normalization.
Independent native runs at 2% and 10% contrast on C, the close C-C pair,
glycine and 1TIU give the same relative residuals as 5%, checking transfer of
the fixed weight mapping rather than refitting at each contrast.
The native/pySEB shell counts are 24 on one atom, 32 on the close pair,
344 on 1TIU, and 190 on 1CRN, with matching shell radii of gyration.

Complete Pepsi curves still include atomic-group and excluded-solvent
differences. Their I(0)-normalized relative RMS differences are 0.64% for
1TIU and 1.86% for 1CRN, using the fixed 5% contrast convention. These are
not experimental goodness-of-fit values. The [FoXS full-curve solvent
discrepancy](atomic-saxs-foxs.md) also remains separate from hydration.

## Reproduce and regression coverage

```bash
.venv-sympy/bin/python Examples/validate_hydration_ladder.py \
  --pepsi /path/to/Pepsi-SAXS --imp-python /path/to/IMP/python \
  --output /tmp/hydration-ladder --proteins
```

Add `--extra-protein /path/to/1CRN.pdb` for the crambin control. Coordinates
for 1TIU come from the existing tutorial dataset; 1CRN comes from
[RCSB](https://files.rcsb.org/download/1CRN.pdb). Results, environment hashes,
raw component curves, parameter mappings and input hashes are preserved in
the output directory. Summary snapshots are in `hydration-reconciliation/`.
Native numerical fixtures in `Tests/data/hydration_reference/` support
regression tests without installing external executables.

The added unit tests independently check water-factor formulas, Guinier
moments, the 24-point shell geometry, and two-sphere solvent accessibility.
Native-program regression tests cover all ten setups without a network.
The tests retain explicit tolerances for different atomic factors and
fast-mode approximations; they do not fit reference curves.

Final verification: 96 Python tests passed with one skipped; the rebuilt C++
NumericalSubunitTest and PDBParserTest passed. `git diff --check` was clean.
