# Native hydration regression data

Generated on 2026-10-15 by `Examples/validate_hydration_ladder.py` using
Pepsi-SAXS 3.0 (May 2020) and IMP/FoXS 2.25.0. The JSON files contain the
exact coordinate inputs, q grids, native hydration-only intensities,
native FoXS exposure fractions/water forward factors, and Pepsi point counts.
Executables are not distributed. The generation environment and checksums
are recorded in `docs/design/hydration-reconciliation/environment.json`.

Pepsi's hydration component is its `Ihs` output column, converted from cm^-1
using its reported intensity scaling. FoXS's quadratic hydration term is
isolated as `(I(c2=1)+I(c2=-1))/2-I(c2=0)` at c1=0. No fitted molecule-specific
scale, contrast, or background is involved.

C/O/CH3, pairs, and the triad are synthetic diagnostic coordinate inputs.
Glycine comes from `Tests/data/minimal_glycine.pdb`. 1TIU is model 1 of the
existing tutorial heavy-atom dataset; see
`docs/tutorials/data/TITIN-DATA-NOTICE.txt`. 1CRN was downloaded from
https://files.rcsb.org/download/1CRN.pdb; only the same ATOM selection used
in the comparison is relevant to the tests.

Regenerate with the documented command in
`docs/design/hydration-reconciliation.md`. External tools are required only
to regenerate the data; tests use these numerical fixtures without a network.
