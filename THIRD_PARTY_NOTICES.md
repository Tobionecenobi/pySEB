# Third-party notices

## fkYAML

pySEB vendors the fkYAML v0.4.1 single-header YAML parser.

- Upstream: https://github.com/fktn-k/fkYAML
- Header SHA-256: `95424af46a4c599c085a859ea1dbcf819ed0558feab4396919b1ba2c7dae39f8`
- License: MIT; the complete text is in `third_party/fkYAML/LICENSE.txt`.

## Atomic SAXS reference data

`AtomicSAXSParameterProfile` contains neutral-atom Cromer--Mann coefficient
values for H, C, N, O, P, and S. The coefficient convention and tabulation are
identified by the IUCr core CIF dictionary, which cites International Tables
for X-ray Crystallography Vol. IV (1974), Table 2.2B, and International Tables
for Crystallography Vol. C (1991), Tables 6.1.1.4--6.1.1.5:

- https://www.iucr.org/__data/iucr/cifdic_html/3_orig/CORE_DIC/Iatom_type_scat.Cromer_Mann_a1.html

The explicit legacy volume control is calculated from the Bondi van der Waals
radii for H, C, N, O, P, and S. The default protein profile uses the displaced
atomic volumes tabulated by Fraser, MacRae, and Suzuki, while the alternative
protein-core profile uses the chemical-group means reported by Pontius,
Richelle, and Wodak (with documented Fraser fallbacks for groups absent from
the compact Pontius table):

- Fraser, MacRae, and Suzuki (1978), *J. Appl. Cryst.* 11, 693--694,
  https://doi.org/10.1107/S0021889878014296
- Pontius, Richelle, and Wodak (1996), *J. Mol. Biol.* 264, 121--136,
  https://doi.org/10.1006/jmbi.1996.0628
- FoXS method and the Fraser displaced-solvent convention:
  https://doi.org/10.1093/nar/gkq461

## FreeSASA

pySEB vendors the C-library subset of FreeSASA 2.1.3 for automatic
solvent-accessible surface calculations.

- Upstream: https://github.com/mittinatten/freesasa/releases/tag/2.1.3
- Release archive SHA-256:
  `dc4fe377352569299d69329d9ded2f141e527da325960b20d2e30b5335c2b3ff`
- License: MIT; the complete text is in
  `third_party/freesasa-2.1.3/LICENSE`.
- Mitternacht (2016), *F1000Research* 5:189,
  https://doi.org/10.12688/f1000research.7931.1

The bundled build excludes the FreeSASA CLI, language bindings, Gemmi,
JSON/XML output, and threaded execution. FoXS, Pepsi-SAXS, and CRYSOL remain
optional external validation tools and are not bundled or called at runtime.
