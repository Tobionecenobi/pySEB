# Scattering derivations

[Spheroid derivation (Quarto source)](spheroid.qmd) ·
[Rendered HTML](saved/spheroid.html)

The spheroid tutorial derives the homogeneous volume amplitude, orientationally
averaged form factor, reference-point amplitudes and phase factors, surface-area
weights, and small-q size moments with SymPy. It also discusses sphere and
collapsed limits and gives numerical illustrations.

The source, rendered HTML, and supporting assets were imported from the local
`derivations` branch (tutorial commit `9f83c57b`). The saved HTML is the render from that branch, with local asset paths adjusted
for its new directory. Site rendering uses `spheroid.qmd`. Source adaptations update repository discovery for the new location and make the optional
`wave_interference_demo` visual helper optional: that helper was not present on
the branch. It is not needed for the mathematical derivation.
The ThinDisk comparison also includes the current model's explicit zero-argument
limit in its Bessel integrand.

## Relation to the implemented model

The implemented model is [`models/Spheroid.pyseb.yaml`](../../../models/Spheroid.pyseb.yaml).
Both use equatorial semi-axis `a`, polar semi-axis `c`, and physical surface-area
weighting. The tutorial develops a larger reference-point catalogue:

| Implemented reference | Tutorial name |
| --- | --- |
| `center` | `center` |
| `north`, legacy alias `pole` | `pole_c_plus` |
| `south` | `pole_c_minus` |
| `surface` | `surface` |

The tutorial's equatorial axis tips and distributed `equator` reference are not
currently exposed by the bundled YAML model. Its C++ snippets illustrate the
formulas; the YAML file defines the current API.

## Working with the source

Open the HTML directly to read the saved derivation without running anything.
To regenerate it, install Quarto with its Python/Jupyter support and a Python
environment containing SymPy, NumPy, IPython, Plotly, and pySEB, then run from
the repository root:

```sh
quarto render docs/derivations/spheroid/spheroid.qmd
```

A full render performs symbolic calculations and numerical illustrations; it is
not part of the ordinary unit-test run.
