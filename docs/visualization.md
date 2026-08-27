# YAML-driven visualization and OpenUSD

Schema version 1 model files may include a `visualization` block.  Geometry is
declared as named `curve`, `surface`, or `random_walk` patches.  Curves use a
bounded `u` domain and XYZ expressions; surfaces use bounded `u` and `v`
domains.  Random walks accept `closure: open` or `closure: bridge` and are
scaled to `target_rg`.  References select fixed positions, curve fractions, or
distributed geometry locations.  Expressions may use model parameters,
visualization definitions, domain variables (`u`, `v`, `t`), constants, and
the standard mathematical functions.

Links constrain positions only.  A labelled distributed reference is sampled
deterministically from the export seed and instance path, and the linked child
receives one uniformly random 3D rotation.  This free orientation is only a
representative realization: it does not alter pySEB's orientationally averaged
scattering calculation.  Geometry overlap is permitted and changing the seed
produces a different realization.

The exporter writes a deterministic representative realization as text USDA;
no OpenUSD runtime is required.  Use `.usda` for explicit text or `.usd` when
an interchange tool expects that extension:

```python
options = pyseb.USDExportOptions(pyseb.LengthUnit.Angstrom)
options.seed = 42
world.export_usd("structure", "structure.usda", parameters, options)
```

The C++ API is equivalent:

```cpp
USDExportOptions options(LengthUnit::Angstrom);
options.seed = 42;
world.ExportUSD("structure", "structure.usda", parameters, options);
```

Curve and surface resolutions default to 96 and 48 samples respectively.
Display color/opacity can be supplied in YAML style defaults or overridden by
full instance path in `USDExportOptions`.  The resulting hierarchy, curves,
meshes, point instancers, and metadata can be opened by Blender or `usdview`.

For example, a cone needs no C++ registration: define a `surface` patch with
`coordinates: ["(1-u) * R * cos(v)", "(1-u) * R * sin(v)", "u * H"]`,
`u: {lower: 0, upper: 1}`, and `v: {lower: 0, upper: "2*pi"}` in YAML.

The complete sphere-and-rods example is available in
`pyseb/examples/python/usd_sphere_rods.py`.
The corresponding sphere-and-polymers example is in
`pyseb/examples/python/usd_sphere_polymers.py`.
A contour-linked binary rod dendrimer is demonstrated in
`pyseb/examples/python/usd_rod_dendrimer.py`.
The contour-to-contour variant is in
`pyseb/examples/python/usd_contour_linked_rod_dendrimer.py`.
