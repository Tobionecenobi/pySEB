# YAML-driven visualization and OpenUSD

Schema version 1 model files may include a `visualization` block.  Geometry is
declared as named `curve`, `surface`, or `random_walk` patches.  Curves use a
bounded `u` domain and XYZ expressions; surfaces use bounded `u` and `v`
domains.  Random walks accept `closure: open` or `closure: bridge` and are
scaled to `target_rg`.  References select fixed positions, curve fractions, or
distributed geometry locations.  Expressions may use model parameters,
visualization definitions, domain variables (`u`, `v`, `t`), constants, and
the standard mathematical functions.

A distributed reference can expose named visualization-only `variants`.  A
reference such as `ends#top#joint7` selects the `top` variant and uses `joint7`
as the deterministic point identity.  The base scattering reference remains
`ends`; the variant affects only the representative visualization.  The legacy
one-label form `ends#joint7` continues to use the reference's normal sampling
rule.

Links constrain positions only.  A labelled distributed reference is sampled
deterministically from the export seed and instance path.  The default
`Random` layout gives each linked child one uniformly random 3D rotation;
geometry overlap is permitted.  This free orientation is only a representative
realization and does not alter pySEB's orientationally averaged scattering
calculation.  Changing the seed produces a different realization.

The optional `Readable` layout evaluates several seeded rotations for each
new child and selects the one with the best clearance from geometry already
placed.  It then performs deterministic relaxation sweeps that rotate complete
descendant subtrees around their external junctions.  These rigid subtree moves
preserve the external junction and every internal link exactly.  Directly
connected geometry participates in collision scoring after masking a small
neighborhood around its intended contact point.

Readable layout never rejects an export when overlap cannot be avoided.
Candidate scoring uses coarse geometry proxies, so this remains deterministic
best-effort layout rather than a collision-free guarantee.  Its USD metadata
identifies the orientations as optimized for readability rather than uniformly
sampled.

The exporter writes a deterministic representative realization as text USDA;
no OpenUSD runtime is required.  Use `.usda` for explicit text or `.usd` when
an interchange tool expects that extension:

```python
options = pyseb.USDExportOptions(pyseb.LengthUnit.Angstrom)
options.seed = 42
options.layout_mode = pyseb.LayoutMode.Readable
options.orientation_trials = 64
options.relaxation_sweeps = 4
options.minimum_clearance = 0.05  # In the selected model length unit.
world.export_usd("structure", "structure.usda", parameters, options)
```

The C++ API is equivalent:

```cpp
USDExportOptions options(LengthUnit::Angstrom);
options.seed = 42;
options.layoutMode = USDLayoutMode::Readable;
options.orientationTrials = 64;
options.relaxationSweeps = 4;
options.minimumClearance = 0.05;
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
A 100-circle center-to-contour chain is demonstrated in
`pyseb/examples/python/usd_circle_chain.py`.
A 100-sphere center-to-surface chain is demonstrated in
`pyseb/examples/python/usd_sphere_chain.py`.
A 20-cylinder end-to-end chain, selecting opposite cap variants on every
interior cylinder, is demonstrated in
`pyseb/examples/python/usd_cylinder_chain.py`.
