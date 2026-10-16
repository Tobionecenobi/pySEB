# Investigating a sasmodels adapter for pySEB

Run `python Examples/sasmodels_pyseb_adapter.py` in a pySEB environment with
`sasmodels==1.1.0`. This is a proof of concept, not a public SASAssembly API.

The prototype builds the AOT tutorial's linked `SolidSphere` and
`SolidSphericalShell` in pySEB. It obtains both the form factor and signed
amplitude from that graph, removes their beta normalization symbolically,
and exposes the particle as a vectorized Python `Fq` model to sasmodels.
No replacement sasmodels particle model is used. CPU Python evaluation needs
no OpenCL installation.

sasmodels then performs the two-parameter population averaging and volume
normalization. Checks at all three AOT contrasts verify the signed first
moment, second moment, mean particle volume, absolute-scale intensity, and
zero-dispersity limit. Explicit averaging on the same library-generated nodes
checks the adapter independently of differences in integration grids. This
validates the interface, not convergence of an experimental fit. A separate
comparison against sasmodels' built-in `core_shell_sphere` checks the same
polydisperse intensity at all three contrasts to relative tolerance `1e-8`.

Important boundaries:

- Amplitudes are multiplied by `1e4` and squared amplitudes by `1e8` to convert
  pySEB's excess scattering lengths in angstroms into sasmodels' intensity
  convention. `form_volume` returns the entire core-shell particle volume;
  `scale` therefore represents particle volume fraction.
- sasmodels' lognormal center is the median and its relative width is the
  standard deviation in log space. For tutorial mean `m` and coefficient of
  variation `c`, use `s = sqrt(log(1 + c*c))`, center `m*exp(-s*s/2)`, and
  width `s`. The Gaussian shell center remains the underlying untruncated
  normal center. Positive bounds and finite sampling also affect the average.
- This Python prototype uses the library's built-in finite linear grids.
  The earlier tutorial used Gauss-Hermite and Gauss-Legendre quadrature;
  the updated tutorial uses the smooth Gaussian extension described below.
  Different grids need separate convergence checks before comparing fitted
  parameters. Do not silently reinterpret existing fitted widths.
- Squaring the center amplitude is valid for this spherical particle and is
  checked against its pySEB form factor symbolically. For anisotropic or
  flexible particles, the second moment must come from their form factor,
  not the square of a reference-point amplitude. General interaction adapters
  require an appropriate common amplitude reference and approximation.
- Hydrated atomic/analytical hybrid models, correlated parameter distributions,
  and automatic parameter mapping have not been tested by this prototype.

The AOT tutorial now embeds a compiled adapter generated from the same pySEB
expression. It delegates population averaging, volume normalization, the
hard-sphere product with beta correction, and resolution to sasmodels; SciPy
remains the fitter. Its custom `PositiveGaussian` dispersion supplies smooth
positive-thickness quadrature nodes through sasmodels' extension interface.
The built-in Gaussian grid drops negative nodes as its width changes, which
introduced jumps in the objective during the initial fit experiment. The
small extension avoids that problem without restoring manual averaging loops.

The tutorial explicitly converts lognormal conventions, checks the particle
against native pySEB and the independent sasmodels model at all contrasts,
checks the zero-dispersity and zero-interaction limits, and increases both
population and resolution grids. The particle/structure-factor product uses
`structure_factor_mode=1`. To retain independent intensity volume fraction
`phi` and interaction packing fraction `eta`, it uses product `volfraction=eta`
and `scale=phi/eta`, with a separate zero-eta branch using the particle alone.
The effective radius is obtained from sasmodels' mean particle volume rather
than changing the earlier radius convention. These are tutorial-specific
choices, not a general public adapter API.

References: [Python model interface](https://www.sasview.org/docs/user/qtgui/Perspectives/Fitting/plugin.html),
[resolution API](https://www.sasview.org/docs/dev/sasmodels-api/resolution.html).
