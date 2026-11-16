# Historical sampled reference curves

These `.q` files were moved unchanged from `Mathematica/Sampled/`. They contain
saved scattering curves produced by the geometric sampling programs now in
`Tests/sampling/`. Original directory names and file contents are preserved.

Current CTest sampling wrappers run the drivers in their test working directory;
they do not overwrite these saved curves. The legacy expression-to-reference
comparisons use the separate `.dat` files in `Examples/Validation/`.
