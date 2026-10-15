# Vendored FreeSASA core

pySEB vendors the C-library portion of FreeSASA 2.1.3 from
https://github.com/mittinatten/freesasa/releases/tag/2.1.3.

The downloaded release archive had SHA-256
`dc4fe377352569299d69329d9ded2f141e527da325960b20d2e30b5335c2b3ff`.
Only the C API, coordinate algorithms, classifiers, and their required parser
support are retained. The command-line program, Python bindings, Gemmi,
JSON, XML, and threaded calculation are not built.

FreeSASA is distributed under the MIT license; see `LICENSE` in this
directory.

Local patch: `src/sasa_sr.c` returns the full probe-expanded sphere area for
an atom with no neighbours. Upstream 2.1.3 otherwise dereferences the empty
neighbour list in this valid single/isolated-atom case.
