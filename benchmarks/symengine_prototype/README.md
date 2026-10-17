# SymEngine feasibility experiment

This document records an exploratory prototype. The prototype source, generated
outputs and third-party builds are not included in this commit. These are local
measurements, not a reproducible benchmark suite or an installed pySEB feature.
Production code, dependencies and defaults were unchanged.

## Scope and validation

The experiment tested native C++ SymEngine 0.14.0 with GMP, C++14, Python 3.11
and SymPy 1.14.0 in two positions:

- Transfer an existing portable expression into SymEngine, simplify it, then
  export its tree to SymPy.
- Construct expressions directly using an experimental C++ symbolic factory,
  alongside the portable factory, then simplify and export to SymPy.

Both routes used the same tree exporter and the existing World traversal,
without topology-specific algebra or whole-expression text parsing. Unsupported
special functions were retained as explicit tags and restored during export.

All 27 checks of intensity, unnormalized amplitude and phase passed across
chains, sphere-stars and nested structures. Outputs were compared with native
portable numerical evaluation at q = 0.03, 0.1 and 0.3 using three parameter
sets, including signed contrasts, changed sizes and total contrast matching.
Ordinary expressions and CSE outputs were checked with tolerances
`1e-8 + 2e-8*abs(reference)`. Ten special-function forms were also checked for
symbolic preservation, including definite integrals. CSE across bound-variable
operations was explicitly rejected.

The prototype did not validate clouds/hybrids, arbitrary models, general native
numerical evaluation, series, coefficient extraction, packaging or complete
backend integration. It did not compare performance against GiNaC.

Timings below are medians of three fresh processes per route, with rotating
route order and one BLAS thread. Model construction was excluded. First
expression time includes generation, transfer where applicable, native
simplification and SymPy export. No builds or other benchmarks ran concurrently.

## Findings

Both positions work, but simply adding native simplification did not establish
a faster ordinary export. Medians of three fresh-process measurements per
route, with rotating route order, on this machine:

| Intensity case | Current path | Portable -> SymEngine -> SymPy | SymEngine backend -> SymPy |
| --- | ---: | ---: | ---: |
| 16-polymer chain | 96.9 ms | 250.8 ms | 253.9 ms |
| Sphere + 16 arms | 91.5 ms | 89.4 ms | 87.4 ms |
| 8 nested blocks / 16 polymers | 74.5 ms | 101.9 ms | 105.2 ms |

Small star differences are not established speedup claims. Most extra time is
in SymPy reconstruction, not portable-to-SymEngine transfer. Ordinary chain
text decreases only from 24,501 to 22,649 bytes; the star decreases from 18,614
to 18,516 bytes, and the nested case from 11,269 to 10,873 bytes.

Native SymEngine's tested simplifier does **not** factor `a*b+a*c` into
`a*(b+c)`. In contrast, the Python wheel's `.simplify()` method invokes SymPy's
simplifier. That method therefore cannot demonstrate independent pre-SymPy
compaction. This delegation was verified locally with the SymEngine 0.14.1 Python wheel.

CSE does produce much smaller named-definition output, but **SymPy already
provides comparable compaction**:

| Intensity case | Full current expression bytes | SymPy CSE bytes | SymEngine CSE bytes |
| --- | ---: | ---: | ---: |
| 16-polymer chain | 24,501 | 5,464 | 5,722 |
| Sphere + 16 arms | 18,614 | 3,839 | 3,763 |
| 8 nested blocks / 16 polymers | 11,269 | 3,473 | 3,412 |

CSE bytes include all definitions plus the reduced expression, using the same
Python string printer. These outputs have a different contract from one Expr:
definitions must remain available to the printer/evaluator. Re-expanding them
would undo that compaction. CSE and DAG-export timings are measured **after**
ordinary export and are warmed; they must not be compared with cold first
conversion timings to claim an end-to-end speedup. A future experiment should
measure a dedicated native-CSE-first route before any ordinary SymPy export.

The prototype therefore supports the proposed architecture, but does not
justify a mandatory dependency or a default export switch yet. General
factoring policy and/or a named-intermediate output format need a separate
decision and validation. Raw results, third-party libraries and generated
binaries stay under `/tmp`.


## Decision

Do not add SymEngine as an intermediate compactor or change the default export
path on this evidence. Its more useful role is an optional C++ symbolic backend
that users could select instead of GiNaC. SymPy export would remain a separate
capability. Such a backend still requires implementation and acceptance tests;
this experiment does not establish that it is faster than GiNaC.
