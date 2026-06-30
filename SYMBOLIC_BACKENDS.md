# Symbolic Backend Direction

The goal is for SEB to build scattering expressions against a small symbolic
interface instead of depending directly on one symbolic algebra library.

## Core idea

SEB should construct expressions through `SymbolicFactory` and manipulate them
through `SymbolicExpression`. A backend can then wrap GiNaC, SymPy, Sage,
Mathematica/Wolfram, or a lightweight string/AST representation as long as it
implements that contract.

The core construction path should only require:

- symbol and constant creation
- arithmetic operations
- elementary functions
- scattering special functions
- substitution
- string export

More advanced behavior should be optional and discovered through
`SymbolicCapabilities`:

- symbolic simplification
- numeric evaluation
- series expansion and coefficient extraction
- LaTeX output
- Python output
- C code output

This matters because not every useful backend can do every operation inside
C++. For example, the current C++ `SymPySymbolic` backend can build Python
expression strings, but actual SymPy simplification/evaluation happens in
Python.

## Backend roles

Current backends:

- `portable`: always-available C++ expression tree used as the neutral core.
  It can build expressions, substitute symbols, evaluate fully numeric
  expressions for supported elementary functions, and export SymPy/Python text.
- `ginac`: native C++ symbolic algebra with numeric evaluation, series support,
  and printers.
- `sympy`: Python adapter exposed through pySEB. It selects the portable C++
  backend and converts exported expression text into real `sympy.Expr` objects.

Future backends should be added as peer implementations, not by adding
library-specific behavior to `Expression`, `World`, or subunit classes.

## Practical next steps

1. Move GiNaC-specific validation paths behind capability checks.
2. Replace `USE_GINAC`/`USE_SYMPY` branching in binding code with selected
   backend registration where possible.
3. Add backend tests that build the same simple expression with each available
   backend and compare exported forms or numeric values when supported.
4. Keep `Expression` as the stable facade used by subunits and structures.
5. Treat Python/SymPy evaluation as a Python backend feature rather than a C++
   numeric-evaluation feature.
