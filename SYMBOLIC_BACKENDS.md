# Symbolic Backends

SEB builds symbolic expressions through the reusable `seb-symbolic` layer
instead of depending directly on one symbolic algebra library. The current
connectors are:

- `portable`: always available C++ expression tree.
- `ginac`: native C++ symbolic backend when GiNaC is found.
- `sympy`: Python adapter that converts exported SEB expressions to real
  `sympy.Expr` objects.

## Core idea

SEB code should use `sebsym::Expression` and the backend registry. The older
global `Expression` name remains as a compatibility alias.

```cpp
#include "Symbolic.hpp"

sebsym::initialize();
sebsym::set_backend("portable");

auto x = sebsym::symbol("x");
auto y = sebsym::symbol("y");
auto expr = (x * 2.0 + y).sin()
          + (x / 3.0).exp()
          + (y + 4.0).log()
          + (x + sebsym::e()).sqrt();

double value = expr.subs("x", 1.25).subs("y", 2.5).eval();
```

Python users can build expressions with pySEB and convert them to SymPy:

```python
import pyseb

pyseb.set_backend("portable")
x = pyseb.symbol("x")
y = pyseb.symbol("y")
expr = pyseb.sin(x * 2 + y) + pyseb.sqrt(x + pyseb.e())

sympy_expr = pyseb.to_sympy(expr)
value = float(sympy_expr.subs({"x": 1.25, "y": 2.5}).evalf())
```

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

## Backend roles

- `portable`: always-available C++ expression tree used as the neutral core.
  It can build expressions, substitute symbols, evaluate fully numeric
  expressions for supported elementary functions, and export SymPy/Python text.
- `ginac`: native C++ symbolic algebra with numeric evaluation, series support,
  and printers.
- `sympy`: Python adapter exposed through pySEB. It selects the portable C++
  backend and converts exported expression text into real `sympy.Expr` objects.

Future backends should be added as peer implementations, not by adding
library-specific behavior to `Expression`, `World`, or subunit classes.

## Capabilities

Backends advertise optional behavior through `SymbolicCapabilities`.
`sebsym::Expression` checks capabilities before operations such as numeric
evaluation, series expansion, polynomial conversion, and coefficient
extraction. For example, `portable` reports numeric evaluation and Python
export, but not series expansion.

The common tested numeric subset across `portable`, `ginac`, and SymPy
conversion includes arithmetic, explicit substitution, `sin`, `exp`, `log`,
`sqrt`, `pi`, and `e`.

The portable backend also exports special functions such as `erf`, `erfc`, and
Bessel J0/J1 to SymPy syntax. Exact numeric parity for special functions should
be added once every native backend implements the same functions rather than
approximations.

## Practical next steps

1. Replace remaining `USE_GINAC`/`USE_SYMPY` branching in binding code with selected
   backend registration where possible.
2. Add parity tests for more special functions once native backend
   implementations are exact enough to compare numerically.
3. Keep `Expression` as the stable compatibility alias for old user code.
4. Treat Python/SymPy evaluation as a Python backend feature rather than a C++
   numeric-evaluation feature.
