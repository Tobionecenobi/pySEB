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
#include <Symbolic.hpp>

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

Installed C++ consumers can link the symbolic layer through CMake:

```cmake
find_package(SEBSymbolic CONFIG REQUIRED)
target_link_libraries(my_app PRIVATE SEB::symbolic)
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
C++. For example, the portable C++ backend can export SymPy-compatible Python
expression strings, while actual SymPy simplification/evaluation happens in
the Python adapter.

## Backend roles

- `portable`: always-available C++ expression tree used as the neutral core.
  It can build expressions, substitute symbols, evaluate fully numeric
  expressions for supported elementary and special functions, and export
  SymPy/Python text.
- `ginac`: native C++ symbolic algebra with numeric evaluation, series support,
  printers, and GSL-backed numeric special functions.
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

Special-function parity is also tested for `erf`, `erfc`, Bessel J0/J1, and
Dawson. The C++ backends use GSL for numeric special-function evaluation, while
the Python adapter converts portable expressions to real SymPy expressions.

## Practical next steps

1. Keep `Expression` as the stable compatibility alias for old user code.
2. Treat Python/SymPy evaluation as a Python backend feature rather than a C++
   numeric-evaluation feature.
3. Add more special functions as peer backend operations when SEB needs them.
