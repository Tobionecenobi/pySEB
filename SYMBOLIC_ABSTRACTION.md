# Making SEB Independent of GiNaC

This document outlines the approach for making the Scattering Equation Builder (SEB) library independent of the GiNaC symbolic manipulation library, allowing it to work with both GiNaC in C++ and SymPy in Python.

## Overview

The key idea is to create an abstraction layer that can work with either GiNaC (in C++) or SymPy (in Python). This allows SEB to be used in Python without requiring GiNaC to be installed, which simplifies the installation process for Python users.

## Implementation Details

### 1. Symbolic Expression Interface

We've created an abstract interface for symbolic expressions in `SymbolicInterface.hpp`. This interface defines all the operations that SEB needs to perform on symbolic expressions, such as:

- Basic arithmetic operations (add, subtract, multiply, divide)
- Mathematical functions (sin, cos, exp, log, etc.)
- Special functions needed for scattering (Bessel functions, error functions, etc.)
- Substitution and evaluation
- String representation (to_string, to_latex, to_python, to_cform)

### 2. GiNaC Implementation

We've implemented this interface for GiNaC in `GiNaCSymbolic.hpp`. This implementation wraps GiNaC expressions and provides the functionality defined in the interface.

### 3. SymPy Implementation

We've implemented the interface for SymPy in `pyseb/symbolic.py`. This implementation wraps SymPy expressions and provides the same functionality as the GiNaC implementation.

### 4. Factory Pattern

We've used the factory pattern to create symbolic expressions. The `SymbolicFactory` class is an abstract factory that can create symbolic expressions. We've implemented this factory for both GiNaC and SymPy.

### 5. Python Bindings

We've created Python bindings for the symbolic interface in `src/symbolic_bindings.cpp`. This allows Python code to interact with the C++ symbolic interface.

## Usage

### C++ Usage

In C++, SEB can be used with GiNaC as before:

```cpp
#include "SEB.hpp"
#include "GiNaCSymbolic.hpp"

int main() {
    // Set up the GiNaC factory as the default symbolic engine
    SymbolicFactory::setInstance(new GiNaCFactory());

    // Use SEB as before
    World w;
    GraphID g = w.Add("GaussianPolymer", "poly1");
    w.Link("GaussianPolymer", "poly2.end1", "poly1.end2");
    GraphID g2 = w.Add(g, "diblockcopolymer");
    SymExprPtr formFactor = w.FormFactor("diblockcopolymer");

    // Print out form factor expression
    std::cout << formFactor->to_string() << std::endl;

    return 0;
}
```

### Python Usage

In Python, SEB can be used with SymPy:

```python
import pyseb
import numpy as np
import matplotlib.pyplot as plt

# The SymPy factory is set up as the default symbolic engine in __init__.py

# Create a world
world = pyseb.World()

# Add a single polymer, named 'poly1'
g = world.Add("GaussianPolymer", "poly1")

# Add a second polymer, named 'poly2'. Poly2's end1 is linked to poly1's end2 forming a diblock copolymer
world.Link("GaussianPolymer", "poly2.end1", "poly1.end2")

# Wrap the g structure naming it "diblockcopolymer"
g2 = world.Add(g, "diblockcopolymer")

# Request the symbolic expression for the form factor of this structure
form_factor = world.FormFactor("diblockcopolymer")

# Print out form factor expression
print(form_factor)

# Print out LaTeX expression for form factor
print(pyseb.to_latex(form_factor))

# Evaluate the form factor for a range of q values
params = {
    "beta_poly1": 1.0,
    "beta_poly2": 1.0,
    "Rg_poly1": 10.0,
    "Rg_poly2": 20.0
}

q_values = np.logspace(-3, 0, 100)
intensity = pyseb.evaluate_expression(form_factor, params, q_values)

# Plot the form factor
plt.figure(figsize=(10, 6))
plt.loglog(q_values, intensity)
plt.xlabel('q')
plt.ylabel('I(q)')
plt.title('Form Factor of Diblock Copolymer')
plt.grid(True, which="both", ls="-")
plt.show()
```

## Benefits

This approach has several benefits:

1. **Simplified Installation**: Python users don't need to install GiNaC, which can be difficult to install on some platforms.
2. **Flexibility**: SEB can be used with either GiNaC or SymPy, depending on the user's preference.
3. **Maintainability**: The symbolic operations are isolated in a well-defined interface, making it easier to maintain and extend.
4. **Performance**: SymPy can be faster than GiNaC for some operations, especially when using NumPy for numerical evaluation.

## Limitations

There are some limitations to this approach:

1. **Feature Parity**: SymPy and GiNaC have different feature sets, so some operations may not be available in both.
2. **Performance**: SymPy may be slower than GiNaC for some operations, especially for complex symbolic manipulations.
3. **Compatibility**: The interface may need to be updated if either GiNaC or SymPy changes significantly.

## Future Work

Future work could include:

1. **More Symbolic Engines**: Implement the interface for other symbolic engines, such as Mathematica or Maple.
2. **Optimizations**: Optimize the interface for better performance.
3. **Extended Functionality**: Add more functionality to the interface as needed.

## Conclusion

By creating an abstraction layer for symbolic expressions, we've made SEB independent of GiNaC, allowing it to work with both GiNaC in C++ and SymPy in Python. This simplifies the installation process for Python users and provides more flexibility for all users.
