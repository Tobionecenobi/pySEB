# Making SEB Independent of GiNaC: Implementation Summary

This document summarizes the approach for making the Scattering Equation Builder (SEB) library independent of the GiNaC symbolic manipulation library, allowing it to work with both GiNaC in C++ and SymPy in Python.

## Overview

We've created an abstraction layer that allows SEB to work with different symbolic engines. This makes it possible to use SEB in Python without requiring GiNaC to be installed, which significantly simplifies the installation process for Python users.

## Key Components

### 1. Symbolic Expression Interface

We've defined an abstract interface (`SymbolicInterface.hpp`) that specifies all the operations SEB needs to perform on symbolic expressions:

```cpp
class SymbolicExpression {
public:
    virtual ~SymbolicExpression() = default;

    // Basic operations
    virtual SymExprPtr add(const SymExprPtr& other) const = 0;
    virtual SymExprPtr subtract(const SymExprPtr& other) const = 0;
    // ... other operations ...

    // String representation
    virtual std::string to_string() const = 0;
    virtual std::string to_latex() const = 0;
    virtual std::string to_python() const = 0;
    virtual std::string to_cform() const = 0;
};
```

### 2. Factory Pattern

We've used the factory pattern to create symbolic expressions:

```cpp
class SymbolicFactory {
public:
    virtual ~SymbolicFactory() = default;
    
    virtual SymExprPtr createSymbol(const std::string& name) = 0;
    virtual SymExprPtr createConstant(double value) = 0;
    virtual SymExprPtr createPi() = 0;
    virtual SymExprPtr createE() = 0;
    virtual SymExprPtr createI() = 0;
    
    // Get the singleton instance
    static SymbolicFactory* instance();
    
    // Set the factory implementation
    static void setInstance(SymbolicFactory* factory);
};
```

### 3. GiNaC Implementation

We've implemented the interface for GiNaC in `GiNaCSymbolic.hpp` and `GiNaCSymbolic.cpp`:

```cpp
class GiNaCExpression : public SymbolicExpression {
public:
    GiNaCExpression(const GiNaC::ex& expr) : _expr(expr) {}
    
    // Implementation of interface methods
    // ...
    
private:
    GiNaC::ex _expr;
};

class GiNaCFactory : public SymbolicFactory {
public:
    // Implementation of factory methods
    // ...
};
```

### 4. SymPy Implementation

We've implemented the interface for SymPy in `pyseb/symbolic.py`:

```python
class SymPyExpression:
    def __init__(self, expr):
        self._expr = expr
    
    # Implementation of interface methods
    # ...

class SymPyFactory:
    # Implementation of factory methods
    # ...
```

### 5. Modified SymbolInterface

We've updated the `SymbolInterface` class to use our abstraction layer instead of directly using GiNaC:

```cpp
class SymbolInterface {
private:
    static SymbolInterface* myInstance;
    map<string, SymExprPtr> symbolDirectory;
    
    // ...
    
public:
    // ...
    
    const SymExprPtr& get(const string& s);
    const SymExprPtr& getSymbol(const string& s);
    // ...
};
```

### 6. Python Bindings

The Python bindings expose the symbolic facade through the main pySEB module.
Common expression construction and `pyseb.to_sympy` are registered in
`src/bindingsTypes.cpp`, while shared world-expression conversion lives in
`src/bindingsSymbolic.cpp`. Backend selection is exposed with
`pyseb.available_backends()`, `pyseb.get_backend()`, and `pyseb.set_backend()`.

### 7. CMake Configuration

We've updated the CMake configuration to support building with or without GiNaC:

```cmake
option(USE_GINAC "Use GiNaC for symbolic manipulation" ON)
option(BUILD_PYTHON "Build Python bindings" OFF)

# Add GiNaC implementation if enabled
if(USE_GINAC)
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(GINAC REQUIRED ginac)
    pkg_check_modules(CLN REQUIRED cln)
    
    list(APPEND SEB_SOURCES
        SEB/GiNaCSymbolic.hpp
        SEB/GiNaCSymbolic.cpp
    )
    
    add_definitions(-DUSE_GINAC)
endif()
```

## Usage

### C++ Usage with GiNaC

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

### Python Usage with SymPy

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

1. **Simplified Installation**: Python users don't need to install GiNaC, which can be difficult to install on some platforms.
2. **Flexibility**: SEB can be used with either GiNaC or SymPy, depending on the user's preference.
3. **Maintainability**: The symbolic operations are isolated in a well-defined interface, making it easier to maintain and extend.
4. **Performance**: SymPy can be faster than GiNaC for some operations, especially when using NumPy for numerical evaluation.

## Next Steps

1. **Testing**: Thoroughly test both the GiNaC and SymPy implementations to ensure they produce the same results.
2. **Documentation**: Update the documentation to explain how to use SEB with either GiNaC or SymPy.
3. **Performance Optimization**: Optimize the symbolic operations for better performance.
4. **Additional Symbolic Engines**: Consider adding support for other symbolic engines, such as Mathematica or Maple.
