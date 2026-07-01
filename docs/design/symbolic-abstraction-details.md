# Symbolic Expression Abstraction Layer for SEB

This document provides a detailed explanation of the implementation of the symbolic expression abstraction layer for SEB, which allows the library to work with both GiNaC in C++ and SymPy in Python.

## Overview

The goal of this project was to make SEB independent of GiNaC by creating an abstraction layer that supports both GiNaC in C++ and SymPy in Python. This abstraction layer consists of several key components:

1. **SymbolicInterface**: An abstract interface that defines the operations SEB needs to perform on symbolic expressions.
2. **GiNaCSymbolic**: An implementation of the interface for GiNaC.
3. **SymPyExpression**: An implementation of the interface for SymPy.
4. **Expression**: A wrapper class that provides operator overloading to make the code compatible with the existing GiNaC-based code.

## Implementation Details

### 1. SymbolicInterface (seb-symbolic/SymbolicInterface.hpp, seb-symbolic/SymbolicInterface.cpp)

The `SymbolicInterface` class defines the abstract interface for symbolic expressions. It includes methods for:

- Basic arithmetic operations (add, subtract, multiply, divide)
- Mathematical functions (exp, log, sin, cos, etc.)
- Substitution and evaluation
- Series expansion and coefficient extraction
- String representation

This interface allows SEB to work with different symbolic expression libraries without being tied to a specific implementation.

### 2. GiNaCSymbolic (seb-symbolic/GiNaCSymbolic.hpp, seb-symbolic/GiNaCSymbolic.cpp)

The `GiNaCSymbolic` class implements the `SymbolicInterface` for GiNaC. It wraps GiNaC's `ex` type and provides implementations for all the methods defined in the interface. This allows SEB to continue using GiNaC as the symbolic engine while adhering to the abstraction layer.

Key features:
- Wraps GiNaC's `ex` type in a `GiNaCExpression` class
- Implements all methods defined in the `SymbolicInterface`
- Provides factory methods for creating symbols, constants, and special values (pi, e, i)

### 3. Expression (seb-symbolic/Expression.hpp)

The `Expression` class is a wrapper that provides operator overloading to make the code compatible with the existing GiNaC-based code. It delegates all operations to the underlying `SymExprPtr` (a shared pointer to a `SymbolicInterface` implementation).

Key features:
- Overloaded operators (+, -, *, /, etc.)
- Methods for mathematical functions (exp, log, sin, cos, etc.)
- Methods for substitution and evaluation
- Methods for series expansion and coefficient extraction
- Methods for string representation

### 4. SymPy Implementation (pyseb/symbolic.py)

The SymPy implementation provides a Python-based alternative to GiNaC. It implements the same interface as the C++ version but uses SymPy as the symbolic engine.

Key features:
- `SymPyExpression` class that wraps SymPy's expressions
- Implementation of all methods defined in the C++ interface
- Factory methods for creating symbols, constants, and special values

## Updates to Existing Code

To make SEB work with the abstraction layer, we had to update several existing classes:

### 1. World Class (seb/World.hpp, seb/World.cpp)

The `World` class was updated to use the `Expression` class instead of GiNaC's `ex` type. This involved:
- Changing the type of member variables from `ex` to `Expression`
- Updating method signatures to use `Expression` instead of `ex`
- Modifying the implementation to work with the `Expression` class

### 2. Subunit Class (seb/Subunit.hpp, seb/Subunit.cpp)

The `Subunit` class was updated to use the `Expression` class instead of GiNaC's `ex` type. This involved:
- Changing the type of member variables from `ex` to `Expression`
- Updating method signatures to use `Expression` instead of `ex`
- Modifying the implementation to work with the `Expression` class

### 3. Specific Subunit Classes (seb/Subunits/*.hpp)

All specific subunit classes (GaussianPolymer, SolidSphere, SolidCylinder, etc.) were updated to use the `Expression` class instead of GiNaC's `ex` type. This involved:
- Changing the type of local variables from `ex` to `Expression`
- Updating method implementations to work with the `Expression` class

## Series Expansion and Coefficient Extraction

One of the key features we added to the abstraction layer was support for series expansion and coefficient extraction, which is used in the `ValidateFunctionSymbolically`, `ValidateGraph`, and `ValidateExpressionFile` methods in the `Subunit` class.

### 1. C++ Implementation

In the C++ implementation, we added the following methods to the `SymbolicInterface`:
- `series`: Computes the series expansion of an expression around a point
- `series_to_poly`: Converts a series expansion to a polynomial
- `coeff`: Extracts the coefficient of a variable raised to a power

These methods are implemented in the `GiNaCExpression` class using GiNaC's corresponding functions.

### 2. Python Implementation

In the Python implementation, we added the same methods to the `SymPyExpression` class:
- `series`: Uses SymPy's `series` method
- `series_to_poly`: Uses SymPy's `removeO` method to remove the order term
- `coeff`: Uses SymPy's `coeff` method

## Testing

To ensure that the abstraction layer works correctly, we created a test program
(`seb/test_expression.cpp`) that exercises the key functionality of the
`Expression` class:
- Basic arithmetic operations
- Mathematical functions
- Substitution and evaluation
- Series expansion and coefficient extraction
- Numeric evaluation and conversion

## Conclusion

The implementation of the symbolic expression abstraction layer for SEB allows the library to work with both GiNaC in C++ and SymPy in Python. This makes SEB more flexible and accessible to users who prefer Python over C++.

The abstraction layer is designed to be extensible, so additional symbolic expression libraries could be supported in the future by implementing the `SymbolicInterface`.
