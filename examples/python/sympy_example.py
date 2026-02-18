#!/usr/bin/env python3
"""
Example demonstrating the use of pySEB with SymPy.
"""

import pyseb
import sympy

def main():
    # Create a world
    world = pyseb.World("SymPyWorld")
    
    # Add a diblock copolymer
    world.Add("diblock", "diblock")
    
    # Get the form factor expression
    form_factor = world.FormFactor("diblock")
    
    print("Form Factor Expression:")
    print(form_factor)
    
    # Convert to LaTeX
    latex_expr = sympy.latex(form_factor.expr)
    print("\nLaTeX Expression:")
    print(latex_expr)
    
    # Convert to C code
    c_code = sympy.ccode(form_factor.expr)
    print("\nC Code:")
    print(c_code)
    
    # Evaluate the form factor at q=0.1
    q_value = 0.1
    # Substitute q with the value
    q_symbol = sympy.Symbol('q')
    evaluated = form_factor.expr.subs(q_symbol, q_value)
    print(f"\nEvaluated at q={q_value}:")
    print(evaluated)

if __name__ == "__main__":
    main()
