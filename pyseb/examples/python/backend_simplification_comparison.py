#!/usr/bin/env python3
"""Compare SymPy simplification of expressions from available C++ backends."""

from pathlib import Path
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

import numpy as np
import pyseb
import sympy


def build_figure13_expression():
    world = pyseb.World()

    diblock = world.Add("GaussianPolymer", "A")
    world.Link("GaussianPolymer", "B.end1", "A.end2")

    star = world.Add(diblock, "diblock1")
    world.Link(diblock, "diblock2:A.end1", "diblock1:A.end1")
    world.Link(diblock, "diblock3:A.end1", "diblock1:A.end1")
    world.Link(diblock, "diblock4:A.end1", "diblock1:A.end1")

    chain = world.Add(star, "star1")
    world.Link(star, "star2:diblock1:B.end2", "star1:diblock3:B.end2")
    world.Link(star, "star3:diblock1:B.end2", "star2:diblock3:B.end2")
    world.Link(star, "star4:diblock1:B.end2", "star3:diblock3:B.end2")
    world.Link(star, "star5:diblock1:B.end2", "star4:diblock3:B.end2")
    world.Add(chain, "chain")

    return pyseb.to_sympy(world.FormFactor("chain"))


def simplify_with_stats(expr):
    start = time.perf_counter()
    simplified = sympy.factor(sympy.cancel(expr))
    elapsed = time.perf_counter() - start
    return {
        "expr": simplified,
        "elapsed": elapsed,
        "original_ops": sympy.count_ops(expr),
        "simplified_ops": sympy.count_ops(simplified),
        "original_chars": len(str(expr)),
        "simplified_chars": len(str(simplified)),
    }


def evaluate(expr, q_values):
    params = {
        "Rg_A": 1.0,
        "Rg_B": 1.0,
        "beta_A": 1.0,
        "beta_B": 1.0,
    }
    q = sympy.Symbol("q")
    fn = sympy.lambdify(q, expr.subs(params), modules=["numpy"])
    return np.asarray(fn(q_values), dtype=float)


def main():
    print("Backend simplification comparison")
    print("Available backends:", ", ".join(pyseb.available_backends()))

    pyseb.set_backend("portable")
    portable_expr = build_figure13_expression()
    portable_stats = simplify_with_stats(portable_expr)

    print("\nportable -> SymPy")
    print(f"operations: {portable_stats['original_ops']} -> {portable_stats['simplified_ops']}")
    print(f"characters: {portable_stats['original_chars']} -> {portable_stats['simplified_chars']}")
    print(f"simplification seconds: {portable_stats['elapsed']:.3f}")

    if "ginac" not in pyseb.available_backends():
        print("\nGiNaC backend is not available in this build; skipping backend comparison.")
        return

    pyseb.set_backend("ginac")
    ginac_expr = build_figure13_expression()
    ginac_stats = simplify_with_stats(ginac_expr)

    print("\nginac -> SymPy")
    print(f"operations: {ginac_stats['original_ops']} -> {ginac_stats['simplified_ops']}")
    print(f"characters: {ginac_stats['original_chars']} -> {ginac_stats['simplified_chars']}")
    print(f"simplification seconds: {ginac_stats['elapsed']:.3f}")

    symbolic_difference = sympy.simplify(portable_stats["expr"] - ginac_stats["expr"])
    print("\nSymbolic equivalence after simplification:", symbolic_difference == 0)

    q_values = np.logspace(-2, 1, 12)
    portable_values = evaluate(portable_stats["expr"], q_values)
    ginac_values = evaluate(ginac_stats["expr"], q_values)
    max_abs_delta = float(np.max(np.abs(portable_values - ginac_values)))
    print(f"max |portable - ginac| over q grid: {max_abs_delta:.3e}")


if __name__ == "__main__":
    main()
