from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


class TestExamples(unittest.TestCase):
    def test_sympy_example_runs_when_launched_by_absolute_path(self):
        repo_root = Path(__file__).resolve().parents[3]
        script = repo_root / "pyseb" / "examples" / "python" / "sympy_example.py"

        with tempfile.TemporaryDirectory() as tmpdir:
            result = subprocess.run(
                [sys.executable, str(script)],
                cwd=tmpdir,
                text=True,
                capture_output=True,
                timeout=30,
            )

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("Form Factor Expression:", result.stdout)

    def test_python_examples_run_when_launched_by_absolute_path(self):
        repo_root = Path(__file__).resolve().parents[3]
        examples_dir = repo_root / "pyseb" / "examples" / "python"

        expected_output = {
            "backend_simplification_comparison.py": "Backend simplification comparison",
            "diblock_evaluation.py": "Diblock copolymer form factor",
            "figure13_star_chain.py": "Figure 13 star-chain form factor",
            "micelle_evaluation.py": "Micelle form factor",
            "nested_structures.py": "Nested diblock-star chain form factor",
            "symbolic_integral.py": "Solid cylinder symbolic integral",
            "symbolic_backend.py": "Available backends:",
        }

        with tempfile.TemporaryDirectory() as tmpdir:
            for filename, expected in expected_output.items():
                result = subprocess.run(
                    [sys.executable, str(examples_dir / filename)],
                    cwd=tmpdir,
                    text=True,
                    capture_output=True,
                    timeout=30,
                )

                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertIn(expected, result.stdout)


if __name__ == "__main__":
    unittest.main()
