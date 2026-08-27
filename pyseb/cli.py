"""Command-line tools for file-defined pySEB subunits."""

from __future__ import annotations

import argparse
import sys


def subunit_main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="pyseb-subunit")
    commands = parser.add_subparsers(dest="command", required=True)
    validate = commands.add_parser("validate", help="validate a .pyseb.yaml model")
    validate.add_argument("file")
    args = parser.parse_args(argv)

    if args.command == "validate":
        try:
            from . import validate_subunit_file

            report = validate_subunit_file(args.file)
        except Exception as error:
            print(f"error: {error}", file=sys.stderr)
            return 2

        for warning in report.warnings:
            print(f"warning: {warning}")
        for failure in report.failures:
            print(f"failed: {failure.case_name}: {failure.message}", file=sys.stderr)
        if not report.ok:
            return 1
        print(f"valid: {report.model_id} ({report.case_count} reference case(s))")
        return 0

    return 2


if __name__ == "__main__":
    raise SystemExit(subunit_main())
