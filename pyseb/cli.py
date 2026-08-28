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


def world_main(argv: list[str] | None = None) -> int:
    """Validate or inspect a .pyseb-world.yaml file."""
    parser = argparse.ArgumentParser(prog="pyseb-world")
    commands = parser.add_subparsers(dest="command", required=True)
    for name, help_text in (
        ("validate", "validate and construct a world"),
        ("inspect", "print the canonical world representation"),
    ):
        command = commands.add_parser(name, help=help_text)
        command.add_argument("file")
        command.add_argument(
            "--model-file", action="append", default=[],
            help="register an external .pyseb.yaml model (repeatable)",
        )

    args = parser.parse_args(argv)
    try:
        import pyseb

        world = pyseb.load_world(args.file, model_files=args.model_file)
        if args.command == "inspect":
            print(world.to_yaml(), end="")
        else:
            print(f"valid: {args.file}")
        return 0
    except Exception as error:
        print(f"error: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(subunit_main())
