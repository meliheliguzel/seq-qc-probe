import argparse
import sys

from seq_qc_probe.checker import FASTQChecker
from seq_qc_probe.reporter import QCReporter


def build_parser() -> argparse.ArgumentParser:
    """Build and configure the top-level argument parser and subcommands."""
    parser = argparse.ArgumentParser(
        prog="seq-qc-probe",
        description="Fast pre-flight FASTQ structural integrity validator.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    check_parser = subparsers.add_parser(
        "check",
        help="Validate FASTQ and enforce QC thresholds",
    )
    check_parser.add_argument(
        "-i",
        "--input",
        dest="file",
        type=str,
        required=True,
        help="Path to the FASTQ file to validate (or '-' for standard input).",
    )
    check_parser.add_argument(
        "--min-quality",
        "--min-q",
        dest="min_quality",
        type=float,
        default=0.0,
        help="Minimum acceptable Phred quality score (default: 0.0).",
    )
    check_parser.add_argument(
        "--max-errors",
        type=int,
        default=10,
        help="Maximum number of validation errors to display before aborting (default: 10).",
    )
    check_parser.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="Suppress all standard output; only return exit code.",
    )
    check_parser.add_argument(
        "--json",
        action="store_true",
        help="Output validation results and metrics in JSON format to stdout.",
    )
    check_parser.add_argument(
        "-o",
        "--output-json",
        dest="output_json",
        type=str,
        default=None,
        help="Write validation report in JSON format to the specified file path.",
    )
    return parser


def handle_check(args: argparse.Namespace) -> int:
    """Delegates check command execution to FASTQChecker and QCReporter."""
    checker = FASTQChecker(
        input_path=args.file,
        min_quality=args.min_quality,
        max_errors=args.max_errors,
        quiet=args.quiet,
    )
    result = checker.run()

    if result.exit_code == 2 and result.error_message:
        if not args.quiet:
            sys.stderr.write(f"{result.error_message}\n")
        return 2

    reporter = QCReporter(result)

    if args.output_json:
        reporter.write_json(args.output_json)

    if args.json:
        reporter.print_json()
    elif not args.quiet:
        reporter.print_terminal_summary()

    return result.exit_code


def main(argv: list[str] | None = None) -> None:
    """Main CLI entrypoint."""
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command == "check":
        sys.exit(handle_check(args))
    sys.exit(0)


if __name__ == "__main__":
    main()
