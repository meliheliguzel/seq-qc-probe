import argparse
import sys
from pathlib import Path

from seq_qc_probe.checker import FASTQChecker
from seq_qc_probe.reporter import QCReporter, VCFReporter
from seq_qc_probe.vcf_validator import VCFValidator


def build_parser() -> argparse.ArgumentParser:
    """Build and configure the top-level argument parser and subcommands."""
    parser = argparse.ArgumentParser(
        prog="seq-qc-probe",
        description="Fast pre-flight FASTQ and VCF structural integrity validator.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # 1. FASTQ check subcommand
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

    # 2. VCF vcf-check subcommand
    vcf_parser = subparsers.add_parser(
        "vcf-check",
        help="Validate VCF 4.2 headers, columns, and records",
    )
    vcf_parser.add_argument(
        "-i",
        "--input",
        dest="file",
        type=str,
        required=True,
        help="Path to the VCF file to validate (or '-' for standard input).",
    )
    vcf_parser.add_argument(
        "--max-errors",
        type=int,
        default=10,
        help="Maximum number of validation errors to display before aborting (default: 10).",
    )
    vcf_parser.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="Suppress all standard output; only return exit code.",
    )
    vcf_parser.add_argument(
        "--json",
        action="store_true",
        help="Output validation results and summary in JSON format to stdout.",
    )
    vcf_parser.add_argument(
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


def handle_vcf_check(args: argparse.Namespace) -> int:
    """Delegates vcf-check command execution to VCFValidator and VCFReporter."""
    if args.file != "-":
        input_path = Path(args.file)
        if not input_path.is_file():
            if not args.quiet:
                sys.stderr.write(f"Error: File not found: {args.file}\n")
            return 2

    validator = VCFValidator(
        file_path=args.file,
        max_errors=args.max_errors,
    )
    summary = validator.validate()

    reporter = VCFReporter(
        summary=summary,
        input_file=args.file,
        max_errors=args.max_errors,
    )

    if args.output_json:
        reporter.write_json(args.output_json)

    if args.json:
        reporter.print_json()
    elif not args.quiet:
        reporter.print_terminal_summary()

    return 0 if summary.is_valid else 1


def main(argv: list[str] | None = None) -> None:
    """Main CLI entrypoint."""
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command == "check":
        sys.exit(handle_check(args))
    elif args.command == "vcf-check":
        sys.exit(handle_vcf_check(args))
    sys.exit(0)


if __name__ == "__main__":
    main()
