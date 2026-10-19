"""Read every VCF record parse_vcf gives for a file and print how many there were: the
workload of the benchmarks of the VCF reader, timed as a command.

    uv run --no-sync python benchmarks/run_parse_vcf.py VCF FASTA [SPAN_LIMIT] [--threads N]

SPAN_LIMIT is a whole number of 1 or more, or `none` for no limit; without it,
parse_vcf's default, DEFAULT_SPAN_LIMIT. N is parse_vcf's num_threads; without it,
parse_vcf's default, DEFAULT_NUM_THREADS.
"""

import argparse
from pathlib import Path

from pvc_eval.fasta import IndexedFasta
from pvc_eval.vcf import DEFAULT_NUM_THREADS, DEFAULT_SPAN_LIMIT, parse_vcf


def main() -> None:
    args = _parse_args()
    with IndexedFasta(args.fasta) as reference:
        num_records = _count_records(args.vcf, reference, args.span_limit, args.threads)
    print(num_records)


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("vcf", type=Path, help="the VCF, plain or bgzipped")
    parser.add_argument("fasta", type=Path, help="the reference, with its .fai beside it")
    parser.add_argument(
        "span_limit",
        nargs="?",
        type=_parse_span_limit,
        default=DEFAULT_SPAN_LIMIT,
        help=f"a whole number, or none (default {DEFAULT_SPAN_LIMIT})",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=DEFAULT_NUM_THREADS,
        help=f"the threads that parse the samples (default {DEFAULT_NUM_THREADS})",
    )
    return parser.parse_args()


def _parse_span_limit(text: str) -> int | None:
    return None if text == "none" else int(text)


def _count_records(
    vcf: Path, reference: IndexedFasta, span_limit: int | None, num_threads: int
) -> int:
    return sum(1 for _ in parse_vcf(vcf, reference, span_limit=span_limit, num_threads=num_threads))


if __name__ == "__main__":
    main()
