"""Print the number of VCF records parse_vcf gives for a file and a sha256 of the repr of
every record, in the order given: the output an experiment compares with its baseline's,
since the count alone would not change if a span or a genotype did.

    uv run --no-sync python benchmarks/digest_parse_vcf.py VCF FASTA [SPAN_LIMIT] [--threads N]

SPAN_LIMIT is a whole number of 1 or more, or `none` for no limit; without it,
parse_vcf's default, DEFAULT_SPAN_LIMIT. N is parse_vcf's num_threads; without it,
parse_vcf's default, DEFAULT_NUM_THREADS. Not timed: the hashing adds to the run.
"""

import argparse
import hashlib
from pathlib import Path

from pvc_eval.fasta import IndexedFasta
from pvc_eval.vcf import DEFAULT_NUM_THREADS, DEFAULT_SPAN_LIMIT, parse_vcf


def main() -> None:
    args = _parse_args()
    digest = hashlib.sha256()
    num_records = 0
    with IndexedFasta(args.fasta) as reference:
        for vcf_record in parse_vcf(
            args.vcf, reference, span_limit=args.span_limit, num_threads=args.threads
        ):
            # No separator: every repr starts with "VcfRecord(", so none is needed.
            digest.update(repr(vcf_record).encode("utf-8"))
            num_records += 1
    print(num_records, digest.hexdigest())


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


if __name__ == "__main__":
    main()
