"""Where the peak memory of the workload of run_parse_vcf.py comes from: the peak bytes
tracemalloc traced, and the source lines that held the most when the bytes held between
two records were at their highest.

    uv run --no-sync python benchmarks/peak_parse_vcf.py VCF FASTA [SPAN_LIMIT]

The peak inside a call, while a record is made, can be above the highest bytes held
between two records; both are printed. tracemalloc slows the run, so its time is not a
measurement.
"""

import argparse
import tracemalloc
from pathlib import Path

from pvc_eval.fasta import IndexedFasta
from pvc_eval.vcf import DEFAULT_SPAN_LIMIT, parse_vcf

_NUM_TOP_LINES = 15
_NUM_FRAMES = 10
# A new snapshot is taken only when the bytes held grow by this share over the last one,
# so that the snapshots, which are slow, are few.
_SNAPSHOT_GROWTH = 1.02
_MIB = 1024 * 1024


def main() -> None:
    args = _parse_args()
    tracemalloc.start(_NUM_FRAMES)
    with IndexedFasta(args.fasta) as reference:
        num_records, snapshot, snapshot_bytes = _find_largest_snapshot(
            args.vcf, reference, args.span_limit
        )
    _, peak = tracemalloc.get_traced_memory()
    print(f"records: {num_records}")
    print(f"peak traced: {peak / _MIB:.1f} MiB")
    print(f"highest held between two records: {snapshot_bytes / _MIB:.1f} MiB")
    print(f"top {_NUM_TOP_LINES} lines then:")
    for statistic in snapshot.statistics("lineno")[:_NUM_TOP_LINES]:
        print(statistic)


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
    return parser.parse_args()


def _parse_span_limit(text: str) -> int | None:
    return None if text == "none" else int(text)


def _find_largest_snapshot(
    vcf: Path, reference: IndexedFasta, span_limit: int | None
) -> tuple[int, tracemalloc.Snapshot, int]:
    """The number of records, and a snapshot taken when the bytes held between two
    records were within _SNAPSHOT_GROWTH of their highest, with those bytes."""
    snapshot = tracemalloc.take_snapshot()
    snapshot_bytes, _ = tracemalloc.get_traced_memory()
    num_records = 0
    for _ in parse_vcf(vcf, reference, span_limit=span_limit):
        num_records += 1
        held, _ = tracemalloc.get_traced_memory()
        if held > snapshot_bytes * _SNAPSHOT_GROWTH:
            snapshot = tracemalloc.take_snapshot()
            snapshot_bytes = held
    return num_records, snapshot, snapshot_bytes


if __name__ == "__main__":
    main()
