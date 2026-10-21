"""Write stretches of a reference as the contigs of a new FASTA, with its .fai, for the
workloads of the STR finder (benchmarks/README.md).

    uv run python benchmarks/make_str_slices.py REFERENCE OUTPUT CONTIG:START-END...

START and END are 0-based and half-open, as in a BED. Each stretch becomes a contig named
CONTIG_START_END, its bases on one line, as tests/test_whole_port.py writes its slices.
"""

import argparse
from pathlib import Path

from pvc_eval.fasta import IndexedFasta


def _parse_stretch(text: str) -> tuple[str, int, int]:
    contig, _, span = text.rpartition(":")
    start, _, end = span.replace(",", "").partition("-")
    return contig, int(start), int(end)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("stretches", nargs="+", type=_parse_stretch)
    arguments = parser.parse_args()
    fai_lines = []
    offset = 0
    with IndexedFasta(arguments.reference) as reference, arguments.output.open("w") as fasta:
        for contig, start, end in arguments.stretches:
            name = f"{contig}_{start}_{end}"
            bases = reference.get_bases(contig, start, end)
            header = f">{name}\n"
            offset += len(header)
            fai_lines.append(f"{name}\t{len(bases)}\t{offset}\t{len(bases)}\t{len(bases) + 1}\n")
            fasta.write(header)
            fasta.write(f"{bases}\n")
            offset += len(bases) + 1
    Path(f"{arguments.output}.fai").write_text("".join(fai_lines), encoding="utf-8")


if __name__ == "__main__":
    main()
