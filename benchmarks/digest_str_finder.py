"""The number and a sha256 of every segment segment_contig gives on each contig of a FASTA,
generic regions included, at the default parameters, so that an experiment can show it
changed no segment (the BED of str-regions leaves the generic regions out).

    uv run --no-sync python benchmarks/digest_str_finder.py FASTA
"""

import hashlib
import sys
from pathlib import Path

from pvc_eval.fasta import IndexedFasta
from pvc_eval.segments import segment_contig


def main() -> None:
    digest = hashlib.sha256()
    num_segments = 0
    with IndexedFasta(Path(sys.argv[1])) as reference:
        for contig in reference.contigs:
            for segment in segment_contig(reference, contig):
                digest.update(repr(segment).encode())
                num_segments += 1
    print(f"{num_segments} segments {digest.hexdigest()}")


if __name__ == "__main__":
    main()
