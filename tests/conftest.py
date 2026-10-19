"""What the tests of reconciliation share: CaseFiles, which writes a case's reference and
VCF lines to a directory and reads them back with parse_vcf, its fixture, and the helpers
that build a case's sample pairs and read its locus's bases.
"""

from collections.abc import Iterator, Sequence
from contextlib import ExitStack
from pathlib import Path

import pytest

from pvc_eval.fasta import IndexedFasta
from pvc_eval.locus import Locus
from pvc_eval.samples import SamplePair
from pvc_eval.vcf import DEFAULT_SPAN_LIMIT, VcfRecord, parse_vcf

CONTIG = "chr1"
# The contig of every test: N, which equals no base of a VCF record, so no record slides
# into it, around the bases a case gives.
CONTIG_LENGTH = 80


class CaseFiles:
    """Writes a case's reference and VCF files to a directory, and parses the VCF files
    back with parse_vcf."""

    def __init__(self, directory: Path, stack: ExitStack) -> None:
        self._directory = directory
        self._stack = stack
        self._num_files = 0

    def write_reference(self, ref_start_1based: int, bases: str) -> IndexedFasta:
        """A FASTA of CONTIG, N but for bases from ref_start_1based on, with its .fai."""
        contig = "N" * (ref_start_1based - 1) + bases
        contig += "N" * (CONTIG_LENGTH - len(contig))
        path = self._make_path("fa")
        path.write_text(f">{CONTIG}\n{contig}\n", encoding="utf-8")
        fai = (
            f"{CONTIG}\t{CONTIG_LENGTH}\t{len(CONTIG) + 2}\t{CONTIG_LENGTH}\t{CONTIG_LENGTH + 1}\n"
        )
        Path(f"{path}.fai").write_text(fai, encoding="utf-8")
        return self._stack.enter_context(IndexedFasta(path))

    def write_and_parse_records(
        self,
        reference: IndexedFasta,
        lines: Sequence[str],
        *,
        num_samples: int = 1,
        ploidy: int = 2,
        span_limit: int = DEFAULT_SPAN_LIMIT,
    ) -> tuple[VcfRecord, ...]:
        """The records parse_vcf reads from lines `POS REF ALT GT...`, one GT per sample,
        each written `GT` or `GT:PS`."""
        samples = "\t".join(f"S{index}" for index in range(num_samples))
        text = "##fileformat=VCFv4.3\n"
        text += f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{samples}\n"
        for line in lines:
            pos, ref, alt, *gts = line.split()
            sample_fields = [gt if ":" in gt else f"{gt}:." for gt in gts]
            text += (
                f"{CONTIG}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\tGT:PS\t{'\t'.join(sample_fields)}\n"
            )
        path = self._make_path("vcf")
        path.write_text(text, encoding="utf-8")
        return tuple(parse_vcf(path, reference, span_limit=span_limit, ploidy=ploidy))

    def _make_path(self, suffix: str) -> Path:
        self._num_files += 1
        return self._directory / f"file{self._num_files}.{suffix}"


@pytest.fixture
def case_files(tmp_path: Path) -> Iterator[CaseFiles]:
    with ExitStack() as stack:
        yield CaseFiles(tmp_path, stack)


def get_locus_bases(reference: IndexedFasta, locus: Locus) -> str:
    return reference.get_bases(locus.contig, locus.start, locus.end)


def pair_samples(num_samples: int) -> tuple[SamplePair, ...]:
    return tuple(
        SamplePair(name_a=f"S{index}", index_a=index, name_b=f"S{index}", index_b=index)
        for index in range(num_samples)
    )
