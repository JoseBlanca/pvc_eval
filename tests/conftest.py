"""What the tests of reconciliation and of the re-spelling program share: CaseFiles, which
writes a case's reference and VCF lines to a directory and reads them back with parse_vcf,
its fixture; write_case_vcf, which writes a case's lines as a VCF file, with
make_chrom_line and FILEFORMAT_LINE; and the helpers that build a case's sample pairs and
read its locus's bases; and EnvironmentPath and get_paths_from_environment, which give
the paths a test reads from the environment, or skip it.
"""

import os
from collections.abc import Iterator, Mapping, Sequence
from contextlib import ExitStack
from dataclasses import dataclass
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
        return self.write_contigs({CONTIG: contig})

    def write_contigs(self, contigs: Mapping[str, str]) -> IndexedFasta:
        """A FASTA of these contigs, by name, each on one line, with its .fai."""
        path = self._make_path("fa")
        fasta_lines = []
        fai_lines = []
        offset = 0
        for name, bases in contigs.items():
            offset += len(name) + 2  # `>name` and its newline
            fai_lines.append(f"{name}\t{len(bases)}\t{offset}\t{len(bases)}\t{len(bases) + 1}\n")
            fasta_lines.append(f">{name}\n{bases}\n")
            offset += len(bases) + 1
        path.write_text("".join(fasta_lines), encoding="utf-8")
        Path(f"{path}.fai").write_text("".join(fai_lines), encoding="utf-8")
        return self._stack.enter_context(IndexedFasta(path))

    def write_and_parse_records(
        self,
        reference: IndexedFasta,
        lines: Sequence[str],
        *,
        num_samples: int = 1,
        ploidy: int = 2,
        span_limit: int = DEFAULT_SPAN_LIMIT,
        has_chrom: bool = False,
    ) -> tuple[VcfRecord, ...]:
        """The records parse_vcf reads from lines `POS REF ALT GT...`, one GT per sample,
        each written `GT` or `GT:PS`, on CONTIG; with has_chrom, the lines are
        `CHROM POS REF ALT GT...`."""
        path = self.write_vcf(lines, num_samples=num_samples, has_chrom=has_chrom)
        return tuple(parse_vcf(path, reference, span_limit=span_limit, ploidy=ploidy))

    def write_vcf(
        self, lines: Sequence[str], *, num_samples: int = 1, has_chrom: bool = False
    ) -> Path:
        """A VCF file of lines `POS REF ALT GT...`, as write_and_parse_records reads them."""
        return write_case_vcf(
            self._make_path("vcf"), lines, num_samples=num_samples, has_chrom=has_chrom
        )

    def _make_path(self, suffix: str) -> Path:
        self._num_files += 1
        return self._directory / f"file{self._num_files}.{suffix}"


@pytest.fixture
def case_files(tmp_path: Path) -> Iterator[CaseFiles]:
    with ExitStack() as stack:
        yield CaseFiles(tmp_path, stack)


FILEFORMAT_LINE = "##fileformat=VCFv4.3"


def make_chrom_line(num_samples: int) -> str:
    """The #CHROM line of a VCF with samples S0, S1..., and no FORMAT column with none."""
    columns = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
    if num_samples:
        columns += ["FORMAT", *(f"S{index}" for index in range(num_samples))]
    return "\t".join(columns)


def write_case_vcf(
    path: Path,
    lines: Sequence[str],
    *,
    num_samples: int | None = None,
    has_chrom: bool = False,
    header_lines: Sequence[str] = (FILEFORMAT_LINE,),
    fills_ps: bool = True,
) -> Path:
    """A VCF file at path of a case's lines, `POS REF ALT sample...`, on CONTIG, or
    `CHROM POS REF ALT sample...` with has_chrom; ID, QUAL, FILTER and INFO are `.`.

    num_samples: the samples of the #CHROM line, by default as many as the first line has.
    fills_ps: FORMAT is GT:PS and a sample field with no PS gets `.`; otherwise FORMAT is
        GT:PS when some field of the line has a PS and GT when none has, and the fields
        are written as given, as the cases of doc/specs/respell.md are.
    """
    if num_samples is None:
        num_samples = len(lines[0].split()) - (4 if has_chrom else 3)
    text_lines = [*header_lines, make_chrom_line(num_samples)]
    for line in lines:
        fields = line.split()
        chrom = fields.pop(0) if has_chrom else CONTIG
        pos, ref, alt, *sample_fields = fields
        columns = [chrom, pos, ".", ref, alt, ".", ".", "."]
        if num_samples:
            has_ps = fills_ps or any(":" in sample_field for sample_field in sample_fields)
            if fills_ps:
                sample_fields = [field if ":" in field else f"{field}:." for field in sample_fields]
            columns += ["GT:PS" if has_ps else "GT", *sample_fields]
        text_lines.append("\t".join(columns))
    path.write_text("".join(f"{text_line}\n" for text_line in text_lines), encoding="utf-8")
    return path


def get_locus_bases(reference: IndexedFasta, locus: Locus) -> str:
    return reference.get_bases(locus.contig, locus.start, locus.end)


def pair_samples(num_samples: int) -> tuple[SamplePair, ...]:
    return tuple(
        SamplePair(name_a=f"S{index}", index_a=index, name_b=f"S{index}", index_b=index)
        for index in range(num_samples)
    )


@dataclass(frozen=True, slots=True)
class EnvironmentPath:
    """A path the test reads from the environment: its variable, what it names, and the
    suffixes of the files that must be beside it."""

    variable: str
    what: str
    beside: tuple[str, ...] = ()


def get_paths_from_environment(*needed: EnvironmentPath) -> list[Path]:
    """The paths the variables name, in their order. The test is skipped, naming each that
    is missing, when a variable is unset, or its file or a file beside it does not
    exist."""
    paths = []
    missing = []
    for environment_path in needed:
        variable, what = environment_path.variable, environment_path.what
        value = os.environ.get(variable)
        if not value:
            missing.append(f"{variable} is not set; it names {what}")
            continue
        path = Path(value)
        absent = [
            f"{path.name}{suffix}"
            for suffix in environment_path.beside
            if not Path(f"{path}{suffix}").is_file()
        ]
        if not path.is_file():
            missing.append(f"{variable} names {path}, which does not exist; it names {what}")
        elif absent:
            missing.append(f"{variable} names {path}, which has no {', '.join(absent)} beside it")
        else:
            paths.append(path)
    if missing:
        pytest.skip(". ".join(missing))
    return paths
