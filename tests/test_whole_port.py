"""The whole-port check of doc/specs/str_finder.md: pvc_eval's STR regions equal those of
`pvc dev type-regions` on the same sequences, on 600 generated contigs and on three slices
of GRCh38, at pvc's defaults and with --min-period 2 --max-period 5.

It needs two paths from the environment, and a part whose path is unset or missing is
skipped with a message that names it:

- PVC_EVAL_PVC, a pvc binary built at bbcc9f4, the commit the port follows. The pvc of
  pvc's target/release lacks `pvc dev`, so it is built from an export of that commit,
  which leaves pvc's repository untouched:

      git -C <pvc repository> archive bbcc9f4 | tar -x -C <scratch directory>
      cd <scratch directory> && cargo build --release --bin pvc

  and PVC_EVAL_PVC is <scratch directory>/target/release/pvc.
- PVC_EVAL_GRCH38, the GRCh38 FASTA, with its .fai beside it, for the slices alone. The
  slices are written into the test's own directory; nothing is written beside it.

For each input the test writes a FASTA and its .fai, builds pvc's catalogue beside it with
`pvc repeat-catalog`, runs `pvc dev type-regions` and `pvc_eval str-regions --min-copies
8,6,6,6,5,4`, pvc's calling floors, and compares their rows that are not generic: contig,
start, end, type, and the motif of an STR segment.
"""

import random
import subprocess
from collections import Counter
from collections.abc import Sequence
from pathlib import Path
from typing import NamedTuple

import pytest
from tests.conftest import EnvironmentPath, get_paths_from_environment

from pvc_eval.cli import main
from pvc_eval.fasta import IndexedFasta
from pvc_eval.segments import (
    StrFinderParameters,
    _drop_near_contig_ends,
    _find_candidates,
)

pytestmark = pytest.mark.slow

# pvc's calling floors, the defaults of `pvc dev type-regions` (doc/specs/str_finder.md,
# "The whole-port check").
_PVC_COPY_FLOORS = (8, 6, 6, 6, 5, 4)
_SETTINGS = {
    "defaults": (),
    "periods_2_to_5": ("--min-period", "2", "--max-period", "5"),
}
# pvc's region types and pvc_eval's (doc/specs/str_finder.md, "The whole-port check").
_PVC_TYPES = {"ssr_locus": "str", "ssr_bundle": "str-bundle", "satellite": "satellite"}

_GENERATED_SEED = 20261019
_GENERATED_CONTIGS = 600

# The stretches of GRCh38 written as contigs of their own, 0-based and half-open. The
# first holds chr1's first 10,000 bases of N, its telomere of (CCCTAA) right after them,
# and three gaps of 50,000 N, so that runs of N meet repeats; the other two, stretches of
# ordinary sequence picked by position, hold no N. On 19 October 2026 each slice held each
# of the three types of STR region.
_GRCH38_SLICES = (
    ("chr1", 0, 1_000_000),
    ("chr1", 10_000_000, 11_000_000),
    ("chr20", 1_000_000, 2_000_000),
)


_PVC = EnvironmentPath(
    "PVC_EVAL_PVC", "a pvc binary built at bbcc9f4, as this module's docstring says"
)
_GRCH38 = EnvironmentPath(
    "PVC_EVAL_GRCH38", "the GRCh38 FASTA, with its .fai beside it", beside=(".fai",)
)


def _write_fasta(path: Path, contigs: dict[str, str]) -> Path:
    """A FASTA of the contigs, one line each, and its .fai, as samtools faidx writes it."""
    fasta = ""
    fai = ""
    for name, bases in contigs.items():
        fasta += f">{name}\n"
        fai += f"{name}\t{len(bases)}\t{len(fasta)}\t{len(bases)}\t{len(bases) + 1}\n"
        fasta += f"{bases}\n"
    path.write_text(fasta, encoding="utf-8")
    Path(f"{path}.fai").write_text(fai, encoding="utf-8")
    return path


def _make_random_bases(generator: random.Random, length: int) -> str:
    return "".join(generator.choices("ACGT", k=length))


def _make_planted_repeat(generator: random.Random) -> str:
    """A repeat of period 1 to 6 with a random motif: one in five longer than 100 bases, so
    that it is a satellite; two in five impure, each base drawn again at random with a
    chance of 0.08."""
    period = generator.randint(1, 6)
    motif = _make_random_bases(generator, period)
    if generator.random() < 0.2:
        length = generator.randint(101, 300)
    else:
        length = generator.randint(2 * period, 60)
    bases = list((motif * (length // period + 1))[:length])
    if generator.random() < 0.4:
        for position in range(length):
            if generator.random() < 0.08:
                bases[position] = generator.choice("ACGT")
    return "".join(bases)


def _make_spacer(generator: random.Random, short_chance: float, short_most: int) -> str:
    """Random bases: 0 to short_most of them with a chance of short_chance, and
    short_most + 1 to 300 otherwise."""
    if generator.random() < short_chance:
        length = generator.randint(0, short_most)
    else:
        length = generator.randint(short_most + 1, 300)
    return _make_random_bases(generator, length)


def _make_generated_contig(generator: random.Random) -> str:
    """One to five planted repeats, with random bases between them and at both ends. Each
    gap between two repeats is 0 to 25 bases with a chance of 0.5, so that the repeats
    bundle, and 26 to 300 otherwise. The flank at each end is 0 to 30 bases with a chance
    of 0.3, so that a repeat comes near that end of the contig, and 31 to 300 otherwise."""
    pieces = [_make_spacer(generator, short_chance=0.3, short_most=30)]
    for index in range(generator.randint(1, 5)):
        if index:
            pieces.append(_make_spacer(generator, short_chance=0.5, short_most=25))
        pieces.append(_make_planted_repeat(generator))
    pieces.append(_make_spacer(generator, short_chance=0.3, short_most=30))
    return "".join(pieces)


def _make_generated_contigs() -> dict[str, str]:
    generator = random.Random(_GENERATED_SEED)
    return {f"g{index}": _make_generated_contig(generator) for index in range(_GENERATED_CONTIGS)}


def _make_grch38_slices(grch38: Path) -> dict[str, str]:
    with IndexedFasta(grch38) as reference:
        return {
            f"{contig}_{start}_{end}": reference.get_bases(contig, start, end)
            for contig, start, end in _GRCH38_SLICES
        }


def _run_command(command: Sequence[str | Path]) -> None:
    """Run a program; the test fails with the command, its exit status and what it wrote
    to stderr when the status is not 0."""
    finished = subprocess.run(command, capture_output=True, text=True, check=False)
    if finished.returncode != 0:
        shown = " ".join(map(str, command))
        pytest.fail(f"{shown} exited with status {finished.returncode}:\n{finished.stderr}")


def _make_catalogue(pvc: Path, fasta: Path) -> None:
    """pvc's repeat catalogue, beside the FASTA, at its defaults."""
    _run_command([pvc, "repeat-catalog", "--reference", fasta, "--quiet"])


def _find_pvc_str_regions(pvc: Path, fasta: Path, options: Sequence[str]) -> list[str]:
    """The STR regions `pvc dev type-regions` gives, as pvc_eval writes them: pvc's types
    renamed, and the motif of an STR segment alone."""
    output = fasta.with_suffix(".pvc.bed")
    _run_command(
        [pvc, "dev", "type-regions", "--reference", fasta, "--output", output, "--quiet", *options]
    )
    rows = []
    for line in output.read_text(encoding="utf-8").splitlines():
        if line.startswith("#"):
            continue
        contig, start, end, pvc_type, motif = line.split("\t")[:5]
        if pvc_type == "generic":
            continue
        region_type = _PVC_TYPES[pvc_type]
        shown_motif = motif if region_type == "str" else "."
        rows.append(f"{contig}\t{start}\t{end}\t{region_type}\t{shown_motif}")
    return rows


def _find_pvc_eval_str_regions(fasta: Path, options: Sequence[str]) -> list[str]:
    """The STR regions `pvc_eval str-regions` gives with pvc's calling floors, as rows of
    its BED."""
    output = fasta.with_suffix(".pvc_eval.bed")
    copy_floors = ",".join(map(str, _PVC_COPY_FLOORS))
    main(
        [
            "str-regions",
            "--reference",
            str(fasta),
            "--min-copies",
            copy_floors,
            "--output",
            str(output),
            *options,
        ]
    )
    lines = output.read_text(encoding="utf-8").splitlines()
    return [line for line in lines if not line.startswith("#")]


def _count_region_types(rows: Sequence[str]) -> Counter[str]:
    return Counter(row.split("\t")[3] for row in rows)


class _Input(NamedTuple):
    """A comparison input: the pvc binary, and a FASTA with pvc's catalogue beside it."""

    pvc: Path
    fasta: Path


@pytest.fixture(scope="module")
def generated_input(tmp_path_factory: pytest.TempPathFactory) -> _Input:
    (pvc,) = get_paths_from_environment(_PVC)
    directory = tmp_path_factory.mktemp("generated")
    fasta = _write_fasta(directory / "generated.fa", _make_generated_contigs())
    _make_catalogue(pvc, fasta)
    return _Input(pvc=pvc, fasta=fasta)


@pytest.fixture(scope="module")
def grch38_input(tmp_path_factory: pytest.TempPathFactory) -> _Input:
    pvc, grch38 = get_paths_from_environment(_PVC, _GRCH38)
    directory = tmp_path_factory.mktemp("grch38")
    fasta = _write_fasta(directory / "slices.fa", _make_grch38_slices(grch38))
    _make_catalogue(pvc, fasta)
    return _Input(pvc=pvc, fasta=fasta)


def test_generated_contigs_have_candidates_near_each_end_of_a_contig() -> None:
    # The plan's "What could go wrong": repeats that never come near a contig end would
    # leave step 2 untested, and step 2 tests each end on its own. A dropped candidate that
    # starts in the first half of its contig is counted as near the start, any other as
    # near the end. On 19 October 2026, 76 were dropped near the start and 41 near the end.
    parameters = StrFinderParameters(copy_floors=_PVC_COPY_FLOORS)
    dropped_near_start = 0
    dropped_near_end = 0
    for bases in _make_generated_contigs().values():
        candidates = _find_candidates(bases, parameters)
        kept = set(_drop_near_contig_ends(candidates, bases, parameters))
        for candidate in candidates:
            if candidate in kept:
                continue
            if candidate.start < len(bases) // 2:
                dropped_near_start += 1
            else:
                dropped_near_end += 1
    assert dropped_near_start >= 30
    assert dropped_near_end >= 30


# The STR regions compared on 19 October 2026 (STR segments, STR bundles, satellites): on
# the generated contigs, 726, 62 and 309 at the defaults, 457, 20 and 294 with periods 2 to
# 5; on the slices of GRCh38, 1,928, 138 and 23, and 291, 16 and 23.
@pytest.mark.parametrize("setting", _SETTINGS)
def test_generated_contigs_give_the_str_regions_of_pvc(
    generated_input: _Input, setting: str
) -> None:
    options = _SETTINGS[setting]
    pvc_rows = _find_pvc_str_regions(generated_input.pvc, generated_input.fasta, options)
    assert _find_pvc_eval_str_regions(generated_input.fasta, options) == pvc_rows
    counts = _count_region_types(pvc_rows)
    assert counts["str"] > 0
    assert counts["str-bundle"] > 0
    assert counts["satellite"] > 0


@pytest.mark.parametrize("setting", _SETTINGS)
def test_grch38_slices_give_the_str_regions_of_pvc(grch38_input: _Input, setting: str) -> None:
    options = _SETTINGS[setting]
    pvc_rows = _find_pvc_str_regions(grch38_input.pvc, grch38_input.fasta, options)
    assert _find_pvc_eval_str_regions(grch38_input.fasta, options) == pvc_rows
    counts = _count_region_types(pvc_rows)
    assert counts["str"] > 0
    assert counts["str-bundle"] > 0
    assert counts["satellite"] > 0
