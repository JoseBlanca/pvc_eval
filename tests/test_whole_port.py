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

import os
import random
import subprocess
from collections import Counter
from collections.abc import Sequence
from pathlib import Path

import pytest

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
_PVC_COPY_FLOORS = "8,6,6,6,5,4"
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
# ordinary sequence picked by position, hold no N. Run on 19
# October 2026, the three together held 1,928 STR segments, 138 STR bundles and 23
# satellites at the defaults, and each slice held each of the three types.
_GRCH38_SLICES = (
    ("chr1", 0, 1_000_000),
    ("chr1", 10_000_000, 11_000_000),
    ("chr20", 1_000_000, 2_000_000),
)


_PVC = ("PVC_EVAL_PVC", "a pvc binary built at bbcc9f4, as this module's docstring says")
_GRCH38 = ("PVC_EVAL_GRCH38", "the GRCh38 FASTA, with its .fai beside it")


def _get_paths_from_environment(*needed: tuple[str, str]) -> list[Path]:
    """The paths the variables name, in their order. Each is given as its name and what it
    names; the test is skipped, naming each that is missing, when a variable is unset or
    its file, or a FASTA's .fai, does not exist."""
    paths = []
    missing = []
    for variable, what in needed:
        value = os.environ.get(variable)
        path = Path(value) if value else None
        if path is None:
            missing.append(f"{variable} is not set; it names {what}")
        elif not path.is_file():
            missing.append(f"{variable} names {path}, which does not exist; it names {what}")
        elif variable == _GRCH38[0] and not Path(f"{path}.fai").is_file():
            missing.append(f"{variable} names {path}, which has no {path.name}.fai beside it")
        else:
            paths.append(path)
    if missing:
        pytest.skip(". ".join(missing))
    return paths


def _write_fasta(path: Path, contigs: dict[str, str]) -> Path:
    """A FASTA of the contigs, one line each, and its .fai, as samtools faidx writes it."""
    fasta = ""
    fai = ""
    for name, bases in contigs.items():
        fasta += f">{name}\n"
        fai += f"{name}\t{len(bases)}\t{len(fasta)}\t{len(bases)}\t{len(bases) + 1}\n"
        fasta += f"{bases}\n"
    path.write_text(fasta, encoding="ascii")
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


def _make_generated_contig(generator: random.Random) -> str:
    """Random bases with one to five planted repeats. Half the gaps between repeats are
    under 26 bases, so that the repeats bundle; three in ten of the flanks at the ends are
    under 31 bases, so that repeats come near a contig end."""
    pieces = []
    for _ in range(generator.randint(1, 5)):
        gap = generator.randint(0, 25) if generator.random() < 0.5 else generator.randint(20, 300)
        pieces += [_make_random_bases(generator, gap), _make_planted_repeat(generator)]
    for at in (0, len(pieces)):
        near = generator.random() < 0.3
        flank = generator.randint(0, 30) if near else generator.randint(30, 300)
        pieces.insert(at, _make_random_bases(generator, flank))
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


def _build_catalogue(pvc: Path, fasta: Path) -> None:
    """pvc's repeat catalogue, beside the FASTA, at its defaults."""
    subprocess.run(
        [pvc, "repeat-catalog", "--reference", fasta, "--quiet"], check=True, capture_output=True
    )


def _read_pvc_rows(pvc: Path, fasta: Path, output: Path, options: Sequence[str]) -> list[str]:
    """The STR regions of `pvc dev type-regions`, as pvc_eval writes them: pvc's types
    renamed, and the motif of an STR segment alone."""
    subprocess.run(
        [pvc, "dev", "type-regions", "--reference", fasta, "--output", output, "--quiet", *options],
        check=True,
        capture_output=True,
    )
    rows = []
    for line in output.read_text(encoding="utf-8").splitlines():
        if line.startswith("#"):
            continue
        contig, start, end, kind, motif = line.split("\t")[:5]
        if kind == "generic":
            continue
        region_type = _PVC_TYPES[kind]
        shown_motif = motif if region_type == "str" else "."
        rows.append(f"{contig}\t{start}\t{end}\t{region_type}\t{shown_motif}")
    return rows


def _read_pvc_eval_rows(fasta: Path, output: Path, options: Sequence[str]) -> list[str]:
    main(
        [
            "str-regions",
            "--reference",
            str(fasta),
            "--min-copies",
            _PVC_COPY_FLOORS,
            "--output",
            str(output),
            *options,
        ]
    )
    lines = output.read_text(encoding="utf-8").splitlines()
    return [line for line in lines if not line.startswith("#")]


def _compare(pvc: Path, fasta: Path, options: Sequence[str]) -> Counter[str]:
    """Check that both programs give the same STR regions of the FASTA, and count them by
    type."""
    pvc_rows = _read_pvc_rows(pvc, fasta, fasta.with_suffix(".pvc.bed"), options)
    pvc_eval_rows = _read_pvc_eval_rows(fasta, fasta.with_suffix(".pvc_eval.bed"), options)
    assert pvc_eval_rows == pvc_rows
    return Counter(row.split("\t")[3] for row in pvc_rows)


# A comparison input: the pvc binary, and a FASTA with pvc's catalogue beside it.
type _Input = tuple[Path, Path]


@pytest.fixture(scope="module")
def generated_input(tmp_path_factory: pytest.TempPathFactory) -> _Input:
    (pvc,) = _get_paths_from_environment(_PVC)
    directory = tmp_path_factory.mktemp("generated")
    fasta = _write_fasta(directory / "generated.fa", _make_generated_contigs())
    _build_catalogue(pvc, fasta)
    return pvc, fasta


@pytest.fixture(scope="module")
def grch38_input(tmp_path_factory: pytest.TempPathFactory) -> _Input:
    pvc, grch38 = _get_paths_from_environment(_PVC, _GRCH38)
    directory = tmp_path_factory.mktemp("grch38")
    fasta = _write_fasta(directory / "slices.fa", _make_grch38_slices(grch38))
    _build_catalogue(pvc, fasta)
    return pvc, fasta


def test_generated_contigs_have_candidates_near_a_contig_end() -> None:
    # The plan's "What could go wrong": repeats that never come near a contig end would
    # leave step 2 untested. 486 candidates that reach pvc's calling floors were dropped at
    # step 2 on 19 October 2026. That the repeats bundle and make satellites is asserted
    # where the STR regions are compared.
    parameters = StrFinderParameters(copy_floors=(8, 6, 6, 6, 5, 4))
    dropped_near_an_end = 0
    for bases in _make_generated_contigs().values():
        candidates = _find_candidates(bases, parameters)
        kept = list(_drop_near_contig_ends(candidates, bases, parameters))
        dropped_near_an_end += len(candidates) - len(kept)
    assert dropped_near_an_end > 100


# The STR regions compared on 19 October 2026 (STR segments, STR bundles, satellites): on
# the generated contigs, 568, 41 and 224 at the defaults, 338, 12 and 213 with periods 2 to
# 5; on the slices of GRCh38, 1,928, 138 and 23, and 291, 16 and 23.
@pytest.mark.parametrize("setting", _SETTINGS)
def test_generated_contigs_give_the_str_regions_of_pvc(
    generated_input: _Input, setting: str
) -> None:
    counts = _compare(*generated_input, _SETTINGS[setting])
    assert counts["str"] > 0
    assert counts["str-bundle"] > 0
    assert counts["satellite"] > 0


@pytest.mark.parametrize("setting", _SETTINGS)
def test_grch38_slices_give_the_str_regions_of_pvc(grch38_input: _Input, setting: str) -> None:
    counts = _compare(*grch38_input, _SETTINGS[setting])
    assert counts["str"] > 0
    assert counts["str-bundle"] > 0
    assert counts["satellite"] > 0
