"""The time of each step of the STR finder on every contig of a FASTA, as _find_segments
runs them, and the counts that drive each: the candidates, the bundle members, the STR
segments, the satellites before and after absorption, and the passes of absorption.

    uv run --no-sync python benchmarks/stages_str_finder.py FASTA

It runs the finder's private functions in the order _find_segments calls them, at the
default parameters, and prints one line per contig and step.
"""

import sys
import time
from collections.abc import Sequence
from pathlib import Path

from pvc_eval import segments
from pvc_eval.fasta import IndexedFasta
from pvc_eval.segments import StrFinderParameters


def _count_passes() -> list[int]:
    """Wrap _make_satellites_absorbing to count its calls; the list holds the count."""
    calls = [0]
    original = segments._make_satellites_absorbing

    def counted(
        satellites: Sequence[tuple[int, int]], span: tuple[int, int], bundle_distance: int
    ) -> list[tuple[int, int]] | None:
        calls[0] += 1
        return original(satellites, span, bundle_distance)

    segments._make_satellites_absorbing = counted  # ty: ignore[invalid-assignment]  # a counting wrapper of the same signature
    return calls


def main() -> None:
    parameters = StrFinderParameters()
    calls = _count_passes()
    with IndexedFasta(Path(sys.argv[1])) as reference:
        for contig in reference.contigs:
            started = time.perf_counter()
            bases = reference.get_bases(contig, 0, reference.get_length(contig))
            read = time.perf_counter()
            found = segments._find_candidates(bases, parameters)
            step1 = time.perf_counter()
            candidates = list(segments._drop_near_contig_ends(found, bases, parameters))
            step2 = time.perf_counter()
            classification = segments._find_str_segments(candidates, bases, parameters)
            steps4to6 = time.perf_counter()
            runs = [
                run
                for run in segments._make_repeat_runs(candidates)
                if run[1] - run[0] > parameters.max_str_len
            ]
            features = segments._find_features(candidates, classification, bases, parameters)
            step7 = time.perf_counter()
            segments._make_segments(features, len(bases))
            step8 = time.perf_counter()
            num_features = len(classification.bundle_members) + len(classification.str_segments)
            num_satellites = sum(
                1 for feature in features if feature.region_type.value == "satellite"
            )
            print(f"{contig}: {len(bases):,} bases")
            print(f"  read the bases        {read - started:8.2f} s")
            print(f"  step 1 and 3          {step1 - read:8.2f} s  {len(found):,} candidates")
            print(f"  step 2                {step2 - step1:8.2f} s  {len(candidates):,} left")
            print(
                f"  steps 4 to 6          {steps4to6 - step2:8.2f} s"
                f"  {len(classification.str_segments):,} STR segments,"
                f" {len(classification.bundle_members):,} bundle members"
            )
            print(
                f"  step 7                {step7 - steps4to6:8.2f} s"
                f"  {len(runs):,} satellites from repeat runs, {num_satellites:,} after;"
                f" {calls[0]:,} absorption calls for {num_features:,} features"
            )
            print(f"  step 8                {step8 - step7:8.2f} s")
            calls[0] = 0


if __name__ == "__main__":
    main()
