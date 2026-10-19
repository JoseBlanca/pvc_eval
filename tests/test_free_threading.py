import sys

from pvc_eval import fasta, segments, vcf  # noqa: F401 -- imported for the GIL check, not used


def test_importing_pvc_eval_leaves_the_gil_off() -> None:
    assert not sys._is_gil_enabled()
