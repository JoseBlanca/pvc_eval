import sys

# Imported for the GIL check, not used: an import that turns the GIL on fails the test.
import pvc_eval.fasta
import pvc_eval.vcf  # noqa: F401 -- see the comment above


def test_importing_pvc_eval_leaves_the_gil_off() -> None:
    assert not sys._is_gil_enabled()
