import sys

import pvc_eval.vcf


def test_importing_pvc_eval_leaves_the_gil_off() -> None:
    assert pvc_eval.vcf is not None
    assert not sys._is_gil_enabled()
