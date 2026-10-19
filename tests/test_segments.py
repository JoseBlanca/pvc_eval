import pytest

from pvc_eval.segments import GenomeSegment, RegionType


def test_a_segment_holds_its_contig_span_type_and_motif() -> None:
    segment = GenomeSegment("chr1", 58, 80, RegionType.STR, "CA")
    assert (segment.contig, segment.start, segment.end) == ("chr1", 58, 80)
    assert segment.region_type is RegionType.STR
    assert segment.motif == "CA"


def test_region_types_are_the_four_of_the_contract() -> None:
    assert [region_type.value for region_type in RegionType] == [
        "generic",
        "str",
        "str-bundle",
        "satellite",
    ]


@pytest.mark.parametrize(("start", "end"), [(10, 10), (10, 9)], ids=["empty", "reversed"])
def test_a_segment_whose_end_is_not_after_its_start_raises(start: int, end: int) -> None:
    with pytest.raises(ValueError, match="end"):
        GenomeSegment("chr1", start, end, RegionType.GENERIC, None)
