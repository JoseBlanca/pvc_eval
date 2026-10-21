"""The 50 cases of doc/specs/hard_cases.md §4.5, S1 to O5, one HardCase each, in HARD_CASES.

Every value comes from §4.5, which the owner approved on 21 October 2026: the excerpt and
its offset, the sample, the records of each file, the locus with its kind and STR regions,
the locus alleles, the possible genotypes of each file and the row. The line of §4.5 that
gives the excerpt, each record, the locus, the possible genotypes and the row stands
beside it as a comment; the records, the locus and the STR regions are written in excerpt
positions, the position minus the offset.

The 301 bases of each excerpt were read on 21 October 2026 from the stretch §4.5 gives,
with pvc_eval.fasta.IndexedFasta, in upper case: E12 to E16 from SL4.0,
S_lycopersicum_chromosomes.4.00.fa, and the other 45 from GRCh38,
GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna. The excerpt tests of
tests/test_hard_cases.py check them against those references.

Public: HardCase, Reference, StrRegion1Based, HARD_CASES.
"""

from dataclasses import dataclass
from enum import Enum

from pvc_eval.kind import Kind
from pvc_eval.outcome import InvalidReason, LocusGenotype, Outcome
from pvc_eval.segments import RegionType


class Reference(Enum):
    """The reference an excerpt was read from."""

    GRCH38 = "GRCh38"
    SL4_0 = "SL4.0"


# An STR region in excerpt positions, 1-based and inclusive, with its type and motif.
type StrRegion1Based = tuple[int, int, RegionType, str | None]


@dataclass(frozen=True, slots=True, kw_only=True)
class HardCase:
    """One case of doc/specs/hard_cases.md §4.5.

    name: the case, S1 to O5.
    reference, contig, offset: the reference the excerpt was read from, and its contig;
        an excerpt position, 1-based, is the position in the contig minus offset.
    sample: the sample of both files.
    bases: the excerpt, 301 bases.
    lines_a, lines_b: the records of each file, `POS REF ALT GT`, a GT perhaps with its
        `:PS`, POS in excerpt positions.
    locus_1based: the locus in excerpt positions, 1-based and inclusive.
    kind, str_regions_1based: the locus's kind and its STR regions.
    alleles: the locus alleles, REF first.
    possible_genotypes_a, possible_genotypes_b: the possible genotypes of each file over
        the locus, each as indices into alleles, sorted; None for a sample that is
        invalid.
    outcome, genotype_a, genotype_b, possible_a, possible_b, records_a, records_b,
        invalid_reason_a, invalid_reason_b: the row, in the fields of OutcomeRow.
    """

    name: str
    reference: Reference
    contig: str
    offset: int
    sample: str
    bases: str
    lines_a: tuple[str, ...]
    lines_b: tuple[str, ...]
    locus_1based: tuple[int, int]
    kind: Kind
    str_regions_1based: tuple[StrRegion1Based, ...]
    alleles: tuple[str, ...]
    possible_genotypes_a: frozenset[LocusGenotype] | None
    possible_genotypes_b: frozenset[LocusGenotype] | None
    outcome: Outcome
    genotype_a: LocusGenotype | None
    genotype_b: LocusGenotype | None
    possible_a: int
    possible_b: int
    records_a: int
    records_b: int
    invalid_reason_a: InvalidReason | None
    invalid_reason_b: InvalidReason | None


HARD_CASES = (
    HardCase(
        name="S1",
        reference=Reference.GRCH38,
        # excerpt  chr1:4416100-4416400, 301 bases; excerpt position = position - 4416099; sample HG002
        contig="chr1",
        offset=4416099,
        sample="HG002",
        bases=(
            "CCACTGATTTGGGAAAAGATGATCAATGAATTCCCATCTTCAAGTGCACTTATTAGCTGCTGCCTGACATGACAGA"
            "CTGAAAAACTTAGAGATGAAATTAGTGGTGTTGGTGTGTCATGACCCAGCTCAGTCTGAGTGGATATGCGTGTGTG"
            "TGTGTGTGTGTGAAAGAGAGAGAGCGAGCATTTTCTTAAATAACTAAAATTCTCTCGCAGGTGTCTCCACATTAGA"
            "AGCATCCTATGTAGCTGTTCTGCTTGCAGCTGCCCAAGTCTGCCTGGGCCTCCTGCATGAGGTCATGAGTCTT"
        ),
        lines_a=(
            "145 C CGTGT 0|1",  # A  TR:180  4416244 C CGTGT 0|1
            "145 C CGTGTGT 1|0",  # A  TR:181  4416244 C CGTGTGT 1|0
        ),
        lines_b=(
            "146 GTGTGTGTGTGTGTGTGT GTGTGTGTGTGTGTGTGTGTGT,GTGTGTGTGTGTGTGTGTGTGTGT 1/2",  # B  pvc-ng30:2677  4416245 GTGTGTGTGTGTGTGTGT GTGTGTGTGTGTGTGTGTGTGT,GTGTGTGTGTGTGTGTGTGTGTGT 1/2
        ),
        # locus    chr1:4416244-4416275 (excerpt 145-176), kind STR; STR regions: 4416245-4416275 str-bundle
        locus_1based=(145, 176),
        kind=Kind.STR,
        str_regions_1based=((146, 176, RegionType.STR_BUNDLE, None),),
        alleles=(
            "CGTGTGTGTGTGTGTGTGTGAAAGAGAGAGAG",
            "CGTGTGTGTGTGTGTGTGTGTGTGTGAAAGAGAGAGAG",
            "CGTGTGTGTGTGTGTGTGTGTGTGAAAGAGAGAGAG",
        ),
        # possible genotypes  A: 1/2   B: 1/2
        possible_genotypes_a=frozenset({(1, 2)}),
        possible_genotypes_b=frozenset({(1, 2)}),
        # row  A 1/2  B 1/2  possible 1, 1  records 2, 1  -> same genotype
        outcome=Outcome.SAME_GENOTYPE,
        genotype_a=(1, 2),
        genotype_b=(1, 2),
        possible_a=1,
        possible_b=1,
        records_a=2,
        records_b=1,
        invalid_reason_a=None,
        invalid_reason_b=None,
    ),
    HardCase(
        name="S2",
        reference=Reference.GRCH38,
        # excerpt  chr1:150328890-150329190, 301 bases; excerpt position = position - 150328889; sample HG002
        contig="chr1",
        offset=150328889,
        sample="HG002",
        bases=(
            "CCATGTTGGGCAGGCTGGTCTTCAACTTCTGACCTCAGGTCATCTGCCCACCTTCGCCTCCCAGAGTGCTGGGATT"
            "ATAGGCGTGAGCCAGCGTGCCTCGCCTTATTGTGATCTTGATTTGAATTAGTGTTTAGCATTTGGGGTAAACTTTT"
            "TTTTTTGAGATAGAGTTTCACTCTGTTGCCCATGCTGGAGTGTAGTGGCGTGATCTCAGCTCACGGCAACCTCCAC"
            "CTCCCAGGTTCAAGCGATTCTCCTGCCTAAGCCTCCGGAGTAGCTGGGACTACATGCGAACACCACTACGCCT"
        ),
        lines_a=(
            "144 T A 1|1",  # A  TR:1583  150329033 T A 1|1
            "146 A ACTTTTT 1|1",  # A  TR:1584  150329035 A ACTTTTT 1|1
            "147 A T 1|1",  # A  TR:1585  150329036 A T 1|1
            "148 C T 1|1",  # A  TR:1586  150329037 C T 1|1
        ),
        lines_b=(
            "143 GT G 1/1",  # B  pvc-ng30:3583  150329032 GT G 1/1
            "149 TTTTTTTTTT TTTTTTTTTTTTTTTTT 1/1",  # B  pvc-ng30:3584  150329038 TTTTTTTTTT TTTTTTTTTTTTTTTTT 1/1
        ),
        # locus    chr1:150329032-150329048 (excerpt 143-159), kind STR; STR regions: 150329038-150329047 str T
        locus_1based=(143, 159),
        kind=Kind.STR,
        str_regions_1based=((149, 158, RegionType.STR, "T"),),
        alleles=(
            "GTAAACTTTTTTTTTTG",
            "GAAACTTTTTTTTTTTTTTTTTG",
        ),
        # possible genotypes  A: 1/1   B: 1/1
        possible_genotypes_a=frozenset({(1, 1)}),
        possible_genotypes_b=frozenset({(1, 1)}),
        # row  A 1/1  B 1/1  possible 1, 1  records 4, 2  -> same genotype
        outcome=Outcome.SAME_GENOTYPE,
        genotype_a=(1, 1),
        genotype_b=(1, 1),
        possible_a=1,
        possible_b=1,
        records_a=4,
        records_b=2,
        invalid_reason_a=None,
        invalid_reason_b=None,
    ),
    HardCase(
        name="S3",
        reference=Reference.GRCH38,
        # excerpt  chr1:9955260-9955560, 301 bases; excerpt position = position - 9955259; sample HG002
        contig="chr1",
        offset=9955259,
        sample="HG002",
        bases=(
            "TGGGCGTGGTGGTATGCACCTGTAATCCCAGCTACTCGGGAGGCTGAGGCAGGAGAATCGCTTGAACCCGGGAGGC"
            "GGAGGTTGCAGTGAGCTGGGATCGCGCCATTGCACTCCAGCCTGGGGACAAGAGCGAGACTTTGTCTCAAAAAAAA"
            "AAAAAAAAAAAGAAAGAAACCCCTTCTCTACTAAAAATACAAAATTAGCTGGGTGTGGTGGCGCATGCCTGTAATT"
            "GCAGCTACTAGGGAGGCTGAGGCAGGAGAATCACTTGAACTTGGGAGGCAGAGGTTGCAGTGAGCCGAGATTG"
        ),
        lines_a=(
            "144 C CAAA 1|1",  # A  TR:283  9955403 C CAAA 1|1
        ),
        lines_b=(
            "145 AAAAAAAAAAAAAAAAAAA AAAAAAAAAAAAAAAAAAAAAA 1/1",  # B  pvc-ng30:2747  9955404 AAAAAAAAAAAAAAAAAAA AAAAAAAAAAAAAAAAAAAAAA 1/1
            "164 G A 0/1",  # B  pvc-ng30:2748  9955423 G A 0/1
        ),
        # locus    chr1:9955403-9955423 (excerpt 144-164), kind STR; STR regions: 9955404-9955422 str A
        locus_1based=(144, 164),
        kind=Kind.STR,
        str_regions_1based=((145, 163, RegionType.STR, "A"),),
        alleles=(
            "CAAAAAAAAAAAAAAAAAAAG",
            "CAAAAAAAAAAAAAAAAAAAAAAG",
            "CAAAAAAAAAAAAAAAAAAAAAAA",
        ),
        # possible genotypes  A: 1/1   B: 1/2
        possible_genotypes_a=frozenset({(1, 1)}),
        possible_genotypes_b=frozenset({(1, 2)}),
        # row  A 1/1  B 1/2  possible 1, 1  records 1, 2  -> different
        outcome=Outcome.DIFFERENT,
        genotype_a=(1, 1),
        genotype_b=(1, 2),
        possible_a=1,
        possible_b=1,
        records_a=1,
        records_b=2,
        invalid_reason_a=None,
        invalid_reason_b=None,
    ),
    HardCase(
        name="S4",
        reference=Reference.GRCH38,
        # excerpt  chr1:14722006-14722306, 301 bases; excerpt position = position - 14722005; sample HG002
        contig="chr1",
        offset=14722005,
        sample="HG002",
        bases=(
            "TTAGCCAGGCATGGTGGTGCATGCCTGTAATCTCAGCTACTTGGGAGGCTGAGGTAGAATTGCTTGAACCCAGGAG"
            "GCGGAGGTTGCAGTGAGCTGAGACTGCACCACTGCACTCAGCCTAAGCGTCAGAGTGAGACTCTGTCTCAAAAAAA"
            "AAAAAGAATAAAAAATAAAAAATGGTGAGGATTGATCCACATGGACATGGAGAGATGACCAGGAAACAGTTACATT"
            "TAAGCAAATTTCCATGTGCGTTTTAGTATGTGATGATCCATTTTGTAGAATAAAATTATATATCCATCACAGG"
        ),
        lines_a=(
            "153 A T 0|1",  # A  TR:317  14722158 A T 0|1
            "156 AAG A 0|1",  # A  TR:318  14722161 AAG A 0|1
        ),
        lines_b=(
            "146 AAAAAAAAAAAA AAAAAAATAA 0/1",  # B  pvc-ng30:2776  14722151 AAAAAAAAAAAA AAAAAAATAA 0/1
        ),
        # locus    chr1:14722151-14722179 (excerpt 146-174), kind STR; STR regions: 14722151-14722179 str-bundle
        locus_1based=(146, 174),
        kind=Kind.STR,
        str_regions_1based=((146, 174, RegionType.STR_BUNDLE, None),),
        alleles=(
            "AAAAAAAAAAAAGAATAAAAAATAAAAAA",
            "AAAAAAATAAAAATAAAAAATAAAAAA",
            "AAAAAAATAAGAATAAAAAATAAAAAA",
        ),
        # possible genotypes  A: 0/1   B: 0/2
        possible_genotypes_a=frozenset({(0, 1)}),
        possible_genotypes_b=frozenset({(0, 2)}),
        # row  A 0/1  B 0/2  possible 1, 1  records 2, 1  -> different
        outcome=Outcome.DIFFERENT,
        genotype_a=(0, 1),
        genotype_b=(0, 2),
        possible_a=1,
        possible_b=1,
        records_a=2,
        records_b=1,
        invalid_reason_a=None,
        invalid_reason_b=None,
    ),
    HardCase(
        name="S5",
        reference=Reference.GRCH38,
        # excerpt  chr1:14722006-14722306, 301 bases; excerpt position = position - 14722005; sample HG002
        contig="chr1",
        offset=14722005,
        sample="HG002",
        bases=(
            "TTAGCCAGGCATGGTGGTGCATGCCTGTAATCTCAGCTACTTGGGAGGCTGAGGTAGAATTGCTTGAACCCAGGAG"
            "GCGGAGGTTGCAGTGAGCTGAGACTGCACCACTGCACTCAGCCTAAGCGTCAGAGTGAGACTCTGTCTCAAAAAAA"
            "AAAAAGAATAAAAAATAAAAAATGGTGAGGATTGATCCACATGGACATGGAGAGATGACCAGGAAACAGTTACATT"
            "TAAGCAAATTTCCATGTGCGTTTTAGTATGTGATGATCCATTTTGTAGAATAAAATTATATATCCATCACAGG"
        ),
        lines_a=(
            "153 A T 0|1",  # A  TR:317  14722158 A T 0|1
            "156 AAG A 0|1",  # A  TR:318  14722161 AAG A 0|1
        ),
        lines_b=(
            "146 AAAAAAAAAAAAG AAAAAAATAAA 0/1",  # B  pvc-main30:2770  14722151 AAAAAAAAAAAAG AAAAAAATAAA 0/1
        ),
        # locus    chr1:14722151-14722179 (excerpt 146-174), kind STR; STR regions: 14722151-14722179 str-bundle
        locus_1based=(146, 174),
        kind=Kind.STR,
        str_regions_1based=((146, 174, RegionType.STR_BUNDLE, None),),
        alleles=(
            "AAAAAAAAAAAAGAATAAAAAATAAAAAA",
            "AAAAAAATAAAAATAAAAAATAAAAAA",
        ),
        # possible genotypes  A: 0/1   B: 0/1
        possible_genotypes_a=frozenset({(0, 1)}),
        possible_genotypes_b=frozenset({(0, 1)}),
        # row  A 0/1  B 0/1  possible 1, 1  records 2, 1  -> same genotype
        outcome=Outcome.SAME_GENOTYPE,
        genotype_a=(0, 1),
        genotype_b=(0, 1),
        possible_a=1,
        possible_b=1,
        records_a=2,
        records_b=1,
        invalid_reason_a=None,
        invalid_reason_b=None,
    ),
    HardCase(
        name="S6",
        reference=Reference.GRCH38,
        # excerpt  chr1:69233295-69233595, 301 bases; excerpt position = position - 69233294; sample HG002
        contig="chr1",
        offset=69233294,
        sample="HG002",
        bases=(
            "CAGTGTCCCAAGGACAAGTCTGGTTTTGGTGGACATACCACTCTATTAACTTAGTCATCTAAATTAATCCTCCAGC"
            "CTGAATTGGCTTAAAACAATCTTGGCTGCTTTCACATTAGCTAATATGTTATATTATGCTATAATAATAATAATAA"
            "TAATAATAATAATAATAATGGTTTTCATCTATTGATAATTGTGCATTATCTCTATACCTAAAAGCATTCCCAAAGG"
            "GTGAAAATTTATTCCCTTCATTTTGCAAATGAAAACTTTTAGGCACCATTAGATTAGAAGGCATATCCAAGCC"
        ),
        lines_a=(
            "136 TATAATAATA T 1|1",  # A  TR:966  69233430 TATAATAATA T 1|1
        ),
        lines_b=(
            "137 ATAATAATAATAATAATAATAATAATAATAATA ATAATAATAATAATAATAATAATA 1/1",  # B  pvc-ng30:3258  69233431 ATAATAATAATAATAATAATAATAATAATAATA ATAATAATAATAATAATAATAATA 1/1
        ),
        # locus    chr1:69233430-69233465 (excerpt 136-171), kind STR; STR regions: 69233431-69233463 str ATA
        locus_1based=(136, 171),
        kind=Kind.STR,
        str_regions_1based=((137, 169, RegionType.STR, "ATA"),),
        alleles=(
            "TATAATAATAATAATAATAATAATAATAATAATAAT",
            "TATAATAATAATAATAATAATAATAAT",
        ),
        # possible genotypes  A: 1/1   B: 1/1
        possible_genotypes_a=frozenset({(1, 1)}),
        possible_genotypes_b=frozenset({(1, 1)}),
        # row  A 1/1  B 1/1  possible 1, 1  records 1, 1  -> same genotype
        outcome=Outcome.SAME_GENOTYPE,
        genotype_a=(1, 1),
        genotype_b=(1, 1),
        possible_a=1,
        possible_b=1,
        records_a=1,
        records_b=1,
        invalid_reason_a=None,
        invalid_reason_b=None,
    ),
    HardCase(
        name="S7",
        reference=Reference.GRCH38,
        # excerpt  chr11:123865020-123865320, 301 bases; excerpt position = position - 123865019; sample HG003
        contig="chr11",
        offset=123865019,
        sample="HG003",
        bases=(
            "TTTGATAAAACAGCGTTTTTGAAATAAGGCTTTTTTTGTTTGTTTGTTTTGTTTTGTTTTAATCATCTCTAACCCC"
            "ACAGCTGGCAGTGATCTAAAAAGAAAGAAAGAATTTATGGTATTAAAAATTCTGTAACCGTGTGTGTGTGTGTGTG"
            "TGTGTGTGTGTGTGTGTGTGTAACAAGCTCAACAAGACATTGTTCACAAGGATGATTGATTAAAGACATAAATAGC"
            "TTAATTTTATATGACTTTTATTGGTAGTTACTAAGGTATGTCTCTTCTAAAGTATTTTCCATGGTTCCTATTT"
        ),
        lines_a=(
            "135 CGTGTGTGT CGTGTGT,C 2/1",  # A  GIAB-HG003:2756918  123865154 CGTGTGTGT CGTGTGT,C 2/1
            "152 G A 0/1",  # A  GIAB-HG003:2756919  123865171 G A 0/1
        ),
        lines_b=(
            "136 GTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGT GTGTGTGTGTGTGTATGTGTGTGTGTGTGTGTGTGT,GTGTGTGTGTGTGTGTGTGTGTGTGTGTGT 1/2",  # B  pvc-300x-HG003:3112  123865155 GTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGT GTGTGTGTGTGTGTATGTGTGTGTGTGTGTGTGTGT,GTGTGTGTGTGTGTGTGTGTGTGTGTGTGT 1/2
        ),
        # locus    chr11:123865140-123865192 (excerpt 121-173), kind STR; STR regions: 123865140-123865192 str-bundle
        locus_1based=(121, 173),
        kind=Kind.STR,
        str_regions_1based=((121, 173, RegionType.STR_BUNDLE, None),),
        alleles=(
            "AAAAATTCTGTAACCGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGT",
            "AAAAATTCTGTAACCGTGTGTGTGTGTGTGTGTGTGTGTGTGTGT",
            "AAAAATTCTGTAACCGTGTGTGTGTGTGTATGTGTGTGTGTGTGTGTGTGT",
            "AAAAATTCTGTAACCGTGTGTGTATGTGTGTGTGTGTGTGTGTGT",
            "AAAAATTCTGTAACCGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGT",
        ),
        # possible genotypes  A: 1/2, 3/4   B: 1/2
        possible_genotypes_a=frozenset({(1, 2), (3, 4)}),
        possible_genotypes_b=frozenset({(1, 2)}),
        # row  A 1/2  B 1/2  possible 2, 1  records 2, 1  -> same genotype
        outcome=Outcome.SAME_GENOTYPE,
        genotype_a=(1, 2),
        genotype_b=(1, 2),
        possible_a=2,
        possible_b=1,
        records_a=2,
        records_b=1,
        invalid_reason_a=None,
        invalid_reason_b=None,
    ),
    HardCase(
        name="S8",
        reference=Reference.GRCH38,
        # excerpt  chr14:90749262-90749562, 301 bases; excerpt position = position - 90749261; sample HG004
        contig="chr14",
        offset=90749261,
        sample="HG004",
        bases=(
            "AAAAACATCCAGCCATTTTATCTTCAAAAATGTTTATGTCACCCTCCTCCTTCCTTAGGGGACTTTACGCACATAT"
            "TAGGCCATGTGAAGTTGTCCTACAGATCACTTATGGTCTTTCCACTTAAAATATAATAATTTTTCTTTTTTTTTTT"
            "TTTTTTTTTGAGACAGAGTTTCGCTCTTGTTGCCCAGGCTGAAGTGCAATGGCGCAATCTCAGCTCACTGCAACCT"
            "CCACCTCCCAGGTTCAAGTGATTCTCCTGCTTCAGCCTCCCGAGTAGCTGGGATTGCAGGCATGCACCACCAC"
        ),
        lines_a=(
            "141 CTT CT,C 2/1",  # A  GIAB-HG004:3265346  90749402 CTT CT,C 2/1
            "161 T G 0/1",  # A  GIAB-HG004:3265347  90749422 T G 0/1
        ),
        lines_b=(
            "142 TTTTTTTTTTTTTTTTTTTT TTTTTTTTTTTTTTTTTG,TTTTTTTTTTTTTTTTTTT 1/2",  # B  pvc-300x-HG004:3106  90749403 TTTTTTTTTTTTTTTTTTTT TTTTTTTTTTTTTTTTTG,TTTTTTTTTTTTTTTTTTT 1/2
        ),
        # locus    chr14:90749397-90749422 (excerpt 136-161), kind STR; STR regions: 90749397-90749422 str-bundle
        locus_1based=(136, 161),
        kind=Kind.STR,
        str_regions_1based=((136, 161, RegionType.STR_BUNDLE, None),),
        alleles=(
            "TTTTTCTTTTTTTTTTTTTTTTTTTT",
            "TTTTTCTTTTTTTTTTTTTTTTTT",
            "TTTTTCTTTTTTTTTTTTTTTTTTG",
            "TTTTTCTTTTTTTTTTTTTTTTTG",
            "TTTTTCTTTTTTTTTTTTTTTTTTT",
        ),
        # possible genotypes  A: 1/2, 3/4   B: 3/4
        possible_genotypes_a=frozenset({(1, 2), (3, 4)}),
        possible_genotypes_b=frozenset({(3, 4)}),
        # row  A 3/4  B 3/4  possible 2, 1  records 2, 1  -> same genotype
        outcome=Outcome.SAME_GENOTYPE,
        genotype_a=(3, 4),
        genotype_b=(3, 4),
        possible_a=2,
        possible_b=1,
        records_a=2,
        records_b=1,
        invalid_reason_a=None,
        invalid_reason_b=None,
    ),
    HardCase(
        name="S9",
        reference=Reference.GRCH38,
        # excerpt  chr10:4853300-4853600, 301 bases; excerpt position = position - 4853299; sample HG004
        contig="chr10",
        offset=4853299,
        sample="HG004",
        bases=(
            "CCCCAGCCTCATAAATTAGTAGTTGGAAAAAGGAGGCATATTTAGAGCTTTTTCAGATAATTGTATTTCTTTGATA"
            "CATTAGACTGGACACACAGTAGTTTGTTTAAGGTTAATTGCAATATTGCAATGAGGAACTTGTTGCCAAAAAAAAA"
            "AAAAAAACCAAACTCTGTAAAATATCTGAAGAGATTGATTCTAAGCCAAATATGAGAACCATGACCCATGACACAG"
            "CCTCAGGAGATCCTGAGGACATGTGCCCATTGTGGTCTGGTTACAGCTTAGTTTTATATGTTTTAGGGAGACA"
        ),
        lines_a=(
            "143 CA C 0/1",  # A  GIAB-HG004:2407486  4853442 CA C 0/1
            "159 A C 1/1",  # A  GIAB-HG004:2407487  4853458 A C 1/1
        ),
        lines_b=(
            "144 AAAAAAAAAAAAAAAA AAAAAAAAAAAAAAAC,AAAAAAAAAAAAAAC 1/2",  # B  pvc-300x-HG004:3028  4853443 AAAAAAAAAAAAAAAA AAAAAAAAAAAAAAAC,AAAAAAAAAAAAAAC 1/2
        ),
        # locus    chr10:4853442-4853458 (excerpt 143-159), kind STR; STR regions: 4853443-4853458 str A
        locus_1based=(143, 159),
        kind=Kind.STR,
        str_regions_1based=((144, 159, RegionType.STR, "A"),),
        alleles=(
            "CAAAAAAAAAAAAAAAA",
            "CAAAAAAAAAAAAAAAC",
            "CAAAAAAAAAAAAAAC",
        ),
        # possible genotypes  A: 1/2   B: 1/2
        possible_genotypes_a=frozenset({(1, 2)}),
        possible_genotypes_b=frozenset({(1, 2)}),
        # row  A 1/2  B 1/2  possible 1, 1  records 2, 1  -> same genotype
        outcome=Outcome.SAME_GENOTYPE,
        genotype_a=(1, 2),
        genotype_b=(1, 2),
        possible_a=1,
        possible_b=1,
        records_a=2,
        records_b=1,
        invalid_reason_a=None,
        invalid_reason_b=None,
    ),
    HardCase(
        name="S10",
        reference=Reference.GRCH38,
        # excerpt  chr20:48207270-48207570, 301 bases; excerpt position = position - 48207269; sample HG004
        contig="chr20",
        offset=48207269,
        sample="HG004",
        bases=(
            "CCAGATGACTAGCCCACATGGCTTTCTCCCTCAACAAATCCAGATCTTTTCTCAAATGTCACCTTATCAGTGAGCC"
            "CTTCCCAATTACTCAATATATAATTTAGACCCACACAATCACACACATGCATGTGTGTCATGCACACACACACACA"
            "CACACACACACACATCTTATTCCTCTCACCTGCATTCTTTTTCTCCTTTGTAGTTACCACCACTCAACATGTTATT"
            "TATTTACTTAATTACTCTGTTTACTGCTGGCTCCACTGATAGAATATCAGCCGCATGAGGTCAACAGCTGTTT"
        ),
        lines_a=(
            "138 GCACA GCACACA,G 2/1",  # A  GIAB-HG004:3904085  48207407 GCACA GCACACA,G 2/1
            "166 A G 0/1",  # A  GIAB-HG004:3904086  48207435 A G 0/1
        ),
        lines_b=(
            "139 CACACACACACACACACACACACACACA CACACACACACACACACACACACACACACA,CACACACACACACACACACACACG 1/2",  # B  pvc-300x-HG004:3235  48207408 CACACACACACACACACACACACACACA CACACACACACACACACACACACACACACA,CACACACACACACACACACACACG 1/2
        ),
        # locus    chr20:48207407-48207436 (excerpt 138-167), kind STR; STR regions: 48207408-48207435 str CA
        locus_1based=(138, 167),
        kind=Kind.STR,
        str_regions_1based=((139, 166, RegionType.STR, "CA"),),
        alleles=(
            "GCACACACACACACACACACACACACACAT",
            "GCACACACACACACACACACACACAT",
            "GCACACACACACACACACACACACACACACGT",
            "GCACACACACACACACACACACACGT",
            "GCACACACACACACACACACACACACACACAT",
        ),
        # possible genotypes  A: 1/2, 3/4   B: 3/4
        possible_genotypes_a=frozenset({(1, 2), (3, 4)}),
        possible_genotypes_b=frozenset({(3, 4)}),
        # row  A 3/4  B 3/4  possible 2, 1  records 2, 1  -> same genotype
        outcome=Outcome.SAME_GENOTYPE,
        genotype_a=(3, 4),
        genotype_b=(3, 4),
        possible_a=2,
        possible_b=1,
        records_a=2,
        records_b=1,
        invalid_reason_a=None,
        invalid_reason_b=None,
    ),
    HardCase(
        name="S11",
        reference=Reference.GRCH38,
        # excerpt  chr8:24760150-24760450, 301 bases; excerpt position = position - 24760149; sample HG004
        contig="chr8",
        offset=24760149,
        sample="HG004",
        bases=(
            "TGAGACATAAGAGAATAAAAATGGATAATTCAATAAAATCAGGAAAATCATTCAATATCTGAAGGAAATATTCAAT"
            "AGGTATTTTAAAAAGAACCAAACAGAAATCTTTGAGCTAAGGAAGTAAGCTAATGAAATTTAAAAATATACACATA"
            "TATATATATATATATAAGAGCTTCAGCAGCAGACTAGATCAAGTAGAAGAAAGAATCTGTGAACGTGAAGATAGAT"
            "AGTAATTTTGAAATTATTTAGTCAGAGGGAAAAAAAAGAATAAAAAGAACAAAGACAGACTGTGGGACTTATT"
        ),
        lines_a=(
            "145 TAC T 0/1",  # A  GIAB-HG004:2066390  24760294 TAC T 0/1
            "147 CACAT C 0/1",  # A  GIAB-HG004:2066391  24760296 CACAT C 0/1
        ),
        lines_b=(
            "143 TATACAC TATAC,* 1/2",  # B  pvc-300x-HG004:2886  24760292 TATACAC TATAC,* 1/2
            "148 ACATATATATATATATATAT ATATATATATATATAT 0/1",  # B  pvc-300x-HG004:2887  24760297 ACATATATATATATATATAT ATATATATATATATAT 0/1
        ),
        # locus    chr8:24760287-24760317 (excerpt 138-168), kind STR; STR regions: 24760287-24760317 str-bundle
        locus_1based=(138, 168),
        kind=Kind.STR,
        str_regions_1based=((138, 168, RegionType.STR_BUNDLE, None),),
        alleles=(
            "AAAAATATACACATATATATATATATATATA",
            "AAAAATATATATATATATATATATA",
            "AAAAATATACATATATATATATATATA",
            "AAAAATATACATATATATATATATATATA",
        ),
        # possible genotypes  A: 0/1, 2/3   B: 2/3
        possible_genotypes_a=frozenset({(0, 1), (2, 3)}),
        possible_genotypes_b=frozenset({(2, 3)}),
        # row  A 2/3  B 2/3  possible 2, 1  records 2, 2  -> same genotype
        outcome=Outcome.SAME_GENOTYPE,
        genotype_a=(2, 3),
        genotype_b=(2, 3),
        possible_a=2,
        possible_b=1,
        records_a=2,
        records_b=2,
        invalid_reason_a=None,
        invalid_reason_b=None,
    ),
    HardCase(
        name="S12",
        reference=Reference.GRCH38,
        # excerpt  chr4:28574875-28575175, 301 bases; excerpt position = position - 28574874; sample HG002
        contig="chr4",
        offset=28574874,
        sample="HG002",
        bases=(
            "TAGATATTACCATGAAACTTAATTTGACACAGTCATAATTTATAATACCCAATAACCATGGCTTGTATATCACTTG"
            "TCTTGTGCATATTAAGGTTTGCTTTGGGGTGAAAGCATATCTCTGTGTATATATTTACACACCCCCCCCCCCACAC"
            "ACACACACACACAGAAAGAGAGTGAGAGAGACTCTGTATCTCTAATTCTATACCCTGCTTTTGATTCCGTGCTTGT"
            "CCCTTGGAAGTCAGAATCTATCTCTTTTTGAGTCATCATCATCATTGTCTTAGTCCTTTCAGGCTGCTATAAC"
        ),
        lines_a=(
            "137 ACC A 1|0",  # A  TR:7780  28575011 ACC A 1|0
            "143 CCCCCCA C 0|1",  # A  TR:7781  28575017 CCCCCCA C 0|1
            "166 G C 1|0",  # A  TR:7782  28575040 G C 1|0
        ),
        lines_b=(
            "133 ACACACCCCCCCCCCCACACACACACACACACAG ACACACCCCCCCACACACACACACACAG,ACACACCCCCCCCCACACACACACACACACAC 1/2",  # B  pvc-main30:7107  28575007 ACACACCCCCCCCCCCACACACACACACACACAG ACACACCCCCCCACACACACACACACAG,ACACACCCCCCCCCACACACACACACACACAC 1/2
        ),
        # locus    chr4:28575007-28575040 (excerpt 133-166), kind STR; STR regions: 28575007-28575039 str-bundle
        locus_1based=(133, 166),
        kind=Kind.STR,
        str_regions_1based=((133, 165, RegionType.STR_BUNDLE, None),),
        alleles=(
            "ACACACCCCCCCCCCCACACACACACACACACAG",
            "ACACACCCCCCCCCACACACACACACACACAC",
            "ACACACCCCCCCACACACACACACACAG",
        ),
        # possible genotypes  A: 1/2   B: 1/2
        possible_genotypes_a=frozenset({(1, 2)}),
        possible_genotypes_b=frozenset({(1, 2)}),
        # row  A 1/2  B 1/2  possible 1, 1  records 3, 1  -> same genotype
        outcome=Outcome.SAME_GENOTYPE,
        genotype_a=(1, 2),
        genotype_b=(1, 2),
        possible_a=1,
        possible_b=1,
        records_a=3,
        records_b=1,
        invalid_reason_a=None,
        invalid_reason_b=None,
    ),
    HardCase(
        name="S13",
        reference=Reference.GRCH38,
        # excerpt  chr9:11024574-11024874, 301 bases; excerpt position = position - 11024573; sample HG004
        contig="chr9",
        offset=11024573,
        sample="HG004",
        bases=(
            "ATTATATTACTCAGCGTGGGGTTTATTTACATTCATGTGTATGTATTTCCCCAAGGTTATTCTGCCCTTCTTTGTC"
            "CTACGTATCATTTTCTTCAAGAGAAAACAATTTTATACACGCAGTTGTATGTACACACACACACACACACACACAC"
            "AGAGAGAGAGAGAGAGGTCAAGCACACATAAAACCCAAGAGAGAAGTCAAGCACACATAAAACCCAATCTAAAACA"
            "AAACACATTTAGGGAGTATTTTTTGTTTCTCTTATTAGTAGCGTTCAGTTGTACTTTGTATCCACAAAGTTAA"
        ),
        lines_a=(
            "150 C CAG 1/1",  # A  GIAB-HG004:2245326  11024723 C CAG 1/1
            "152 C G 1/1",  # A  GIAB-HG004:2245327  11024725 C G 1/1
        ),
        lines_b=(
            "150 CAC CAGAG 1/1",  # B  pvc-explain30-HG004:2917  11024723 CAC CAGAG 1/1
        ),
        # locus    chr9:11024723-11024725 (excerpt 150-152), kind complex; STR regions: none
        locus_1based=(150, 152),
        kind=Kind.COMPLEX,
        str_regions_1based=(),
        alleles=(
            "CAC",
            "CAGAG",
        ),
        # possible genotypes  A: 1/1   B: 1/1
        possible_genotypes_a=frozenset({(1, 1)}),
        possible_genotypes_b=frozenset({(1, 1)}),
        # row  A 1/1  B 1/1  possible 1, 1  records 2, 1  -> same genotype
        outcome=Outcome.SAME_GENOTYPE,
        genotype_a=(1, 1),
        genotype_b=(1, 1),
        possible_a=1,
        possible_b=1,
        records_a=2,
        records_b=1,
        invalid_reason_a=None,
        invalid_reason_b=None,
    ),
    HardCase(
        name="S14",
        reference=Reference.GRCH38,
        # excerpt  chr9:11024574-11024874, 301 bases; excerpt position = position - 11024573; sample HG004
        contig="chr9",
        offset=11024573,
        sample="HG004",
        bases=(
            "ATTATATTACTCAGCGTGGGGTTTATTTACATTCATGTGTATGTATTTCCCCAAGGTTATTCTGCCCTTCTTTGTC"
            "CTACGTATCATTTTCTTCAAGAGAAAACAATTTTATACACGCAGTTGTATGTACACACACACACACACACACACAC"
            "AGAGAGAGAGAGAGAGGTCAAGCACACATAAAACCCAAGAGAGAAGTCAAGCACACATAAAACCCAATCTAAAACA"
            "AAACACATTTAGGGAGTATTTTTTGTTTCTCTTATTAGTAGCGTTCAGTTGTACTTTGTATCCACAAAGTTAA"
        ),
        lines_a=(
            "150 C CAG 1/1",  # A  GIAB-HG004:2245326  11024723 C CAG 1/1
            "152 C G 1/1",  # A  GIAB-HG004:2245327  11024725 C G 1/1
        ),
        lines_b=(
            "150 C CAG,G 0/1",  # B  pvc-ng30-HG004:2916  11024723 C CAG,G 0/1
            "152 C G 1/1",  # B  pvc-ng30-HG004:2917  11024725 C G 1/1
        ),
        # locus    chr9:11024723-11024725 (excerpt 150-152), kind complex; STR regions: none
        locus_1based=(150, 152),
        kind=Kind.COMPLEX,
        str_regions_1based=(),
        alleles=(
            "CAC",
            "CAGAG",
            "CAG",
        ),
        # possible genotypes  A: 1/1   B: 1/2
        possible_genotypes_a=frozenset({(1, 1)}),
        possible_genotypes_b=frozenset({(1, 2)}),
        # row  A 1/1  B 1/2  possible 1, 1  records 2, 2  -> different
        outcome=Outcome.DIFFERENT,
        genotype_a=(1, 1),
        genotype_b=(1, 2),
        possible_a=1,
        possible_b=1,
        records_a=2,
        records_b=2,
        invalid_reason_a=None,
        invalid_reason_b=None,
    ),
    HardCase(
        name="E1",
        reference=Reference.GRCH38,
        # excerpt  chr5:7613015-7613315, 301 bases; excerpt position = position - 7613014; sample HG002
        contig="chr5",
        offset=7613014,
        sample="HG002",
        bases=(
            "GAGTGAGACTCTGACTCAAAAAAATATATATATACCTAATGTAAATGATGAGTTAATGGGTGCAGCACACCAACAT"
            "GGCACATGTATACATATGTAACAAACCTGCACATTGTACACATGTACCCTAAAACTTAAAGTATAATAATAAAAAA"
            "AAAAAGATTGCACACAAGAGCTGTTAACTCATAAGCTGAGCTCTGAGACAGGCACACTGGAGAAGCCCTAGGAATC"
            "CCAGGAAACCAGAGATGGAGAAGGCACAAATAGCTCAACTACAGAAATGCAGGCTGATATGGTGTGTCCTATC"
        ),
        lines_a=(
            "144 AAT A 1|0",  # A  TR:10175  7613158 AAT A 1|0
        ),
        lines_b=(
            "145 ATAAAAAAAAAAA AAAAAAAAAAA 0/1",  # B  pvc-touch30:8329  7613159 ATAAAAAAAAAAA AAAAAAAAAAA 0/1
        ),
        # locus    chr5:7613158-7613171 (excerpt 144-157), kind STR; STR regions: 7613161-7613171 str A
        locus_1based=(144, 157),
        kind=Kind.STR,
        str_regions_1based=((147, 157, RegionType.STR, "A"),),
        alleles=(
            "AATAAAAAAAAAAA",
            "AAAAAAAAAAAA",
        ),
        # possible genotypes  A: 0/1   B: 0/1
        possible_genotypes_a=frozenset({(0, 1)}),
        possible_genotypes_b=frozenset({(0, 1)}),
        # row  A 0/1  B 0/1  possible 1, 1  records 1, 1  -> same genotype
        outcome=Outcome.SAME_GENOTYPE,
        genotype_a=(0, 1),
        genotype_b=(0, 1),
        possible_a=1,
        possible_b=1,
        records_a=1,
        records_b=1,
        invalid_reason_a=None,
        invalid_reason_b=None,
    ),
    HardCase(
        name="E2",
        reference=Reference.GRCH38,
        # excerpt  chr13:31685085-31685385, 301 bases; excerpt position = position - 31685084; sample HG002
        contig="chr13",
        offset=31685084,
        sample="HG002",
        bases=(
            "AGTTCGGAGGACCAGCCTCTCATTTATGTTACTGACAGCCAATGCTCCTCACGGACAAATCTTTCCCAGACCTTTC"
            "TTATTCTGAATAAAGAAAATTTGGAGCGCAGTTTGCTCAAGGTCTACTTTTCATTCAGAAAAACTATATATATATA"
            "TATATATATACTGACTTCTTAGTTTTGCATCTCTGATATATCTTGTTTTTCCTTTAGTGTCATTAGGCTTCCTCTT"
            "GTAACAGTCTGTGCCCACATCCAGTTTTCTCATTACAATATCTTGCATTTATATTATGTTTCACAATTTTGCA"
        ),
        lines_a=(
            "139 ACTAT A 1|1",  # A  TR:26438  31685223 ACTAT A 1|1
            "144 A G 1|1",  # A  TR:26439  31685228 A G 1|1
        ),
        lines_b=(
            "139 AC AG,* 1/2",  # B  pvc-touch30:17007  31685223 AC AG,* 1/2
            "140 CTATATATATATATATATATATA GTATATATATATATATATA,CTATATATATATATATATA 1/2",  # B  pvc-touch30:17008  31685224 CTATATATATATATATATATATA GTATATATATATATATATA,CTATATATATATATATATA 1/2
        ),
        # locus    chr13:31685219-31685246 (excerpt 135-162), kind STR; STR regions: 31685219-31685246 str-bundle
        locus_1based=(135, 162),
        kind=Kind.STR,
        str_regions_1based=((135, 162, RegionType.STR_BUNDLE, None),),
        alleles=(
            "AAAAACTATATATATATATATATATATA",
            "AAAAAGTATATATATATATATATA",
        ),
        # possible genotypes  A: 1/1   B: 1/1
        possible_genotypes_a=frozenset({(1, 1)}),
        possible_genotypes_b=frozenset({(1, 1)}),
        # row  A 1/1  B 1/1  possible 1, 1  records 2, 2  -> same genotype
        outcome=Outcome.SAME_GENOTYPE,
        genotype_a=(1, 1),
        genotype_b=(1, 1),
        possible_a=1,
        possible_b=1,
        records_a=2,
        records_b=2,
        invalid_reason_a=None,
        invalid_reason_b=None,
    ),
    HardCase(
        name="E3",
        reference=Reference.GRCH38,
        # excerpt  chr20:3028985-3029285, 301 bases; excerpt position = position - 3028984; sample HG002
        contig="chr20",
        offset=3028984,
        sample="HG002",
        bases=(
            "AACATTTAAAATACAAAGTATGTTTTACATTAGTTGTCTCAAGTCACTAGTTGTTCCCCATCTGTACTTGTGTTAA"
            "CTTTGGTACAAGGCTTTACACATATTTTCTCACCCTGATGTTTCATTACATCAGGATTTTTTTTTTTTTTTTGGTT"
            "TTGTTTGTTTGTTTGTTTGTTTGAGAGAGGGTCTCACTTTGTCACCCAGGCTGGAGTGCAGTGGCACGATCTTGGC"
            "TCACTGCAGCCTCAACCTCCCAGGTTCAAGTGATCCTCCCACATAGCCCTCCAAGTAGCTAGGACTAAAGGCA"
        ),
        lines_a=(
            "132 A G 1|1",  # A  TR:34237  3029116 A G 1|1
            "148 T G 1|1",  # A  TR:34238  3029132 T G 1|1
            "150 G T 1|1",  # A  TR:34239  3029134 G T 1|1
            "151 T TTTG 1|1",  # A  TR:34240  3029135 T TTTG 1|1
        ),
        lines_b=(
            "131 GA GG,* 1/2",  # B  pvc-touch30:21736  3029115 GA GG,* 1/2
            "132 ATTTTTTTTTTTTTTTTGGTTTTGTTTGTTTGTTTGTTTGTTTG GTTTTTTTTTTTTTTTGGTTTTGTTTGTTTGTTTGTTTGTTTGTTTG,ATTTTTTTTTTTTTTTGGTTTTGTTTGTTTGTTTGTTTGTTTGTTTG 1/2",  # B  pvc-touch30:21737  3029116 ATTTTTTTTTTTTTTTTGGTTTTGTTTGTTTGTTTGTTTGTTTG GTTTTTTTTTTTTTTTGGTTTTGTTTGTTTGTTTGTTTGTTTGTTTG,ATTTTTTTTTTTTTTTGGTTTTGTTTGTTTGTTTGTTTGTTTGTTTG 1/2
        ),
        # locus    chr20:3029115-3029159 (excerpt 131-175), kind STR; STR regions: 3029117-3029159 str-bundle
        locus_1based=(131, 175),
        kind=Kind.STR,
        str_regions_1based=((133, 175, RegionType.STR_BUNDLE, None),),
        alleles=(
            "GATTTTTTTTTTTTTTTTGGTTTTGTTTGTTTGTTTGTTTGTTTG",
            "GGTTTTTTTTTTTTTTTGGTTTTGTTTGTTTGTTTGTTTGTTTGTTTG",
        ),
        # possible genotypes  A: 1/1   B: 1/1
        possible_genotypes_a=frozenset({(1, 1)}),
        possible_genotypes_b=frozenset({(1, 1)}),
        # row  A 1/1  B 1/1  possible 1, 1  records 4, 2  -> same genotype
        outcome=Outcome.SAME_GENOTYPE,
        genotype_a=(1, 1),
        genotype_b=(1, 1),
        possible_a=1,
        possible_b=1,
        records_a=4,
        records_b=2,
        invalid_reason_a=None,
        invalid_reason_b=None,
    ),
    HardCase(
        name="E4",
        reference=Reference.GRCH38,
        # excerpt  chr20:3028985-3029285, 301 bases; excerpt position = position - 3028984; sample HG002
        contig="chr20",
        offset=3028984,
        sample="HG002",
        bases=(
            "AACATTTAAAATACAAAGTATGTTTTACATTAGTTGTCTCAAGTCACTAGTTGTTCCCCATCTGTACTTGTGTTAA"
            "CTTTGGTACAAGGCTTTACACATATTTTCTCACCCTGATGTTTCATTACATCAGGATTTTTTTTTTTTTTTTGGTT"
            "TTGTTTGTTTGTTTGTTTGTTTGAGAGAGGGTCTCACTTTGTCACCCAGGCTGGAGTGCAGTGGCACGATCTTGGC"
            "TCACTGCAGCCTCAACCTCCCAGGTTCAAGTGATCCTCCCACATAGCCCTCCAAGTAGCTAGGACTAAAGGCA"
        ),
        lines_a=(
            "132 A G 1|1",  # A  TR:34237  3029116 A G 1|1
            "148 T G 1|1",  # A  TR:34238  3029132 T G 1|1
            "150 G T 1|1",  # A  TR:34239  3029134 G T 1|1
            "151 T TTTG 1|1",  # A  TR:34240  3029135 T TTTG 1|1
        ),
        lines_b=(
            "131 GA G,GG 1/2",  # B  pvc-base30:21766  3029115 GA G,GG 1/2
            "133 TTTTTTTTTTTTTTTTGGTTTTGTTTGTTTGTTTGTTTGTTTG GTTTTTTTTTTTTTTTGGTTTTGTTTGTTTGTTTGTTTGTTTGTTTG,TTTTTTTTTTTTTTTGGTTTTGTTTGTTTGTTTGTTTGTTTGTTTG 1/2",  # B  pvc-base30:21767  3029117 TTTTTTTTTTTTTTTTGGTTTTGTTTGTTTGTTTGTTTGTTTG GTTTTTTTTTTTTTTTGGTTTTGTTTGTTTGTTTGTTTGTTTGTTTG,TTTTTTTTTTTTTTTGGTTTTGTTTGTTTGTTTGTTTGTTTGTTTG 1/2
        ),
        # locus    chr20:3029115-3029159 (excerpt 131-175), kind STR; STR regions: 3029117-3029159 str-bundle
        locus_1based=(131, 175),
        kind=Kind.STR,
        str_regions_1based=((133, 175, RegionType.STR_BUNDLE, None),),
        alleles=(
            "GATTTTTTTTTTTTTTTTGGTTTTGTTTGTTTGTTTGTTTGTTTG",
            "GGTTTTTTTTTTTTTTTGGTTTTGTTTGTTTGTTTGTTTGTTTGTTTG",
            "GTTTTTTTTTTTTTTTGGTTTTGTTTGTTTGTTTGTTTGTTTGTTTG",
            "GGGTTTTTTTTTTTTTTTGGTTTTGTTTGTTTGTTTGTTTGTTTGTTTG",
        ),
        # possible genotypes  A: 1/1   B: 1/1, 2/3
        possible_genotypes_a=frozenset({(1, 1)}),
        possible_genotypes_b=frozenset({(1, 1), (2, 3)}),
        # row  A 1/1  B 1/1  possible 1, 2  records 4, 2  -> same genotype
        outcome=Outcome.SAME_GENOTYPE,
        genotype_a=(1, 1),
        genotype_b=(1, 1),
        possible_a=1,
        possible_b=2,
        records_a=4,
        records_b=2,
        invalid_reason_a=None,
        invalid_reason_b=None,
    ),
    HardCase(
        name="E5",
        reference=Reference.GRCH38,
        # excerpt  chr8:104588825-104589125, 301 bases; excerpt position = position - 104588824; sample HG002
        contig="chr8",
        offset=104588824,
        sample="HG002",
        bases=(
            "CCCCAGCGAGGAAAAGCAAGAGCAACGCAGACCTCCACCGCGGAGACTCTTTTGTGCTCCAGCGACAGGCCATAAC"
            "CACAGCAGATGGAGAGAGAGAGGAGGAGACGGAGGAGGAGGGAGGAGAAGCTGGAGGTAGACGACGCCGACGCCGC"
            "CGCCGCCGCCGCCGCCGCCGCCGAGCCACCGGCTGCTCCCTGCGCTCTCCGCGGCTGCGGGAGGGGGAAGGGAGGG"
            "GCCGCCGCCGCCCGCGCGCGCTCCCTCCTCCCTCCTCCCTCCGGCTCGCCTGCTCCCACCCCGGGCGGTGCGA"
        ),
        lines_a=(
            "137 A ACGACGC 0|1",  # A  TR:17773  104588961 A ACGACGC 0|1
            "140 ACGCCGACGCCGCCGCCGC A 1|0",  # A  TR:17774  104588964 ACGCCGACGCCGCCGCCGC A 1|0
        ),
        lines_b=(
            "137 ACGACGCCGA ACGACGCCGACGCCGA,* 1/2",  # B  pvc-touch30:12607  104588961 ACGACGCCGA ACGACGCCGACGCCGA,* 1/2
            "141 CGCCGACGCCGCCGCCGCCGCCGCCGCCGCCGC CGCCGCCGCCGCCGC 0/1",  # B  pvc-touch30:12608  104588965 CGCCGACGCCGCCGCCGCCGCCGCCGCCGCCGC CGCCGCCGCCGCCGC 0/1
        ),
        # locus    chr8:104588961-104588997 (excerpt 137-173), kind STR; STR regions: 104588971-104588997 str CGC
        locus_1based=(137, 173),
        kind=Kind.STR,
        str_regions_1based=((147, 173, RegionType.STR, "CGC"),),
        alleles=(
            "ACGACGCCGACGCCGCCGCCGCCGCCGCCGCCGCCGC",
            "ACGACGCCGCCGCCGCCGC",
            "ACGACGCCGACGCCGACGCCGCCGCCGCCGCCGCCGCCGCCGC",
        ),
        # possible genotypes  A: 1/2   B: 1/2
        possible_genotypes_a=frozenset({(1, 2)}),
        possible_genotypes_b=frozenset({(1, 2)}),
        # row  A 1/2  B 1/2  possible 1, 1  records 2, 2  -> same genotype
        outcome=Outcome.SAME_GENOTYPE,
        genotype_a=(1, 2),
        genotype_b=(1, 2),
        possible_a=1,
        possible_b=1,
        records_a=2,
        records_b=2,
        invalid_reason_a=None,
        invalid_reason_b=None,
    ),
    HardCase(
        name="E6",
        reference=Reference.GRCH38,
        # excerpt  chr12:13880855-13881155, 301 bases; excerpt position = position - 13880854; sample HG002
        contig="chr12",
        offset=13880854,
        sample="HG002",
        bases=(
            "TGCCTTCCTGTTTACTGTGTTCACTCTCATGCTGTGCAAAGGTTGTTCCACAGCCTATTTGAGACCTCGCATCTGG"
            "AGGCAGCTTTGTTTCAGTTTTCAAATTCTCTTGCCTTGCTTTGGCTGACAAAGGGGTATAAGGTTGTGTGTGTGTG"
            "TGTGTGTGCGCGTGCGCGCACGCATGTGCGCAGGTGTGCATGCACATGTGCACACGAGAGTATGTATAAAAAAGAT"
            "AACTTTATTTCTATTCATTTTGCCTCTCCTACTCCATCCTCTTTGTTCTATTCCTCACCTTGTGTCTTAAATC"
        ),
        lines_a=(
            "140 T TTGTGTG 0|1",  # A  TR:23634  13880994 T TTGTGTG 0|1
            "155 TGTGTGC T 1|0",  # A  TR:23635  13881009 TGTGTGC T 1|0
            "161 C T 0|1",  # A  TR:23636  13881015 C T 0|1
        ),
        lines_b=(
            "141 TGTGTGTGTGTGTGTGTGTGC TGTGTGTGTGTGTGT,TGTGTGTGTGTGTGTGTGTGTGTGTGC 1/2",  # B  pvc-touch30:15740  13880995 TGTGTGTGTGTGTGTGTGTGC TGTGTGTGTGTGTGT,TGTGTGTGTGTGTGTGTGTGTGTGTGC 1/2
            "161 C T,* 1/2",  # B  pvc-touch30:15741  13881015 C T,* 1/2
        ),
        # locus    chr12:13880994-13881016 (excerpt 140-162), kind STR; STR regions: 13880995-13881014 str TG
        locus_1based=(140, 162),
        kind=Kind.STR,
        str_regions_1based=((141, 160, RegionType.STR, "TG"),),
        alleles=(
            "TTGTGTGTGTGTGTGTGTGTGCG",
            "TTGTGTGTGTGTGTGTG",
            "TTGTGTGTGTGTGTGTGTGTGTGTGTGTG",
        ),
        # possible genotypes  A: 1/2   B: 1/2
        possible_genotypes_a=frozenset({(1, 2)}),
        possible_genotypes_b=frozenset({(1, 2)}),
        # row  A 1/2  B 1/2  possible 1, 1  records 3, 2  -> same genotype
        outcome=Outcome.SAME_GENOTYPE,
        genotype_a=(1, 2),
        genotype_b=(1, 2),
        possible_a=1,
        possible_b=1,
        records_a=3,
        records_b=2,
        invalid_reason_a=None,
        invalid_reason_b=None,
    ),
    HardCase(
        name="E7",
        reference=Reference.GRCH38,
        # excerpt  chr12:13880855-13881155, 301 bases; excerpt position = position - 13880854; sample HG002
        contig="chr12",
        offset=13880854,
        sample="HG002",
        bases=(
            "TGCCTTCCTGTTTACTGTGTTCACTCTCATGCTGTGCAAAGGTTGTTCCACAGCCTATTTGAGACCTCGCATCTGG"
            "AGGCAGCTTTGTTTCAGTTTTCAAATTCTCTTGCCTTGCTTTGGCTGACAAAGGGGTATAAGGTTGTGTGTGTGTG"
            "TGTGTGTGCGCGTGCGCGCACGCATGTGCGCAGGTGTGCATGCACATGTGCACACGAGAGTATGTATAAAAAAGAT"
            "AACTTTATTTCTATTCATTTTGCCTCTCCTACTCCATCCTCTTTGTTCTATTCCTCACCTTGTGTCTTAAATC"
        ),
        lines_a=(
            "140 T TTGTGTG 0/1",  # A  GIAB-HG002:2834424  13880994 T TTGTGTG 0/1
            "155 TGTGTGC T 0/1",  # A  GIAB-HG002:2834425  13881009 TGTGTGC T 0/1
        ),
        lines_b=(
            "141 TGTGTGTGTGTGTGTGTGTGC TGTGTGTGTGTGTGT,TGTGTGTGTGTGTGTGTGTGTGTGTGC 1/2",  # B  pvc-touch30:15740  13880995 TGTGTGTGTGTGTGTGTGTGC TGTGTGTGTGTGTGT,TGTGTGTGTGTGTGTGTGTGTGTGTGC 1/2
            "161 C T,* 1/2",  # B  pvc-touch30:15741  13881015 C T,* 1/2
        ),
        # locus    chr12:13880994-13881016 (excerpt 140-162), kind STR; STR regions: 13880995-13881014 str TG
        locus_1based=(140, 162),
        kind=Kind.STR,
        str_regions_1based=((141, 160, RegionType.STR, "TG"),),
        alleles=(
            "TTGTGTGTGTGTGTGTGTGTGCG",
            "TTGTGTGTGTGTGTGTGTGTGTG",
            "TTGTGTGTGTGTGTGTG",
            "TTGTGTGTGTGTGTGTGTGTGTGTGTGCG",
            "TTGTGTGTGTGTGTGTGTGTGTGTGTGTG",
        ),
        # possible genotypes  A: 0/1, 2/3   B: 2/4
        possible_genotypes_a=frozenset({(0, 1), (2, 3)}),
        possible_genotypes_b=frozenset({(2, 4)}),
        # row  A 0/1  B 2/4  possible 2, 1  records 2, 2  -> different
        outcome=Outcome.DIFFERENT,
        genotype_a=(0, 1),
        genotype_b=(2, 4),
        possible_a=2,
        possible_b=1,
        records_a=2,
        records_b=2,
        invalid_reason_a=None,
        invalid_reason_b=None,
    ),
    HardCase(
        name="E8",
        reference=Reference.GRCH38,
        # excerpt  chr4:98053390-98053690, 301 bases; excerpt position = position - 98053389; sample HG002
        contig="chr4",
        offset=98053389,
        sample="HG002",
        bases=(
            "TAGAGTATAGGCATAGTCATTCCATATCCTTTAATATTCTGGATTCTATAGACACTTCCAGCAGAGACAAACTCTA"
            "TTAGAGCCATTAACTATAGCTATATAAAACCTTTAAGAATCAATCTCTAAATACTTATAAAAGCAAATAACAAAAA"
            "AAAAATCTAAATTTGAGTAAAAGTAAACACAGTTTCTCTGACTTTAGGCTTTTATTTAGTGCCTGTCTTCAAATGT"
            "ACTTCCTCCTTCCATTTTCTTCTGTACATCAGTTTAGATGACCAATAATTTTATCTTTGTCCTTCACAAGAAT"
        ),
        lines_a=(
            "144 TAAC T 1|0",  # A  TR:8408  98053533 TAAC T 1|0
            "150 A T 1|0",  # A  TR:8409  98053539 A T 1|0
        ),
        lines_b=(
            "144 TAAC T,* 1/2",  # B  pvc-touch30:7554  98053533 TAAC T,* 1/2
            "147 CAAAAAAAAAA TAAAAAAA 0/1",  # B  pvc-touch30:7555  98053536 CAAAAAAAAAA TAAAAAAA 0/1
        ),
        # locus    chr4:98053533-98053546 (excerpt 144-157), kind STR; STR regions: 98053537-98053546 str A
        locus_1based=(144, 157),
        kind=Kind.STR,
        str_regions_1based=((148, 157, RegionType.STR, "A"),),
        alleles=(
            "TAACAAAAAAAAAA",
            "TAATAAAAAAA",
            "TAAAAAAAAAA",
        ),
        # possible genotypes  A: 0/1   B: 1/2
        possible_genotypes_a=frozenset({(0, 1)}),
        possible_genotypes_b=frozenset({(1, 2)}),
        # row  A 0/1  B 1/2  possible 1, 1  records 2, 2  -> different
        outcome=Outcome.DIFFERENT,
        genotype_a=(0, 1),
        genotype_b=(1, 2),
        possible_a=1,
        possible_b=1,
        records_a=2,
        records_b=2,
        invalid_reason_a=None,
        invalid_reason_b=None,
    ),
    HardCase(
        name="E9",
        reference=Reference.GRCH38,
        # excerpt  chr4:98053390-98053690, 301 bases; excerpt position = position - 98053389; sample HG002
        contig="chr4",
        offset=98053389,
        sample="HG002",
        bases=(
            "TAGAGTATAGGCATAGTCATTCCATATCCTTTAATATTCTGGATTCTATAGACACTTCCAGCAGAGACAAACTCTA"
            "TTAGAGCCATTAACTATAGCTATATAAAACCTTTAAGAATCAATCTCTAAATACTTATAAAAGCAAATAACAAAAA"
            "AAAAATCTAAATTTGAGTAAAAGTAAACACAGTTTCTCTGACTTTAGGCTTTTATTTAGTGCCTGTCTTCAAATGT"
            "ACTTCCTCCTTCCATTTTCTTCTGTACATCAGTTTAGATGACCAATAATTTTATCTTTGTCCTTCACAAGAAT"
        ),
        lines_a=(
            "144 TAAC T 1|0",  # A  TR:8408  98053533 TAAC T 1|0
            "150 A T 1|0",  # A  TR:8409  98053539 A T 1|0
        ),
        lines_b=(
            "144 TAAC T 0/1",  # B  pvc-base30:7565  98053533 TAAC T 0/1
            "148 AAAAAAAAAA AATAAAAAAA 0/1",  # B  pvc-base30:7566  98053537 AAAAAAAAAA AATAAAAAAA 0/1
        ),
        # locus    chr4:98053533-98053546 (excerpt 144-157), kind STR; STR regions: 98053537-98053546 str A
        locus_1based=(144, 157),
        kind=Kind.STR,
        str_regions_1based=((148, 157, RegionType.STR, "A"),),
        alleles=(
            "TAACAAAAAAAAAA",
            "TAATAAAAAAA",
            "TAACAATAAAAAAA",
            "TAAAAAAAAAA",
        ),
        # possible genotypes  A: 0/1   B: 0/1, 2/3
        possible_genotypes_a=frozenset({(0, 1)}),
        possible_genotypes_b=frozenset({(0, 1), (2, 3)}),
        # row  A 0/1  B 0/1  possible 1, 2  records 2, 2  -> same genotype
        outcome=Outcome.SAME_GENOTYPE,
        genotype_a=(0, 1),
        genotype_b=(0, 1),
        possible_a=1,
        possible_b=2,
        records_a=2,
        records_b=2,
        invalid_reason_a=None,
        invalid_reason_b=None,
    ),
    HardCase(
        name="E10",
        reference=Reference.GRCH38,
        # excerpt  chr16:5228700-5229000, 301 bases; excerpt position = position - 5228699; sample HG002
        contig="chr16",
        offset=5228699,
        sample="HG002",
        bases=(
            "GGAGTAAGAGAGAGGGAGGGAGGGAAACGGAGGAAGGAAGAGAGAGAAACAGAGGGAGAGAGAGAGAAAGGGAGGG"
            "AGAGAGAAACAAAGATGCAGAAGCAGGGAGAGAGAGAAAGAGACAAGGGAGGGAGGAAGAGAGCGAGAGAGAGAGA"
            "GAGAGAGAGAGAGAGAGAGAGAGAGAGAGAAAGGAGAGAAGGAGAGAGACAGAGAGGGGCTGATAGAATGCTAACC"
            "CCTGCCTAGATTAGGGGTCAAATGCCAATGCCACAGGTTCATGTTGGGGAGACACACTTTTTTTTTTGCCAGC"
        ),
        lines_a=(
            "133 AAGAGAGCGAGAGAGAGAGAG A 1|1",  # A  TR:29998  5228832 AAGAGAGCGAGAGAGAGAGAG A 1|1
            "160 A T 1|1",  # A  TR:29999  5228859 A T 1|1
        ),
        lines_b=(
            "133 AAGAGAGC A,* 1/2",  # B  pvc-touch30:18978  5228832 AAGAGAGC A,* 1/2
            "140 CGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGA CGAGAGAGAGAGAGAGAGAGAGA,TGAGAGAGAGAGAGAGAGAGAGA 1/2",  # B  pvc-touch30:18979  5228839 CGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGA CGAGAGAGAGAGAGAGAGAGAGA,TGAGAGAGAGAGAGAGAGAGAGA 1/2
        ),
        # locus    chr16:5228832-5228881 (excerpt 133-182), kind STR; STR regions: 5228840-5228881 str GA
        locus_1based=(133, 182),
        kind=Kind.STR,
        str_regions_1based=((141, 182, RegionType.STR, "GA"),),
        alleles=(
            "AAGAGAGCGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGA",
            "AAGAGAGTGAGAGAGAGAGAGAGAGAGAGA",
            "AGAGAGAGAGAGAGAGAGAGAGA",
        ),
        # possible genotypes  A: 1/1   B: 1/2
        possible_genotypes_a=frozenset({(1, 1)}),
        possible_genotypes_b=frozenset({(1, 2)}),
        # row  A 1/1  B 1/2  possible 1, 1  records 2, 2  -> different
        outcome=Outcome.DIFFERENT,
        genotype_a=(1, 1),
        genotype_b=(1, 2),
        possible_a=1,
        possible_b=1,
        records_a=2,
        records_b=2,
        invalid_reason_a=None,
        invalid_reason_b=None,
    ),
    HardCase(
        name="E11",
        reference=Reference.GRCH38,
        # excerpt  chr1:64355880-64356180, 301 bases; excerpt position = position - 64355879; sample HG002
        contig="chr1",
        offset=64355879,
        sample="HG002",
        bases=(
            "TCTCTCCAAAAGTTACAAAAACTTCTCCCTACTACAAATAATAGGAATTTAAAAAGCTTCCATTTATTAACTATCT"
            "ACTATGTACCAGATACTTTTCACTCTTCCATCATTCATTCACTAATTCATTCCATATATATATATATATATATATA"
            "TATATATATATATATATATATAAAAGTTCCATGTGTTATGTCAGTATACGGTGGATTAAACAAAGTTAGAAGCCAT"
            "CCCTTCTCTCCTAGAGATGGCTCTTTTATAAGGGAGGCAAACAGTAAAGAAGAAAACAGAGACCCATATTAGG"
        ),
        lines_a=(
            "151 TATATATATATATATATATATATAA T 0|1",  # A  TR:932  64356030 TATATATATATATATATATATATAA T 0|1
            "161 TATATATATATATAA T 1|0",  # A  TR:933  64356040 TATATATATATATAA T 1|0
        ),
        lines_b=(
            "130 ATATATATATATATATATATATATATATATATATATATATATATAA ATATATATATATATATATATATAA,ATATATATATATATATATATATATATATATAT 1/2",  # B  pvc-touch30:3169  64356009 ATATATATATATATATATATATATATATATATATATATATATATAA ATATATATATATATATATATATAA,ATATATATATATATATATATATATATATATAT 1/2
            "173 TAA T,* 1/2",  # B  pvc-touch30:3170  64356052 TAA T,* 1/2
        ),
        # locus    chr1:64356008-64356056 (excerpt 129-177), kind STR; STR regions: 64356009-64356052 str AT
        locus_1based=(129, 177),
        kind=Kind.STR,
        str_regions_1based=((130, 173, RegionType.STR, "AT"),),
        alleles=(
            "CATATATATATATATATATATATATATATATATATATATATATATAAAA",
            "CATATATATATATATATATATATATATATATATAA",
            "CATATATATATATATATATATATAA",
        ),
        # possible genotypes  A: 1/2   B: none
        possible_genotypes_a=frozenset({(1, 2)}),
        possible_genotypes_b=None,
        # row  A 1/2  B .  possible 1, 0  records 2, 2  -> invalid in B, contradicting records
        outcome=Outcome.INVALID_IN_B,
        genotype_a=(1, 2),
        genotype_b=None,
        possible_a=1,
        possible_b=0,
        records_a=2,
        records_b=2,
        invalid_reason_a=None,
        invalid_reason_b=InvalidReason.CONTRADICTING_RECORDS,
    ),
    HardCase(
        name="E12",
        reference=Reference.SL4_0,
        # excerpt  SL4.0ch05:576000-576300, 301 bases; excerpt position = position - 575999; sample SRS3394689
        contig="SL4.0ch05",
        offset=575999,
        sample="SRS3394689",
        bases=(
            "TTGTTCTTTTCTTTTTCTCTAGATTCACCAAGCATCATTCTGTTTGCCATCATTCGCCAAGAGGGTGATCGAAACC"
            "AAAGGAAAATAATTTCGAAGGAGTGCTGAGTTAGTTGTACTACGCCCTATAGTTATCGCAGTTCTTTTTTTTTTTT"
            "ATAAATTAAAGACAATGTTTTTAGGTTATTTAAGCCAAATGGGTTGAATGATGGATTATAATTTCATCTTTGATAT"
            "AATATTCTAGTTTGTTTATCTTGATGATTTTCATTCACTTTACATTATAGAGTTGGACAACATTTTGGTATTT"
        ),
        lines_a=(
            "152 TA T 1/1",  # A  GATK-tomato:65098  576151 TA T 1/1
            "153 A T,* 2/2",  # A  GATK-tomato:65099  576152 A T,* 2/2
        ),
        lines_b=(
            "141 TTTTTTTTTTTTA TTTTTTTTTTT 0/0",  # B  pvc-tomato:71841  576140 TTTTTTTTTTTTA TTTTTTTTTTT 0/0
            "153 ATA TA,TTA,* 1/1",  # B  pvc-tomato:71842  576152 ATA TA,TTA,* 1/1
        ),
        # locus    SL4.0ch05:576140-576154 (excerpt 141-155), kind STR; STR regions: 576140-576151 str T
        locus_1based=(141, 155),
        kind=Kind.STR,
        str_regions_1based=((141, 152, RegionType.STR, "T"),),
        alleles=(
            "TTTTTTTTTTTTATA",
            "TTTTTTTTTTTTTA",
        ),
        # possible genotypes  A: 1/1   B: 1/1
        possible_genotypes_a=frozenset({(1, 1)}),
        possible_genotypes_b=frozenset({(1, 1)}),
        # row  A 1/1  B 1/1  possible 1, 1  records 2, 2  -> same genotype
        outcome=Outcome.SAME_GENOTYPE,
        genotype_a=(1, 1),
        genotype_b=(1, 1),
        possible_a=1,
        possible_b=1,
        records_a=2,
        records_b=2,
        invalid_reason_a=None,
        invalid_reason_b=None,
    ),
    HardCase(
        name="E13",
        reference=Reference.SL4_0,
        # excerpt  SL4.0ch05:576000-576300, 301 bases; excerpt position = position - 575999; sample SRS3394663
        contig="SL4.0ch05",
        offset=575999,
        sample="SRS3394663",
        bases=(
            "TTGTTCTTTTCTTTTTCTCTAGATTCACCAAGCATCATTCTGTTTGCCATCATTCGCCAAGAGGGTGATCGAAACC"
            "AAAGGAAAATAATTTCGAAGGAGTGCTGAGTTAGTTGTACTACGCCCTATAGTTATCGCAGTTCTTTTTTTTTTTT"
            "ATAAATTAAAGACAATGTTTTTAGGTTATTTAAGCCAAATGGGTTGAATGATGGATTATAATTTCATCTTTGATAT"
            "AATATTCTAGTTTGTTTATCTTGATGATTTTCATTCACTTTACATTATAGAGTTGGACAACATTTTGGTATTT"
        ),
        lines_a=(
            "152 TA T 1/1",  # A  GATK-tomato:65098  576151 TA T 1/1
            "153 A T,* 2/2",  # A  GATK-tomato:65099  576152 A T,* 2/2
        ),
        lines_b=(
            "141 TTTTTTTTTTTTA TTTTTTTTTTT 0/0",  # B  pvc-tomato:71841  576140 TTTTTTTTTTTTA TTTTTTTTTTT 0/0
            "153 ATA TA,TTA,* 0/1",  # B  pvc-tomato:71842  576152 ATA TA,TTA,* 0/1
        ),
        # locus    SL4.0ch05:576140-576154 (excerpt 141-155), kind STR; STR regions: 576140-576151 str T
        locus_1based=(141, 155),
        kind=Kind.STR,
        str_regions_1based=((141, 152, RegionType.STR, "T"),),
        alleles=(
            "TTTTTTTTTTTTATA",
            "TTTTTTTTTTTTTA",
        ),
        # possible genotypes  A: 1/1   B: 0/1
        possible_genotypes_a=frozenset({(1, 1)}),
        possible_genotypes_b=frozenset({(0, 1)}),
        # row  A 1/1  B 0/1  possible 1, 1  records 2, 2  -> same alleles, other genotype
        outcome=Outcome.SAME_ALLELES,
        genotype_a=(1, 1),
        genotype_b=(0, 1),
        possible_a=1,
        possible_b=1,
        records_a=2,
        records_b=2,
        invalid_reason_a=None,
        invalid_reason_b=None,
    ),
    HardCase(
        name="E14",
        reference=Reference.SL4_0,
        # excerpt  SL4.0ch05:576000-576300, 301 bases; excerpt position = position - 575999; sample SRS3394701
        contig="SL4.0ch05",
        offset=575999,
        sample="SRS3394701",
        bases=(
            "TTGTTCTTTTCTTTTTCTCTAGATTCACCAAGCATCATTCTGTTTGCCATCATTCGCCAAGAGGGTGATCGAAACC"
            "AAAGGAAAATAATTTCGAAGGAGTGCTGAGTTAGTTGTACTACGCCCTATAGTTATCGCAGTTCTTTTTTTTTTTT"
            "ATAAATTAAAGACAATGTTTTTAGGTTATTTAAGCCAAATGGGTTGAATGATGGATTATAATTTCATCTTTGATAT"
            "AATATTCTAGTTTGTTTATCTTGATGATTTTCATTCACTTTACATTATAGAGTTGGACAACATTTTGGTATTT"
        ),
        lines_a=(
            "152 TA T 0/0",  # A  GATK-tomato:65098  576151 TA T 0/0
            "153 A T,* 1/1",  # A  GATK-tomato:65099  576152 A T,* 1/1
        ),
        lines_b=(
            "141 TTTTTTTTTTTTA TTTTTTTTTTT 0/0",  # B  pvc-tomato:71841  576140 TTTTTTTTTTTTA TTTTTTTTTTT 0/0
            "153 ATA TA,TTA,* 0/2",  # B  pvc-tomato:71842  576152 ATA TA,TTA,* 0/2
        ),
        # locus    SL4.0ch05:576140-576154 (excerpt 141-155), kind STR; STR regions: 576140-576151 str T
        locus_1based=(141, 155),
        kind=Kind.STR,
        str_regions_1based=((141, 152, RegionType.STR, "T"),),
        alleles=(
            "TTTTTTTTTTTTATA",
            "TTTTTTTTTTTTTTA",
        ),
        # possible genotypes  A: 1/1   B: 0/1
        possible_genotypes_a=frozenset({(1, 1)}),
        possible_genotypes_b=frozenset({(0, 1)}),
        # row  A 1/1  B 0/1  possible 1, 1  records 2, 2  -> same alleles, other genotype
        outcome=Outcome.SAME_ALLELES,
        genotype_a=(1, 1),
        genotype_b=(0, 1),
        possible_a=1,
        possible_b=1,
        records_a=2,
        records_b=2,
        invalid_reason_a=None,
        invalid_reason_b=None,
    ),
    HardCase(
        name="E15",
        reference=Reference.SL4_0,
        # excerpt  SL4.0ch05:58429278-58429578, 301 bases; excerpt position = position - 58429277; sample SRS3394663
        contig="SL4.0ch05",
        offset=58429277,
        sample="SRS3394663",
        bases=(
            "AATTTCAAAATGGTGCATCAATTTTTTTTAAAACCAAAAAGGTAAAATCACTCCTGACAGAGTGATTTTTTTCTGA"
            "CAAAATTGATGTGTTAATTTTTTTTAAAAAAATACAAAAATCGCTGCTGTGACAGCGATTTCATTCAAAAAAAAAT"
            "AAATAAATTAACAAAGCTGGCCCTGCCTCTGTTGCTGGCTTGCTGCAAAGAGACGAAGACCAAAATGCCTGCCTCT"
            "ATTAACTTTTTTAAGGAAAAGCAAATCGCTACAATATTTTAATTTGATTTTTTTTTTTTTGAAAAACGCTTCA"
        ),
        lines_a=(
            "151 ATAAAT A 1/1",  # A  GATK-tomato:74913  58429428 ATAAAT A 1/1
            "152 T A,* 2/2",  # A  GATK-tomato:74914  58429429 T A,* 2/2
            "156 T A,* 2/2",  # A  GATK-tomato:74915  58429433 T A,* 2/2
        ),
        lines_b=(
            "143 AAAAAAAAAT AAAAAAAA,AAAAAT,AAAAAAA 0/2",  # B  pvc-tomato:83029  58429420 AAAAAAAAAT AAAAAAAA,AAAAAT,AAAAAAA 0/2
            "152 TAAATAAA AAAAAAA,AAA,* 2/2",  # B  pvc-tomato:83030  58429429 TAAATAAA AAAAAAA,AAA,* 2/2
        ),
        # locus    SL4.0ch05:58429419-58429436 (excerpt 142-159), kind STR; STR regions: 58429420-58429428 str A
        locus_1based=(142, 159),
        kind=Kind.STR,
        str_regions_1based=((143, 151, RegionType.STR, "A"),),
        alleles=(
            "CAAAAAAAAATAAATAAA",
            "CAAAAAAAAAAAA",
            "CAAAAAAAA",
        ),
        # possible genotypes  A: 1/1   B: 1/2
        possible_genotypes_a=frozenset({(1, 1)}),
        possible_genotypes_b=frozenset({(1, 2)}),
        # row  A 1/1  B 1/2  possible 1, 1  records 3, 2  -> different
        outcome=Outcome.DIFFERENT,
        genotype_a=(1, 1),
        genotype_b=(1, 2),
        possible_a=1,
        possible_b=1,
        records_a=3,
        records_b=2,
        invalid_reason_a=None,
        invalid_reason_b=None,
    ),
    HardCase(
        name="E16",
        reference=Reference.SL4_0,
        # excerpt  SL4.0ch05:58429278-58429578, 301 bases; excerpt position = position - 58429277; sample SRS3394698
        contig="SL4.0ch05",
        offset=58429277,
        sample="SRS3394698",
        bases=(
            "AATTTCAAAATGGTGCATCAATTTTTTTTAAAACCAAAAAGGTAAAATCACTCCTGACAGAGTGATTTTTTTCTGA"
            "CAAAATTGATGTGTTAATTTTTTTTAAAAAAATACAAAAATCGCTGCTGTGACAGCGATTTCATTCAAAAAAAAAT"
            "AAATAAATTAACAAAGCTGGCCCTGCCTCTGTTGCTGGCTTGCTGCAAAGAGACGAAGACCAAAATGCCTGCCTCT"
            "ATTAACTTTTTTAAGGAAAAGCAAATCGCTACAATATTTTAATTTGATTTTTTTTTTTTTGAAAAACGCTTCA"
        ),
        lines_a=(
            "151 ATAAAT A 0/0",  # A  GATK-tomato:74913  58429428 ATAAAT A 0/0
            "152 T A,* 1|1:58429429",  # A  GATK-tomato:74914  58429429 T A,* 1|1:58429429
            "156 T A,* 1|1:58429429",  # A  GATK-tomato:74915  58429433 T A,* 1|1:58429429
        ),
        lines_b=(
            "143 AAAAAAAAAT AAAAAAAA,AAAAAT,AAAAAAA 0/0",  # B  pvc-tomato:83029  58429420 AAAAAAAAAT AAAAAAAA,AAAAAT,AAAAAAA 0/0
            "152 TAAATAAA AAAAAAA,AAA,* 2/2",  # B  pvc-tomato:83030  58429429 TAAATAAA AAAAAAA,AAA,* 2/2
        ),
        # locus    SL4.0ch05:58429419-58429436 (excerpt 142-159), kind STR; STR regions: 58429420-58429428 str A
        locus_1based=(142, 159),
        kind=Kind.STR,
        str_regions_1based=((143, 151, RegionType.STR, "A"),),
        alleles=(
            "CAAAAAAAAATAAATAAA",
            "CAAAAAAAAAAAAAAAAA",
            "CAAAAAAAAAAAA",
        ),
        # possible genotypes  A: 1/1   B: 2/2
        possible_genotypes_a=frozenset({(1, 1)}),
        possible_genotypes_b=frozenset({(2, 2)}),
        # row  A 1/1  B 2/2  possible 1, 1  records 3, 2  -> different
        outcome=Outcome.DIFFERENT,
        genotype_a=(1, 1),
        genotype_b=(2, 2),
        possible_a=1,
        possible_b=1,
        records_a=3,
        records_b=2,
        invalid_reason_a=None,
        invalid_reason_b=None,
    ),
    HardCase(
        name="H1",
        reference=Reference.GRCH38,
        # excerpt  chr1:113339083-113339383, 301 bases; excerpt position = position - 113339082; sample HG002
        contig="chr1",
        offset=113339082,
        sample="HG002",
        bases=(
            "ACAATAGTTATTGTGTTTGACATTATGGATCAGGGTTTCTGAATTTCAGTGAAGCATTTTACAATCTTTTATATTT"
            "GTTACTAAAGCACATAAATATGGATTCCTTGGAGGAAGACACAGGAGGGTGTGTGTATGTGTATGTGTGTGTGTGT"
            "GTGTGTGTGTGTGTGTGTGTGTGTGTTTACATCATCTGTATCATAAGCACCAAGCATAGTGCCTGGCACATGATAG"
            "ATGACAAGTAAAAGTTTAATGAATAAATACCTGGGTTGGTTGAACAGGTATATGTTTTTTAAAAGTGTGCACT"
        ),
        lines_a=(
            "139 ATGTGTGTG A 1|0",  # A  TR:1458  113339221 ATGTGTGTG A 1|0
        ),
        lines_b=(
            "125 GTGTGTGTATGTGTATGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGT GTGTGTGTATGTGTATGTGTGTGTGTGTGTGTGTGTGTGTGTGTGT 0|1",  # B  HipSTR30:3213  113339207 GTGTGTGTATGTGTATGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGT GTGTGTGTATGTGTATGTGTGTGTGTGTGTGTGTGTGTGTGTGTGT 0|1
        ),
        # locus    chr1:113339221-113339260 (excerpt 139-178), kind STR; STR regions: 113339222-113339259 str TG
        locus_1based=(139, 178),
        kind=Kind.STR,
        str_regions_1based=((140, 177, RegionType.STR, "TG"),),
        alleles=(
            "ATGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGT",
            "ATGTGTGTGTGTGTGTGTGTGTGTGTGTGTGT",
        ),
        # possible genotypes  A: 0/1   B: 0/1
        possible_genotypes_a=frozenset({(0, 1)}),
        possible_genotypes_b=frozenset({(0, 1)}),
        # row  A 0/1  B 0/1  possible 1, 1  records 1, 1  -> same genotype
        outcome=Outcome.SAME_GENOTYPE,
        genotype_a=(0, 1),
        genotype_b=(0, 1),
        possible_a=1,
        possible_b=1,
        records_a=1,
        records_b=1,
        invalid_reason_a=None,
        invalid_reason_b=None,
    ),
    HardCase(
        name="H2",
        reference=Reference.GRCH38,
        # excerpt  chr1:113339083-113339383, 301 bases; excerpt position = position - 113339082; sample HG002
        contig="chr1",
        offset=113339082,
        sample="HG002",
        bases=(
            "ACAATAGTTATTGTGTTTGACATTATGGATCAGGGTTTCTGAATTTCAGTGAAGCATTTTACAATCTTTTATATTT"
            "GTTACTAAAGCACATAAATATGGATTCCTTGGAGGAAGACACAGGAGGGTGTGTGTATGTGTATGTGTGTGTGTGT"
            "GTGTGTGTGTGTGTGTGTGTGTGTGTTTACATCATCTGTATCATAAGCACCAAGCATAGTGCCTGGCACATGATAG"
            "ATGACAAGTAAAAGTTTAATGAATAAATACCTGGGTTGGTTGAACAGGTATATGTTTTTTAAAAGTGTGCACT"
        ),
        lines_a=(
            "125 GTGTGTGTATGTGTATGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGT GTGTGTGTATGTGTATGTGTGTGTGTGTGTGTGTGTGTGTGTGTGT 0|1",  # A  HipSTR30:3213  113339207 GTGTGTGTATGTGTATGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGT GTGTGTGTATGTGTATGTGTGTGTGTGTGTGTGTGTGTGTGTGTGT 0|1
        ),
        lines_b=(
            "140 TGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTG TGTGTGTGTGTGTGTGTGTGTGTGTGTGTG 0/1",  # B  pvc-m07-30:3467  113339222 TGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTG TGTGTGTGTGTGTGTGTGTGTGTGTGTGTG 0/1
        ),
        # locus    chr1:113339221-113339260 (excerpt 139-178), kind STR; STR regions: 113339222-113339259 str TG
        locus_1based=(139, 178),
        kind=Kind.STR,
        str_regions_1based=((140, 177, RegionType.STR, "TG"),),
        alleles=(
            "ATGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGT",
            "ATGTGTGTGTGTGTGTGTGTGTGTGTGTGTGT",
        ),
        # possible genotypes  A: 0/1   B: 0/1
        possible_genotypes_a=frozenset({(0, 1)}),
        possible_genotypes_b=frozenset({(0, 1)}),
        # row  A 0/1  B 0/1  possible 1, 1  records 1, 1  -> same genotype
        outcome=Outcome.SAME_GENOTYPE,
        genotype_a=(0, 1),
        genotype_b=(0, 1),
        possible_a=1,
        possible_b=1,
        records_a=1,
        records_b=1,
        invalid_reason_a=None,
        invalid_reason_b=None,
    ),
    HardCase(
        name="H3",
        reference=Reference.GRCH38,
        # excerpt  chr1:72941135-72941435, 301 bases; excerpt position = position - 72941134; sample HG002
        contig="chr1",
        offset=72941134,
        sample="HG002",
        bases=(
            "TGTTATTTACCAATAATAATGTGGTTCACAAAAATGTGAATAATCTAGACTCTCGAAGATAAATTTTATTCATCAT"
            "GTATATGTCTCACCTAGTGGAGCTGTTGTACTTGAAAAACAGTCTATCATCTATCTATCTATCTATCTATCTATCT"
            "ATCTATCTATCCATCTATCTATCTATCATTTTTCATCCATCCATACACCTATCATTCCATCCATCACTTAATTTAT"
            "CCATGATATTCTTTAGTGAGTAAAGCTAATACCTGTGCCATTCTTGGATTTTTCACAGCTTATTTTCACTATT"
        ),
        lines_a=(
            "160 T TATCTATCC 1|1",  # A  TR:1010  72941294 T TATCTATCC 1|1
            "164 C T 1|1",  # A  TR:1011  72941298 C T 1|1
        ),
        lines_b=(
            "126 TCTATCTATCTATCTATCTATCTATCTATCTATCTATCCATCTATCTATCTA TCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCCATCTATCTATCTATCTA 1|1",  # B  HipSTR30:3016  72941260 TCTATCTATCTATCTATCTATCTATCTATCTATCTATCCATCTATCTATCTA TCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCCATCTATCTATCTATCTA 1|1
        ),
        # locus    chr1:72941245-72941319 (excerpt 111-185), kind STR; STR regions: 72941245-72941319 str-bundle
        locus_1based=(111, 185),
        kind=Kind.STR,
        str_regions_1based=((111, 185, RegionType.STR_BUNDLE, None),),
        alleles=(
            "AAAAACAGTCTATCATCTATCTATCTATCTATCTATCTATCTATCTATCTATCCATCTATCTATCTATCATTTTT",
            "AAAAACAGTCTATCATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCCATCTATCTATCTATCTATCATTTTT",
        ),
        # possible genotypes  A: 1/1   B: 1/1
        possible_genotypes_a=frozenset({(1, 1)}),
        possible_genotypes_b=frozenset({(1, 1)}),
        # row  A 1/1  B 1/1  possible 1, 1  records 2, 1  -> same genotype
        outcome=Outcome.SAME_GENOTYPE,
        genotype_a=(1, 1),
        genotype_b=(1, 1),
        possible_a=1,
        possible_b=1,
        records_a=2,
        records_b=1,
        invalid_reason_a=None,
        invalid_reason_b=None,
    ),
    HardCase(
        name="H4",
        reference=Reference.GRCH38,
        # excerpt  chr1:72941135-72941435, 301 bases; excerpt position = position - 72941134; sample HG002
        contig="chr1",
        offset=72941134,
        sample="HG002",
        bases=(
            "TGTTATTTACCAATAATAATGTGGTTCACAAAAATGTGAATAATCTAGACTCTCGAAGATAAATTTTATTCATCAT"
            "GTATATGTCTCACCTAGTGGAGCTGTTGTACTTGAAAAACAGTCTATCATCTATCTATCTATCTATCTATCTATCT"
            "ATCTATCTATCCATCTATCTATCTATCATTTTTCATCCATCCATACACCTATCATTCCATCCATCACTTAATTTAT"
            "CCATGATATTCTTTAGTGAGTAAAGCTAATACCTGTGCCATTCTTGGATTTTTCACAGCTTATTTTCACTATT"
        ),
        lines_a=(
            "160 T TATCTATCC 1|1",  # A  TR:1010  72941294 T TATCTATCC 1|1
            "164 C T 1|1",  # A  TR:1011  72941298 C T 1|1
        ),
        lines_b=(
            "163 C CTATCCAT 1/1",  # B  GIAB-HG002:98215  72941297 C CTATCCAT 1/1
            "164 C CT 1/1",  # B  GIAB-HG002:98216  72941298 C CT 1/1
        ),
        # locus    chr1:72941245-72941319 (excerpt 111-185), kind STR; STR regions: 72941245-72941319 str-bundle
        locus_1based=(111, 185),
        kind=Kind.STR,
        str_regions_1based=((111, 185, RegionType.STR_BUNDLE, None),),
        alleles=(
            "AAAAACAGTCTATCATCTATCTATCTATCTATCTATCTATCTATCTATCTATCCATCTATCTATCTATCATTTTT",
            "AAAAACAGTCTATCATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCCATCTATCTATCTATCTATCATTTTT",
        ),
        # possible genotypes  A: 1/1   B: 1/1
        possible_genotypes_a=frozenset({(1, 1)}),
        possible_genotypes_b=frozenset({(1, 1)}),
        # row  A 1/1  B 1/1  possible 1, 1  records 2, 2  -> same genotype
        outcome=Outcome.SAME_GENOTYPE,
        genotype_a=(1, 1),
        genotype_b=(1, 1),
        possible_a=1,
        possible_b=1,
        records_a=2,
        records_b=2,
        invalid_reason_a=None,
        invalid_reason_b=None,
    ),
    HardCase(
        name="H5",
        reference=Reference.GRCH38,
        # excerpt  chr1:80541386-80541686, 301 bases; excerpt position = position - 80541385; sample HG002
        contig="chr1",
        offset=80541385,
        sample="HG002",
        bases=(
            "GTTTGAGGCTGAAATGCCAGAATTATAATATTTAAAACAAAAGATAATCATTTTATTGTAGTTAAATAATACATTT"
            "AGTTCTAAATATAAGACACACACATACACACACAGACAGATAGATGGATGGATGGATGGATAGACAGATAGATAGA"
            "TAGATAGATAGATAGATAGGTAGATAGATAGATAGATATTCTCCTGGAATATTTTTAGGCAAGAAAAATAAGTATA"
            "ACTAAAATACCTGAAACTCTGAAAGCTTGCTTTGCCCAGAAAAAGGAAGACGATGGTAGAAAGGATCAGGAAA"
        ),
        lines_a=(
            "141 C CAGATAGATAGATAGATAGATAGATAGAT 1|0",  # A  TR:1070  80541526 C CAGATAGATAGATAGATAGATAGATAGAT 1|0
            "172 G A 1|0",  # A  TR:1071  80541557 G A 1|0
        ),
        lines_b=(
            "114 AGATAGATGGATGGATGGATGGATAGACAGATAGATAGATAGATAGATAGATAGATAGGTAGATAGATAGATAGAT AGATAGATGGATGGATGGATGGATAGACAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGAT 0|1",  # B  HipSTR30:3058  80541499 AGATAGATGGATGGATGGATGGATAGACAGATAGATAGATAGATAGATAGATAGATAGGTAGATAGATAGATAGAT AGATAGATGGATGGATGGATGGATAGACAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGAT 0|1
        ),
        # locus    chr1:80541526-80541557 (excerpt 141-172), kind complex; STR regions: none
        locus_1based=(141, 172),
        kind=Kind.COMPLEX,
        str_regions_1based=(),
        alleles=(
            "CAGATAGATAGATAGATAGATAGATAGATAGG",
            "CAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGA",
        ),
        # possible genotypes  A: 0/1   B: 0/1
        possible_genotypes_a=frozenset({(0, 1)}),
        possible_genotypes_b=frozenset({(0, 1)}),
        # row  A 0/1  B 0/1  possible 1, 1  records 2, 1  -> same genotype
        outcome=Outcome.SAME_GENOTYPE,
        genotype_a=(0, 1),
        genotype_b=(0, 1),
        possible_a=1,
        possible_b=1,
        records_a=2,
        records_b=1,
        invalid_reason_a=None,
        invalid_reason_b=None,
    ),
    HardCase(
        name="H6",
        reference=Reference.GRCH38,
        # excerpt  chr1:80541386-80541686, 301 bases; excerpt position = position - 80541385; sample HG002
        contig="chr1",
        offset=80541385,
        sample="HG002",
        bases=(
            "GTTTGAGGCTGAAATGCCAGAATTATAATATTTAAAACAAAAGATAATCATTTTATTGTAGTTAAATAATACATTT"
            "AGTTCTAAATATAAGACACACACATACACACACAGACAGATAGATGGATGGATGGATGGATAGACAGATAGATAGA"
            "TAGATAGATAGATAGATAGGTAGATAGATAGATAGATATTCTCCTGGAATATTTTTAGGCAAGAAAAATAAGTATA"
            "ACTAAAATACCTGAAACTCTGAAAGCTTGCTTTGCCCAGAAAAAGGAAGACGATGGTAGAAAGGATCAGGAAA"
        ),
        lines_a=(
            "141 C CAGATAGATAGATAGATAGATAGATAGAT 0/1",  # A  GIAB-HG002:109203  80541526 C CAGATAGATAGATAGATAGATAGATAGAT 0/1
            "172 G A 0/1",  # A  GIAB-HG002:109204  80541557 G A 0/1
        ),
        lines_b=(
            "114 AGATAGATGGATGGATGGATGGATAGACAGATAGATAGATAGATAGATAGATAGATAGGTAGATAGATAGATAGAT AGATAGATGGATGGATGGATGGATAGACAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGAT 0|1",  # B  HipSTR30:3058  80541499 AGATAGATGGATGGATGGATGGATAGACAGATAGATAGATAGATAGATAGATAGATAGGTAGATAGATAGATAGAT AGATAGATGGATGGATGGATGGATAGACAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGAT 0|1
        ),
        # locus    chr1:80541526-80541557 (excerpt 141-172), kind complex; STR regions: none
        locus_1based=(141, 172),
        kind=Kind.COMPLEX,
        str_regions_1based=(),
        alleles=(
            "CAGATAGATAGATAGATAGATAGATAGATAGG",
            "CAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGA",
            "CAGATAGATAGATAGATAGATAGATAGATAGA",
            "CAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGG",
        ),
        # possible genotypes  A: 0/1, 2/3   B: 0/1
        possible_genotypes_a=frozenset({(0, 1), (2, 3)}),
        possible_genotypes_b=frozenset({(0, 1)}),
        # row  A 0/1  B 0/1  possible 2, 1  records 2, 1  -> same genotype
        outcome=Outcome.SAME_GENOTYPE,
        genotype_a=(0, 1),
        genotype_b=(0, 1),
        possible_a=2,
        possible_b=1,
        records_a=2,
        records_b=1,
        invalid_reason_a=None,
        invalid_reason_b=None,
    ),
    HardCase(
        name="H7",
        reference=Reference.GRCH38,
        # excerpt  chr1:14451792-14452092, 301 bases; excerpt position = position - 14451791; sample HG002
        contig="chr1",
        offset=14451791,
        sample="HG002",
        bases=(
            "CCAAATCCAAACCCAAGGACCAGAGAAGGCTAATGTTCCACCTTGAAGACAGTCAGGTAGACTGAGTTCTCCCTTC"
            "TCAGCCTTTTTGTTCTGTTCAGGCACTTGGCAGATTGGACGCTGTTTTGGGATGAGATTAACTCTTTTTTGTTTTG"
            "TTTTGTTTTGTTTTAAGATGGAGTCTCGCTCTGTCACCCAGGCTGGAGTGCACTGGCTCCGTCCCGGCCCACTGCA"
            "GCCTCCGTCTCCTGGATTCTAGCGATTCTCCTGCCTCAGCCTCCCAGGTAGCTGGGATTACAGACACGCCCCA"
        ),
        lines_a=(
            "140 C T 1|0",  # A  TR:311  14451931 C T 1|0
            "146 T G 1|0",  # A  TR:312  14451937 T G 1|0
            "147 G T 1|0",  # A  TR:313  14451938 G T 1|0
            "149 T G 1|0",  # A  TR:314  14451940 T G 1|0
        ),
        lines_b=(
            "139 TCTTTTTTGTTTTGTTTTGTTTTG TTTTTTTGTTGTTGTTTTGTTTTG 0|1",  # B  HipSTR30:2683  14451930 TCTTTTTTGTTTTGTTTTGTTTTG TTTTTTTGTTGTTGTTTTGTTTTG 0|1
        ),
        # locus    chr1:14451931-14451957 (excerpt 140-166), kind STR; STR regions: 14451932-14451957 str-bundle
        locus_1based=(140, 166),
        kind=Kind.STR,
        str_regions_1based=((141, 166, RegionType.STR_BUNDLE, None),),
        alleles=(
            "CTTTTTTGTTTTGTTTTGTTTTGTTTT",
            "TTTTTTGTTGTTGTTTTGTTTTGTTTT",
        ),
        # possible genotypes  A: 0/1   B: 0/1
        possible_genotypes_a=frozenset({(0, 1)}),
        possible_genotypes_b=frozenset({(0, 1)}),
        # row  A 0/1  B 0/1  possible 1, 1  records 4, 1  -> same genotype
        outcome=Outcome.SAME_GENOTYPE,
        genotype_a=(0, 1),
        genotype_b=(0, 1),
        possible_a=1,
        possible_b=1,
        records_a=4,
        records_b=1,
        invalid_reason_a=None,
        invalid_reason_b=None,
    ),
    HardCase(
        name="H8",
        reference=Reference.GRCH38,
        # excerpt  chr1:14451792-14452092, 301 bases; excerpt position = position - 14451791; sample HG002
        contig="chr1",
        offset=14451791,
        sample="HG002",
        bases=(
            "CCAAATCCAAACCCAAGGACCAGAGAAGGCTAATGTTCCACCTTGAAGACAGTCAGGTAGACTGAGTTCTCCCTTC"
            "TCAGCCTTTTTGTTCTGTTCAGGCACTTGGCAGATTGGACGCTGTTTTGGGATGAGATTAACTCTTTTTTGTTTTG"
            "TTTTGTTTTGTTTTAAGATGGAGTCTCGCTCTGTCACCCAGGCTGGAGTGCACTGGCTCCGTCCCGGCCCACTGCA"
            "GCCTCCGTCTCCTGGATTCTAGCGATTCTCCTGCCTCAGCCTCCCAGGTAGCTGGGATTACAGACACGCCCCA"
        ),
        lines_a=(
            "139 TC T 0/1",  # A  GIAB-HG002:19368  14451930 TC T 0/1
            "149 T TG 0/1",  # A  GIAB-HG002:19369  14451940 T TG 0/1
        ),
        lines_b=(
            "139 TCTTTTTTGTTTTGTTTTGTTTTG TTTTTTTGTTGTTGTTTTGTTTTG 0|1",  # B  HipSTR30:2683  14451930 TCTTTTTTGTTTTGTTTTGTTTTG TTTTTTTGTTGTTGTTTTGTTTTG 0|1
        ),
        # locus    chr1:14451930-14451957 (excerpt 139-166), kind STR; STR regions: 14451932-14451957 str-bundle
        locus_1based=(139, 166),
        kind=Kind.STR,
        str_regions_1based=((141, 166, RegionType.STR_BUNDLE, None),),
        alleles=(
            "TCTTTTTTGTTTTGTTTTGTTTTGTTTT",
            "TTTTTTTGTTGTTGTTTTGTTTTGTTTT",
            "TCTTTTTTGTTGTTGTTTTGTTTTGTTTT",
            "TTTTTTTGTTTTGTTTTGTTTTGTTTT",
        ),
        # possible genotypes  A: 0/1, 2/3   B: 0/1
        possible_genotypes_a=frozenset({(0, 1), (2, 3)}),
        possible_genotypes_b=frozenset({(0, 1)}),
        # row  A 0/1  B 0/1  possible 2, 1  records 2, 1  -> same genotype
        outcome=Outcome.SAME_GENOTYPE,
        genotype_a=(0, 1),
        genotype_b=(0, 1),
        possible_a=2,
        possible_b=1,
        records_a=2,
        records_b=1,
        invalid_reason_a=None,
        invalid_reason_b=None,
    ),
    HardCase(
        name="H9",
        reference=Reference.GRCH38,
        # excerpt  chr1:61522820-61523120, 301 bases; excerpt position = position - 61522819; sample HG002
        contig="chr1",
        offset=61522819,
        sample="HG002",
        bases=(
            "CACCACTTCAACCAGAACATTTTCGTACTTATCTGATTTAGAAAATGAATTTCTGTATAAGATTTCATCTGAGGAA"
            "AGGTTTCTGCTGCTATCTATCTACCTATCTGATTATATATTATATATAAATACATATACGTGTGTGTGTATATATA"
            "TATATATGTTTGAGTTCAATAGTTCCCAAGAAATGAACTATTCTCTCCGTAAACTCAAAAATCTGCTCCTAGGTAA"
            "GCACCCAACCTCTGTCCTCCAAAGGTAAGCAAGGGTGGGGGTCTTCTTGGTAGAATGTGGCCAGATGTTTTCA"
        ),
        lines_a=(
            "144 G A 1|0",  # A  TR:910  61522963 G A 1|0
            "156 A ATATATATATATATATATATATATATATATG 1|0",  # A  TR:911  61522975 A ATATATATATATATATATATATATATATATG 1|0
        ),
        lines_b=(
            "143 TGTATATATATATATA TATATATATATATATATATATATATATATATATATATATATATGTA 0|1",  # B  HipSTR30:2961  61522962 TGTATATATATATATA TATATATATATATATATATATATATATATATATATATATATATGTA 0|1
        ),
        # locus    chr1:61522963-61522979 (excerpt 144-160), kind complex; STR regions: none
        locus_1based=(144, 160),
        kind=Kind.COMPLEX,
        str_regions_1based=(),
        alleles=(
            "GTATATATATATATATG",
            "ATATATATATATATATATATATATATATATATATATATATATGTATG",
        ),
        # possible genotypes  A: 0/1   B: 0/1
        possible_genotypes_a=frozenset({(0, 1)}),
        possible_genotypes_b=frozenset({(0, 1)}),
        # row  A 0/1  B 0/1  possible 1, 1  records 2, 1  -> same genotype
        outcome=Outcome.SAME_GENOTYPE,
        genotype_a=(0, 1),
        genotype_b=(0, 1),
        possible_a=1,
        possible_b=1,
        records_a=2,
        records_b=1,
        invalid_reason_a=None,
        invalid_reason_b=None,
    ),
    HardCase(
        name="H10",
        reference=Reference.GRCH38,
        # excerpt  chr1:102475645-102475945, 301 bases; excerpt position = position - 102475644; sample HG002
        contig="chr1",
        offset=102475644,
        sample="HG002",
        bases=(
            "TGATTAATTGCTTTTTAATAAAAGTGCTAAAACAATTTAATGGGGAAATATTTTTTTTTACAAGAAATGGAGCTGA"
            "AACAACTGGATATTCACCTGTCGAATAAATGAATTTGAACCTTGCCTCACATTACACACACACACACACACACACA"
            "CACACACACACACACACACATATACTTCTGTTAAAATGAGACAGACCTAAAAGCAAGCTGAAACTATAAACATATT"
            "TTAGGGAAATGGCATAAGATACATATTCACCATCTTGTGGTAAGGAAAAGACTTTTAGTATATAAAAATCAGT"
        ),
        lines_a=(
            "129 T TAC 0|1",  # A  TR:1332  102475773 T TAC 0|1
            "129 TACAC T 1|0",  # A  TR:1333  102475773 TACAC T 1|0
            "173 T C 0|1",  # A  TR:1334  102475817 T C 0|1
        ),
        lines_b=(
            "130 ACACACACACACACACACACACACACACACACACACACACACAT ACACACACACACACACACACACACACACACACACACACAT,ACACACACACACACACACACACACACACACACACACACACACACAC 1|2",  # B  HipSTR30:3163  102475774 ACACACACACACACACACACACACACACACACACACACACACAT ACACACACACACACACACACACACACACACACACACACAT,ACACACACACACACACACACACACACACACACACACACACACACAC 1|2
        ),
        # locus    chr1:102475773-102475817 (excerpt 129-173), kind STR; STR regions: 102475774-102475815 str AC
        locus_1based=(129, 173),
        kind=Kind.STR,
        str_regions_1based=((130, 171, RegionType.STR, "AC"),),
        alleles=(
            "TACACACACACACACACACACACACACACACACACACACACACAT",
            "TACACACACACACACACACACACACACACACACACACACAT",
            "TACACACACACACACACACACACACACACACACACACACACACACAC",
        ),
        # possible genotypes  A: 1/2   B: 1/2
        possible_genotypes_a=frozenset({(1, 2)}),
        possible_genotypes_b=frozenset({(1, 2)}),
        # row  A 1/2  B 1/2  possible 1, 1  records 3, 1  -> same genotype
        outcome=Outcome.SAME_GENOTYPE,
        genotype_a=(1, 2),
        genotype_b=(1, 2),
        possible_a=1,
        possible_b=1,
        records_a=3,
        records_b=1,
        invalid_reason_a=None,
        invalid_reason_b=None,
    ),
    HardCase(
        name="H11",
        reference=Reference.GRCH38,
        # excerpt  chr1:102475645-102475945, 301 bases; excerpt position = position - 102475644; sample HG002
        contig="chr1",
        offset=102475644,
        sample="HG002",
        bases=(
            "TGATTAATTGCTTTTTAATAAAAGTGCTAAAACAATTTAATGGGGAAATATTTTTTTTTACAAGAAATGGAGCTGA"
            "AACAACTGGATATTCACCTGTCGAATAAATGAATTTGAACCTTGCCTCACATTACACACACACACACACACACACA"
            "CACACACACACACACACACATATACTTCTGTTAAAATGAGACAGACCTAAAAGCAAGCTGAAACTATAAACATATT"
            "TTAGGGAAATGGCATAAGATACATATTCACCATCTTGTGGTAAGGAAAAGACTTTTAGTATATAAAAATCAGT"
        ),
        lines_a=(
            "130 ACACACACACACACACACACACACACACACACACACACACACAT ACACACACACACACACACACACACACACACACACACACAT,ACACACACACACACACACACACACACACACACACACACACACACAC 1|2",  # A  HipSTR30:3163  102475774 ACACACACACACACACACACACACACACACACACACACACACAT ACACACACACACACACACACACACACACACACACACACAT,ACACACACACACACACACACACACACACACACACACACACACACAC 1|2
        ),
        lines_b=(
            "130 ACACACACACACACACACACACACACACACACACACACACAC ACACACACACACACACACACACACACACACACACACAC,ACACACACACACACACACACACACACACACACACACACACACAC 1/2",  # B  pvc-m07-30:3413  102475774 ACACACACACACACACACACACACACACACACACACACACAC ACACACACACACACACACACACACACACACACACACAC,ACACACACACACACACACACACACACACACACACACACACACAC 1/2
            "172 AT AC 0/1",  # B  pvc-m07-30:3414  102475816 AT AC 0/1
        ),
        # locus    chr1:102475773-102475817 (excerpt 129-173), kind STR; STR regions: 102475774-102475815 str AC
        locus_1based=(129, 173),
        kind=Kind.STR,
        str_regions_1based=((130, 171, RegionType.STR, "AC"),),
        alleles=(
            "TACACACACACACACACACACACACACACACACACACACACACAT",
            "TACACACACACACACACACACACACACACACACACACACAT",
            "TACACACACACACACACACACACACACACACACACACACACACACAC",
            "TACACACACACACACACACACACACACACACACACACACAC",
            "TACACACACACACACACACACACACACACACACACACACACACACAT",
        ),
        # possible genotypes  A: 1/2   B: 1/2, 3/4
        possible_genotypes_a=frozenset({(1, 2)}),
        possible_genotypes_b=frozenset({(1, 2), (3, 4)}),
        # row  A 1/2  B 1/2  possible 1, 2  records 1, 2  -> same genotype
        outcome=Outcome.SAME_GENOTYPE,
        genotype_a=(1, 2),
        genotype_b=(1, 2),
        possible_a=1,
        possible_b=2,
        records_a=1,
        records_b=2,
        invalid_reason_a=None,
        invalid_reason_b=None,
    ),
    HardCase(
        name="H12",
        reference=Reference.GRCH38,
        # excerpt  chr1:172778235-172778535, 301 bases; excerpt position = position - 172778234; sample HG002
        contig="chr1",
        offset=172778234,
        sample="HG002",
        bases=(
            "AACAGGGGGAGTCAGTGCCTCCGCCAGTTTTAGCCTCTTGCTTATGCCATTGCAGCATGTGATTAAAAGCTGAATC"
            "CTGTCATTTTTACTGGTTGTTTTAATAGGCTGCGCTGATACCTGGCCACTCAGCTCTCTCTCTCTCTCTCTCTCTC"
            "TCTCTCTCTCTCTCTCTCTCAGCTGAGTTCCTGACAGTCCCAGCTTTGCCTCTTACTGGTTGTGTGACCTTAGGCA"
            "AATCATTTAATGTTTCTCAGCTTCATTTTATCATATGTAAATTTAGGACAATAATACTCGATGAAGTTCTTAT"
        ),
        lines_a=(
            "129 G GCT 0|1",  # A  TR:1784  172778363 G GCT 0|1
            "129 G GCTCTCTCTCTCTCTCTCT 1|0",  # A  TR:1785  172778363 G GCTCTCTCTCTCTCTCTCT 1|0
        ),
        lines_b=(
            "128 AGCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCT AGCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCT,ACCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCT 1|2",  # B  HipSTR30:3357  172778362 AGCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCT AGCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCT,ACCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCT 1|2
        ),
        # locus    chr1:172778363-172778407 (excerpt 129-173), kind STR; STR regions: 172778364-172778405 str CT
        locus_1based=(129, 173),
        kind=Kind.STR,
        str_regions_1based=((130, 171, RegionType.STR, "CT"),),
        alleles=(
            "GCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCA",
            "GCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCA",
            "GCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCA",
            "CCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCA",
        ),
        # possible genotypes  A: 1/2   B: 2/3
        possible_genotypes_a=frozenset({(1, 2)}),
        possible_genotypes_b=frozenset({(2, 3)}),
        # row  A 1/2  B 2/3  possible 1, 1  records 2, 1  -> different
        outcome=Outcome.DIFFERENT,
        genotype_a=(1, 2),
        genotype_b=(2, 3),
        possible_a=1,
        possible_b=1,
        records_a=2,
        records_b=1,
        invalid_reason_a=None,
        invalid_reason_b=None,
    ),
    HardCase(
        name="H13",
        reference=Reference.GRCH38,
        # excerpt  chr10:12391846-12392146, 301 bases; excerpt position = position - 12391845; sample HG002
        contig="chr10",
        offset=12391845,
        sample="HG002",
        bases=(
            "TAAATTTATGGGATAGTTACTGTCATTAACAGCTCATATAACTGAAAATGGTGAAATGATTTGAAAATCAGATAAA"
            "TGGTAATTGCACCACTGCACTCCAGCCTGGATAGCAGAGCGAGACTTGTCTTAAAACACACACACACACACACACA"
            "CACACACACACACACACACAAACAGTCTGGGTATGGTGGCTCACACCTATAATCCTAGCACTTTGGGAGGCCAAGG"
            "TGGTGGATCACCTGTGGTCAGGAGTTCGAGACCAGCCTGACCAACATGGTGAAACTCCGTTTCTACTAAACAT"
        ),
        lines_a=(
            "131 AAC A 1|0",  # A  TR:20498  12391976 AAC A 1|0
        ),
        lines_b=(
            "132 ACACACACACACACACACACACACACACACACACACACAC ACACACACACACACACACACACACACACACACACACAC 1|1",  # B  HipSTR30:3816  12391977 ACACACACACACACACACACACACACACACACACACACAC ACACACACACACACACACACACACACACACACACACAC 1|1
        ),
        # locus    chr10:12391976-12392017 (excerpt 131-172), kind STR; STR regions: 12391977-12392016 str AC
        locus_1based=(131, 172),
        kind=Kind.STR,
        str_regions_1based=((132, 171, RegionType.STR, "AC"),),
        alleles=(
            "AACACACACACACACACACACACACACACACACACACACACA",
            "AACACACACACACACACACACACACACACACACACACACA",
        ),
        # possible genotypes  A: 0/1   B: 1/1
        possible_genotypes_a=frozenset({(0, 1)}),
        possible_genotypes_b=frozenset({(1, 1)}),
        # row  A 0/1  B 1/1  possible 1, 1  records 1, 1  -> same alleles, other genotype
        outcome=Outcome.SAME_ALLELES,
        genotype_a=(0, 1),
        genotype_b=(1, 1),
        possible_a=1,
        possible_b=1,
        records_a=1,
        records_b=1,
        invalid_reason_a=None,
        invalid_reason_b=None,
    ),
    HardCase(
        name="H14",
        reference=Reference.GRCH38,
        # excerpt  chr1:93319126-93319426, 301 bases; excerpt position = position - 93319125; sample HG002
        contig="chr1",
        offset=93319125,
        sample="HG002",
        bases=(
            "ACAGAATACTGCTAAATGTTGGAGTATGCAGAACTCTGGTCCTCGGCTCCTCTCTTTCTATCCAAGTGTAGGGAAT"
            "GCTATGCATTTACAAATTCGTAAACTTTCTTAAAACATTATGAGATATTGTGTGTGTGTGTGTGTGTGTGTGTGTG"
            "TGTGTGTGTGTGTGTGTGTGTGTGTGTTCTATCCCAGAGTAGGGAATGCTATACATTTATAAATTCGTAAACGTTC"
            "TTAAAACATTATGAGATATTTTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTAGCTCATCAGCTATTGT"
        ),
        lines_a=(
            "123 AT A 1|0",  # A  TR:1161  93319248 AT A 1|0
        ),
        lines_b=(
            "125 TGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTG TGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTG 1|1",  # B  HipSTR30:3117  93319250 TGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTG TGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTG 1|1
        ),
        # locus    chr1:93319248-93319304 (excerpt 123-179), kind STR; STR regions: 93319250-93319303 str TG
        locus_1based=(123, 179),
        kind=Kind.STR,
        str_regions_1based=((125, 178, RegionType.STR, "TG"),),
        alleles=(
            "ATTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGT",
            "ATGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGT",
            "ATTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGT",
        ),
        # possible genotypes  A: 0/1   B: 2/2
        possible_genotypes_a=frozenset({(0, 1)}),
        possible_genotypes_b=frozenset({(2, 2)}),
        # row  A 0/1  B 2/2  possible 1, 1  records 1, 1  -> different
        outcome=Outcome.DIFFERENT,
        genotype_a=(0, 1),
        genotype_b=(2, 2),
        possible_a=1,
        possible_b=1,
        records_a=1,
        records_b=1,
        invalid_reason_a=None,
        invalid_reason_b=None,
    ),
    HardCase(
        name="H15",
        reference=Reference.GRCH38,
        # excerpt  chr1:93319126-93319426, 301 bases; excerpt position = position - 93319125; sample HG002
        contig="chr1",
        offset=93319125,
        sample="HG002",
        bases=(
            "ACAGAATACTGCTAAATGTTGGAGTATGCAGAACTCTGGTCCTCGGCTCCTCTCTTTCTATCCAAGTGTAGGGAAT"
            "GCTATGCATTTACAAATTCGTAAACTTTCTTAAAACATTATGAGATATTGTGTGTGTGTGTGTGTGTGTGTGTGTG"
            "TGTGTGTGTGTGTGTGTGTGTGTGTGTTCTATCCCAGAGTAGGGAATGCTATACATTTATAAATTCGTAAACGTTC"
            "TTAAAACATTATGAGATATTTTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTAGCTCATCAGCTATTGT"
        ),
        lines_a=(
            "124 TTGTG T 0/1",  # A  GIAB-HG002:128803  93319249 TTGTG T 0/1
        ),
        lines_b=(
            "125 TGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTG TGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTG 1|1",  # B  HipSTR30:3117  93319250 TGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTG TGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTG 1|1
        ),
        # locus    chr1:93319249-93319304 (excerpt 124-179), kind STR; STR regions: 93319250-93319303 str TG
        locus_1based=(124, 179),
        kind=Kind.STR,
        str_regions_1based=((125, 178, RegionType.STR, "TG"),),
        alleles=(
            "TTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGT",
            "TTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGT",
        ),
        # possible genotypes  A: 0/1   B: 1/1
        possible_genotypes_a=frozenset({(0, 1)}),
        possible_genotypes_b=frozenset({(1, 1)}),
        # row  A 0/1  B 1/1  possible 1, 1  records 1, 1  -> same alleles, other genotype
        outcome=Outcome.SAME_ALLELES,
        genotype_a=(0, 1),
        genotype_b=(1, 1),
        possible_a=1,
        possible_b=1,
        records_a=1,
        records_b=1,
        invalid_reason_a=None,
        invalid_reason_b=None,
    ),
    HardCase(
        name="O1",
        reference=Reference.GRCH38,
        # excerpt  chr6:41399460-41399760, 301 bases; excerpt position = position - 41399459; sample HG002
        contig="chr6",
        offset=41399459,
        sample="HG002",
        bases=(
            "ATTCCCCTCCCCTTGTGTGTGGGTTGGACTTAGTGACTTGCCTTTAACAGATAAATTATGTTGGAAGTCACAAGTA"
            "TGTGACTTGTAAGACCAAGTCATAAAAGGCATTGTGGCTGCCTCCACCGCCTCCTCCTCCTCCTCCTCCTCCTCCT"
            "CCTCCTCTTCTTCTTCCTCCTCCTCCTCCTCTTCCTCCTCTCTCTTTCTCTCTTTCTCACTCACTCGCTCGCTCAC"
            "AGAAACCAGCTGCCACGAGGTGAGGACTTAAAGTAACCCTGTGGAAAGGCCCACATGGGAGAAACGGAGGCCT"
        ),
        lines_a=(
            "125 GCCT G 0|1",  # A  TR:12602  41399584 GCCT G 0|1
            "154 CTCCTCT C 1|0",  # A  TR:12603  41399613 CTCCTCT C 1|0
            "166 TTCC T 0|1",  # A  TR:12604  41399625 TTCC T 0|1
            "169 C T 1|0",  # A  TR:12605  41399628 C T 1|0
        ),
        lines_b=(
            "126 CCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCTTCTTCTTCCTCCTCCTCCTCCT CCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCTTCTTCTTCCTCCTCCTCCT 1/1",  # B  pvc-touch30:9746  41399585 CCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCTTCTTCTTCCTCCTCCTCCTCCT CCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCTTCTTCTTCCTCCTCCTCCT 1/1
        ),
        # locus    chr6:41399584-41399642 (excerpt 125-183), kind STR; STR regions: 41399585-41399641 str CCT
        locus_1based=(125, 183),
        kind=Kind.STR,
        str_regions_1based=((126, 182, RegionType.STR, "CCT"),),
        alleles=(
            "GCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCTTCTTCTTCCTCCTCCTCCTCCTC",
            "GCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCTTCTTCTTCCTCCTCCTCCTC",
        ),
        # possible genotypes  A: 1/1   B: 1/1
        possible_genotypes_a=frozenset({(1, 1)}),
        possible_genotypes_b=frozenset({(1, 1)}),
        # row  A 1/1  B 1/1  possible 1, 1  records 4, 1  -> same genotype
        outcome=Outcome.SAME_GENOTYPE,
        genotype_a=(1, 1),
        genotype_b=(1, 1),
        possible_a=1,
        possible_b=1,
        records_a=4,
        records_b=1,
        invalid_reason_a=None,
        invalid_reason_b=None,
    ),
    HardCase(
        name="O2",
        reference=Reference.GRCH38,
        # excerpt  chr6:41399460-41399760, 301 bases; excerpt position = position - 41399459; sample HG002
        contig="chr6",
        offset=41399459,
        sample="HG002",
        bases=(
            "ATTCCCCTCCCCTTGTGTGTGGGTTGGACTTAGTGACTTGCCTTTAACAGATAAATTATGTTGGAAGTCACAAGTA"
            "TGTGACTTGTAAGACCAAGTCATAAAAGGCATTGTGGCTGCCTCCACCGCCTCCTCCTCCTCCTCCTCCTCCTCCT"
            "CCTCCTCTTCTTCTTCCTCCTCCTCCTCCTCTTCCTCCTCTCTCTTTCTCTCTTTCTCACTCACTCGCTCGCTCAC"
            "AGAAACCAGCTGCCACGAGGTGAGGACTTAAAGTAACCCTGTGGAAAGGCCCACATGGGAGAAACGGAGGCCT"
        ),
        lines_a=(
            "125 GCCT G 0|1",  # A  TR:12602  41399584 GCCT G 0|1
            "154 CTCCTCT C 1|0",  # A  TR:12603  41399613 CTCCTCT C 1|0
            "166 TTCC T 0|1",  # A  TR:12604  41399625 TTCC T 0|1
            "169 C T 1|0",  # A  TR:12605  41399628 C T 1|0
        ),
        lines_b=(
            "154 CTCCTCTT CT,CTCCTCCT 1/2",  # B  freebayes30:31443  41399613 CTCCTCTT CT,CTCCTCCT 1/2
            "166 T C 0/1",  # B  freebayes30:31444  41399625 T C 0/1
            "169 C T 0/1",  # B  freebayes30:31445  41399628 C T 0/1
        ),
        # locus    chr6:41399584-41399642 (excerpt 125-183), kind STR; STR regions: 41399585-41399641 str CCT
        locus_1based=(125, 183),
        kind=Kind.STR,
        str_regions_1based=((126, 182, RegionType.STR, "CCT"),),
        alleles=(
            "GCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCTTCTTCTTCCTCCTCCTCCTCCTC",
            "GCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCTTCTTCTTCCTCCTCCTCCTC",
            "GCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCTTCTTCCTCCTCCTCCTCCTC",
            "GCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCTTCCTCTTCCTCCTCCTCCTC",
            "GCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCTTCCTCCTCCTCCTCCTCCTC",
            "GCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCTTCTTCTTCCTCCTCCTCCTC",
            "GCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCTTCCTCCTCCTCCTCCTCCTC",
            "GCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCTTCCTCTTCCTCCTCCTCCTC",
            "GCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCCTCTTCTTCCTCCTCCTCCTCCTC",
        ),
        # possible genotypes  A: 1/1   B: 1/6, 2/3, 4/5, 7/8
        possible_genotypes_a=frozenset({(1, 1)}),
        possible_genotypes_b=frozenset({(1, 6), (2, 3), (4, 5), (7, 8)}),
        # row  A 1/1  B 1/6  possible 1, 4  records 4, 3  -> different
        outcome=Outcome.DIFFERENT,
        genotype_a=(1, 1),
        genotype_b=(1, 6),
        possible_a=1,
        possible_b=4,
        records_a=4,
        records_b=3,
        invalid_reason_a=None,
        invalid_reason_b=None,
    ),
    HardCase(
        name="O3",
        reference=Reference.GRCH38,
        # excerpt  chr3:37126720-37127020, 301 bases; excerpt position = position - 37126719; sample HG002
        contig="chr3",
        offset=37126719,
        sample="HG002",
        bases=(
            "ATGGTGGCGGGCGCCTGTAATCCCAGCTACTTGAGATGCTGAGGCAGGAGAATCGCTTGAACCCAGGAGACGGAGG"
            "TTGCAGTGAGCTGAGATTGTGCCACTGAACTCCAGCCTGGGCAACAAGCGCGAAACTCCATCTCAAAAAAAAAAAA"
            "GAAAGAAAGAAAGAAAGAAAGTTGCTTGTTTAGCAGCAATCAGACAGTGGTTCTTAACCCTTTAGGAATTACAGAC"
            "TACTTTGAGAATCTGACAATAAAAATACCTGTGGGCATTTCCCATAAAAATAATTACATTTGAAAGTCTGCAG"
        ),
        lines_a=(
            "149 A AAAAAG 1|1",  # A  TR:5893  37126868 A AAAAAG 1|1
        ),
        lines_b=(
            "141 AAAAAAAAAAAA AAAAAAAAAAAAAGAAA 1/1",  # B  pvc-touch30:5857  37126860 AAAAAAAAAAAA AAAAAAAAAAAAAGAAA 1/1
        ),
        # locus    chr3:37126860-37126872 (excerpt 141-153), kind STR; STR regions: 37126860-37126871 str A
        locus_1based=(141, 153),
        kind=Kind.STR,
        str_regions_1based=((141, 152, RegionType.STR, "A"),),
        alleles=(
            "AAAAAAAAAAAAG",
            "AAAAAAAAAAAAAGAAAG",
        ),
        # possible genotypes  A: 1/1   B: 1/1
        possible_genotypes_a=frozenset({(1, 1)}),
        possible_genotypes_b=frozenset({(1, 1)}),
        # row  A 1/1  B 1/1  possible 1, 1  records 1, 1  -> same genotype
        outcome=Outcome.SAME_GENOTYPE,
        genotype_a=(1, 1),
        genotype_b=(1, 1),
        possible_a=1,
        possible_b=1,
        records_a=1,
        records_b=1,
        invalid_reason_a=None,
        invalid_reason_b=None,
    ),
    HardCase(
        name="O4",
        reference=Reference.GRCH38,
        # excerpt  chr3:37126720-37127020, 301 bases; excerpt position = position - 37126719; sample HG002
        contig="chr3",
        offset=37126719,
        sample="HG002",
        bases=(
            "ATGGTGGCGGGCGCCTGTAATCCCAGCTACTTGAGATGCTGAGGCAGGAGAATCGCTTGAACCCAGGAGACGGAGG"
            "TTGCAGTGAGCTGAGATTGTGCCACTGAACTCCAGCCTGGGCAACAAGCGCGAAACTCCATCTCAAAAAAAAAAAA"
            "GAAAGAAAGAAAGAAAGAAAGTTGCTTGTTTAGCAGCAATCAGACAGTGGTTCTTAACCCTTTAGGAATTACAGAC"
            "TACTTTGAGAATCTGACAATAAAAATACCTGTGGGCATTTCCCATAAAAATAATTACATTTGAAAGTCTGCAG"
        ),
        lines_a=(
            "149 A AAAAAG 1|1",  # A  TR:5893  37126868 A AAAAAG 1|1
        ),
        lines_b=(
            "150 AAAGAAAGAAAGAAAGAAAGAAAG AAAAGAAAGAAAGAAAGAAAGAAAGAAAG 1|1",  # B  HipSTR30:10361  37126869 AAAGAAAGAAAGAAAGAAAGAAAG AAAAGAAAGAAAGAAAGAAAGAAAGAAAG 1|1
        ),
        # locus    chr3:37126860-37126872 (excerpt 141-153), kind STR; STR regions: 37126860-37126871 str A
        locus_1based=(141, 153),
        kind=Kind.STR,
        str_regions_1based=((141, 152, RegionType.STR, "A"),),
        alleles=(
            "AAAAAAAAAAAAG",
            "AAAAAAAAAAAAAGAAAG",
        ),
        # possible genotypes  A: 1/1   B: 1/1
        possible_genotypes_a=frozenset({(1, 1)}),
        possible_genotypes_b=frozenset({(1, 1)}),
        # row  A 1/1  B 1/1  possible 1, 1  records 1, 1  -> same genotype
        outcome=Outcome.SAME_GENOTYPE,
        genotype_a=(1, 1),
        genotype_b=(1, 1),
        possible_a=1,
        possible_b=1,
        records_a=1,
        records_b=1,
        invalid_reason_a=None,
        invalid_reason_b=None,
    ),
    HardCase(
        name="O5",
        reference=Reference.GRCH38,
        # excerpt  chr3:37126720-37127020, 301 bases; excerpt position = position - 37126719; sample HG002
        contig="chr3",
        offset=37126719,
        sample="HG002",
        bases=(
            "ATGGTGGCGGGCGCCTGTAATCCCAGCTACTTGAGATGCTGAGGCAGGAGAATCGCTTGAACCCAGGAGACGGAGG"
            "TTGCAGTGAGCTGAGATTGTGCCACTGAACTCCAGCCTGGGCAACAAGCGCGAAACTCCATCTCAAAAAAAAAAAA"
            "GAAAGAAAGAAAGAAAGAAAGTTGCTTGTTTAGCAGCAATCAGACAGTGGTTCTTAACCCTTTAGGAATTACAGAC"
            "TACTTTGAGAATCTGACAATAAAAATACCTGTGGGCATTTCCCATAAAAATAATTACATTTGAAAGTCTGCAG"
        ),
        lines_a=(
            "149 A AAAAAG 1|1",  # A  TR:5893  37126868 A AAAAAG 1|1
        ),
        lines_b=(
            "140 CAAAAAAAAAAAAGAAAGAAAGAAAGAAAGAAAGT CAAAAAAAAAAAAAGAAAGAAAGAAAGAAAGAAAGAAAGT 1/1",  # B  freebayes30:16445  37126859 CAAAAAAAAAAAAGAAAGAAAGAAAGAAAGAAAGT CAAAAAAAAAAAAAGAAAGAAAGAAAGAAAGAAAGAAAGT 1/1
        ),
        # locus    chr3:37126860-37126872 (excerpt 141-153), kind STR; STR regions: 37126860-37126871 str A
        locus_1based=(141, 153),
        kind=Kind.STR,
        str_regions_1based=((141, 152, RegionType.STR, "A"),),
        alleles=(
            "AAAAAAAAAAAAG",
            "AAAAAAAAAAAAAGAAAG",
        ),
        # possible genotypes  A: 1/1   B: 1/1
        possible_genotypes_a=frozenset({(1, 1)}),
        possible_genotypes_b=frozenset({(1, 1)}),
        # row  A 1/1  B 1/1  possible 1, 1  records 1, 1  -> same genotype
        outcome=Outcome.SAME_GENOTYPE,
        genotype_a=(1, 1),
        genotype_b=(1, 1),
        possible_a=1,
        possible_b=1,
        records_a=1,
        records_b=1,
        invalid_reason_a=None,
        invalid_reason_b=None,
    ),
)
