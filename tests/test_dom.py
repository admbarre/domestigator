import pytest
from domestigator import REnzyme
from domestigator import WildSite

# Cut sites
bsai = REnzyme("BsaI","GGTCTC")
bsmbi = REnzyme("BsmbI","CGTCTC")
paqci = REnzyme("PaqcI","CACCTGC")

# Courtesy of Mr. GPT
# NOTE: MR GPT SOMETIMES LIES MANUALLY DOUBLE CHECK
tests = [
    # BsaI (GGTCTC / GAGACC)
    (bsai, "ATCGATCGATCG", []),  # No site
    (bsai, "AAAAGGTCTCTTTT", [("GGTCTC", 4)]),  # One forward
    (bsai, "TTTGAGACCGGG", [("GAGACC", 3)]),  # One reverse
    (bsai, "GGTCTCaaaGGTCTC", [("GGTCTC", 0), ("GGTCTC", 9)]),  # Multiple forward
    (bsai, "GAGACCaaaGAGACC", [("GAGACC", 0), ("GAGACC", 9)]),  # Multiple reverse

    # BsmBI (CGTCTC / GAGACG)
    (bsmbi, "ATCGATCGATCG", []),  # No site
    (bsmbi, "TTTCGTCTCGGG", [("CGTCTC", 3)]),  # One forward
    (bsmbi, "TTTGAGACGGGG", [("GAGACG", 3)]),  # One reverse
    (bsmbi, "CGTCTCxxxCGTCTC", [("CGTCTC", 0), ("CGTCTC", 9)]),  # Multiple forward
    (bsmbi, "GAGACGyyyGAGACG", [("GAGACG", 0), ("GAGACG", 9)]),  # Multiple reverse

    # PaqCI (CACCTGC / GCAGGTG)
    (paqci, "ATCGATCGATCG", []),  # No site
    (paqci, "TTTCACCTGCGGG", [("CACCTGC", 3)]),  # One forward
    (paqci, "TTTGCAGGTGGGG", [("GCAGGTG", 3)]),  # One reverse
    (paqci, "CACCTGCaaaCACCTGC", [("CACCTGC", 0), ("CACCTGC", 10)]),  # Multiple forward
    (paqci, "GCAGGTGcccGCAGGTG", [("GCAGGTG", 0), ("GCAGGTG", 10)]),  # Multiple reverse
]

ws_tests = [
    # BsaI (forward: GGTCTC, reverse: GAGACC)
    ("GGTCTCaaa", "GGTCTC", 0, "GGTCTC"),        # BsaI forward, in frame
    ("aGGTCTCaaa", "GGTCTC", 1, "aGGTCTCaa"),     # BsaI forward, +1
    ("aaGGTCTCaaa", "GGTCTC", 2, "aaGGTCTCa"),    # BsaI forward, +2
    ("GAGACCaaa", "GAGACC", 0, "GAGACC"),        # BsaI reverse, in frame
    ("aGAGACCaaa", "GAGACC", 1, "aGAGACCaa"),     # BsaI reverse, +1
    ("aaGAGACCaaa", "GAGACC", 2, "aaGAGACCa"),    # BsaI reverse, +2

    # BsmBI (forward: CGTCTC, reverse: GAGACG)
    ("CGTCTCaaa", "CGTCTC", 0, "CGTCTC"),        # BsmBI forward, in frame
    ("aCGTCTCaaa", "CGTCTC", 1, "aCGTCTCaa"),     # BsmBI forward, +1
    ("aaCGTCTCaaa", "CGTCTC", 2, "aaCGTCTCa"),    # BsmBI forward, +2
    ("GAGACGaaa", "GAGACG", 0, "GAGACG"),        # BsmBI reverse, in frame
    ("aGAGACGaaa", "GAGACG", 1, "aGAGACGaa"),     # BsmBI reverse, +1
    ("aaGAGACGaaa", "GAGACG", 2, "aaGAGACGa"),    # BsmBI reverse, +2

    # PaqCI (forward: CACCTGC, reverse: GCAGGTG)
    ("CACCTGCaaa", "CACCTGC", 0, "CACCTGCaa"),        # PaqCI forward, in frame
    ("aCACCTGCaaa", "CACCTGC", 1, "aCACCTGCa"),      # PaqCI forward, +1
    ("aaCACCTGCaaa", "CACCTGC", 2, "aaCACCTGC"),     # PaqCI forward, +2
    ("GCAGGTGaaa", "GCAGGTG", 0, "GCAGGTGaa"),        # PaqCI reverse, in frame
    ("aGCAGGTGaaa", "GCAGGTG", 1, "aGCAGGTGa"),      # PaqCI reverse, +1
    ("aaGCAGGTGaaa", "GCAGGTG", 2, "aaGCAGGTG"),     # PaqCI reverse, +2
]
check_tests = [
    ("x12345GGTCTC12345yyy", "GGTCTC", 6, ("12345","12345")), # +0 L/R full length
    ("xx12345GGTCTC12345yyy", "GGTCTC", 7, ("x1234","345yy")), # +1 L/R full length
    ("xxx12345GGTCTC12345yyy", "GGTCTC", 8, ("xx123","2345y")), # +2 L/R full length

    ("345GGTCTC12345yyy", "GGTCTC", 3, ("345","12345")), # +0 R full length
    ("2345GGTCTC12345yyy", "GGTCTC", 4, ("234","345yy")), # +1 R full length
    ("12345GGTCTC12345yyy", "GGTCTC", 5, ("123","2345y")), # +2 R full length

    ("x12345GGTCTC123", "GGTCTC", 6, ("12345","123")), # +0 L full length
    ("xx12345GGTCTC1234", "GGTCTC", 7, ("x1234","34")), # +1 L full length
    ("xxx12345GGTCTC123", "GGTCTC", 8, ("xx123","23")), # +2 L full length

    ("345GGTCTC123", "GGTCTC", 3, ("345","123")), # +0
    ("2345GGTCTC1234", "GGTCTC", 4, ("234","34")), # +1
    ("12345GGTCTC123", "GGTCTC", 5, ("123","23")), # +2
]
@pytest.mark.parametrize("enz,seq,expected",tests)
def test_find_cuts(enz,seq,expected):
    assert enz.find_cuts(seq) == expected

@pytest.mark.parametrize("seq,cut,index,expected",ws_tests)
def test_wild_site(seq,cut,index,expected):
    ws = WildSite(seq,cut,index)
    assert ws.orf_seq == expected

@pytest.mark.parametrize("seq,cut,index,expected",check_tests)
def test_checks(seq,cut,index,expected):
    left,right = expected
    ws = WildSite(seq,cut,index)
    assert ws.left == left
    assert ws.right == right

        
