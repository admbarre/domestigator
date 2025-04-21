import pytest
from domestigator import REnzyme

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

@pytest.mark.parametrize("enz,seq,expected",tests)
def test_find_cuts(enz,seq,expected):
    assert enz.find_cuts(seq) == expected
