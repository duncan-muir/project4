# Importing Dependencies
import pytest
from align import *
import numpy as np

def test_nw_alignment():
    """
    TODO: Write your unit test for NW alignment
    using test_seq1.fa and test_seq2.fa by
    asserting that you have correctly filled out
    the your 3 alignment matrices.
    Use the BLOSUM62 matrix and a gap open penalty
    of -10 and a gap extension penalty of -1.
    """
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")

    nw = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", -10, -1)

    score, _, _ = nw.align(seq1, seq2)
    assert score == 4

    align_mat = nw.get_align_matrix()

    assert np.array_equal(align_mat[GAP_A_IDX],
                          np.array([[-10, -11, -12, -13],
                                    [float("-inf"), -22, -6, -7],
                                    [float("-inf"), -23, -17, -7],
                                    [float("-inf"), -24, -18, -12],
                                    [float("-inf"), -25, -19, -17]]))



def test_nw_backtrace():
    """
    TODO: Write your unit test for NW backtracing
    using test_seq3.fa and test_seq4.fa by
    asserting that the backtrace is correct.
    Use the BLOSUM62 matrix. Use a gap open
    penalty of -10 and a gap extension penalty of -1.
    """
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")
    pass




