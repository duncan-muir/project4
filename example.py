from align import *


nw = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", -10, -1)
seq1, _ = read_fasta("./data/test_seq1.fa")
seq2, _ = read_fasta("./data/test_seq2.fa")

print(nw.align(seq1, seq2))
align_mat = nw.get_align_matrix()
