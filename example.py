from align import NeedlemanWunsch, read_fasta


nw = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", -10, -1)
seq1, _ = read_fasta("./data/test_seq3.fa")
seq2, _ = read_fasta("./data/test_seq4.fa")

print(nw.align(seq1, seq2))
