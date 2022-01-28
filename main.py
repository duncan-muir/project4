# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch

def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    """
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    gg_seq, gg_header = read_fasta("./data/Gallus_gallus_BRD2.fa")
    mm_seq, mm_header = read_fasta("./data/Mus_musculus_BRD2.fa")
    br_seq, br_header = read_fasta("./data/Balaeniceps_rex_BRD2.fa")
    tt_seq, tt_header = read_fasta("./data/tursiops_truncatus_BRD2.fa")

    # TODO Align all species to humans and print species in order of most similar to human BRD
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    pass

    # TODO print all of the alignment score between each species BRD2 and human BRD2
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    nw = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", -10, -1)

    gg_score, gg_seq_a, gg_seq_b = nw.align(hs_seq, gg_seq)
    print(f"Chicken (Gallus): {gg_score}")
    mm_score, mm_seq_a, mm_seq_b = nw.align(hs_seq, mm_seq)
    print(f"Mouse (Mus): {mm_score}")
    br_score, br_seq_a, br_seq_b = nw.align(hs_seq, br_seq)
    print(f"Stork (Balaeniceps): {br_score}")
    tt_score, tt_seq_a, tt_seq_b = nw.align(hs_seq, tt_seq)
    print(f"Dolphin (Tursiops): {tt_score}")

if __name__ == "__main__":
    main()
