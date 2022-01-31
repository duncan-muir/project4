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

    # get species name (gross string parsing)
    species = [item.split("=")[1].rstrip(" OX") for item in [gg_header, mm_header, br_header, tt_header]]

    nw = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", -10, -1)

    gg_alignment = nw.align(hs_seq, gg_seq)
    mm_alignment = nw.align(hs_seq, mm_seq)
    br_alignment = nw.align(hs_seq, br_seq)
    tt_alignment = nw.align(hs_seq, tt_seq)

    # TODO Align all species to humans and print species in order of most similar to human BRD
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    pass

    species_alignment_zip = list(zip([gg_alignment, mm_alignment, br_alignment, tt_alignment], species))
    print("Species in order of most to least similar to human:")
    for _, spec in sorted(species_alignment_zip, reverse=True):
        print("\t", spec)

    print("Alignment scores:")
    for alignment, spec in species_alignment_zip:
        score, _, _ = alignment
        print(f"\t{spec} alignment score: {score}")


if __name__ == "__main__":
    main()
