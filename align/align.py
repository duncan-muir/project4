# Importing Dependencies
import numpy as np
from typing import Tuple

GAP_A_IDX = 0
ALIGN_IDX = 1
GAP_B_IDX = 2


# Defining class for Needleman-Wunsch Algorithm for Global pairwise alignment
class NeedlemanWunsch:
    """ Class for NeedlemanWunsch Alignment

    Parameters:
        sub_matrix_file: str
            Path/filename of substitution matrix
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty

    Attributes:
        seqA_align: str
            seqA alignment
        seqB_align: str
            seqB alignment
        alignment_score: float
            Score of alignment from algorithm
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty
    """
    def __init__(self, sub_matrix_file: str, gap_open: float, gap_extend: float):
        # Init alignment and gap matrices
        self._align_matrix = None
        self._gapA_matrix = None
        self._gapB_matrix = None

        # Init matrices for backtrace procedure
        self._back = None
        self._back_A = None
        self._back_B = None

        # Init alignment_score
        self.alignment_score = 0

        # Init empty alignment attributes
        self.seqA_align = ""
        self.seqB_align = ""

        # Init empty sequences
        self._seqA = ""
        self._seqB = ""

        # Setting gap open and gap extension penalties
        self.gap_open = gap_open
        assert gap_open < 0, "Gap opening penalty must be negative."
        self.gap_extend = gap_extend
        assert gap_extend < 0, "Gap extension penalty must be negative."

        # Generating substitution matrix
        self.sub_dict = self._read_sub_matrix(sub_matrix_file) # substitution dictionary

    def _read_sub_matrix(self, sub_matrix_file):
        """
        DO NOT MODIFY THIS METHOD! IT IS ALREADY COMPLETE!

        This function reads in a scoring matrix from any matrix like file.
        Where there is a line of the residues followed by substitution matrix.
        This file also saves the alphabet list attribute.

        Parameters:
            sub_matrix_file: str
                Name (and associated path if not in current working directory)
                of the matrix file that contains the scoring matrix.

        Returns:
            dict_sub: dict
                Substitution matrix dictionary with tuple of the two residues as
                the key and score as value e.g. {('A', 'A'): 4} or {('A', 'D'): -8}
        """
        with open(sub_matrix_file, 'r') as f:
            dict_sub = {}  # Dictionary for storing scores from sub matrix
            residue_list = []  # For storing residue list
            start = False  # trigger for reading in score values
            res_2 = 0  # used for generating substitution matrix
            # reading file line by line
            for line_num, line in enumerate(f):
                # Reading in residue list
                if '#' not in line.strip() and start is False:
                    residue_list = [k for k in line.strip().upper().split(' ') if k != '']
                    start = True
                # Generating substitution scoring dictionary
                elif start is True and res_2 < len(residue_list):
                    line = [k for k in line.strip().split(' ') if k != '']
                    # reading in line by line to create substitution dictionary
                    assert len(residue_list) == len(line), "Score line should be same length as residue list"
                    for res_1 in range(len(line)):
                        dict_sub[(residue_list[res_1], residue_list[res_2])] = float(line[res_1])
                    res_2 += 1
                elif start is True and res_2 == len(residue_list):
                    break
        return dict_sub

    def align(self, seqA: str, seqB: str) -> Tuple[float, str, str]:
        """
        # TODO: Fill in the Needleman-Wunsch Algorithm below
        to perform global sequence alignment of seqA and seqB
        and return a tuple with the following format
        (alignment score, seqA alignment, seqB alignment)
        Also, write up a docstring for this function using the
        _read_sub_matrix as an example.
        Don't forget to comment your code!
        """
        # Initialize 6 matrix private attributes for use in alignment
        # create matrices for alignment scores and gaps
        self._align_matrix = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf
        self._gapA_matrix = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf
        self._gapB_matrix = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf

        self._align_3d = np.ones((3, len(seqA) + 1, len(seqB) + 1)) * -np.inf
        # create matrices for pointers used in backtrace procedure
        self._back = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf
        self._back_A = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf
        self._back_B = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf

        # TODO Nan???
        self._back_3d = np.zeros((3, len(seqA) + 1, len(seqB) + 1, 3), dtype=int)
        # Resetting alignment in case method is called more than once
        self.seqA_align = ""
        self.seqB_align = ""

        # Resetting alignment score in case method is called more than once
        self.alignment_score = 0

        # Initializing sequences for use in backtrace method
        self._seqA = seqA
        self._seqB = seqB

        self._align_matrix[0, 0] = 0

        self._gapA_matrix[0] = (self.gap_open + self.gap_extend
                                   * np.array(range(len(self._seqB) + 1)))

        self._gapB_matrix[:, 0] = (self.gap_open + self.gap_extend
                                * np.array(range(len(self._seqA) + 1)))

        self._align_3d[ALIGN_IDX, 0, 0] = 0
        self._align_3d[GAP_A_IDX, 0] = (self.gap_open + self.gap_extend
                                   * np.array(range(len(self._seqB) + 1)))
        self._align_3d[GAP_B_IDX, :, 0] = (self.gap_open + self.gap_extend
                                * np.array(range(len(self._seqA) + 1)))

        # SET SOME INITIAL FRICKIN POINTERS!!!!!!
        self._back_3d[GAP_A_IDX][:,0][:,1] = list(range(len(self._seqA) + 1))
        self._back_3d[GAP_B_IDX][0][:,2] = list(range(len(self._seqB) + 1))
        # print(self._back_3d[GAP_B_IDX][0])

        # Going to have following key:
        # 0 : gapA,
        # 1 : aligned,
        # 2 : gapB
        for i in range(1, len(self._seqA) + 1):
            for j in range(1, len(self._seqB) + 1):

                # trying to implement highroad????????
                curr_idx = (i, j)

                # start with match matrix
                from_i = i - 1
                from_j = j - 1
                score = self.sub_dict[(self._seqA[i - 1], self._seqB[j - 1])]
                aligned = self._align_matrix[from_i, from_j] + score
                gapped_a = self._gapA_matrix[from_i, from_j] + score
                gapped_b = self._gapB_matrix[from_i, from_j] + score

                if (gapped_a >= aligned) and (gapped_a >= gapped_b):
                    self._align_matrix[curr_idx] = gapped_a
                    self._back_3d[ALIGN_IDX, i, j] = (GAP_A_IDX, from_i, from_j)

                elif aligned >= gapped_b:
                    self._align_matrix[curr_idx] = aligned
                    self._back_3d[ALIGN_IDX, i, j] = (ALIGN_IDX, from_i, from_j)

                else:
                    self._align_matrix[curr_idx] = gapped_b
                    self._back_3d[ALIGN_IDX, i, j] = (GAP_B_IDX, from_i, from_j)

                # next do gapA matrix
                from_i = i
                from_j = j - 1
                aligned = self._align_matrix[from_i, from_j] + self.gap_open + self.gap_extend
                gapped_a = self._gapA_matrix[from_i, from_j] + self.gap_extend
                gapped_b = self._gapB_matrix[from_i, from_j] + self.gap_open + self.gap_extend

                if (gapped_a >= aligned) and (gapped_a >= gapped_b):
                    self._gapA_matrix[curr_idx] = gapped_a
                    self._back_3d[GAP_A_IDX, i, j] = (GAP_A_IDX, from_i, from_j)

                elif aligned >= gapped_b:
                    self._gapA_matrix[curr_idx] = aligned
                    self._back_3d[GAP_A_IDX, i, j] = (ALIGN_IDX, from_i, from_j)

                else:
                    self._gapA_matrix[curr_idx] = gapped_b
                    self._back_3d[GAP_A_IDX, i, j] = (GAP_B_IDX, from_i, from_j)

                # next do gapB matrix
                from_i = i - 1
                from_j = j

                aligned = self._align_matrix[from_i, from_j] + self.gap_open + self.gap_extend
                gapped_a = self._gapA_matrix[from_i, from_j] + self.gap_open + self.gap_extend
                gapped_b = self._gapB_matrix[from_i, from_j] + self.gap_extend

                if (gapped_a >= aligned) and (gapped_a >= gapped_b):
                    self._gapB_matrix[curr_idx] = gapped_a
                    self._back_3d[GAP_B_IDX, i, j] = (GAP_A_IDX, from_i, from_j)

                elif aligned >= gapped_b:
                    self._gapB_matrix[curr_idx] = aligned
                    self._back_3d[GAP_B_IDX, i, j] = (ALIGN_IDX, from_i, from_j)

                else:
                    self._gapB_matrix[curr_idx] = gapped_b
                    self._back_3d[GAP_B_IDX, i, j] = (GAP_B_IDX, from_i, from_j)

        return self._backtrace()

    def _backtrace(self) -> Tuple[float, str, str]:
        """
        # TODO Implement the traceback procedure method below
        based on the heuristic you implement in the align method.
        The traceback method should return a tuple of the alignment
        score, the seqA alignment and the seqB alignment respectively.
        """
        # Implement this method based upon the heuristic chosen in the align method above.
        i, j = len(self._seqA), len(self._seqB)

        opt_strat = None

        align_last = self._align_matrix[i, j]
        gap_a_last = self._gapA_matrix[i, j]
        gap_b_last = self._gapB_matrix[i, j]

        if (gap_a_last >= align_last) and (gap_a_last >= gap_b_last):
            opt_strat = GAP_A_IDX
            self.alignment_score = gap_a_last

        elif align_last >= gap_b_last:
            opt_strat = ALIGN_IDX
            self.alignment_score = align_last
        else:
            opt_strat = GAP_B_IDX
            self.alignment_score = gap_b_last

        while (i != 0) and (j != 0):
            if opt_strat == GAP_A_IDX:
                self.seqA_align += "-"
                self.seqB_align += self._seqB[j - 1]
            if opt_strat == ALIGN_IDX:
                self.seqA_align += self._seqA[i - 1]
                self.seqB_align += self._seqB[j - 1]
            if opt_strat == GAP_B_IDX:
                self.seqA_align += self._seqA[i - 1]
                self.seqB_align += "-"

            next_idx = self._back_3d[opt_strat, i, j]
            opt_strat, i, j = tuple(next_idx)

        temp_seq_a = self.seqA_align
        temp_seq_b = self.seqB_align
        self.seqA_align = "".join([temp_seq_a[i] for i in range(len(temp_seq_a) - 1, -1, -1)])
        self.seqB_align = "".join([temp_seq_b[i] for i in range(len(temp_seq_b) - 1, -1, -1)])

        return self.alignment_score, self.seqA_align, self.seqB_align



def read_fasta(fasta_file: str) -> Tuple[str, str]:
    """
    DO NOT MODIFY THIS FUNCTION! IT IS ALREADY COMPLETE!

    This function reads in a FASTA file and returns the associated
    string of characters (residues or nucleotides) and the header.
    This function assumes a single protein or nucleotide sequence
    per fasta file and will only read in the first sequence in the
    file if multiple are provided.

    Parameters:
        fasta_file: str
            name (and associated path if not in current working directory)
            of the Fasta file.

    Returns:
        seq: str
            String of characters from FASTA file
        header: str
            Fasta header
    """
    assert fasta_file.endswith(".fa"), "Fasta file must be a fasta file with the suffix .fa"
    with open(fasta_file) as f:
        seq = ""  # initializing sequence
        first_header = True
        for line in f:
            is_header = line.strip().startswith(">")
            # Reading in the first header
            if is_header and first_header:
                header = line.strip()  # reading in fasta header
                first_header = False
            # Reading in the sequence line by line
            elif not is_header:
                seq += line.strip().upper()  # generating full sequence
            # Breaking if more than one header is provided in the fasta file
            elif is_header and not first_header:
                break
    return seq, header
