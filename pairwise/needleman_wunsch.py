# !/usr/bin/python

# Global Alignment
# Needleman-Wunsch Algorithm
# Biological Sequence Analysis - Page 19

import sys

import numpy as np

from Bio import SeqIO
from Bio.SubsMat import MatrixInfo

def main(seqs, d):

    with open(seqs[0], "rU") as handle:
        seq1_record = SeqIO.parse(handle, "fasta").__next__()
    with open(seqs[1], "rU") as handle:
        seq2_record = SeqIO.parse(handle, "fasta").__next__()

    seq1 = seq1_record.seq
    seq2 = seq2_record.seq

    mtx = build_alignment_matrix(seq1, seq2, d)

    alignment = traceback(seq1, seq2, mtx)

    print(alignment[0])
    print(alignment[1])

def build_alignment_matrix(x, y, d):

    # Initialize a matrix with all 0 values
    mtx = np.array([[[0, 0, 0]] * len(x)] * len(y))

    # Initialize the initial condition
    mtx[0][0] = (0, -1, -1)

    # Let every value in the top row be equal to the gap penalty (theoretical max value)
    for i in range(1, len(x)):
        mtx[0][i] = (-i * d, i - 1, 0)

    # Let every value in the left column be equal to the gap penalty (theoretical max value)
    for j in range(1, len(y)):
        mtx[j][0] = (-j * d, 0, j - 1)

    # For the rest of the cells, use the recursive function
    # We choose either a gap or a pairing
    # Pairings log-odds ratios are taken from the substitution matrix (Blosum50)
    # Gap scores are chosen using the gap-odds penalty
    for n in range(1, len(x)):
        for m in range(1, len(y)):
            pair = mtx[m-1][n-1][0] + (MatrixInfo.blosum50[(x[n], y[m])] if MatrixInfo.blosum50.keys().__contains__((x[n], y[m])) else MatrixInfo.blosum50[(y[m], x[n])])
            gap1 = mtx[m][n-1][0] - d
            gap2 = mtx[m-1][n][0] - d
            maxScore = max(pair, gap1, gap2)
            if maxScore == pair:
                mtx[m][n] = (maxScore, n - 1, m - 1)
            elif maxScore == gap1:
                mtx[m][n] = (maxScore, n - 1, m)
            else:
                mtx[m][n] = (maxScore, n, m - 1)

    return mtx

def traceback(x, y, mtx):

    aligns = [ "", "" ]
    (n, m) = (len(mtx[0]) - 1, len(mtx) - 1)
    while (n, m) != (-1, -1):

        cell = mtx[m][n]
        # Pair
        if cell[1] == n - 1 and cell[2] == m - 1:
            aligns[0] = x[n] + aligns[0]
            aligns[1] = y[m] + aligns[1]
        # Gap 1
        if cell[1] == n - 1 and cell[2] == m:
            aligns[0] = x[n] + aligns[0]
            aligns[1] = "-" + aligns[1]
        # Gap 2
        if cell[1] == n and cell[2] == m - 1:
            aligns[0] = "-" + aligns[0]
            aligns[1] = y[m] + aligns[1]
        n = cell[1]
        m = cell[2]

    return aligns

if __name__ == "__main__":

    seqs = []
    d = 0
    for arg in sys.argv:
        if arg.find("=") != -1:
            argName = arg.split("=")[0]
            argVal = arg.split("=")[1]
            if argName == "--seq":
                seqs.append(arg.split("=")[1])
            if argName == "--go-penalty":
                d = float(argVal)

    main(seqs, d)