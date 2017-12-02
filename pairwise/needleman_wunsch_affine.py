# !/usr/bin/python

# Global Alignment with Affine Gaps
# Modified Needleman-Wunsch Algorithm
# Biological Sequence Analysis - Page 29

import sys

import numpy as np

from Bio import SeqIO
from Bio.SubsMat import MatrixInfo

import common as cmn
from common import AffineElement

def main(seqs, d, e):

    with open(seqs[0], "rU") as handle:
        seq1_record = SeqIO.parse(handle, "fasta").__next__()
    with open(seqs[1], "rU") as handle:
        seq2_record = SeqIO.parse(handle, "fasta").__next__()

    print("Aligning Sequences with FASTA ID %s and %s \n" % (seq1_record.id, seq2_record.id))

    seq1 = seq1_record.seq
    seq2 = seq2_record.seq

    # Build alignment matrix
    mtx = build_alignment_matrix(seq1, seq2, d, e)

    # Get alignments by tracing back through matrix
    alignment = traceback(seq1, seq2, mtx)

    # Color aligns according to matches
    cmn.color_print_alignments(seq1_record, seq2_record, alignment, 150)

def build_alignment_matrix(x, y, d, e):

    # Initialize a matrix with all 0 values
    mtx = np.array([ [ AffineElement() for i in range(len(x)) ] for j in range(len(y)) ])

    # Initialize the initial condition
    mtx[0][0] = AffineElement(0, 0, -1, -1)

    # Let every value in the top row be equal to the gap penalty (theoretical max value)
    for i in range(1, len(x)):
        mtx[0][i] = AffineElement(0, -d - (i - 1) * e, i - 1, 0)

    # Let every value in the left column be equal to the gap penalty (theoretical max value)
    for j in range(1, len(y)):
        mtx[j][0] = AffineElement(0, -d - (j - 1) * e, 0, j - 1)

    # For the rest of the cells, use the recursive function
    # We choose either a gap or a pairing
    # Pairings log-odds ratios are taken from the substitution matrix (Blosum50)
    # Gap scores are chosen using the gap-odds penalty
    for n in range(1, len(x)):

        for m in range(1, len(y)):

            pair_log_odds = (MatrixInfo.blosum50[(x[n], y[m])] if MatrixInfo.blosum50.keys().__contains__((x[n], y[m])) else MatrixInfo.blosum50[(y[m], x[n])])
            pair = max(
                mtx[m-1][n-1].m + pair_log_odds,
                mtx[m-1][n-1].i + pair_log_odds
            )
            gap1 = max(
                mtx[m][n-1].m - d,
                mtx[m][n-1].i - e
            )
            gap2 = max(
                mtx[m-1][n].m - d,
                mtx[m-1][n].i - e
            )

            mtx[m][n].m = pair

            maxScore = max(gap1, gap2)
            mtx[m][n].i = maxScore

            if maxScore < pair:
                mtx[m][n].x = n - 1
                mtx[m][n].y = m - 1
            elif maxScore == gap1:
                mtx[m][n].x = n - 1
                mtx[m][n].y = m
            else:
                mtx[m][n].x = n
                mtx[m][n].y = m - 1

    return mtx

def traceback(x, y, mtx):

    aligns = [ "", "" ]
    (n, m) = (len(mtx[0]) - 1, len(mtx) - 1)
    while (n, m) != (-1, -1):

        cell = mtx[m][n]
        # Pair
        if cell.x == n - 1 and cell.y == m - 1:
            aligns[0] = x[n] + aligns[0]
            aligns[1] = y[m] + aligns[1]
        # Gap 1
        if cell.x == n - 1 and cell.y == m:
            aligns[0] = x[n] + aligns[0]
            aligns[1] = "-" + aligns[1]
        # Gap 2
        if cell.x == n and cell.y == m - 1:
            aligns[0] = "-" + aligns[0]
            aligns[1] = y[m] + aligns[1]
        n = cell.x
        m = cell.y

    return aligns

if __name__ == "__main__":

    seqs = []
    d = 0
    e = 0
    for arg in sys.argv:
        if arg.find("=") != -1:
            argName = arg.split("=")[0]
            argVal = arg.split("=")[1]
            if argName == "--seq":
                seqs.append(arg.split("=")[1])
            if argName == "--go-penalty":
                d = float(argVal)
            if argName == "--ge-penalty":
                e = float(argVal)

    main(seqs, d, e)