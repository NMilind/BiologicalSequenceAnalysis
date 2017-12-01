# !/usr/bin/python

# Local Alignment
# Smith-Waterman Algorithm
# Biological Sequence Analysis - Page 22

import sys

import numpy as np

from Bio import SeqIO
from Bio.SubsMat import MatrixInfo

import common as cmn

def main(seqs, d):

    with open(seqs[0], "rU") as handle:
        seq1_record = SeqIO.parse(handle, "fasta").__next__()
    with open(seqs[1], "rU") as handle:
        seq2_record = SeqIO.parse(handle, "fasta").__next__()

    seq1 = seq1_record.seq
    seq2 = seq2_record.seq

    # Build alignment matrix
    mtx = build_alignment_matrix(seq1, seq2, d)

    # Get alignments by tracing back through matrix
    alignment = traceback(seq1, seq2, mtx)

    cmn.color_print_alignments(alignment, 150)

def build_alignment_matrix(x, y, d):

    # Initialize a matrix with all 0 values
    mtx = np.array([[[0, -1, -1]] * len(x)] * len(y))

    # We choose either a gap or a pairing
    # Pairings log-odds ratios are taken from the substitution matrix (Blosum50)
    # Gap scores are chosen using the gap-odds penalty
    # We only take the max score if it is greater than 0
    for n in range(1, len(x)):
        for m in range(1, len(y)):
            pair = mtx[m-1][n-1][0] + (MatrixInfo.blosum50[(x[n], y[m])] if MatrixInfo.blosum50.keys().__contains__((x[n], y[m])) else MatrixInfo.blosum50[(y[m], x[n])])
            gap1 = mtx[m][n-1][0] - d
            gap2 = mtx[m-1][n][0] - d
            maxScore = max(pair, gap1, gap2, 0)
            if maxScore == pair:
                mtx[m][n] = (maxScore, n - 1, m - 1)
            elif maxScore == gap1:
                mtx[m][n] = (maxScore, n - 1, m)
            elif maxScore == gap2:
                mtx[m][n] = (maxScore, n, m - 1)

    return mtx

def traceback(x, y, mtx):

    aligns = [ "", "" ]
    (n, m) = find_max(mtx)
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

def find_max(mtx):

    maxVal = mtx[0][0][0]
    (x, y) = (0, 0)
    for n in range(len(mtx[0])):
        for m in range(len(mtx)):
            if mtx[m][n][0] > maxVal:
                (x, y) = (n, m)
                maxVal = mtx[m][n][0]
    return (x, y)

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
