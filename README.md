# Biological Sequence Analysis

## Pairwise

This section deals with rudimentary pairwise alignment algorithms for two sequences. I have currently implemented the following algorithms, and have included sample FASTA protein files to test with the algorithms:

### Needleman-Wunsch Algorithm (Global Alignment)

This force both sequences to be aligned completely and globally. It builds an alignment matrix using the Blosum50 log-odds scoring matrix, and traces back the sequence with the highest log-odds ratios. A gaps-odd penalty is included as an argument when running the program, and two sequences are taken in by the program to align properly. The alignment is printed out with gaps, and colored according to similarity (Red = gap align, Green = perfect align, Yellow = misalign, Cyan = biochemically close amino acids).

```
python3 needleman_wunsch.py --seq=fasta/homo_sapiens_hemoglobin.fasta --seq=fasta/homo_sapiens_myoglobin.py --go-penalty=6
```

### Smith-Waterman Algorithm (Local Alignment)

Similar method used as the Needleman-Wunsch algorithm, but the highest score in the matrix is chosen rather than the end of the sequence. No negative scores are allowed, so that only the positive-scoring alignments are kept. Tends to reveal domains or motiffs that are common/conserved in two protein sequences.

```
smith_waterman.py --seq=fasta/homo_sapiens_hemoglobin.fasta --seq=fasta/homo_sapiens_leukotriene_c4_synthase.fasta --go-penalty=8
```

## Citations

Much of the work in this repository is influenced by "Biological Sequence Analysis - Probabilistic Models of Proteins and Nucleic Acids," by R. Durbin, S. Eddy, A. Krogh, and G. Mitchison.

