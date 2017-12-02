# !/usr/bin/python

# Common alignment tools

# The map of similar amino acids
similarities = {

    # Positively Charged Amino Acids
    "R": [ "R", "H", "K" ],
    "H": [ "R", "H", "K" ],
    "K": [ "R", "H", "K" ],

    # Negatively Charged Amino Acids
    "D": [ "D", "E" ],
    "E": [ "D", "E" ],

    # Small Polar Amino Acids
    "S": [ "S", "T", "C" ],
    "T": [ "S", "T", "C" ],
    "C": [ "S", "T", "C" ],

    # Large Polar Amino Acids
    "Y": [ "Y", "N", "Q" ],
    "N": [ "Y", "N", "Q" ],
    "Q": [ "Y", "N", "Q" ],

    # Small Non-Polar Amino Acids
    "G": [ "G", "A", "V", "L", "I" ],
    "A": [ "G", "A", "V", "L", "I" ],
    "V": [ "G", "A", "V", "L", "I" ],
    "L": [ "G", "A", "V", "L", "I" ],
    "I": [ "G", "A", "V", "L", "I" ],

    # Large Non-Polar Amino Acids
    "M": [ "M", "F", "W", "P" ],
    "F": [ "M", "F", "W", "P" ],
    "W": [ "M", "F", "W", "P" ],
    "P": [ "M", "F", "W", "P" ]
}
# Determines if amino acids are similar in a biochemical fashion
def aa_are_similar(x, y):

    global similarities

    if x in similarities:
        return y in similarities[x]
    else:
        return False

# Colors the alignments
# Red for gap alignments
# Green for perfect matches
# Cyan for biochemically similar amino acids
# Yellow for mismatches
def color_print_alignments(seq1, seq2, alignment, n):

    for i in range(0, len(alignment[0]), n):
        aligns = ["", ""]
        for x in range(i, i + n if i + n <= len(alignment[0]) else len(alignment[0])):
            if alignment[0][x] == alignment[1][x]:
                aligns[0] += "\033[1;30;42m" + alignment[0][x]
                aligns[1] += "\033[1;30;42m" + alignment[1][x]
            elif aa_are_similar(alignment[0][x], alignment[1][x]):
                aligns[0] += "\033[1;30;46m" + alignment[0][x]
                aligns[1] += "\033[1;30;46m" + alignment[1][x]
            elif alignment[0][x] == "-" or alignment[1][x] == "-":
                aligns[0] += "\033[1;30;41m" + alignment[0][x]
                aligns[1] += "\033[1;30;41m" + alignment[1][x]
            else:
                aligns[0] += "\033[1;30;43m" + alignment[0][x]
                aligns[1] += "\033[1;30;43m" + alignment[1][x]
        print("\033[0;30;0m" + seq1.id + "\t" + str(i) + "\t" + aligns[0])
        print("\033[0;30;0m" + seq2.id + "\t" + str(i) + "\t" + aligns[1])
        print("\033[0;30;0m")

class AffineElement(object):

    def __init__(self, m=0, i=0, x=0, y=0):
        self.m = m
        self.i = i
        self.x = x
        self.y = y

    def __str__(self):
        return "[M=%s, I=%s] -> (%s, %s)" % (self.m, self.i, self.x, self.y)

    def __repr__(self):
        return self.__str__()