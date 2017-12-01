# !/usr/bin/python

# Common alignment tools

def color_alignments(alignment):

    aligns = ["", ""]

    for x in range(len(alignment[0])):
        if alignment[0][x] == alignment[1][x]:
            aligns[0] += "\033[1;30;42m" + alignment[0][x]
            aligns[1] += "\033[1;30;42m" + alignment[1][x]
        elif alignment[0][x] == "-" or alignment[1][x] == "-":
            aligns[0] += "\033[1;30;41m" + alignment[0][x]
            aligns[1] += "\033[1;30;41m" + alignment[1][x]
        else:
            aligns[0] += "\033[1;30;43m" + alignment[0][x]
            aligns[1] += "\033[1;30;43m" + alignment[1][x]

    return aligns