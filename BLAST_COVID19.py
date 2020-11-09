import os
import numpy as np

SARS_N_FILENAME = "SARS_N.txt"
COVID19_N_FILENAME = "Covid19_N.txt"


def fileToString(filename):
    inputFile = open(filename, "r")
    inputlines = inputFile.readlines()
    CompleteString = ""
    for line in inputlines:
        workingString = line.replace(" ", "")
        workingString = workingString.replace("\n", "")
        CompleteString = CompleteString + workingString
    return CompleteString


sars_n_string = fileToString(SARS_N_FILENAME)
covid_n_string = fileToString(COVID19_N_FILENAME)

# returns alignment Strings
# T = Top
# L = Left
# D = Diagonal
# F = Finished
LOW_SCORE = -9999999999


def recursiveAlignmentLine(DTLF, lengthItr, widthItr):
    global sequence1Alignment
    global sequence2Alignment
    # base case
    if (widthItr == 0 and lengthItr == 0):
        return

    # check which number is higher, top, left, or diagonal
    diagonalScore = None
    if lengthItr == 0:
        diagonalScore = LOW_SCORE
        topScore = LOW_SCORE
    else:
        topScore = allignment_matrix[lengthItr - 1][widthItr]

    if widthItr == 0:
        diagonalScore = LOW_SCORE
        leftScore = LOW_SCORE
    else:
        leftScore = allignment_matrix[lengthItr][widthItr - 1]

    if diagonalScore is None:
        diagonalScore = allignment_matrix[lengthItr - 1][widthItr - 1]

    maxScore = max(diagonalScore, leftScore)
    maxScore = max(maxScore, topScore)

    currentScore = allignment_matrix[lengthItr][widthItr]

    if currentScore == diagonalScore + MATCH_SCORE or currentScore == diagonalScore + MISMATCH_SCORE:
        bestScore = "Dia"
    if (currentScore == topScore + GAP_SCORE or leftScore == LOW_SCORE) and topScore > diagonalScore:
        bestScore = "Top"
    if (
            currentScore == leftScore + GAP_SCORE or topScore == LOW_SCORE) and leftScore > topScore and leftScore > diagonalScore:
        bestScore = "Left"

    # go diagonal
    if bestScore == "Dia":
        print("Dia")
        print(lengthItr, widthItr)
        recursiveAlignmentLine("D", lengthItr - 1, widthItr - 1)
        sequence1Alignment = sequence1Alignment + sequence1[lengthItr - 1]
        sequence2Alignment = sequence2Alignment + sequence2[widthItr - 1]

    # go top
    elif bestScore == "Top":
        print("Top")
        print(lengthItr, widthItr)
        recursiveAlignmentLine("T", lengthItr - 1, widthItr)
        sequence1Alignment = sequence1Alignment + sequence1[lengthItr - 1]
        sequence2Alignment = sequence2Alignment + "_"
    # go left

    elif bestScore == "Left":
        print("Left")
        print(lengthItr, widthItr)
        recursiveAlignmentLine("L", lengthItr, widthItr - 1)
        sequence1Alignment = sequence1Alignment + "_"
        sequence2Alignment = sequence2Alignment + sequence2[widthItr - 1]


MISMATCH_SCORE = -1
GAP_SCORE = -2
MATCH_SCORE = 1








codonTable = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
    'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
}

sars_codons = [(sars_n_string[x:x + 3]) for x in range(0, len(sars_n_string), 3)]

for x in range(len(sars_codons)):
    sars_codons[x] = codonTable[sars_codons[x]]

sars_codons_string = ''

for x in sars_codons:
    sars_codons_string += x

covid_codons = [(covid_n_string[x:x + 3]) for x in range(0, len(covid_n_string), 3)]

for x in range(len(covid_codons)):
    covid_codons[x] = codonTable[covid_codons[x]]

covid_codons_string = ''

for x in covid_codons:
    covid_codons_string += x










# sequence2 = "ACGGCTC"
# sequence1 = "ATGGCCTC"
sequence2 = covid_codons_string
sequence1 = sars_codons_string
#sequence2 = sars_n_string
#sequence1 = covid_n_string

matrix_length = len(sequence1) + 1
matrix_width = len(sequence2) + 1
matrix_size = matrix_width * matrix_length
allignment_matrix = np.arange(matrix_size).reshape((matrix_length, matrix_width))

# prepare first row(length) of global alignment with all _
lengthAdd = 0
for lengthItr in range(matrix_length):
    allignment_matrix[lengthItr][0] = lengthAdd
    lengthAdd += GAP_SCORE

# prepare width
widthAdd = 0
for widthItr in range(matrix_width):
    allignment_matrix[0][widthItr] = widthAdd
    widthAdd += GAP_SCORE

# Find all scores
for lengthItr in range(1, matrix_length):
    for widthItr in range(1, matrix_width):
        # get top left and diagonal scores
        diagonalScore = allignment_matrix[lengthItr - 1][widthItr - 1]
        leftScore = allignment_matrix[lengthItr][widthItr - 1]
        topScore = allignment_matrix[lengthItr - 1][widthItr]

        # get matchscore sequence1 = length sequence2 = width
        if sequence1[lengthItr - 1] == sequence2[widthItr - 1]:
            matchScore = MATCH_SCORE
        else:
            matchScore = MISMATCH_SCORE

        # Compare top left and diagonal to see which one is max score
        maxScore = max(diagonalScore + matchScore, leftScore + GAP_SCORE)
        maxScore = max(maxScore, topScore + GAP_SCORE)

        # set matrix place to maxScore and matchScore
        allignment_matrix[lengthItr][widthItr] = maxScore

lengthIndex = matrix_length - 1
widthIndex = matrix_width - 1
sequence1Alignment = ""
sequence2Alignment = ""
recursiveAlignmentLine("F", lengthIndex, widthIndex)

print(allignment_matrix)

print(sequence1Alignment)
print(sequence2Alignment)

# Codons
# Start == "ATG"
# Stop == "TAA" or "TAG" or "TGA"



