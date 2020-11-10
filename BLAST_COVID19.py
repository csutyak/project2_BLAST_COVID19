import os
import numpy as np

SARS_N_FILENAME = "SARS_N.txt"
COVID19_N_FILENAME = "Covid19_N.txt"
BLASTALIGNMENT_CODON_FILENAME1 = "BlastAlignment1.txt"
BLASTALIGNMENT_CODON_FILENAME2 = "BlastAlignment2.txt"
BLASTALIGNMENT_GENE_FILENAME1 = "BlastAlignmentGene1.txt"
BLASTALIGNMENT_GENE_FILENAME2 = "BlastAlignmentGene2.txt"

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

def findIndels(sequence1, sequence2):
    totalInsertions = 0
    totalDeletions = 0
    totalPointMutations = 0
    totalIndels = 0

    indelsDict = dict()

    for x in range(len(sequence1)):
        if sequence1[x] != sequence2[x]:
            if sequence1[x] == '_':
                totalIndels += 1
                totalInsertions += 1
                indelsDict[x] = 1
            elif sequence2[x] == '_':
                totalIndels += 1
                totalDeletions += 1
                indelsDict[x] = 1
            else:
                totalIndels += 1
                totalPointMutations += 1
                indelsDict[x] = 1

    print("Total insertions: ", totalInsertions)
    print("Total deletions: ", totalDeletions)
    print("Total point mutations: ", totalPointMutations)
    print("Total indels: ", totalIndels)
    return indelsDict


def findBaseMutations(indelsDict, sequence1, sequence2):
    synonymousCount = 0
    nonSynonymousCount = 0
    for x in range(len(sequence1)):
        if sequence1[x] != sequence2[x]:
            if (x//3) in indelsDict:
                nonSynonymousCount += 1
            else:
                synonymousCount += 1

    totalMutations = synonymousCount + nonSynonymousCount
    print("Total synonymous mutations: ", synonymousCount)
    print("Total non-synonymous mutations: ", nonSynonymousCount)
    print("Total mutations: ", totalMutations)
    print("Percent identical: ", 1 - totalMutations / len(sequence1))

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

LOW_SCORE = -9999999999
MISMATCH_SCORE = -1
GAP_SCORE = -2
MATCH_SCORE = 1
def createAlignmentMatrix(sequence1, sequence2):

    matrix_length = len(sequence1) + 1
    matrix_width = len(sequence2) + 1
    matrix_size = matrix_width * matrix_length
    alignment_matrix = np.arange(matrix_size).reshape((matrix_length, matrix_width))

    # prepare first row(length) of global alignment with all _
    for lengthItr in range(matrix_length):
        alignment_matrix[lengthItr][0] = lengthItr * GAP_SCORE

    # prepare width
    for widthItr in range(matrix_width):
        alignment_matrix[0][widthItr] = widthItr * GAP_SCORE

    # Find all scores
    for lengthItr in range(1, matrix_length):
        for widthItr in range(1, matrix_width):
            # get top left and diagonal scores
            diagonalScore = alignment_matrix[lengthItr - 1][widthItr - 1]
            leftScore = alignment_matrix[lengthItr][widthItr - 1]
            topScore = alignment_matrix[lengthItr - 1][widthItr]

            # get matchscore sequence1 = length sequence2 = width
            if sequence1[lengthItr - 1] == sequence2[widthItr - 1]:
                matchScore = MATCH_SCORE
            else:
                matchScore = MISMATCH_SCORE

            # Compare top left and diagonal to see which one is max score
            maxScore = max(diagonalScore + matchScore, leftScore + GAP_SCORE)
            maxScore = max(maxScore, topScore + GAP_SCORE)

            # set matrix place to maxScore and matchScore
            alignment_matrix[lengthItr][widthItr] = maxScore

    # CREATE ALIGNMENT LINE
    sequence1Alignment = ''
    sequence2Alignment = ''

    lengthIndex = matrix_length - 1
    widthIndex = matrix_width - 1
    #build working string to save the order of the alignment line
    workingAlignmentString = ''
    loopFlag = True
    #loops until full line is found
    while loopFlag:
        if (widthItr == 0 and lengthItr == 0):
            break
        diagonalScore = None
        if lengthItr == 0:
            diagonalScore = LOW_SCORE
            topScore = LOW_SCORE
        else:
            topScore = alignment_matrix[lengthItr - 1][widthItr]

        if widthItr == 0:
            diagonalScore = LOW_SCORE
            leftScore = LOW_SCORE
        else:
            leftScore = alignment_matrix[lengthItr][widthItr - 1]

        if diagonalScore is None:
            diagonalScore = alignment_matrix[lengthItr - 1][widthItr - 1]

        maxScore = max(diagonalScore, leftScore)
        maxScore = max(maxScore, topScore)

        currentScore = alignment_matrix[lengthItr][widthItr]

        if currentScore == diagonalScore + MATCH_SCORE or currentScore == diagonalScore + MISMATCH_SCORE:
            bestScore = "Dia"
        if (currentScore == topScore + GAP_SCORE or leftScore == LOW_SCORE) and topScore > diagonalScore:
            bestScore = "Top"
        if (currentScore == leftScore + GAP_SCORE or topScore == LOW_SCORE) and leftScore > topScore and leftScore > diagonalScore:
            bestScore = "Left"

        if bestScore == "Dia":
            workingAlignmentString = workingAlignmentString + "D"
            lengthItr = lengthItr - 1
            widthItr = widthItr - 1
        # go top
        elif bestScore == "Top":
            workingAlignmentString = workingAlignmentString + "T"
            lengthItr = lengthItr - 1

        # go left
        elif bestScore == "Left":
            workingAlignmentString = workingAlignmentString + "L"
            widthItr = widthItr - 1

    #build sequences from alignment string
    lengthItr = 0
    widthItr = 0
    for index in range(len(workingAlignmentString) - 1, 0, -1):
        # diagonal
        
        if workingAlignmentString[index] == "D":
            sequence1Alignment = sequence1Alignment + sequence1[lengthItr]
            sequence2Alignment = sequence2Alignment + sequence2[widthItr]
            lengthItr += 1 
            widthItr += 1
        elif workingAlignmentString[index] == "T":
            sequence1Alignment = sequence1Alignment + sequence1[lengthItr]
            sequence2Alignment = sequence2Alignment + "_"
            lengthItr += 1
        elif workingAlignmentString[index] == "L":
            sequence1Alignment = sequence1Alignment + "_"
            sequence2Alignment = sequence2Alignment + sequence2[widthItr]
            widthItr += 1

    return alignment_matrix, sequence1Alignment, sequence2Alignment

def convertToCodons (sequence):
    codon_list = [(sequence[x:x + 3]) for x in range(0, len(sequence), 3)]

    for x in range(len(codon_list)):
        codon_list[x] = codonTable[codon_list[x]]

    codon_string = ''

    for x in codon_list:
        codon_string += x

    return codon_string


sars_codons_string = convertToCodons(sars_n_string)
covid_codons_string = convertToCodons(covid_n_string)

gene_alignment_matrix, gene1Alignment, gene2Alignment = createAlignmentMatrix(sars_n_string, covid_n_string)
codon_alignment_matrix, codon1Alignment, codon2Alignment = createAlignmentMatrix(sars_codons_string, covid_codons_string)

print("Gene data summary:")
indelsDict = findIndels(codon1Alignment, codon2Alignment)
findBaseMutations(indelsDict, gene1Alignment, gene2Alignment)

showGeneAlignment = False
showCodonOutput = False
showBlastOutput = True

print()
print("Gene matrix")
print(gene_alignment_matrix)
if showGeneAlignment:
    print(gene1Alignment)
    print(gene2Alignment)

if showCodonOutput:
    print()
    print("Codon matrix and alignment")
    print(codon_alignment_matrix)
    print(codon1Alignment)
    print(codon2Alignment)

if showBlastOutput:
    print()
    print("Summary of BLAST data:")
    blastAlignment1Codon = fileToString(BLASTALIGNMENT_CODON_FILENAME1)
    blastAlignment2Codon = fileToString(BLASTALIGNMENT_CODON_FILENAME2)

    blastAlignment1Gene = fileToString(BLASTALIGNMENT_GENE_FILENAME1)
    blastAlignment2Gene = fileToString(BLASTALIGNMENT_GENE_FILENAME2)

    blastIndelsDict = findIndels(blastAlignment1Codon, blastAlignment2Codon)
    findBaseMutations(blastIndelsDict, blastAlignment1Gene, blastAlignment2Gene)
