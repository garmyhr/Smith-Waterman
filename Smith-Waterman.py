import numpy as np
import sys

# Defining gap penalties #
openGap = 11
extGap = 1

# Reading the BLOSUM62 matrix to memory #
file = open("BLOSUM62.txt", 'r')
source_info = file.read()
source_lines = source_info.splitlines()
score_letters = source_lines[0].split()
blossum_matrix = np.zeros((len(score_letters), len(score_letters)), dtype=int)

for i in range(1, len(source_lines)):
    scores = source_lines[i].split()
    scores.pop(0)
    scores_int = []
    for item in scores:
        scores_int.append(int(item))
    blossum_matrix[i-1, 0:] = scores_int

# Method for getting the sequence from FASTA-file #
def get_sequence(fasta):
    file = open(fasta, 'r')
    lines = file.readlines()
    sequence = ""
    for i in range(1, len(lines)):
        lines[i] = lines[i].strip('\n')
        sequence += lines[i]

    return list(sequence)

# Gets the index of a letter in proteinalphabet #
def get_letter_index(letter):
    for i in range(0, len(score_letters)):
        if (score_letters[i] == letter):
            return i

# Gets the matchscore from BLOSUM62-matrix #
def match_score(letter_1, letter_2):

    i = get_letter_index(letter_1)
    j = get_letter_index(letter_2)

    return blossum_matrix[i, j]

# The smith-waterman algorithm #
# Takes two sequences and prints the local sequence alignment to the terminal #
# Also prints the raw allignment score with it's location in the matrix #
def smith_waterman(seq_1, seq_2):

    # Initializing variables used in the algorithm #
    m, n = len(seq_1), len(seq_2)
    S = np.zeros((m+1, n+1), dtype=int)
    Ix = np.zeros((m+1, n+1), dtype=int)
    Iy = np.zeros((m+1, n+1), dtype=int)
    align1, middle_string, align2 = '', '', ''
    max_score, max_i, max_j = 0, 0, 0


    # Scoring process #
    for i in range(1, m+1):
        for j in range(1, n+1):

            seq_1_letter = seq_1[i-1]
            seq_2_letter = seq_2[j-1]

            Iy[i][j] = max((S[i-1, j] - openGap - extGap), (Iy[i-1, j]-extGap))
            Ix[i][j] = max((S[i, j-1] - openGap - extGap), (Ix[i, j-1]-extGap))
            score_diagonal = S[i-1,j-1] + match_score(seq_1_letter, seq_2_letter)

            S[i][j] = max(0, score_diagonal, Ix[i,j], Iy[i,j])

            if S[i][j] >= max_score:
                max_score = S[i][j]
                max_i = i; max_j = j


    # Traceback #
    i, j = max_i, max_j
    while S[i][j] != 0:

        if S[i][j] == S[i-1][j-1] + match_score(seq_1[i-1], seq_2[j-1]):
            align1 += seq_1[i-1]
            align2 += seq_2[j-1]
            i -= 1; j -= 1
        elif S[i][j] == Iy[i][j]:
            align1 += seq_1[i-1]
            align2 += '-'
            i -= 1
        elif S[i][j] == Ix[i][j]:
            align1 += '-'
            align2 += seq_2[j-1]
            j -= 1

    # Reversing the strings created from traceback #
    align1 = align1[::-1]
    align2 = align2[::-1]

    # Create middle string to see matches and mismatches with positive score in BLOSUM62 #
    for i in range(0, len(align1)):

        letter_1 = align1[i]
        letter_2 = align2[i]

        if (letter_1 == '-' or letter_2 == '-'):
            middle_string += ' '

        elif (letter_1 == letter_2):
            middle_string += letter_1

        else:
            if match_score(letter_1, letter_2) > 0:
                middle_string += '+'
            else:
                middle_string += ' '


    # Formating and printing all the strings #
    for i in range(0, len(align1), 60):
        print("Seq_1: " + align1[i:i+60])
        print("       " + middle_string[i:i+60])
        print("Seq_2: " + align2[i:i+60])
        print("\n")

    print("Max score: %d" % max_score)
    print("Max score at: (%d, %d)\n" % (max_i, max_j))


seq_1 = get_sequence(sys.argv[1])
seq_2 = get_sequence(sys.argv[2])
smith_waterman(seq_1, seq_2)
