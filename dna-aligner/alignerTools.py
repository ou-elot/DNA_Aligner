import numpy
import sys

def read_sequences(filepath):
    """
    Reads a file containing two DNA sequences separated by a newline.

    :param filepath: Path to the file
    :return: Tuple (seq1, seq2)
    """
    with open(filepath, "r") as file:
        lines = file.read().strip().split("\n")
        if len(lines) != 2:
            raise ValueError("File must contain exactly two lines for two DNA sequences.")
        return (lines[0], lines[1])  # Returns as (sequence1, sequence2)

def global_alignment(seq1: str, seq2: str, match_reward: int, mismatch_penalty: int, gap_opening: int, gap_extension: int):
    """
    Compute the affine alignment of two strings based on match reward, mismatch penalty, 
    gap opening penalty, and gap extension penalty.
    """
    lower = []
    middle = []
    upper = []

    for i in range (len(seq1)+1):
        lower.append([])
        middle.append([])
        upper.append([])
        for j in range (len(seq2)+1):
            lower[i].append(0)
            middle[i].append(0)
            upper[i].append(0)

    for j in range (len(seq2)+1):
        lower[0][j] = float('-inf')
    for i in range (len(seq1)+1):
        lower[1][0] = float('-inf')

    middle[0][0] = 0
    middle[0][1] = -1 * gap_opening
    middle[1][0] = -1 * gap_opening
    lower[1][0] = -1 * gap_opening
    upper[0][1] = -1 * gap_opening

    for j in range (2, len(seq2)+1):
        middle[0][j] = middle[0][j-1] - gap_extension_penalty
        upper[0][j] = upper[0][j-1] - gap_extension_penalty
    for i in range (2, len(seq1)+1):
        middle[i][0] = middle[i-1][0] - gap_extension_penalty
        lower[i][0] = lower[i-1][0]- gap_extension_penalty

    for i in range (1,len(seq1)+1):
        for j in range (1,len(seq2)+1):
            lower[i][j] = max(lower[i-1][j]-gap_extension_penalty, middle[i-1][j]-gap_opening_penalty)
            upper[i][j] = max(upper[i][j-1]-gap_extension_penalty, middle[i][j-1]-gap_opening_penalty)
            if (seq1[i-1] == seq2[j-1]):
                match = match_reward
            else:
                match = -1 * mismatch_penalty
            middle[i][j] = max(lower[i][j], middle[i-1][j-1] + match, upper[i][j])
    alignScore = middle[-1][-1]

    i = len(seq1)
    j = len(seq2)
    curr = 'M'
    align1 = []
    align2 = []
    while i>0 or j>0:
        if curr == 'M':
            if s[i-1] == t[j-1]:
                match = match_reward
            else:
                match = -1*mismatch_penalty
            if i>0 and j>0 and middle[i-1][j-1] + match == middle[i][j]:
                align1.append(s[i-1])
                align2.append(t[j-1])
                i = i-1
                j = j-1
            elif middle[i][j] == lower[i][j]:
                curr = 'L'
            elif middle[i][j] == upper[i][j]:
                curr = 'U'
        elif curr =='L':
            if i>0 and lower[i][j]==lower[i-1][j] - gap_extension_penalty:
                align1.append(s[i-1])
                align2.append("-")
                i = i-1
            else:
                align1.append(s[i-1])
                align2.append("-")
                i = i-1
                curr = 'M'
        elif curr == 'U':
            if j>0 and upper[i][j] == upper[i][j-1] - gap_extension_penalty:
                align1.append("-")
                align2.append(t[j-1])
                j = j-1
            else:
                align1.append("-")
                align2.append(t[j-1])
                j = j-1
                curr = 'M'
    alignment1 = ""
    alignment2 = ""
    for x in reversed(range(len(sAlign))):
        alignment1 = alignment1 + align1[x]
        
    for y in reversed(range(len(tAlign))):
        alignment2 = alignment2 + align2[y]
    return (alignScore,sFinal,tFinal)

def local_alignment(seq1: str, seq2: str, match_reward: int, mismatch_penalty: int, indel_penalty: int):
    """
    Compute the local alignment of two strings based on match reward, mismatch penalty, and indel penalty.
    """
    s = seq1
    t = seq2
    movements = [] # 1 for down, 2 for lateral, 3 for diagonal
    alignments = []
    for i in range (len(s)+1):
        alignments.append([])
        movements.append([])
        for j in range (len(t)+1):
            alignments[i].append(0)
            movements[i].append(0)
    alignments[0][0] = 0 # two empty characters
    for i in range (1,len(alignments[0])):
        alignments[0][i] = 0
        movements[0][i] = 0
    for j in range (1, len(alignments)):
        alignments[j][0] = 0
        movements[j][0] = 0
    for i in range (1, len(s)+1):
        for j in range (1, len(t)+1):
            if s[i-1] == t[j-1]:
                match = alignments[i-1][j-1] + match_reward
            else:
                match = alignments[i-1][j-1] - mismatch_penalty
            score1 = alignments[i-1][j] - indel_penalty
            score2 = alignments[i][j-1] - indel_penalty
            alignments[i][j] = max(0, score1, score2, match)
            if max(0,score1,score2,match) == score1:
                movements[i][j] = 1
            elif max(0,score1,score2,match) == score2:
                movements[i][j] = 2
            elif max(0,score1,score2,match) == match:
                movements[i][j] = 3
            elif max(0,score1,score2,match) == 0:
                movements[i][j] = 4

    maxScore = float("-inf")
    iMax = -1
    jMax = -1
    for x in range(len(alignments)):
        for y in range(len(alignments[x])):
            if alignments[x][y]>maxScore:
                maxScore = alignments[x][y]
                iMax = x
                jMax = y
    i = iMax
    j = jMax
    sFinal = ""
    tFinal = ""
    while movements[i][j] !=0 and movements[i][j] !=4:
        if movements[i][j] == 1:
            sFinal = s[i-1] + sFinal 
            tFinal = "-" + tFinal  
            i = i-1
        elif movements[i][j] == 2:
            sFinal = "-" + sFinal  
            tFinal = t[j-1] + tFinal
            j = j-1
        elif movements[i][j] == 3:
            sFinal = s[i-1] + sFinal 
            tFinal = t[j-1] + tFinal
            i = i-1
            j=j-1
        
    return (maxScore, sFinal, tFinal)


    
