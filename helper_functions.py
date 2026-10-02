from Bio.Align import substitution_matrices

def global_alignment(seq1, seq2, scoring_function):
    """Global sequence alignment using the Needleman–Wunsch algorithm.

    Indels should be denoted with the "-" character.

    Parameters
    ----------
    seq1: str
        First sequence to be aligned.
    seq2: str
        Second sequence to be aligned.
    scoring_function: Callable

    Returns
    -------
    str
        First aligned sequence.
    str
        Second aligned sequence.
    float
        Final score of the alignment.

    Examples
    --------
    >>> global_alignment("abracadabra", "dabarakadara", lambda x, y: [-1, 1][x == y])
    ('-ab-racadabra', 'dabarakada-ra', 5.0)

    Other alignments are not possible.

    """
    #Initialise Matrix and set first column
    len1 = len(seq1)
    len2 = len(seq2)
    matrixBuild = []

    for i in range(len1 + 1):
        row = []
        for j in range(len2 + 1):
            row.append(0)
        matrixBuild.append(row)
    
    for i in range(1, len1 + 1):
        matrixBuild[i][0] = matrixBuild[i-1][0] + scoring_function(seq1[i-1], "-")
    for j in range(1, len2 + 1):
        matrixBuild[0][j] = matrixBuild[0][j-1] + scoring_function("-", seq2[j-1])

    for i in range(1, len1 + 1):
        for j in range(1, len2 + 1):
            char1 = seq1[i-1]
            char2 = seq2[j-1]
            
            # Calculate score if we match/mismatch these chars 
            diagonal_score = matrixBuild[i-1][j-1] + scoring_function(char1, char2)
            
            # Calculate score if we add a gap in seq 2
            down_score = matrixBuild[i-1][j] + scoring_function(char1, "-")
            
            # Calculate score if we add a gap in seq 1 
            right_score = matrixBuild[i][j-1] + scoring_function("-", char2)
            
            # Set highest
            matrixBuild[i][j] = max(diagonal_score, down_score, right_score)

    aligned_seq1 = ""
    aligned_seq2 = ""
    i = len1
    j = len2
    
    while i > 0 or j > 0:
        current_score = matrixBuild[i][j]
        
        # Check if from diag? (Match/Mismatch)
        if i > 0 and j > 0:
            char1 = seq1[i-1]
            char2 = seq2[j-1]
            diagonal_score = matrixBuild[i-1][j-1] + scoring_function(char1, char2)
            
            if current_score == diagonal_score:
                aligned_seq1 += char1
                aligned_seq2 += char2
                i -= 1
                j -= 1
                continue 
                
        # Check if from above? (Gap inserted in seq2)
        if i > 0:
            char1 = seq1[i-1]
            down_score = matrixBuild[i-1][j] + scoring_function(char1, "-")
            
            if current_score == down_score:
                aligned_seq1 += char1
                aligned_seq2 += "-"
                i -= 1
                continue
                
        # If it wasn't diagonal or above, must have come from left (Gap in seq1)
        char2 = seq2[j-1]
        aligned_seq1 += "-"
        aligned_seq2 += char2
        j -= 1
            
    # Because we traced backwards, the strings are backwards as well, reverse this.
    final_score = float(matrixBuild[len1][len2])
    
    return aligned_seq1[::-1], aligned_seq2[::-1], final_score


def local_alignment(seq1, seq2, scoring_function):
    """Local sequence alignment using the Smith-Waterman algorithm.

    Indels should be denoted with the "-" character.

    Parameters
    ----------
    seq1: str
        First sequence to be aligned.
    seq2: str
        Second sequence to be aligned.
    scoring_function: Callable

    Returns
    -------
    str
        First aligned sequence.
    str
        Second aligned sequence.
    float
        Final score of the alignment.

    Examples
    --------
    >>> local_alignment("pending itch", "unending glitch", lambda x, y: [-1, 1][x == y])
    ('ending --itch', 'ending glitch', 9.0)

    Other alignments are not possible.

    """
    
    # Initialise Matrix
    len1 = len(seq1)
    len2 = len(seq2)
    matrixBuild = []

    for i in range(len1 + 1):
        row = []
        for j in range(len2 + 1):
            row.append(0)
        matrixBuild.append(row)
        
    # track the highest score and its location
    best_score = 0
    start_i = 0
    start_j = 0
    
 
    for i in range(1, len1 + 1):
        for j in range(1, len2 + 1):
            char1 = seq1[i-1]
            char2 = seq2[j-1]
            
            # Calculate paths
            diagonal_score = matrixBuild[i-1][j-1] + scoring_function(char1, char2)
            down_score = matrixBuild[i-1][j] + scoring_function(char1, "-")
            right_score = matrixBuild[i][j-1] + scoring_function("-", char2)
            
            # In local alignment, a score cannot be negative. If all are negative, it resets to 0.
            matrixBuild[i][j] = max(0, diagonal_score, down_score, right_score)
            
            # Check and replace if highest
            if matrixBuild[i][j] > best_score:
                best_score = matrixBuild[i][j]
                start_i = i
                start_j = j
                
    aligned_seq1 = ""
    aligned_seq2 = ""
    
    # Start traceback from the highest score found in the matrix
    i = start_i
    j = start_j
    
    # Keep tracing back until we hit the edge OR we hit a cell with a score of 0
    while i > 0 and j > 0 and matrixBuild[i][j] > 0:
        current_score = matrixBuild[i][j]
        
        # Check if from diag? (Match/Mismatch)
        char1 = seq1[i-1]
        char2 = seq2[j-1]
        diagonal_score = matrixBuild[i-1][j-1] + scoring_function(char1, char2)
        
        if current_score == diagonal_score:
            aligned_seq1 += char1
            aligned_seq2 += char2
            i -= 1
            j -= 1
            continue 
            
        # Check if from above? (Gap inserted in seq2)
        down_score = matrixBuild[i-1][j] + scoring_function(char1, "-")
        
        if current_score == down_score:
            aligned_seq1 += char1
            aligned_seq2 += "-"
            i -= 1
            continue
            
        # If it wasn't diagonal or above, must have come from left (Gap in seq1)
        aligned_seq1 += "-"
        aligned_seq2 += char2
        j -= 1
            
    # Reverse the final strings just like in global alignment
    return aligned_seq1[::-1], aligned_seq2[::-1], float(best_score)


## This is an example scoring function, you should implement a version which uses a scoring matrix 
matrix = substitution_matrices.load("BLOSUM62")

def blosum62_scoring(char1, char2):
    gap_score = -4
    
    if char1 == "-" or char2 == "-":
        return gap_score
        
    try:
        score = matrix[char1, char2]
        return score
    except KeyError:
        pass
        
    try:
        score = matrix[char2, char1]
        return score
    except KeyError:
        return -1
