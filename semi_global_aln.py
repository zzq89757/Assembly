import numpy as np

def local_alignment(seq1, seq2, match_score=2, mismatch_score=-2, gap_penalty=-2):
    n = len(seq1)
    m = len(seq2)

    # 初始化得分矩阵和轨迹矩阵
    score_matrix = np.zeros((n + 1, m + 1))
    traceback_matrix = np.zeros((n + 1, m + 1), dtype=int)

    # 方向定义：0=停止，1=对角，2=向上，3=向左
    STOP, DIAG, UP, LEFT = 0, 1, 2, 3

    # 初始化最大得分和对应的位置
    max_score = 0
    max_pos = (0, 0)

    # 填充得分矩阵和轨迹矩阵
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match = score_matrix[i - 1][j - 1] + (match_score if seq1[i - 1] == seq2[j - 1] else mismatch_score)
            delete = score_matrix[i - 1][j] + gap_penalty
            insert = score_matrix[i][j - 1] + gap_penalty
            score = max(0, match, delete, insert)
            score_matrix[i][j] = score

            if score == 0:
                traceback_matrix[i][j] = STOP
            elif score == match:
                traceback_matrix[i][j] = DIAG
            elif score == delete:
                traceback_matrix[i][j] = UP
            elif score == insert:
                traceback_matrix[i][j] = LEFT

            if score > max_score:
                max_score = score
                max_pos = (i, j)

    # 回溯得到比对结果
    align1, align2 = '', ''
    i, j = max_pos
    # print(score_matrix)
    while traceback_matrix[i][j] != STOP:
        if traceback_matrix[i][j] == DIAG:
            align1 = seq1[i - 1] + align1
            align2 = seq2[j - 1] + align2
            i -= 1
            j -= 1
        elif traceback_matrix[i][j] == UP:
            align1 = seq1[i - 1] + align1
            align2 = '-' + align2
            i -= 1
        elif traceback_matrix[i][j] == LEFT:
            align1 = '-' + align1
            align2 = seq2[j - 1] + align2
            j -= 1

    return align1, align2, max_score

# 示例使用
seq1 = "ACCCTTACACAGACGTAGTAG"
seq2 = "TAGTCGCGGTGCTAGACA"
alignment1, alignment2, score = local_alignment(seq1, seq2)
print(f"Alignment 1: {alignment1}")
print(f"Alignment 2: {alignment2}")
print(f"Score: {score}")
