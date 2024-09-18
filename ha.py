import numpy as np

def heuristic_alignment(seq1, seq2, match_score=1, mismatch_penalty=-1, gap_penalty=-1, mode='global'):
    # 初始化得分矩阵
    len1, len2 = len(seq1), len(seq2)
    dp = np.zeros((len1 + 1, len2 + 1))
    
    # 记录局部比对最高分及其位置（仅在局部比对模式下使用）
    max_score = 0
    max_pos = (0, 0)
    
    # 初始化边界条件（局部比对中边界条件为0）
    if mode == 'global':
        for i in range(len1 + 1):
            dp[i, 0] = i * gap_penalty
        for j in range(len2 + 1):
            dp[0, j] = j * gap_penalty

    # 填充得分矩阵
    for i in range(1, len1 + 1):
        for j in range(1, len2 + 1):
            match = dp[i-1, j-1] + (match_score if seq1[i-1] == seq2[j-1] else mismatch_penalty)
            delete = dp[i-1, j] + gap_penalty
            insert = dp[i, j-1] + gap_penalty
            
            if mode == 'global':
                dp[i, j] = max(match, delete, insert)
            else:  # 局部比对时，不允许负分数，最小得分为0
                dp[i, j] = max(0, match, delete, insert)
                if dp[i, j] > max_score:
                    max_score = dp[i, j]
                    max_pos = (i, j)

    # 回溯获得比对结果
    alignment1, alignment2 = "", ""
    
    if mode == 'global':
        i, j = len1, len2
    else:  # 局部比对从得分最高的位置开始回溯
        i, j = max_pos
    
    while (i > 0 or j > 0) and (mode == 'global' or dp[i, j] != 0):
        current_score = dp[i, j]
        if i > 0 and j > 0 and dp[i-1, j-1] + (match_score if seq1[i-1] == seq2[j-1] else mismatch_penalty) == current_score:
            alignment1 = seq1[i-1] + alignment1
            alignment2 = seq2[j-1] + alignment2
            i -= 1
            j -= 1
        elif i > 0 and dp[i-1, j] + gap_penalty == current_score:
            alignment1 = seq1[i-1] + alignment1
            alignment2 = "-" + alignment2
            i -= 1
        else:
            alignment1 = "-" + alignment1
            alignment2 = seq2[j-1] + alignment2
            j -= 1

    # 对局部比对的情况，可能存在不回溯到矩阵的左上角
    if mode == 'global':
        return alignment1, alignment2, dp[len1, len2]
    else:
        return alignment1, alignment2, max_score

# 示例
seq1 = 'GCCAAGTGCCCAGCGGGGCTGCTAAAGCGCATGCTCCAGACTGCCTTGGGAAAAGCGCCTCCCCTACCCGGTAGAATTGGATCCCCAAGGTCGGGCAGGAAGAGGGCCTATTTCCCATGATTCCTTCATATTTGCATATACGATACAAGGCTGTTAGAGAGATAATTAGAATTAATTTGACTGTAAACACAAAGATATTAGTACAAAATACGTGACGTAGAAAGTAATAATTTCTTGGGTAGTTTGCAGTTTTAAAATTATGTTTTAAAATGGACTATCATATGCTTACCGTAACTTGAAAGTATTTCGATTTCTTGGCTTTATATATCTTGTGGAAAGGACGAAACACCGTTCAACCCGGAGGCGCCCGAGTTTAAGAGCTAAGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTTGTTGACAATTAATCATCGGCATAGTATATCGGCATAGTATAATACGACAAGGTGAGGAACTAAACCATGGGTCAAAGTAGCGATGAAGCCAACGCTCCCGTTGCAGGGCAGTTTGCGCTTCCCCTGAGTGCCACCTTTGGCTTAGGGGATCGCGTACGCAAGAAATCTGGTGCCGCTTGGCAGGGTCAAGTCGTCGGTTGGTATTGCACAAAACTCACTCCTGAAGGCTATGCGGTCGAGTCCGAATCCCACCCAGGCTCAGTGCAAATTTATCCTGTGGCTGCAC'
seq2 = 'CTGACTAGGGGAGGAGTAGAAGGTGGCGCGAAGGGGCCACCAAAGAACGGAGCCGGTTGGCGCCTACCGGTGGATGTGGAATGTGTGCGAGGCCAGAGGCCACTTGTGTAGCGCCAAGTGCCCAGCGGGGCTGCTAAAGCGCATGCTCCAGACTGCCTTGGGAAAAGCGCCTCCCCTACCCGGTAGAATTGGATCCCCAAGGTCGGGCAGGAAGAGGGCCTATTTCCCATGATTCCTTCATATTTGCATATACGATACAAGGCTGTTAGAGAGATAATTAGAATTAATTTGACTGTAAACACAAAGATATTAGTACAAAATACGTGACGTAGAAAGTAATAATTTCTTGGGTAGTTTGCAGTTTTAAAATTATGTTTTAAAATGGACTATCATATGCTTACCGTAACTTGAAAGTATTTCGATTTCTTGGCTTTATATATCTTGTGGAAAGGACGAAACACCGtTCAACCCGGAGGCGCCCGAGTTTAAGAGCTAAGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTTGTTGACAATTAATCATCGGCATAGTATATCGGCATAGTATAATACGACAAGGTGAGGAACTAAACCATGGGTCAAAGTAGCGATGAAGCCAACGCTCCCGTTGCAGGGCAGTTTGCGCTTCCCCTGAGTGCCACCTTTGGCTTAGGGGATCGCGTACGCAAGAAATCTGGTGCCGCTTGGCAGGGTCAAGTCGTCGGTTGGTATTGCACAAAACTCACTCCTGAAGGCTATGCGGTCGAGTCCGAATCCCACCCAGGCTCAGTGCAAATTTATCCTGTGGCTGCACTTGAACGTGTGGCCTAA'
 

# 全局比对
# alignment1, alignment2, score = heuristic_alignment(seq1, seq2, mode='global')
# print("Global Alignment:")
# print("Alignment 1:", alignment1)
# print("Alignment 2:", alignment2)
# print("Score:", score)

# 局部比对
alignment1, alignment2, score = heuristic_alignment(seq1, seq2, mode='local')
print("\nLocal Alignment:")
print("Alignment 1:", alignment1)
print("Alignment 2:", alignment2)
print("Score:", score)
