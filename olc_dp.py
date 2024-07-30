


def read_sequences(file_path):
    with open(file_path, 'r') as file:
        sequences = [line.strip() for line in file.readlines() if line.strip()]
    return sequences

def overlap(s1, s2, min_length=3, max_mismatches=0, max_indels=0):
    len1, len2 = len(s1), len(s2)
    max_overlap = 0
    best_i, best_j = -1, -1
    
    # Initialize DP table
    dp = [[0] * (len2 + 1) for _ in range(len1 + 1)]
    
    for i in range(1, len1 + 1):
        for j in range(1, len2 + 1):
            if s1[i - 1] == s2[j - 1]:
                dp[i][j] = dp[i - 1][j - 1] + 1
            else:
                dp[i][j] = max(dp[i - 1][j] - 1, dp[i][j - 1] - 1, dp[i - 1][j - 1] - 1)
            
            if dp[i][j] > max_overlap and dp[i][j] >= min_length:
                mismatches = (i + j - 2 * dp[i][j])
                if mismatches <= max_mismatches:
                    best_i, best_j = i, j
                    max_overlap = dp[i][j]
    
    if max_overlap >= min_length:
        return max_overlap, best_i, best_j
    else:
        return 0, -1, -1

def build_overlap_graph(sequences, min_length=3, max_mismatches=0, max_indels=0):
    graph = {}
    for i, seq1 in enumerate(sequences):
        for j, seq2 in enumerate(sequences):
            if i != j:
                olen, best_i, best_j = overlap(seq1, seq2, min_length, max_mismatches, max_indels)
                if olen > 0:
                    graph[(seq1, seq2)] = (olen, best_i, best_j)
    return graph

def assemble_sequence(graph):
    seqs = list(set([seq for edge in graph.keys() for seq in edge]))
    while len(seqs) > 1:
        max_overlap = 0
        best_pair = None
        best_i, best_j = -1, -1
        for (seq1, seq2), (olen, i, j) in graph.items():
            if olen > max_overlap:
                max_overlap = olen
                best_pair = (seq1, seq2)
                best_i, best_j = i, j
        if best_pair:
            seq1, seq2 = best_pair
            combined_seq = seq1 + seq2[best_j:]
            seqs.remove(seq1)
            seqs.remove(seq2)
            seqs.append(combined_seq)
            graph = {k: v for k, v in graph.items() if seq1 not in k and seq2 not in k}
            for seq in seqs:
                if seq != combined_seq:
                    olen1, i1, j1 = overlap(seq, combined_seq)
                    if olen1 > 0:
                        graph[(seq, combined_seq)] = (olen1, i1, j1)
                    olen2, i2, j2 = overlap(combined_seq, seq)
                    if olen2 > 0:
                        graph[(combined_seq, seq)] = (olen2, i2, j2)
        else:
            break
    return seqs[0] if seqs else ""

# 示例使用
file_path = 'sequences.txt'  # 包含DNA序列的文件路径，每行一个序列
sequences = read_sequences(file_path)
print(sequences)
min_overlap_length = 3  # 设置最小重叠长度
max_mismatches = 2  # 设置最大允许错配数
max_indels = 1  # 设置最大允许的插入/缺失数
overlap_graph = build_overlap_graph(sequences, min_overlap_length, max_mismatches, max_indels)
assembled_sequence = assemble_sequence(overlap_graph)
print("Assembled Sequence:", assembled_sequence)
