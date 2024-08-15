from typing import Callable, List
from collections import defaultdict
import numpy as np

import time
import tracemalloc
from functools import wraps

def measure_time_and_memory(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        # 开始测量时间和内存
        start_time = time.time()
        tracemalloc.start()

        # 执行函数
        result = func(*args, **kwargs)

        # 结束测量时间和内存
        end_time = time.time()
        current, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()

        # 打印耗费时间和内存
        print(f"Time taken: {end_time - start_time:.4f} seconds")
        print(f"Current memory usage: {current / 10**6:.4f} MB")
        print(f"Peak memory usage: {peak / 10**6:.4f} MB")
        
        return result

    return wrapper


def reverse_complement(seq: str) -> str:
    trantab = str.maketrans("ACGTNacgtnRYMKrymkVBHDvbhd", "TGCANtgcanYRKMyrkmBVDHbvdh")
    return seq.translate(trantab)[::-1]

import numpy as np

def aln_by_dp(seq1: str, seq2: str, is_global_aln: bool = True) -> None:
    # 定义方向常量
    UP = 1
    X = 4
    LEFT = 2
    # 罚分规则
    gap_penalty = -1
    mis_penalty = -1
    match_reward = 2
    # 矩阵大小
    col_num = len(seq1) + 1
    row_num = len(seq2) + 1
    # 初始化矩阵
    score_matrix = np.zeros((row_num, col_num), dtype=int)
    trace_matrix = np.zeros((row_num, col_num), dtype=int)
    
    # 构建打分和方向矩阵
    for i in range(row_num):
        for j in range(col_num):
            if i == j == 0:
                continue
            # 初始化首行首列（全局比对）
            if i == 0:
                score_matrix[i, j] = score_matrix[i, j - 1] + gap_penalty if is_global_aln else 0
                trace_matrix[i, j] = LEFT
                continue
            if j == 0:
                score_matrix[i, j] = score_matrix[i - 1, j] + gap_penalty if is_global_aln else 0
                trace_matrix[i, j] = UP
                continue
            # 不同方向来源的最终得分
            score_from_up_left = score_matrix[i-1, j-1] + (match_reward if seq1[j - 1] == seq2[i - 1] else mis_penalty)
            score_from_up = score_matrix[i - 1, j] + gap_penalty
            score_from_left = score_matrix[i, j - 1] + gap_penalty
            # 获取得分最高的来源存入轨迹矩阵 最高分存入打分矩阵
            scores = [score_from_up, score_from_up_left, score_from_left]
            best_score = max(scores)
            
            # 对于局部比对，最小值设为0
            if not is_global_aln:
                best_score = max(0, best_score)
            
            score_matrix[i, j] = best_score
            trace_matrix[i, j] = (UP if scores[0] == best_score else 0) + (X if scores[1] == best_score else 0) + (LEFT if scores[2] == best_score else 0)
    
    # 回溯并生成比对结果
    if is_global_aln:
        x_pos = col_num - 1
        y_pos = row_num - 1
    else:
        # 对于局部比对，找到最大值的位置开始回溯
        y_pos, x_pos = np.unravel_index(np.argmax(score_matrix), score_matrix.shape)

    aln1 = ""
    aln2 = ""
    while (x_pos != 0 or y_pos != 0) and (is_global_aln or score_matrix[y_pos, x_pos] > 0):
        path_score = trace_matrix[y_pos, x_pos]
        if path_score == UP:
            aln1 += "-"
            aln2 += seq2[y_pos - 1]
            y_pos -= 1
        elif path_score == LEFT:
            aln1 += seq1[x_pos - 1]
            aln2 += "-"
            x_pos -= 1
        elif path_score == X:
            aln1 += seq1[x_pos - 1]
            aln2 += seq2[y_pos - 1]
            x_pos -= 1
            y_pos -= 1
        elif path_score == UP + LEFT:
            if score_matrix[y_pos - 1, x_pos] > score_matrix[y_pos, x_pos - 1]:
                aln1 += "-"
                aln2 += seq2[y_pos - 1]
                y_pos -= 1
            else:
                aln1 += seq1[x_pos - 1]
                aln2 += "-"
                x_pos -= 1
        elif path_score == LEFT + X:
            if score_matrix[y_pos - 1, x_pos - 1] > score_matrix[y_pos, x_pos - 1]:
                aln1 += seq1[x_pos - 1]
                aln2 += seq2[y_pos - 1]
                x_pos -= 1
                y_pos -= 1
            else:
                aln1 += seq1[x_pos - 1]
                aln2 += "-"
                x_pos -= 1
        elif path_score == UP + X:
            if score_matrix[y_pos - 1, x_pos - 1] > score_matrix[y_pos - 1, x_pos]:
                aln1 += seq1[x_pos - 1]
                aln2 += seq2[y_pos - 1]
                x_pos -= 1
                y_pos -= 1
            else:
                aln1 += "-"
                aln2 += seq2[y_pos - 1]
                y_pos -= 1
        elif path_score == UP + LEFT + X:
            if score_matrix[y_pos - 1, x_pos - 1] > score_matrix[y_pos - 1, x_pos] and score_matrix[y_pos - 1, x_pos - 1] > score_matrix[y_pos, x_pos - 1]:
                aln1 += seq1[x_pos - 1]
                aln2 += seq2[y_pos - 1]
                x_pos -= 1
                y_pos -= 1
            elif score_matrix[y_pos - 1, x_pos] > score_matrix[y_pos - 1, x_pos - 1] and score_matrix[y_pos - 1, x_pos] > score_matrix[y_pos, x_pos - 1]:
                aln1 += "-"
                aln2 += seq2[y_pos - 1]
                y_pos -= 1
            else:
                aln1 += seq1[x_pos - 1]
                aln2 += "-"
                x_pos -= 1
            
    print("Alignment1 ", end="")
    print("".join(aln1[::-1]))
    print("Alignment2 ", end="")
    print("".join(aln2[::-1]))




def construct_dict(seq_li:List[str]) -> defaultdict:
    seq_dict = {i + 1:v for i,v in enumerate(seq_li)}
    seq_dict.update({-i:reverse_complement(v) for i,v in seq_dict.items()})
    return seq_dict

    

def is_valid_aln(seq1:str, seq2:str) -> None:
    '''
    计算两序列比对结果 overlap长度和identity是否达到阈值
    '''
    print(f"{seq1}->{seq2}")
    aln_by_dp(seq1, seq2, False)


# 序列及其反向互补两两比对 计算重叠区identity
# 使用装饰器
@measure_time_and_memory
def assembly_by_olc(seq_li:List[str]) -> str:
    # 构建序列字典
    seq_dict = construct_dict(seq_li)
    n = len(seq_li)
    # 序列索引为1~n 反向互补为-1~-n 
    index_record = []  # 两两比对时 计算索引乘积是否已存在 存在则跳过
    for i in range(n, -n - 1, -1):
        if not i:continue # skip while index eq 0
        for j in range(i - 1, -n - 1, -1):
            if not j:continue # skip while index eq 0
            if i + j == 0:continue # skip while loop to itself
            if i * j in index_record and i * j > 0:continue # skip while index exists in record
            index_record.append(i * j) # append index to record
            # 两两比对并计算重叠区identity 记录拼接中间结果
            is_valid_aln(seq_dict[i], seq_dict[j])
            
        
if __name__ == "__main__":
    seq_li = ['CGTAGCT','GCTGTAGT', 'ACAGTGTAC']
    assembly_by_olc(seq_li)
    # construct_dict(['CGTAGCT','GCTGTAGT', 'ACAGTGTAC'])