from typing import Callable, List
from collections import defaultdict

def reverse_complement(seq: str) -> str:
    trantab = str.maketrans("ACGTNacgtnRYMKrymkVBHDvbhd", "TGCANtgcanYRKMyrkmBVDHbvdh")
    return seq.translate(trantab)[::-1]


def construct_dict(seq_li:List[str]) -> defaultdict:
    seq_dict = {i + 1:v for i,v in enumerate(seq_li)}
    seq_dict.update({-i:reverse_complement(v) for i,v in seq_dict.items()})
    return seq_dict
    

def is_valid_aln(seq1:str, seq2:str) -> None:
    print(f"{seq1}->{seq2}")


# 序列及其反向互补两两比对 计算重叠区identity

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