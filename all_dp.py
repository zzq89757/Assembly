from typing import Callable, List



# 
    # 序列及其反向互补两两比对 计算重叠区identity

def assembly_by_olc(seq_li:List[str]) -> str:
    # 序列索引为1~n 反向互补为-1~-n 
    index_record = []  # 两两比对时 计算索引乘积是否已存在 存在则跳过
    # 