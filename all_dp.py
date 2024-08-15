from typing import Callable, List



# 
    # 序列及其反向互补两两比对 计算重叠区identity

def assembly_by_olc(seq_li:List[str]) -> str:
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
            # 比对并计算重叠区identity
            print(f"{i}*{j}")
        
if __name__ == "__main__":
    assembly_by_olc(['sds','sd', 'asds'])