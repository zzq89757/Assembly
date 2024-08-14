from numpy import array, ndarray
from collections import Counter, defaultdict
from typing import Optional
from bitarray import bitarray

def bitarray_to_decimal(bits: bitarray) -> int:
    # 将 bitarray 转换为二进制字符串
    bin_str = bits.to01()  # 'to01' 方法将 bitarray 转换为 '0' 和 '1' 字符串
    # 将二进制字符串转换为十进制整数
    return int(bin_str, 2)
a = bitarray()
a.append(0b11)
print(a)
exit()
    
def twobit_converted_by_char(ch):
    # Convert a char T,C,A,G into a two bit num, respectively
    ret = -1
    if ch == 'T' or ch == 'U':
        ret = 0b01
    elif ch == 'C':
        ret = 0b10
    elif ch == 'A':
        ret = 0b00
    elif ch == 'G':
        ret = 0b11
    return ret


def char_converted_by_base4(i):
    # Convert a digit 0,1,2,3 into a char T,C,A,G, respectively
    ch = '-'
    if i == 0:
        ch = 'T'
    elif i == 1:
        ch = 'C'
    elif i == 2:
        ch = 'A'
    elif i == 3:
        ch = 'G'
    return ch


def kmer_count(seq_li:ndarray[str], k:int) -> Optional[defaultdict]:
  '''
  将输入的序列及其反向互补分别切kmer并统计数目
  '''
  print(seq_li)
  print(k)


def main() -> None:
  # 根据设置的重叠区最小长度 确定kmer长度
  kmer_count(seq_li, k)


if __name__ == "__main__":
  pass