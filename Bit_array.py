from bitarray import bitarray
from numpy import ndarray
import numpy as np


def bit_li_converted_by_char(ch):
    # Convert a char T,C,A,G into a two bit num, respectively
    ret:ndarray
    if ch == 'T' or ch == 'U':
        ret = np.array([0,1])
    elif ch == 'C':
        ret = np.array([1,0])
    elif ch == 'A':
        ret = np.array([0,0])
    elif ch == 'G':
        ret = np.array([1,1])
    return ret


class BitArray:
    def __init__(self, seq: str):
        # 初始化 BitArray，bit_str 是二进制字符串
        self.seq = seq
        self.bit_li = []
        self.to_bit_li()
        self.bits = bitarray(self.bit_li.tolist())
        # self.to_base4()
        # self.to_seq()
        self.to_base10()

    def to_bit_li(self) -> ndarray[int]:
        self.bit_li:ndarray[int] = np.array([],int)
        for base in self.seq:
            self.bit_li = np.append(self.bit_li, bit_li_converted_by_char(base))
    
    def append(self, base:str):
        # 追加新的两位二进制数到当前 bitarray
        bits = bit_li_converted_by_char(base)
        # bit li append
        self.bit_li = np.append(self.bit_li, bits)
        self.bits.extend(bits.tolist())

    def remove(self):
        # 去掉前两位
        self.bits = self.bits[2:]
        self.bit_li = self.bit_li[2:]
        self.seq = self.seq[1:]
        self.to_base10()
    
    # def to_base4(self) -> str:
    #     # 将 bitarray 每两位转换为四进制
    #     self.base4_str = ''
    #     # 确保 bitarray 长度为2的倍数
    #     if len(self.bits) % 2 != 0:
    #         self.bits = bitarray('0') + self.bits
        
    #     # 遍历每两位并转换为四进制
    #     for i in range(0, len(self.bits), 2):
    #         # 每两位二进制数转为四进制
    #         binary_pair = self.bits[i:i+2]
    #         base4_digit = int(binary_pair.to01(), 2)
    #         self.base4_str += str(base4_digit)
    
    # def to_seq(self):
    #     trantab = str.maketrans("0123","ATCG")
    #     self.seq =  self.base4_str.translate(trantab)
        
    def to_base10(self):
        bits_str = self.bits.to01()
        self.base10 = int(bits_str,2)
        # print(self.base10)
            

seq = "GTGATCGCCACA"
# print(BitArray(seq).bits)


# 示例使用
_bitarray = BitArray(seq)  # 初始化为 '1101'

# 追加两位二进制数（使用整数表示）
_bitarray.append("A")  # 2 -> '10' 结果为 '110110'
_bitarray.remove()  # 2 -> '10' 结果为 '110110'


