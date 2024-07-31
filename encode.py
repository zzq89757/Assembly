class BitSet:
    def __init__(self, size = 0):
        self.size = size
        self.bitset = 0

    def set_bit(self, pos):
        if pos >= self.size or pos < 0:
            raise ValueError("Position out of range")
        self.bitset |= (1 << pos)

    def clear_bit(self, pos):
        if pos >= self.size or pos < 0:
            raise ValueError("Position out of range")
        self.bitset &= ~(1 << pos)

    def check_bit(self, pos):
        if pos >= self.size or pos < 0:
            raise ValueError("Position out of range")
        return bool(self.bitset & (1 << pos))

    def toggle_bit(self, pos):
        if pos >= self.size or pos < 0:
            raise ValueError("Position out of range")
        self.bitset ^= (1 << pos)

    def get_binary(self):
        return bin(self.bitset)[2:]
    
    def get_base10(self):
        return self.bitset

    def add_decimal_as_binary(self, number):
        binary_str = bin(number)[2:]  # Convert the number to binary string without '0b' prefix
        for digit in binary_str:
            self.bitset <<= 1  # Shift all bits to the left by 1 position
            self.size += 1
            if digit == '1':
                self.bitset |= 1  # Set the last bit to 1

    def __repr__(self):
        return f'BitSet({self.size}): ' + self.get_binary()




def base2int(base:str):
    tran_table = str.maketrans("ATCG",'0123')
    return int(base.upper().translate(tran_table))

def seq2bin(seq:str):
    c = BitSet()
    for base in seq:
        c.add_decimal_as_binary(base2int(base))
    return c.get_binary()






if __name__ == "__main__":
    c = seq2bin("CC")
    print(c)
    # print(bin(33))
    c = BitSet()
    c.add_decimal_as_binary(3)
    c.add_decimal_as_binary(3)
    print(c.get_binary())
    # print(bin(33)