import numpy as np

def bit_li_converted_by_char(ch):
    # Convert a char T,C,A,G into a two bit num, respectively
    ret:np.ndarray
    if ch == 'T' or ch == 'U':
        ret = np.array([0,1])
    elif ch == 'C':
        ret = np.array([1,0])
    elif ch == 'A':
        ret = np.array([0,0])
    elif ch == 'G':
        ret = np.array([1,1])
    return ret

def array_converted_by_seq(seq) -> np.ndarray:
    bit_li:np.ndarray[int] = np.array([],int)
    for base in seq:
        bit_li = np.append(bit_li, bit_li_converted_by_char(base))
    return bit_li



if __name__ == "__main__":
    seq = "CGATGCTAGTCG"
    a = array_converted_by_seq(seq)
    for i in range(0,len(a),2):
        print(i)