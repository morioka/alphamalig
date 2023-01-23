#
# alphabet matrixを作成する
#

def print_alphabet_matrix(alphabets, costs):

    # alphabetが印刷可能でない文字を含むか
    printable_alphabets =  list(range(ord('A'), ord('Z')+1))+ list(range(ord('a'), ord('z')+1))+ list(range(ord('0'), ord('9')+1))+ [ord('-')]
    printable_alphabets = [chr(c) for c in printable_alphabets]

    assert (set(map(type, alphabets)) == set([str, int])) or (set(map(type, alphabets)) == set([str])) or (set(map(type, alphabets)) == set([int])) 
    assert len(set(alphabets)) == len(alphabets)
    assert alphabets[-1] == '-' # Gap

    alphabets_ = [chr(a) if type(a) == int else a for a in alphabets]

    # 以下の文字は含まないこと
    assert chr(0x00) not in alphabets_  # NUL (0x00)
    assert chr(0x3e) not in alphabets_  # > (0x3e)
    assert chr(0x3d) not in alphabets_  # = (0x3d)
    assert chr(0x3c) not in alphabets_  # < (0x3c)
    assert chr(0x20) not in alphabets_  # Space (0x20)
    assert chr(0x0d) not in alphabets_  # Carriage Return (0x0d) 
    assert chr(0x0a) not in alphabets_  # Line Feed (0x0a)

    # Gap文字は必ず含むこと
    assert chr(0x2d) in alphabets_  # - (0x2d)

    # 以下 出力
    print(len(alphabets_))

    if sum([c in printable_alphabets for c in alphabets_ ]) != len(alphabets_):
        print(" ".join([f'{ord(c):x}' for c in alphabets_ ]))
    else:
        print(" ".join(alphabets_))

    # match はアルファベットによらず固定値。
    # mismatch はアルファベットによらず固定値。かつ対称
    # gap_penalty は前後のアルファベットや継続長によらず固定値。
    for i in range(len(alphabets_)):
        c = alphabets_[i]
        s = alphabets_[i:]

        if i == len(alphabets_) - 1:
            print(" ".join(map(str, [costs['gap_penalty']] * len(alphabets_))))
        else:
            print(" ".join(map(str, [costs['mismatch']] * (i) + [costs['match']])))


if __name__ == '__test__':

    result_gold ="""6 
    o p s c n -
    2 
    -1 15
    -2 -2 1
    -2 -2 0 1
    -2 -2 -1 -1 1
    -2 -2 0 0 0 0
    """

    costs = {
        'match': 100.0,     # match はアルファベットによらず固定値。
        'mismatch': -10.0,  # mismatch はアルファベットによらず固定値。かつ対称
        'gap_penalty': 0.0, # gap_penalty は前後のアルファベットや継続長によらず固定値。
    }

    alphabets =  [ 'a', 'b', 'c', '-']

    print_alphabet_matrix(alphabets, costs)

    print('#' * 32)

    alphabets =  list(range(ord('A'), ord('Z')+1))+ list(range(ord('a'), ord('z')+1))+ list(range(ord('0'), ord('9')+1))+ [ord('-')]
    alphabets =  [chr(c) for c in alphabets]

    print_alphabet_matrix(alphabets, costs)

    print('#' * 32)

    alphabets =  list(range(ord('A'), ord('Z')+1)) +  [ord('-')]
    alphabets =  [chr(c) for c in alphabets]
    alphabets =  [0x1a, 0x1b] + alphabets

    print_alphabet_matrix(alphabets, costs)

    print('#' * 32)

if __name__ == '__main__':

    costs = {
        'match': 100.0,     # match はアルファベットによらず固定値。
        'mismatch': -10.0,  # mismatch はアルファベットによらず固定値。かつ対称
        'gap_penalty': 0.0, # gap_penalty は前後のアルファベットや継続長によらず固定値。
    }

    alphabets =  list(range(ord('A'), ord('Z')+1))+ list(range(ord('a'), ord('z')+1))+ list(range(ord('0'), ord('9')+1))+ [ord('-')]
    alphabets =  [chr(c) for c in alphabets]

    print_alphabet_matrix(alphabets, costs)

