def hamming_distance(s, t):
    if len(s) != len(t):
        raise ValueError(f'Two strings should be in equal length')


    if len(s) > 1000:
        raise ValueError(f'Length should not exceed 1kbp')


    if s == '' or t == '':
        raise ValueError(f'String should not be empty')

    index = 0
    s = s.upper()
    t = t.upper()
    for i in range(len(s)):
        if s[i] != t[i]:
            index += 1

    return index


