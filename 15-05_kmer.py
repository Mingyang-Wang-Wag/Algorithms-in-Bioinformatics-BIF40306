from sys import argv

def kmer_library():
    lib = []
    alphabet = 'ACGT'

    for first in alphabet:
        for second in alphabet:
            for third in alphabet:
                for fourth in alphabet:
                    kmer = first + second + third + fourth
                    lib.append(kmer)

    return lib

def reader(filename):
    dna_seq = ''
    with open(filename) as f:
        for line in f:
            line = line.strip()

            if line == '':
                continue

            if line.startswith('>'):
                continue

            dna_seq += line

    return dna_seq

def parser(dna_seq):
    kmer_dict = {}
    for i in range(len(dna_seq)):
        if i+4 <= len(dna_seq):
            kmer = dna_seq[i:i+4]
            if kmer not in kmer_dict:
                kmer_dict[kmer] = 0
            else:
                kmer_dict[kmer] += 1

    return kmer_dict

def output(lib, kmer_dict):
    kmer_order_dict = {}
    for kmer in lib:
        if kmer in kmer_dict:
            freq = kmer_dict[kmer]
            kmer_order_dict[kmer] = freq
        else:
            kmer_order_dict[kmer] = 0

    freq_ls = []
    for value in kmer_order_dict.values():
        freq_ls.append(value)

    print(freq_ls)
    return freq_ls

def main():
    filename = argv[1]
    lib = kmer_library()
    dna_seq = reader(filename)
    kmer_dict = parser(dna_seq)
    output(lib, kmer_dict)

if __name__ == "__main__":
    main()