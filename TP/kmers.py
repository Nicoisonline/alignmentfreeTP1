
def kmer2str(val, k):
    """ Transform a kmer integer into a its string representation
    :param int val: An integer representation of a kmer
    :param int k: The number of nucleotides involved into the kmer.
    :return str: The kmer string formatted
    """
    letters = ['A', 'C', 'T', 'G']
    str_val = []
    for _ in range(k):
        str_val.append(letters[val & 0b11])
        val >>= 2

    str_val.reverse()
    return "".join(str_val)

def encode_kmer(text, k):
    kmer = 0
    kmer_inv = 0
    for letter in text[0:k]:
        kmer_inv<<=2 
        kmer<<=2
        kmer = kmer + encode_nucl(letter)
        kmer_inv = kmer_inv + rev_nuc(encode_nucl(letter))
    return kmer, kmer_inv

def stream_kmers(text, k):
    """enumerate_kmer"""
    text = text[0]
    mask = (1<<(2*(k-1)))-1
    kmer, kmer_inv = encode_kmer(text,k)
    for i in range(len(text)-(k)):
        yield min(kmer, kmer_inv)
        kmer &= mask
        kmer_inv &= mask
        kmer <<= 2
        kmer_inv <<= 2
        kmer = kmer + encode_nucl(text[i+k])
        kmer_inv = kmer_inv + rev_nuc(encode_nucl(text[i+k]))
    yield min(kmer, kmer_inv)

def encode_nucl(letter):
    encoding = {"A" : 0, "C" : 1, "T" : 2, "G" : 3}
    return encoding[letter]

def rev_nuc(n_int) :
    return (n_int + 2) & 3