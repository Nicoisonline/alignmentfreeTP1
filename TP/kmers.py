
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


def stream_kmers(text, k):
    """Aussi appelé enumerate_kmers"""
    # --- To complete ---
    mask = (1 << (2(k-1))) - 1
    kmer = text[0:k-1]
    for i in range(len(text) - (k - 1)):
        yield kmer
        kmer &= mask
        kmer <<= 2
        kmer += encode_nucl(text[i + k])

def encode_nucl(nucl):
    """ Encode a nucleotide into a 2-bits integer"""
    return {'A': 0, 'C': 1, 'T': 2, 'G': 3}[nucl]
