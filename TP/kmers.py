
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
    """ input:  text: 'str' sequence
                k: 'int' taille de kmer
        output: kmer: 'int' encodage du kmer
        encode les kmer de str en int en utilisant les bits"""
    kmer = 0
    for letter in text[0:k]:
        kmer<<=2
        kmer = kmer + encode_nucl(letter)
    return kmer

def stream_kmers(text, k):
    """ input:  text: 'list' liste de sequences
                k: 'int' taille de kmer
        output: Generator[int] kmer produits
        genere des kmer de taille k a partir de la sequence text"""
    text = text[0]
    mask = (1<<(2*(k-1)))-1
    kmer = encode_kmer(text,k)
    kmer_inv = encode_kmer_rev(kmer,k)
    for i in range(len(text)-(k)):
        yield min(kmer, kmer_inv)
        kmer &= mask
        kmer <<= 2
        kmer_inv >>= 2
        kmer = kmer + encode_nucl(text[i+k])
        kmer_inv = (encode_nucl(rev_nuc(text[i+k]))<<(2*(k-1))) + kmer_inv
    yield min(kmer, kmer_inv)

def encode_nucl(letter):
    """ input:   letter: 'str' nucleotide
        output:  'int' codage du nucleotide
        encode le nucleotide en int"""
    return {"A" : 0, "C" : 1, "T" : 2, "G" : 3}[letter]


def rev_nuc(n) :
    """ input:   n: 'str' nucleotide
        output:  'str' complementaire du nucleotide
        encode le complementaire du nucleotide"""
    return {"A" : "T", "G": "C", "T": "A" , "C":"G"}[n]

def encode_kmer_rev(kmer,k) : 
    """ input:  kmer: 'int' kmer encode
                k: 'int' taille du kmer
        output: kmer_inv: 'int' reverse complementaire encode
        encode le reverse complementaire d'un kmer """
    kmer_inv = 0
    text = kmer2str(kmer,k)
    text_rev = text[::-1]
    for letter in text_rev:
        kmer_inv<<=2
        kmer_inv = kmer_inv + encode_nucl(rev_nuc(letter))
    return kmer_inv