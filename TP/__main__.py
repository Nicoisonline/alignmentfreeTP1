from TP.loading import load_directory
from TP.kmers import stream_kmers, kmer2str, min_hash


def compare_kmers(km1,km2):
    """ input:  km1 : 'str' premier kmer
                km2 : 'str' deuxieme kmer 
        output : bool 
        renvoie vrai ou faux dans la comparaison entre les deux kmer"""
    for i in range(len(km1)):
        if km1[i]!=km2[i]:
            return False
    return True

def create_index(file, k):
    """ input:  file : 'str' sequence
                k : 'int' taille du kmer 
        output : index : 'dict' dictionnaire des kmer et nombre d'apparition
        crée un dico contenant les kmer de taille k de la sequence et le nombre d'apparition"""
    index = {}
    for idx,km in enumerate(stream_kmers(file,k)):
        kmer = kmer2str(km,k)
        if kmer in index:
            index[kmer] += 1
        else :
            index[kmer] = 1
    return index

'''
def jaccard(fileA, fileB, k):
    """ input:  fileA : 'str' sequence 1
                fileB : 'str' sequence 2
                k : 'int' taille du kmer
        output: jac: 'int' valeur jaccard
        calcule la valeur jaccard entre deux séquences"""
    jac = 0
    index = create_index(fileA, k)
    intersect = 0
    union = sum(index.values())
    for idx,kme in enumerate(stream_kmers(fileB, k)):
        kmer2 = kmer2str(kme,k)
        if (kmer2 in index) and (index[kmer2]>0):
            intersect += 1
            index[kmer2] -= 1
        else:
            union += 1
    if union == 0 :
        return 0
    jac = intersect/union
    return jac
'''
def jaccard(listefA, listefB, k, s):
    """ input:  fileA : 'str' sequence 1
                fileB : 'str' sequence 2
                k : 'int' taille du kmer
        output: jac: 'int' valeur jaccard
        calcule la valeur jaccard entre deux séquences"""
    jac = 0
    intersect = 0
    union = 0
    cptA, cptB = 0, 0
    while ((cptA < s) and (cptB < s)):
        if listefA[cptA] == listefB[cptB]:
            intersect += 1
            union += 1
            cptA += 1
            cptB += 1
        elif listefA[cptA] > listefB[cptB]:
            union += 1
            cptB += 1
        elif listefA[cptA] < listefB[cptB]:
            union += 1
            cptA += 1
    union += len(listefA[cptA:-1]) + len(listefB[cptB:-1])
    jac = intersect/union
    return jac


if __name__ == "__main__":
    print("Computation of Jaccard similarity between files")

    # Load all the files in a dictionary
    files = load_directory("data_all")
    k = 21

    print("Computing Jaccard similarity for all pairs of samples")
    filenames = list(files.keys())

    for i in range(len(files)):
        listefA = min_hash(files[filenames[i]],k,1000)
        for j in range(i+1, len(files)):
            p = files[filenames[j]]
            if len(p) > 1:
                for l in range (1,len(p),1):
                    m = [p[0]]
                    m = [m[0] + p[l]]
                p = m
            listefB = min_hash(p,k,1000)
            jac = jaccard(listefA, listefB, k, 1000)
            print(filenames[i], filenames[j], jac)