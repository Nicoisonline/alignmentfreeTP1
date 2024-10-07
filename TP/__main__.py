from TP.loading import load_directory
from TP.kmers import stream_kmers, kmer2str


def compare_kmers(km1,km2):
    for i in range(len(km1)):
        if km1[i]!=km2[i]:
            return False
    return True

def create_index(file, k):
    index = {}
    for idx,km in enumerate(stream_kmers(file,k)):
        kmer = kmer2str(km,k)
        if kmer in index:
            index[kmer] += 1
        else :
            index[kmer] = 1
    return index

def jaccard(fileA, fileB, k):
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


if __name__ == "__main__":
    print("Computation of Jaccard similarity between files")

    # Load all the files in a dictionary
    files = load_directory("data")
    k = 21

    print("Computing Jaccard similarity for all pairs of samples")
    filenames = list(files.keys())

    for i in range(len(files)):
        for j in range(i+1, len(files)):
            jac = jaccard(files[filenames[i]], files[filenames[j]], k)
            print(filenames[i], filenames[j], jac)