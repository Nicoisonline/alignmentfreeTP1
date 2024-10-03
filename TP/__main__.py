from TP.loading import load_directory
from TP.kmers import stream_kmers, kmer2str

def create_index(text, k):
    """ Create a kmer index from a text"""
    index = {}
    # --- To complete ---
    for kmer in stream_kmers(text, k):
        if kmer in index:
            index[kmer] += 1
        else:
            index[kmer] = 1
    return index

def compare_sequences(seq1, seq2, k): 
    index = create_index(seq1, k)
    intersection = 0
    union = sum(index.values())
    for kmer in stream_kmers(seq2, k):
        if kmer in index:
            intersection += 1
            union += index[kmer] - 1
        else:
            union += 1

def jaccard(fileA, fileB, k):
    """ Compute the Jaccard similarity between two files"""
    j = 0
    # --- To complete ---
    indexA = create_index(fileA, k)
    indexB = create_index(fileB, k)
    intersection = 0
    union = 0
    for kmer in indexA:
        if kmer in indexB:
            intersection += 1
            union += max(indexA[kmer], indexB[kmer])
        else:
            union += indexA[kmer]
    for kmer in indexB:
        if kmer not in indexA:
            union += indexB[kmer]
    j = intersection / union
    return j



if __name__ == "__main__":
    print("Computation of Jaccard similarity between files")

    # Load all the files in a dictionary
    files = load_directory("data")
    k = 21
    
    print("Computing Jaccard similarity for all pairs of samples")
    filenames = list(files.keys())
    for i in range(len(files)):
        for j in range(i+1, len(files)):
            
            # --- Complete here ---
            
            j = jaccard(files[filenames[i]], files[filenames[j]], k)
            print(filenames[i], filenames[j], j)
