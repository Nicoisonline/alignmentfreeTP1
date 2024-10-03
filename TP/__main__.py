from TP.loading import load_directory
from TP.kmers import stream_kmers, kmer2str



def jaccard(fileA, fileB, k):
    """ Compute the Jaccard similarity between two files"""
    j = 0
    # --- To complete ---
    kmersA = set(stream_kmers(fileA, k))
    kmersB = set(stream_kmers(fileB, k))
    intersection = kmersA.intersection(kmersB)
    union = kmersA.union(kmersB)
    j = len(intersection) / len(union)
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
