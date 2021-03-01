def count_homonucleotide_runs(seq):
    k = 1
    num_runs = []
    while True:
        num_runs.append({"A": 0, "C": 0, "G": 0, "T": 0})
        any_runs_found = False
        for i in range(0, len(seq) - k + 1):
            kmer = seq[i:i + k]
            if len(set(kmer)) == 1:
                num_runs[k - 1][kmer[0]] += 1
                any_runs_found = True
        if not any_runs_found: break
        else: k += 1
    return num_runs

with open("binary_gene.fasta", "r") as f:
    g1 = f.read().strip()
with open("a2_gene.fasta", "r") as f:
    g2 = f.read().strip()

count_homonucleotide_runs(g1)
count_homonucleotide_runs(g2)
