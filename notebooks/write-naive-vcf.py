import pileup
spp = pileup.load()
ls
spp["edge_6104"][27553]
spp["edge_6104"][27552]
spp["edge_6104"][27554]
13 / (4925 + 13)
25 / (4912 + 25)
0.5 / 100
vcf = ""
seq2mp = defaultdict(list)
for seq in SEQS:
    for pos, pcol in enumerate(spp[seq][1:], 1):
        cov = pileup.get_cov(pcol)
        if cov >= 1000:
            if pileup.naively_call_mutation(pcol, 0.005):
                # vcf positions are 1-indexed, same as pileup stuff
                seq2mp[seq].append(pos)
                # NOTE this may wind up having ref == alt; doesn't matter much r/n but can change if needed
                ref_nt = "ACGT"[pcol[1]]
                alt_nt = pileup.get_alt_nt(pcol)
                vcf += f"{seq}\t{pos}\t.\t{ref_nt}\t{alt_nt}\t.\t.\t.\n"
from collections import defaultdict
vcf = ""
seq2mp = defaultdict(list)
for seq in SEQS:
    for pos, pcol in enumerate(spp[seq][1:], 1):
        cov = pileup.get_cov(pcol)
        if cov >= 1000:
            if pileup.naively_call_mutation(pcol, 0.005):
                # vcf positions are 1-indexed, same as pileup stuff
                seq2mp[seq].append(pos)
                # NOTE this may wind up having ref == alt; doesn't matter much r/n but can change if needed
                ref_nt = "ACGT"[pcol[1]]
                alt_nt = pileup.get_alt_nt(pcol)
                vcf += f"{seq}\t{pos}\t.\t{ref_nt}\t{alt_nt}\t.\t.\t.\n"
%run "Header.ipynb"
vcf = ""
seq2mp = defaultdict(list)
for seq in SEQS:
    for pos, pcol in enumerate(spp[seq][1:], 1):
        cov = pileup.get_cov(pcol)
        if cov >= 1000:
            if pileup.naively_call_mutation(pcol, 0.005):
                # vcf positions are 1-indexed, same as pileup stuff
                seq2mp[seq].append(pos)
                # NOTE this may wind up having ref == alt; doesn't matter much r/n but can change if needed
                ref_nt = "ACGT"[pcol[1]]
                alt_nt = pileup.get_alt_nt(pcol)
                vcf += f"{seq}\t{pos}\t.\t{ref_nt}\t{alt_nt}\t.\t.\t.\n"
vcf
vcf.split("\n")[0]
print(_)
spp["edge_6104"][1596]
25 / (25+4842)
with open("../seqs/naive_mincov1000_p0.5pct.vcf", "w") as vf:
    vf.write(vcf)
ls
ls ../seqs
ls
%history -f write-naive-vcf.py
