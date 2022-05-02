import time
import pileup
from collections import defaultdict

# Stuff included elsewhere and duplicated here (sorry!)
SEQS = ["edge_6104", "edge_1671", "edge_2358"]
n2i = {"A": 0, "C": 1, "G": 2, "T": 3}
i2n = "ACGT"

spp = pileup.load()
# Header info gleaned by reading over the VCF 4.2 docs
# (https://samtools.github.io/hts-specs/VCFv4.2.pdf) and copying how LoFreq
# organizes their header
vcf = (
    "##fileformat=VCFv4.2\n"
    f"##fileDate={time.strftime('%Y%m%d')}\n"
    # TODO include version number in the source, as shown in VCF 4.2 docs?
    "##source=strainFlye\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
)
for seq in SEQS:
    for pos, pcol in enumerate(spp[seq][1:], 1):
        cov = pileup.get_cov(pcol)
        if cov >= 1000:
            if pileup.naively_call_mutation(pcol, 0.5):
                # vcf positions are 1-indexed, same as pileup stuff
                alt_nt = pileup.get_alt_nt(pcol)
                # extremely silly hack to get the reference nt
                # this might not be necessary but this should consistently
                # at least set it to something that is not the alt nt
                ctscopy = pcol[0][:]
                ctscopy[n2i[alt_nt]] = float("-inf")
                ref_nt = max("ACGT", key=lambda nt: ctscopy[n2i[nt]])
                vcf += f"{seq}\t{pos}\t.\t{ref_nt}\t{alt_nt}\t.\t.\t.\n"

with open("../seqs/naive_mincov1000_p0.5pct.vcf", "w") as vf:
    vf.write(vcf)
