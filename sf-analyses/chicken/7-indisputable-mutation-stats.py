#! /usr/bin/env python3
# just a quick script to demonstrate that there are 503 indisputable
# p-mutations (p = 5%) in target contig s7.ctg000008c
from strainflye.bcf_utils import parse_sf_bcf
from strainflye.call_utils import call_p_mutation

ni = 0
p = 500
contig = "s7.ctg000008c"

bcf, tt, tm = parse_sf_bcf(
    "/Poppy/mfedarko/chicken-gut-meta/sf/call-p-minp100/naive-calls.bcf"
)
for mut in bcf.fetch(contig):
    if call_p_mutation(mut.info["AAD"][0], mut.info["MDP"], p, 2):
        ni += 1

print(f"There are {ni:,} p-mutations for p = {p:,} in {contig}.")
