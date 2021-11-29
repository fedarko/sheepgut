#! /usr/bin/env python3
import os
import subprocess
PROG = "sfdp"
prefix = "graphs/"
for f in sorted(os.listdir(prefix)):
    if f.endswith(".gv"):
        gf = prefix + f
        pf = gf[:-3] + ".png"
        print(f"Converting {gf} -> {pf} using {PROG}...", end=" ", flush=True)
        with open(pf, "w") as pngfile:
            subprocess.run([PROG, gf, "-Tpng"], stdout=pngfile)
        print("Done!")
