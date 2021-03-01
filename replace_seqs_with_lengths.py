# Removes the actual sequences, which reduces a ton of filesize and makes life
# easier for everyone
out = ""
with open("/Poppy/mkolmogo/sheep_meta/flye_big_2.8/assembly_graph.gfa", "r") as gf:
    for line in gf:
        if line.startswith("S\t"):
            split = line.split()
            seq = split[2]
            name = split[1][5:]
            out += "S\t{}\t*\tLN:i:{}\t{}\n".format(
                name, len(seq), split[3]
            )
        elif line.startswith("L\t"):
            out += line.replace("edge_", "")
        else:
            out += line

with open("20201202_graph_noseq.gfa", "w") as f:
    f.write(out)
