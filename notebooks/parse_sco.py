import pandas as pd

def parse_sco(filepath):
    """Returns a DataFrame representing a SCO file produced by e.g. Prodigal.

    This is intended to sort of mimic the way LST files from GeneMark
    can be parsed with Pandas (...with some effort).
    """
    genes = {}
    with open(filepath, "r") as f:
        for line in f:
            if not line.startswith("#"):
                if line.startswith(">"):
                    parts = line[1:].strip().split("_")
                    gene_num = int(parts[0])
                    left_end = int(parts[1])
                    rght_end = int(parts[2])
                    strand = parts[3]
                    length = rght_end - left_end + 1
                    genes[gene_num] = [left_end, rght_end, length, strand]
                else:
                    # If this line doesn't start with # or > and it isn't just
                    # a blank line, then something's up. Throw an error.
                    if len(line.strip()) > 0:
                        raise ValueError(
                            "Unrecognized line prefix for line {}".format(line)
                        )
    df = pd.DataFrame.from_dict(
        genes, orient="index", columns=["LeftEnd", "RightEnd", "Length", "Strand"]
    )
    return df
