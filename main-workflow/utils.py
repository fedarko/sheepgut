import networkx as nx

def load_gfa(filepath):
    
    # We ignore directionality for right now.
    graph = nx.Graph()

    with open(filepath, "r") as gfafile:
        for line in gfafile:
            # Parse sequences. We convert these to "nodes" in the networkx
            # representation of the graph,
            # but they are really edges in the repeat graph produced by
            # metaFlye.
            if line[0] == "S":
                parts = line.strip().split("\t")
                # Strip off the "edge_" name prefix -- turn "edge_123" into
                # "123", for example
                node_name = parts[1][5:]
                node_len = len(parts[2])
                node_cov = None

                # Parse GFA tags
                extra_data = parts[3:]
                for tag in extra_data:
                    # This is overcautious, probably, but whatevs
                    if tag.startswith("LN:i:"):
                        raise ValueError("Duplicate length for node {}".format(node_name))
                    elif tag.startswith("dp:i:") or tag.startswith("KC:i:"):
                        if node_cov is None:
                            node_cov = int(tag[5:])
                        else:
                            raise ValueError("Duplicate coverage for node {}".format(node_name))

                if node_cov is None:
                    raise ValueError("No coverage tag given for node {}".format(node_name))

                graph.add_node(node_name, length=node_len, cov=node_cov)

            # Parse links between sequences.
            # Each link line looks something like:
            #
            # L    edge_74133    -    edge_71431    -    0M
            #
            # (where the whitespace gaps are each a single tab character)
            elif line[0] == "L":
                parts = line.strip().split("\t")
                # Again, the [5:] slices are needed to remove the "edge_"
                # prefix.
                src = parts[1][5:]
                snk = parts[3][5:]
                # As mentioned, we ignore directionality for this specific
                # application. So, e.g.:
                # A+ -> B+
                # A+ -> B-
                # A- -> B+
                # A- -> B-
                # ... would all get treated as the same edge. "Duplicate" edges
                # are implicitly ignored by networkx.
                graph.add_edge(src, snk)
    return graph
