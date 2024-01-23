import sys
import itertools
import random

with open(sys.argv[1], 'r') as f:
    for line in f:
        lsp = [x for x in line.rstrip().split(",") if x]
        if lsp[0] == "Pan_gene_ID": continue
        loci_id = lsp[0]
        loci = [x.split(";") for x in lsp[2:] if x!= "NA"]
        # at least 2 species
        if len(loci) > 1:
            # if there were more than one gene for a certain species in a group, randomly select a gene for this species.
            selected_loci = [random.choice(x) for x in loci]
            # all possible combinations of every two groups
            comb_gr = [x for x in itertools.combinations(selected_loci, 2)]
            for p in comb_gr:
                if "gmap_ID" in p[0] or "gmap_ID" in p[1]: continue 
                print "\t".join((loci_id,) + p)
