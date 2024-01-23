import sys
import itertools
import random

with open(sys.argv[1], 'r') as f:
    for line in f:
        lsp = [x for x in line.rstrip().split("\t") if x]
        if lsp[0] == "Orthogroup": continue
        ogid = lsp[0]
        og = [x.split(", ") for x in lsp[1:]]
        #num = [len(x) for x in og]
        # at least 2 species
        if len(og) > 1:
            # if there were more than one gene for a certain species in a group, randomly select a gene for this species.
            og2 = [random.choice(x) for x in og]
            # all possible combinations of every two groups
            comb_gr = [x for x in itertools.combinations(og2, 2)]
            for p in comb_gr:
                print "\t".join((ogid,) + p)
