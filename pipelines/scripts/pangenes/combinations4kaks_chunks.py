import sys
import itertools

n = 0
num = 1
with open(sys.argv[1], 'r') as f:
    for line in f:
        lsp = [x for x in line.rstrip().split("\t") if x]
        if lsp[0] == "Orthogroup": continue
        ogid = lsp[0]
        og = [x.split(", ") for x in lsp[1:]]
        #num = [len(x) for x in og]
        # at least 2 species
        if len(og) > 1:
            # all possible combinations of every two groups
            comb_gr = [x for x in itertools.combinations(og, 2)]
            for p in comb_gr:
                a, b = p
                pairs = list(itertools.product(a, b))
                for x in pairs:
                    n += 1
                    if n % 100 ==0:
                        num += 1
                    print "\t".join(("chunk" + str(num), ogid) + x)
