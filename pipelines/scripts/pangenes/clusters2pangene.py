import sys
import itertools

ref = 'QDE'
pangene = {}
n = 0
with open(sys.argv[1], 'r') as f:
    for line in f:
        lsp = [x for x in line.split() if x]
        ogid = lsp[0].rstrip(":")
        og = lsp[1:]
        sp = set([x[:3] for x in og])
        sp_num = len(sp)
        if ref not in sp:
            n += 1
            newid = ref + str(n)
            pangene[newid] = sp_num
        else:
            ref_genes = [x for x in og if x[:3] == ref]
            for g in ref_genes:
                pangene[g] = sp_num

for i in pangene:
    print i, pangene[i]

