import sys
import gzip

with open(sys.argv[1], 'r') as f:
    for line in f:
        print(line.rstrip())

svids = {}
with open(sys.argv[2], 'r') as f:
    for line in f:
        lsp = line.rstrip().split()
        svids[lsp[0]] = lsp[1]

with gzip.open(sys.argv[3], 'rt') as f:
    for line in f:
        if line[:2] == "##":
            continue
        else:
            lsp = line.rstrip().split()
            if lsp[0] == "#CHROM":
                print(line.rstrip())
            else:
                ref, alt = lsp[3], lsp[4]
                id = "sv_" + lsp[0] + "_" + lsp[1]
                id2 = id + "_" + ref + "_" + alt
                if id2 in svids and lsp[6] == "PASS":
                    lsp[2] = id
                    lsp[7] = lsp[7] + ";" + svids[id2]
                    print("\t".join(lsp))
