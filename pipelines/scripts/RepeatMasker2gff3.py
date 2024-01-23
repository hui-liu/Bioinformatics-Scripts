import sys

OUT = open(sys.argv[2], 'w')
n=0
with open(sys.argv[1], 'r') as f:
    for line in f:
        lsp = line.split()
        if not lsp: continue
        if lsp[0] in ["SW", "score"]: continue
        n += 1
        score = str(int(round(float(lsp[1]))))
        strand = lsp[8]
        if strand != "+":
            strand = "-"
        length = int(lsp[6]) - int(lsp[5]) + 1
        ID = "ID=" + lsp[4] + ":hit:" + str(n) + ":repeatmasker"
        ID2 = "ID=" + lsp[4] + ":hsp:" + str(n) + ":repeatmasker"
        Name = "Name=" + lsp[4] + "_" + lsp[5] + "_" + lsp[6] + ".1" + "|genus:" + lsp[10]
        Target = "Target=" + lsp[4] + "_" + lsp[5] + "_" + lsp[6] + ".1" + "|genus:" + lsp[10] + " 1 " + str(length) + " " + strand
        OUT.write("\t".join([lsp[4], "RepeatMasker", "match", lsp[5], lsp[6], score, lsp[8], ".", ID + ";" + Name + ";" + Target]) + "\n")

OUT.close()
