import sys


min_sw = 1000
max_div = 30
min_len = 1000
n = 0
OUT = open(sys.argv[2], 'w')
with open(sys.argv[1], 'r') as f:
    for line in f:
        lsp = line.split()
        if not lsp: continue
        if lsp[0] in ["SW", "score"]: continue
        SW_score, div, chr, element_start, element_end, left_len, strand, TE_ID, TE_class = [lsp [i] for i in 0,1,4,5,6,7,8,9,10]
        if strand != "+":
            strand = "-"
        n += 1
        SW_score = int(round(float(lsp[0])))
        div = int(round(float(lsp[1])))
        length = int(element_end) - int(element_start) + 1
        ID = "ID=" + chr + ":hit:" + str(n) + ":edta"
        ID2 = "ID=" + chr + ":hsp:" + str(n) + ":edta"
        type = "Simple_repeat"
        if SW_score >= min_sw and div <= max_div and length >= min_len:
            type = "TE"
        Name = "Name=" + chr + "_" + element_start + "_" + element_end + ".1" + "|genus:" + type
        Target = "Target=" + chr + "_" + element_start + "_" + element_end + ".1" + "|genus:" + type
        OUT.write("\t".join([chr, "EDTA", "match", element_start, element_end, str(div), strand, ".", ID + ";" + Name + ";" + Target + " 1 " + str(length) + " " + strand]) + "\n")
        OUT.write("\t".join([chr, "EDTA", "match_part", element_start, element_end, str(div), strand, ".", ID2 + ";" + Name + ";" + Target + " 1 " + str(length) + " " + strand]) + "\n")
OUT.close()
