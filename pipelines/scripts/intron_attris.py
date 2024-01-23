import sys

intron = {}
lines = []
source = {}
for line in sys.stdin:
    lsp = line.rstrip().split("\t")
    if lsp[2] == "intron":
        id = lsp[8].split("=")[1]
        intron.setdefault(id, []).append(lsp)
    if lsp[2] == "mRNA":
        id = lsp[8].split(";")[0].split("=")[1]
        source[id] = lsp[1]
    lines.append(lsp)

intron_index = {i: [x for x in range(len(intron[i]))] for i in intron}
temp = {i: 0 for i in intron}
out = open(sys.argv[1], 'w')
for line in lines:
    if line[2] == "intron":
        id = line[8].split("=")[1]
        if line[6] == "+":
            attri = "ID=" + id + ".intron" + str(intron_index[id][temp[id]]+1) + ";Parent=" + id
        else:
            attri = "ID=" + id + ".intron" + str(intron_index[id][::-1][temp[id]]+1) + ";Parent=" + id
        temp[id] += 1
        out.write("\t".join([line[0]] + [source[id]] + line[2:-1] + [attri]) + "\n")
    else:
        out.write("\t".join(line) + "\n")
out.close()

