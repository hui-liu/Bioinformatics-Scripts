import sys

orthologs = set()
for line in sys.stdin:
    temp = line.rstrip().split("\t")
    if temp[0] == "Orthogroup": continue
    for i in temp[1].split(", "):
        for j in temp[2].split(", "):
            pair = (i, j)
            if pair[-1::-1] in orthologs: continue
            orthologs.add((i, j))

for i in orthologs:
    print("\t".join(i))
