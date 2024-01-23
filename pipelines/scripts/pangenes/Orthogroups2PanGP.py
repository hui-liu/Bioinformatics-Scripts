import sys

def parseOrthogr(filename):
    n = 0
    res = []
    with open(filename, 'r') as f:
        for line in f:
            n += 1
            lsp = line.rstrip().split("\t")
            if n == 1:
                header = ["ClusterID"] + lsp[1:]
                codes = map(lambda x: x.upper(), [x.split("_")[0][0]+x.split("_")[1][:2] for x in lsp[1:]])
            else:
                d = {}
                for i in lsp[1:]:
                    s = i[:3]
                    t = i.replace(" ", "")
                    d[s] = t
                lst = []
                for v in codes:
                    if v in d:
                        lst.append(d[v])
                    else:
                        lst.append("-")
            if n > 1:
                res.append(lst)
    return header, res

header, cluster = parseOrthogr(sys.argv[1])
_, cluster2 = parseOrthogr(sys.argv[2])
out = open(sys.argv[3], 'w')
out.write("\t".join(header) + "\n")

num = 0
for i in cluster + cluster2:
    num += 1
    out.write("\t".join([str(num)] + i) + "\n")

