import sys

depth = {}
with open(sys.argv[1], 'r') as f:
    for line in f:
        id, pos, d = line.split()
        depth.setdefault(id, []).append(int(d))

out = open(sys.argv[2], 'w')
for i in depth:
    mean_depth = round(sum(depth[i]) * 1.0 / len(depth[i]), 2)
    out.write("\t".join([i, str(mean_depth)]) + "\n")

out.close()
