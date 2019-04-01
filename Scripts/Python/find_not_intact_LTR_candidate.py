#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys

from collections import deque

def breadth_first_search(graph, root):
    visited = set([root])
    queue = deque([root])
    while queue:
        vertex = queue.popleft()
        yield vertex
        for neighbour in graph[vertex]:
            if neighbour not in visited:
                visited.update([neighbour])
                queue.append(neighbour)

def overlaps(a, b):
    result = True
    if a[1] < b[0] or b[1] < a[0]:
        result = False
    return result

def AllPairs(lis):
    for p1 in range(len(lis)):
        for p2 in range(p1+1, len(lis)):
           yield (lis[p1], lis[p2])

# (1) The hits of query should not the same as the position of query
# (2) pident >= 90, qcovs >= 90, evalue < 1e-10
# (3) find the unique hit of query (if multiple hits, filter out by the bitscore value)

ident = float(90)
cov = float(90)
Evalue = 1e-10

aDict = {}
with open(sys.argv[1], 'r') as f:
    for line in f:
        qseqid, sseqid, pident, qcovs, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = line.split()
        strand = "+"
        if int(sstart) > int(send):
            sstart, send = send, sstart
            strand = "-"
        qchrom = qseqid.split("_chr")
        qgstart, qgend = qseqid.split("_")[:-3:-1]
        if qchrom == sseqid and int(qgstart) <= int(sstart) and int(qgend) >= int(send): continue
        if float(pident) >= ident and float(qcovs) >= cov and float(evalue) < Evalue:
            aDict.setdefault(sseqid, []).append([sseqid, int(sstart), int(send), strand, float(bitscore)])

# graph implemented as a dictionary
graph = {}
xy_set = set()
for i in aDict:
    allpairs = AllPairs(aDict[i])
    for x, y in allpairs:
        # same strand and overlap
        if x[3] == y[3] and overlaps(x[1:3], y[1:3]):
            graph.setdefault(tuple(x), []).append(tuple(y))
            graph.setdefault(tuple(y), []).append(tuple(x))
            xy_set.update([tuple(x), tuple(y)])

# find the unique hit
uniq_lis = []
for i in aDict:
    for j in aDict[i]:
        if tuple(j) not in graph:
            uniq_lis.append(j)

# find the biggest connected component
components = []
seen = set()
for g in xy_set:
    if g in seen:
        continue
    else:
        c = [n for n in breadth_first_search(graph, g)]
        seen.update(c)
        components.append(c)

# chromosomes length
chrom_len = {}
with open(sys.argv[2], 'r') as f:
    for line in f:
        lsp = line.split()
        chrom_len[lsp[0]] = lsp[1]
# out
OUT = open(sys.argv[3], 'w')
for i in uniq_lis:
    s, e = i[1] - 3000, i[2] + 3000
    if s < 1:
        s = 1
    if e > chrom_len[i[0]]:
        e = chrom_len[i[0]]
    OUT.write("\t".join([i[0], str(s), str(e), i[3]]) + "\n")

for j in components:
    j.sort(key = lambda x: -x[4])
    j_best = j[0]
    s, e = j_best[1] - 3000, j_best[2] + 3000
    if s < 1:
        s = 1
    if e > chrom_len[j_best[0]]:
        e = chrom_len[j_best[0]]
    OUT.write("\t".join([j_best[0], str(s), str(e), j_best[3]]) + "\n")

OUT.close()
