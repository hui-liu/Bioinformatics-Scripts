#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import re
from collections import deque

# result 1: 5' LTR or 3' LTR unique hit positision file
# result 2: id mapping file, column 1 is the intact LTR ID, column 2 is the hit ID,
#           column 3 is the ID of hit extended 3k, column is the cluster ID
#python ../bin/Intact_LTR_uniq_hit.py Xs_intact_5_3_ltr.outfmt6 \
#/media/12TB/Xso_genome/genome/genome2.fasta.fai \
#../clustering/Xs_intact_5ltr_clustering.txt Intact_LTR_uniq_hit.txt Idmapping.txt


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

ident = 90.0
cov = 90.0
Evalue = 1e-10

regex1 = re.compile("(\w+)_\d+_\d+_\d+")
regex2 = re.compile("(\w+_\d+_\d+)_\d+")

aDict = {}
with open(sys.argv[1], 'r') as f:
    for line in f:
        qseqid, sseqid, pident, qcovs, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = line.split()
        strand = "+"
        if int(sstart) > int(send):
            sstart, send = send, sstart
            strand = "-"
        qchrom = regex1.findall(qseqid)[0]
        qgstart, qgend = qseqid.split("_")[-3:-1]
        if qchrom == sseqid and int(qgstart) <= int(sstart) and int(qgend) >= int(send): continue
        if float(pident) >= ident and float(qcovs) >= cov and float(evalue) < Evalue:
            aDict.setdefault(sseqid, []).append([sseqid, int(sstart), int(send), strand, float(bitscore), qseqid])

# graph implemented as a dictionary
graph = {}
xy_set = set()
for i in aDict:
    print i
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

# clustering nodes
net = {}
with open(sys.argv[3], 'r') as f:
    for line in f:
        lsp = line.split()
        net[lsp[1]] = lsp[0]
# out
OUT = open(sys.argv[4], 'w')
Idmapping = open(sys.argv[5], 'w')
for i in uniq_lis:
    s, e = i[1] - 3000, i[2] + 3000
    if s < 1:
        s = 1
    if e > chrom_len[i[0]]:
        e = chrom_len[i[0]]
    OUT.write("\t".join([i[0], str(s), str(e), i[3]]) + "\n")
    Intact_id = regex2.findall(i[5])[0]
    Hit_id = i[0] + "_" + str(i[1]) + "_" + str(i[2])
    Extended_hit_id = i[0] + "_" + str(s) + "_" + str(e)
    Idmapping.write("\t".join([Intact_id, Hit_id, Extended_hit_id, net[Intact_id]]) + "\n")

for j in components:
    j.sort(key = lambda x: -x[4])
    j_best = j[0]
    s, e = j_best[1] - 3000, j_best[2] + 3000
    if s < 1:
        s = 1
    if e > chrom_len[j_best[0]]:
        e = chrom_len[j_best[0]]
    OUT.write("\t".join([j_best[0], str(s), str(e), j_best[3]]) + "\n")
    Intact_id = regex2.findall(j_best[5])[0]
    Hit_id = j_best[0] + "_" + str(j_best[1]) + "_" + str(j_best[2])
    Extended_hit_id = i[0] + "_" + str(s) + "_" + str(e)
    Idmapping.write("\t".join([Intact_id, Hit_id, Extended_hit_id, net[Intact_id]]) + "\n")

OUT.close()
Idmapping.close()
