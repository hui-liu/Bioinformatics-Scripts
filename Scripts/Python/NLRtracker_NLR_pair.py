import sys
from collections import deque
import re

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

def parsebed(filename):
    gff_dic = {}
    with open(filename, 'r') as fh:
        for line in fh:
            lsp = line.split()
            gene_id = lsp[3]
            chr, start, end, strand = lsp[0], int(lsp[1])+1, int(lsp[2]), lsp[5]
            gff_dic.setdefault(chr, []).append([gene_id, start, end, strand])
    return gff_dic

def AllPairs(lis):
    for p1 in range(len(lis)):
        for p2 in range(p1+1, len(lis)):
           yield (lis[p1], lis[p2])

# parse the gff file
genes_bed = parsebed(sys.argv[1])
genes_info = {}
for chrom in genes_bed: # chroms
    num = 0
    # sort by start per chrom
    genes = sorted(genes_bed[chrom], key=lambda x: x[1])
    for gene in genes: # genes info list
        num += 1 # the order of gene within the chrom
        genes_info[gene[0]] = [chrom, num] + gene[1:]

#
NLR_genes = {}
with open(sys.argv[2], 'r') as f:
    h = f.readline()
    for line in f:
        lsp = line.rstrip().split("\t")
        gene = lsp[0]
        type = lsp[5]
        NLR_genes.setdefault(gene, []).extend(genes_info[gene] + [type])

# graph implemented as a dictionary
graph = {}
# all genes
genes_set = set()
# get the gene name list
genes = [k for k in NLR_genes]
# the generator contained the combination of every two genes
allpairs = AllPairs(genes)
for g1, g2 in allpairs:
    g1_chrom, g2_chrom = genes_info[g1][0], genes_info[g2][0]
    g1_order, g2_order = int(genes_info[g1][1]), int(genes_info[g2][1])
    #g1_start, g1_end = int(genes_info[g1][2]), int(genes_info[g1][3])
    #g2_start, g2_end = int(genes_info[g2][2]), int(genes_info[g2][3])
    g1_strand, g2_strand = genes_info[g1][4], genes_info[g2][4]
    g1_g2_span = abs(g2_order - g1_order) - 1
    if g1_chrom == g2_chrom:
        if g1_g2_span == 0 and g1_strand != g2_strand:
            graph.setdefault(g1, []).append(g2)
            graph.setdefault(g2, []).append(g1)
            genes_set.update([g1, g2])

# the number of genes
n = len(genes_set)
# get the biggest connected component
components = []
seen = set()
for gene in genes_set:
    if gene in seen:
        continue
    else:
        c = [n for n in breadth_first_search(graph, gene)]
        seen.update(c)
        c.sort()
        components.append(c)

# sort
clusters = [[[g, genes_info[g][1], genes_info[g][0]] for g in component] for component in components]
clusters.sort(key = lambda x: (filter(lambda x: x.isalpha(), x[1][2]), int(re.findall(r'\d+', x[1][2])[0], base = 10), x[1][1]))

# output
OUT = open(sys.argv[3], 'w')
num = 0
for cluster in clusters:
    num += 1
    #cluster.sort(key = lambda x: x[1])
    for lis in cluster:
        ID = "pair_" + str(num)
        OUT.write(lis[0] + "\t" + ID + "\t" + "\t".join(map(str, genes_info[lis[0]])) + "\t" + NLR_genes[lis[0]][5] + "\n")
OUT.close()
