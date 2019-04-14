import sys
from collections import OrderedDict


def CheckDomain(strs):
    gypsy_arr = ('GAG', 'PROT', 'RT', 'RH', 'INT')
    copia_arr = ('GAG', 'PROT', 'INT', 'RT', 'RH')
    gypsy = [strs.find(x) for x in gypsy_arr]
    copia = [strs.find(x) for x in copia_arr]
    # check the domains rank
    gypsy_rank = all(gypsy[y] <= gypsy[y+1] for y in xrange(len(gypsy)-1))
    copia_rank = all(copia[y] <= copia[y+1] for y in xrange(len(copia)-1))
    return gypsy_rank, copia_rank

# parse annotatation
family = OrderedDict()
with open(sys.argv[1]) as f:
    for line in f:
        if line[0] == "#":
            continue
        lsp = line.split()
        chr, source, type, start, end, score, strand, phase, attibutes = lsp
        d_attr = dict([v.split('=') for v in attibutes.strip(';').split(';')])
        a = [d_attr["Name"], strand, d_attr['Final_Classification']]
        family.setdefault(chr, []).append(a)

# find subfamily
subfamily = OrderedDict()
for i in family:
    temp = OrderedDict()
    for j in family[i]:
        N, S, F = j
        temp.setdefault((F, S), []).append(N)
    subfamily[i] = temp

# find intact LTR-RTs
OUT = open(sys.argv[2], 'w')
for i in subfamily:
    for j in subfamily[i]:
        subfam, strand = j
        domains = subfamily[i][j]
        if len(set(domains)) >= 5:
            strings = "/".join(domains)
            out_lst = [i, subfam, strings]
            if strand == "+":
                gypsy_Rank, copia_Rank = CheckDomain(strings)
                if gypsy_Rank and not copia_Rank:
                    OUT.write("\t".join(out_lst) + "\n")
                elif copia_Rank and not gypsy_Rank:
                    OUT.write("\t".join(out_lst) + "\n")
            else:
                strings = "/".join(domains[::-1])
                gypsy_Rank, copia_Rank = CheckDomain(strings)
                if gypsy_Rank and not copia_Rank:
                    OUT.write("\t".join(out_lst) + "\n")
                elif copia_Rank and not gypsy_Rank:
                    OUT.write("\t".join(out_lst) + "\n")

OUT.close()
