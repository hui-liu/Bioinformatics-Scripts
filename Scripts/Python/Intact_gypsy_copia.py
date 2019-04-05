import sys
import re

# python ../bin/Intact_gypsy_copia.py xso_ltrdigest_tabout.csv Xs_intact.csv
# xso_ltrdigest_tabout.csv was gernated from ltrdigest

def Find_domain(strings, domain):
    strings = strings + "/"
    starts = [i for i in range(len(strings)) if strings.startswith(domain, i)]
    ends = [strings.find("/", i) for i in starts]
    domains = [strings[starts[i]: ends[i]] for i in xrange(len(starts))]
    return domains

def Gypsy_Copia_class(strings):
    gypsy_arr = ('GAG', 'AP', 'RT', 'RNaseH', 'INT')
    copia_arr = ('GAG', 'AP', 'INT', 'RT', 'RNaseH')
    regex = re.compile("(\w+?)_(\S+)")
    # find all gypsy or copia domains
    domains = [Find_domain(strings, i) for i in gypsy_arr]
    # cluster domains from the same family
    domains_suffix = []
    for x in domains:
        if x:
            domains_suffix.append([regex.findall(y)[0][1] for y in x])
    # should contain five domains
    if all(set(domains_suffix[i]) & set(domains_suffix[i+1]) for i in xrange(len(domains_suffix)-1)) and len(domains_suffix) == 5:
        s1, s2, s3, s4, s5 = [set(i) for i in domains_suffix]
        overlap = list(set.intersection(s1, s2, s3, s4, s5))
        if overlap:
            if len(overlap) == 1:
                family = overlap[0]
                # gypsy and copia sub family
                gypsy = [strings.find(i + "_" + family) for i in gypsy_arr]
                copia = [strings.find(i + "_" + family) for i in copia_arr]
                # check the domains order
                gypsy_order = all(gypsy[i] <= gypsy[i+1] for i in xrange(len(gypsy)-1))
                copia_order = all(copia[i] <= copia[i+1] for i in xrange(len(copia)-1))
                if gypsy_order:
                    return ["Gypsy", "/".join([i + "_" + family for i in gypsy_arr])]
                elif copia_order:
                    return ["Copia", "/".join([i + "_" + family for i in copia_arr])]
                elif gypsy_order and copia_order:
                    return ["nest", "/".join([i + "_" + family for i in copia_arr]), "/".join([i + "_" + family for i in gypsy_arr])]
            else:
                res = []
                for family in overlap:
                    gypsy = [strings.find(i + "_" + family) for i in gypsy_arr]
                    copia = [strings.find(i + "_" + family) for i in copia_arr]
                    gypsy_order = all(gypsy[i] <= gypsy[i+1] for i in xrange(len(gypsy)-1))
                    copia_order = all(copia[i] <= copia[i+1] for i in xrange(len(copia)-1))
                    if gypsy_order:
                        res.append(["Gypsy"] + [i + "_" + family for i in gypsy_arr])
                    elif copia_order:
                        res.append(["Copia"] + [i + "_" + family for i in copia_arr])
                if len(res) == 1:
                    return res[0]
                elif len(res) > 1:
                    return ["nest", "/".join([i + "_" + family for i in copia_arr]), "/".join([i + "_" + family for i in gypsy_arr])]

OUT = open(sys.argv[2], 'w')
with open(sys.argv[1], 'r') as f:
    for line in f:
        lsp = line.rstrip("\n").split("\t")
        if lsp[0] == "element start": continue
        annotation = lsp[-1]
        if annotation:
            Family = Gypsy_Copia_class(annotation)
            if Family:
                start, end, length, seq = map(str, lsp[:4])
                OUT.write("\t".join(["_".join([seq, start, end]), length] + Family) + "\n")            
OUT.close()
