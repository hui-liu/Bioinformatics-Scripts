import sys

# modified from https://github.com/AnimalGenomicsETH/bovine-graphs/blob/d62519594eee194f98e4fd55109e89c94363f303/scripts/get_bialsv.py

def sv_det(line):
    # SDTS1#Chr1 526029 3 2 8998610,8998611,8998612
    # nodes: number of GFA segments in the bubble including the source and the sink of the bubble
    chromo, pos, nodes, _, nodelist = line.strip().split()
    # remove source and sink
    bubble = nodelist.split(",")[1:-1]
    sourcenode = nodelist.split(",")[0]
    sinknode = nodelist.split(",")[-1]
    if nodes in ["3", "4"]:
        if nodes == "3":
            rrank, nodelen = nodeinf[bubble[0]]
            if rrank > 0:
                svtype = "Insertion"
                reflen = 0
                nonreflen = nodelen
                bubnonref = bubble[0]
                bubref = 0
            else:
                svtype = "Deletion"
                reflen = nodelen
                nonreflen = 0
                bubnonref = 0
                bubref = bubble[0]
        elif nodes == "4":
            # do not consider looping in the bubble
            if bubble[0] == bubble[1]:
                return None
            # need to check?
            rrank1, nodelen1 = nodeinf[bubble[0]]
            rrank2, nodelen2 = nodeinf[bubble[1]]
            if rrank1 == 0 and rrank2 == 0:
                #svtype = "Deletion"
                #reflen = str(nodelen1) + "," + str(nodelen2)
                #nonreflen = 0
                #bubnonref = 0
                #bubref = str(bubble[0]) + "," + str(bubble[1])
                return None
            elif rrank1 > 0 and rrank2 > 0:
                #svtype = "Insertion"
                #reflen = 0
                #nonreflen = str(nodelen1) + "," + str(nodelen2)
                #bubnonref = str(bubble[0]) + "," + str(bubble[1])
                #bubref = 0
                return None
            else:
                #svtype = "Divergent"
                for bub in bubble:
                    rrank, nodelen = nodeinf[bub]
                    if rrank > 0:
                        nonreflen = nodelen
                        bubnonref = bub
                    else:
                        reflen = nodelen
                        bubref = bub
                if reflen == nonreflen:
                    if reflen == 1:
                        svtype = "SNP"
                    else:
                        svtype = "Divergent"
                elif reflen > nonreflen:
                    svtype = "AltDeletion"
                else:
                    svtype = "AltInsertion"
        return chromo, pos, svtype, reflen, nonreflen, sourcenode, bubref, bubnonref, sinknode
    else:
        return None

nodeinf = {}
with open(sys.argv[1], 'r') as infile:
    for line in infile:
        temp = line.strip().split()
        if len(temp) > 2:
            nodeid, nodelen, chromo, pos, rrank = temp
            nodeinf[nodeid] = [int(rrank), int(nodelen)]
        else:
            nodeid, nodelen = temp
            # all the rrank of nonfef node were set to 100
            nodeinf[nodeid] = [100, int(nodelen)]

OUT = open(sys.argv[3], 'w')
OUT.write("\t".join(["chro", "pos", "svtype", "reflen", "nonreflen", "sourcenode", "bubref", "bubnonref", "sinknode"]) + "\n")
with open(sys.argv[2], 'r') as infile:
    for line in infile:
        if sv_det(line):
            OUT.write("\t".join(map(str, sv_det(line))) + "\n")
