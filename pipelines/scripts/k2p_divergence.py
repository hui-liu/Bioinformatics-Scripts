import sys
import os
import random

USAGE = "\nusage: python %s pta_ltrdigest_5ltr.fas pta_ltrdigest_3ltr.fas outfile\n" % sys.argv[0]

if len(sys.argv) != 4:
    print USAGE
    sys.exit()


def parseFasta(filename):
    fas = {}
    id = None
    with open(filename, 'r') as fh:
        for line in fh:
            if line[0] == '>':
                header = line[1:].rstrip()
                id = header.split()[0]
                fas[id] = []
            else:
                fas[id].append(line.rstrip())
        for id, seq in fas.iteritems():
            fas[id] = ''.join(seq)
    return fas
    
def K2Pdistance(seq1,seq2):
    """
    Kimura 2-Parameter distance = -(1/2) * ln(1 - 2P - Q) - (1/4) * ln(1 - 2Q)
    Kimura 2-Parameter distance = -0.5 log( (1 - 2p -q) * sqrt( 1 - 2q ) )
    where:
    p = transition frequency
    q = transversion frequency
    """
    from math import log, sqrt
    transitions = [ "AG", "GA", "CT", "TC"]
    transversions = [ "AC", "CA", "AT", "TA",
                      "GC", "CG", "GT", "TG" ]
    condon_pairs = ["AA", "CC", "GG", "TT"] + transitions + transversions

    # collect valid pairs
    seq1 = seq1.upper()
    seq2 = seq2.upper()

    # delete gap
    pairs = [x+y for (x,y) in zip(seq1,seq2) if x+y in condon_pairs]
    ts_count = 0
    tv_count = 0
    length = float(len(pairs))

    # calculate transition rate (p) and transversion rate (q)
    for i in pairs:
        if i in transitions:
            ts_count += 1
        elif i in transversions:
            tv_count += 1
    p = ts_count / length
    q = tv_count / length

    try: d = -0.5 * log( (1 - 2*p - q) * sqrt( 1 - 2*q ) )
    except ValueError:
        d = "NA"
    return d


fiveltr_dict = parseFasta(sys.argv[1])
threeltr_dict = parseFasta(sys.argv[2])

OUT = open(sys.argv[3], 'w')
for id in fiveltr_dict:
    fiveltrseq = fiveltr_dict[id]
    threeltrseq = threeltr_dict[id]
    idlis = [id + "_5", id + "_3"]
    Rnum = str(random.randint(1,1000000))
    Aln = Rnum + "_" + "aln"
    with open(Rnum, "w") as f:
       f.write(">" + idlis[0] + "\n" + fiveltrseq + "\n" + ">" + idlis[1] + "\n" + threeltrseq + "\n")
    os.system('mafft --auto --quiet --preservecase %s > %s' % (Rnum, Aln))
    tmp_dict = parseFasta(Aln)
    D = str(K2Pdistance(tmp_dict[idlis[0]], tmp_dict[idlis[1]]))
    if "-0.0" == D:
        D = "0.0"
    OUT.write(id + "\t" + D + "\n")
    os.remove(Rnum)
    os.remove(Aln)

OUT.close()
