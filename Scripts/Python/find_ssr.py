
import sys

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
                fas[id].append(line.rstrip().upper())
        for id, seq in fas.iteritems():
            fas[id] = ''.join(seq)
    return fas


def repeats(s, max_length, min_length=1):
    for i in xrange(len(s)):
        for j in xrange(min_length, max_length+1):
            count = 1
            while s[i:i+j] == s[i+j*count:i+j*count+j]: # find repeats
                count += 1
            if count > 1:
                yield s[i:i+j], i, count

USAGE = "\nusage: python  %s sequences.fasta out.txt\n" % sys.argv[0]
if len(sys.argv) != 3:
    print USAGE
    sys.exit()

faSeq = parseFasta(sys.argv[1])
OUT = open(sys.argv[2], 'w')

OUT.write("\t".join(["ID","start", "end", "unit","unit_size","unit_counts", "ssr_seq"]) + "\n")

for i in faSeq:
    pre_ssr_seq = ""
    for unit, position, count in repeats(faSeq[i], 6, 1):
        ssr_seq = faSeq[i][position:position + count*len(unit)]
        if len(unit) * count >= 6 and ssr_seq not in pre_ssr_seq:
            OUT.write("\t".join([i, str(position + 1), str(position + count*len(unit)), unit, str(len(unit)), str(count), ssr_seq]) + "\n")
        else: continue
        pre_ssr_seq = ssr_seq
OUT.close()
