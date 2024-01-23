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

lst_len = []
i = 0
nA, nT, nC, nG, nN = 0, 0, 0, 0, 0

in_fa = sys.argv[1]
seq_dict = parseFasta(in_fa)
for record in seq_dict:
    seq = seq_dict[record]
    lst_len.append(len(seq))
    i += 1
    nA = nA + seq.upper().count('A')
    nT = nT + seq.upper().count('T')
    nC = nC + seq.upper().count('C')
    nG = nG + seq.upper().count('G')
    nN = nN + seq.upper().count('N')

lst_len = sorted(lst_len)
min_len = lst_len[0]
max_len = lst_len[len(lst_len)-1]
median_len = lst_len[int(len(lst_len)/2)]
mean_len = sum(lst_len)/len(lst_len)
up_quatile_len = lst_len[int(len(lst_len)*3/4)]
low_quatile_len = lst_len[int(len(lst_len)/4)]
quatile_90_len = lst_len[int(len(lst_len)*9/10)]
quatile_10_len = lst_len[int(len(lst_len)/10)]
gc = round(100.0*(nG + nC)/(nG + nC + nT + nA),2)
length = sum(lst_len)
others = length-(nG + nC + nT + nA + nN)
sum_len = 0
l = 0

for scaf_len in sorted(lst_len,reverse=1):
    sum_len += scaf_len
    l += 1
    if sum_len > 0.5*length:
        N50 = scaf_len
        L50 = l
        break

sum_len = 0
l = 0
for scaf_len in sorted(lst_len,reverse=1):
    sum_len += scaf_len
    l += 1
    if sum_len > 0.1*length:
        N10 = scaf_len
        L10 = l
        break
sum_len = 0
l = 0
for scaf_len in sorted(lst_len,reverse=1):
    sum_len += scaf_len
    l += 1
    if sum_len > 0.9*length:
        N90 = scaf_len
        L90 = l
        break

#print 'Minimum:\t%s\n10%%-tile:\t%s\n25%%-tile:\t%s\nMedian:\t%s\n75%%-tile:\t%s\n90%%-tile:\t%s\nMaximum:\t%s\nMean:\t%s\n# of seqs:\t%s' % (min_len, quatile_10_len, low_quatile_len, median_len, up_quatile_len, quatile_90_len, max_len, mean_len, i)
print 'genome size\t%s bp\ngenome size without N\t%s bp' % (length, length-nN)

print 'GC content\t%s%%\n  A\t%s\n  T\t%s\n  G\t%s\n  C\t%s\n  N\t%s\n  Others\t%s' % (gc, nA, nT, nC, nG, nN, others)

print '\nLongest Contig\t%s bp\n' % max_len
#for j in range(0,len(lst_len)):
#	print j,lst_len[j]
if i > 1 :
    print 'scaffold numbers\t%s' % i
    print '  Mean\t%d bp\t\tMedian\t%s bp' % (mean_len, median_len)
    print '  N10\t%s bp\t\tL10\t%s\n  N50\t%s bp\t\tL50\t%s\n  N90\t%s bp\t\tL90\t%s' % (N10,L10,N50,L50,N90,L90)
