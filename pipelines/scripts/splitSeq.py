import sys
import os
import shutil
from collections import OrderedDict

def parseFasta(filename):
    fas = OrderedDict()
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


inFa = sys.argv[1]
outDir = sys.argv[2]
try: partNo = int(sys.argv[3])
except: partNo = 16

def mknewDir(Dir):
	if os.path.exists(Dir):
#		shutil.rmtree(Dir)
#		os.makedirs(Dir)
		pass
	else:
		os.makedirs(Dir)

def split_seq(inFa,outDir,split_part):
	i = 0
        seq_dcit = parseFasta(inFa)
	for record in seq_dcit:
		i += 1
		out_file = '%s/%0*d.fa' % (outDir, 3, i/split_part)
		f = open(out_file,'a')
                f.write(">" + record + "\n" + seq_dcit[record] + "\n")
		f.close()
	return i

def seq_num(inFa):
	return sum([1 for record in parseFasta(inFa)])

mknewDir(outDir)

totalSeq = seq_num(inFa)
split_part = totalSeq/partNo + 1
finalSeq = split_seq(inFa,outDir,split_part)
