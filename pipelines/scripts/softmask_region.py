import sys, os
from collections import Counter
from Bio import SeqIO
from Bio.Seq import Seq

def lower(base):
	return base.lower()

def main(inSeq=sys.argv[1], inBed=sys.argv[2], soft=True, suffix=None, outSeq=sys.stdout):
	if soft:
		print >>sys.stderr, 'soft mask'
	if suffix is None:
		suffix = os.path.splitext(inBed)[-1].lstrip('.')
	if suffix in {'bed'}:
		kw = dict(chr_col=0, start_col=1, end_col=2, based=0)
	elif suffix in {'gff3', 'gff', 'gff2', 'gtf'}:
		kw = dict(chr_col=0, start_col=3, end_col=4, based=1)
	else:
		raise ValueError('un-recagnized suffix {}'.format(suffix))
	d_bed = bed2dict(inBed, **kw)
	all_n = 0
	for rc in SeqIO.parse(inSeq, 'fasta'):
		if rc.id not in d_bed:
			SeqIO.write(rc, outSeq, 'fasta')
			continue
		regions = d_bed[rc.id]
#		segments = []
#		left = 0
		seq = list(rc.seq)
		for start, end in sorted(regions):
#			segments += [str(rc.seq[left:start])]
#			segments += ['n'*(end-start)]
#			left = end
			if soft:
				seq[start:end] = map(lower, seq[start:end])
			else:
				seq[start:end] = ['N'] * (end-start)
#		segments += [str(rc.seq[left:])]
#		seq = ''.join(segments)
		d_counter = Counter(list(seq))
		num_N = d_counter['n'] + d_counter['N']
		seq = ''.join(seq)
		seq = Seq(seq)
		if not len(seq) == len(rc.seq):
			print >> sys.stderr, '\tERROR', rc.id, len(seq), len(rc.seq), len(regions)
		print >> sys.stderr, rc.id, len(seq), num_N, 100.0*(len(seq)-num_N)/len(seq)
		if len(seq) == num_N:
			all_n += 1
		assert len(seq) == len(rc.seq)
		rc.seq = seq
		SeqIO.write(rc, outSeq, 'fasta')
	print >> sys.stderr, all_n, 'all N'

def bed2dict(inBed, chr_col=0, start_col=1, end_col=2, based=0):
	d_bed = {}
	for line in open(inBed):
#		if line.startswith('#'):
		temp = line.strip().split()
		try:
			CHR, START, END = temp[chr_col], temp[start_col], temp[end_col]
			START = int(START) - based
			END = int(END)
		except IndexError: continue
		except ValueError: continue
		try: d_bed[CHR] += [(START, END)]
		except KeyError: d_bed[CHR] = [(START, END)]
	return d_bed

if __name__ == '__main__':
	main(inSeq=sys.argv[1], inBed=sys.argv[2])
