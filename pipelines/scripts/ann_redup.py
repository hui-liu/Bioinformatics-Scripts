import sys
from operator import itemgetter, attrgetter

def main(inGff3, outGff3):
#	d_gff = {}
	d_record = {}
	d_parent = {}
	genes = []
	i = 0
	for line in open(inGff3):
		if line.startswith('#'):
			continue
		i += 1
		temp = line.rstrip().split('\t')
		try: chr, source, type, start, end, score, strand, phase, attibutes = temp
		except : print temp; exit
		start, end = int(start), int(end)
		try: d_attr = dict([v.split('=') for v in attibutes.strip(';').split(';')])
		except ValueError: print temp; exit
		if type == 'gene':
			gene_id = d_attr['ID']
			gene_length = end-start+1
			key = (chr,start,end,strand,gene_id,i,gene_length)
			genes.append(key)
			d_record[gene_id] = [line]
		else:
			parent = p_id = d_attr['Parent']
			feat_id = d_attr['ID']
			d_parent[feat_id] = parent
			if parent in d_parent:
				d_parent[feat_id] = d_parent[parent]
			gene_id = d_parent[feat_id]
			d_record[gene_id] += [line]
	
	def _overlap(a, b):
		return max(0, min(a[1], b[1]) - max(a[0], b[0]))
#	genes = sorted(genes, key=itemgetter(3,0,1))
#	print genes[:20]
	def _remove_overlap(genes):
		print len(genes)
		genes = sorted(list(set(genes)), key=itemgetter(3,0,1,5))
		print len(genes)
		k0 = genes[0]
		discards = set([])
		un_idents = []
#		new_genes = []
		for i in range(1, len(genes)):
			k1 = genes[i]
			if not k1[0] == k0[0]:
#				new_genes += [k0,k1]
				k0 = k1
				continue
			r0 = k0[1:3]
			r1 = k1[1:3]
			if _overlap(r0, r1):
#				print k0,k1
				len0, len1 = k0[-1], k1[-1]
				if 100.0*_overlap(r0, r1)/min(len0, len1) < 10:		# too short overlap coverage, not discard
					k0 = k1
					continue
				if r0 == r1:
					pass
				else:
					un_idents.append([k0,k1])
#				print [k0,k1]
				if len0 >= len1:
					print [k0,k1]
					discards.add(k1)
#					new_genes += [k0]
					k0 = k0
				else:
					print [k1,k0]
					discards.add(k0)
#					new_genes += [k1]
					k0 = k1
			else:
#				new_genes += [k0, k1]
				k0 = k1
#		for v in un_idents:
#			print v
		print len(discards)
#		if ('chr_14', 4033441, 4037788, '+', 'augustus-chr_14:3600001-4600000-processed-gene-4.2', 4348) in discards:
#			print ('chr_14', 4033441, 4037788, '+', 'augustus-chr_14:3600001-4600000-processed-gene-4.2', 4348)
#		new_genes = sorted(list(set(new_genes)), key=itemgetter(3,0,1))
		new_genes = list(set(genes) - discards)
		if len(discards) == 0:
			return new_genes
		else:
			return _remove_overlap(new_genes)
#		return new_genes
	genes = _remove_overlap(genes)
	f = open(outGff3,'w')
	for k in genes:
#		if k in discards:
#			continue
		gene_id= k[4]
		lines = ''.join(d_record[gene_id])
		f.write(lines)
	f.close()

if __name__ == '__main__':
	inGff3, outGff3 = sys.argv[1:3]
	main(inGff3, outGff3)
