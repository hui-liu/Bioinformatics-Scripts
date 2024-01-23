import sys

def main(inGff3=sys.argv[1], outGtf=sys.stdout):
	d_rna = {}
	for line in open(inGff3):
		if line.startswith('#'):
			continue
		temp = line.rstrip().split('\t')
		chrom, source, feat_type, start, end, score, strand, frame, info = temp
		d_info = dict([v.split('=') for v in info.split(';')])
		if feat_type == 'gene':
			gene_id = d_info['ID']
			temp[8] = 'gene_id "%s";' % (gene_id, )
		elif 'RNA' in feat_type or 'repeat' in feat_type or 'gene' in feat_type:
			gene_id = d_info['Parent']
			rna_id = d_info['ID']
			d_rna[rna_id] = gene_id
			temp[2] = 'transcript'
			temp[8] = 'gene_id "%s"; transcript_id "%s";' % (gene_id, rna_id)
		elif feat_type == 'exon':
			rna_id = d_info['Parent']
			gene_id = d_rna[rna_id]
			temp[8] = 'gene_id "%s"; transcript_id "%s";' % (gene_id, rna_id)
		else:
			continue
		print >> outGtf, '\t'.join(temp)

if __name__ == '__main__':
	main()
