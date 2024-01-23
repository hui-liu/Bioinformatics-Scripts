import sys
def main(inDS='sample_list', inTab='st_out/%s/gene_exp.tab', outMatrix=sys.stdout, exp_col=9):
	exp_col = exp_col-1 # 0-based
	samples = [line.strip().split()[0] for line in open(inDS)]
	d_exp = {}
	for sample in samples:
		expTab = inTab % sample
                with open(expTab, 'r') as f:
                	header = f.readline()
			for line in f:
				temp = line.strip().split('\t')
				gene_id = temp[0]
				exp_value = float(temp[exp_col])
				try:
					d_exp[gene_id][sample] = exp_value
				except KeyError:
					d_exp[gene_id] = {sample:exp_value}

	Hline  = ['gene_id'] + ["_".join(i.split("_")[:2]) for i in samples]
	print >>outMatrix, '\t'.join(Hline)
	for gene_id, d_values in d_exp.items():
		values = [d_values[sample] for sample in samples]
                if sum(values) > 0:
			line = [gene_id] + map(str, values)
			print >>outMatrix, '\t'.join(line)
if __name__ == '__main__':
	main()
