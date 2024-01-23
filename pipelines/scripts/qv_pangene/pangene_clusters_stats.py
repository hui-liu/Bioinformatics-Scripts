import sys

core_gr = set()
core_gene = {}
shell_gr = set()
shell_gene = {}
acc_spec_gr = set()
acc_spec_gene = {}
total_genes = 0
n = int(sys.argv[3])

sp_gene = {}
with open(sys.argv[1], 'r') as f:
    for line in f:
        lsp = line.split()
        sp_gene[lsp[1]] = lsp[0]

with open(sys.argv[2], 'r') as f:
    for line in f:
        lsp = line.split()
        og = lsp[0].rstrip(":")
        genes = lsp[1:]
        sp_num = len(set([sp_gene[x] for x in genes]))
        if sp_num == n:
            core_gr.add(og)
            for i in genes:
                sp = sp_gene[i]
                core_gene.setdefault(sp, []).append(i)
                total_genes += 1
        elif 2 <= sp_num < n:
            shell_gr.add(og)
            for i in genes:
                sp = sp_gene[i]
                shell_gene.setdefault(sp, []).append(i)
                total_genes += 1
        else:
            acc_spec_gr.add(og)
            for i in genes:
                sp = sp_gene[i]
                acc_spec_gene.setdefault(sp, []).append(i)
                total_genes += 1

total_gr = len(core_gr) + len(shell_gr) + len(acc_spec_gr)
print "total species: %s (%s)" % (n, ",".join(core_gene.keys()))
print "total clusters:", total_gr
print "total genes:", total_genes

print "\ncore clusters (%s spceies):" % (n)
print len(core_gr), str(round(len(core_gr) * 1.0 / total_gr * 100, 2)) + "%"

print "\ncore genes:"
total_core_gene = 0
for i in core_gene:
    print i, len(core_gene[i])
    total_core_gene += len(core_gene[i])
print "total core genes:", total_core_gene, str(round(total_core_gene * 1.0 / total_genes * 100, 2)) + "%"

print "\nshell clusters (2-%s species):" % (n-1)
print len(shell_gr), str(round(len(shell_gr) * 1.0 / total_gr * 100, 2)) + "%"

print "\nshell genes:"
total_shell_gene = 0
for i in shell_gene:
    print i, len(shell_gene[i])
    total_shell_gene += len(shell_gene[i])
print "total shell genes:", total_shell_gene, str(round(total_shell_gene * 1.0 / total_genes * 100, 2)) + "%"

print "\naccession-specific clusters (many clusters contained only 1 gene):"
print len(acc_spec_gr), str(round(len(acc_spec_gr) * 1.0 / total_gr * 100, 2)) + "%"

print "\naccession-specific genes:"
total_acc_spec_gene = 0
for i in acc_spec_gene:
    print i, len(acc_spec_gene[i])
    total_acc_spec_gene += len(acc_spec_gene[i])
print "total accession-specific genes:",total_acc_spec_gene, str(round(total_acc_spec_gene * 1.0 / total_genes * 100, 2)) + "%"
