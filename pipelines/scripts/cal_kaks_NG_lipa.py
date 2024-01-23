import sys
import os

def getGenePair(filename):
    p = []
    with open(filename, 'r') as f:
        for line in f:
            _, g1, g2 = line.split()
            p.append([g1, g2])
    return p

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
        for id, seq in fas.items():
            fas[id] = ''.join(seq)
    return fas

def pal2nal(pal, nuc_seqs):
    nal = {}
    nal_ = ''
    for pid, seq in pal.items():
        if pid not in nuc_seqs:
            print("{0} not in cds file".format(pid))
        else:
            nuc = nuc_seqs[pid]
            nal_ = ''
            j = 0
            for i in range(len(seq)):
                if seq[i] == '-':
                    nal_ += '---'
                else:
                    nal_ += nuc[j:j + 3]
                    j += 3
            nal[pid] = nal_
    return nal

def aln2axt(alndict, rnd_name, ID, out):
    outname = rnd_name + ".axt"
    with open(out + "/" + outname, 'w') as f:
        f.write(ID + "\n")
        for s in alndict:
            f.write(alndict[s] + "\n")
    return outname

def alignment(pep_dict, nuc_dict, pair, rnd_name, out):
    id = "_".join(pair)
    fasfile = rnd_name + ".fasta"
    pepaln = rnd_name  + "." + "aln"
    with open(out + "/" + fasfile, "w") as f:
       for s in pep_dict:
           f.write(">" + s + "\n" + pep_dict[s] + "\n")
    os.system('mafft --localpair --maxiterate 1000 --quiet --preservecase %s/%s > %s/%s' % (out, fasfile, out, pepaln))
    pepaln_dict = parseFasta(out + "/" + pepaln)
    codon_aln = pal2nal(pepaln_dict, nuc_dict)
    axt = aln2axt(codon_aln, rnd_name, id, out)
    return axt


genepair = getGenePair(sys.argv[1])
pep_dict = parseFasta(sys.argv[2])
nuc_dict = parseFasta(sys.argv[3])
tmpdir = sys.argv[4]

out = tmpdir + "/KaKs_Calculator/OUT"
if not os.path.exists(out):
    os.system('mkdir -p %s' % (out))

m = 0
for p in genepair:
    m += 1
    mm = "p" + str(m)
    pseqs_dict = {k:pep_dict[k] for k in p}
    nseqs_dict = {k:nuc_dict[k] for k in p}
    aln_name = alignment(pseqs_dict, nuc_dict, p, mm, out)
    os.system('KaKs -i %s/%s -o %s/%s.kaks -m NG > %s/%s.log' % (out, aln_name, out, aln_name, out, aln_name))
