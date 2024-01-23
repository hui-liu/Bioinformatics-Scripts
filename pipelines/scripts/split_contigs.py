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
                fas[id].append(line.rstrip())
        for id, seq in fas.items():
            fas[id] = ''.join(seq)
    return fas

seq_dict = parseFasta(sys.argv[1])
size_cutoff = int(sys.argv[2])
for id in seq_dict:
    seq = seq_dict[id]
    size = len(seq)
    if size >= size_cutoff:
        print id
        with open(id + ".fa", 'w') as fw:
            fw.write(">" + id + "\n" + seq + "\n")
    else:
        with open("Contigs.fa", 'a') as fw:
            fw.write(">" + id + "\n" + seq + "\n")

print "Contigs"
