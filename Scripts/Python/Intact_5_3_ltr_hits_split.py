import sys
import os

# python /media/40TB/Liuhui/Xsorbifolium/LTR/bin/Intact_LTR_not_overlap_hit.py \
# Xs_intact_5_3_ltr.outfmt6 \
# /media/12TB/Xso_genome/genome/genome2.fasta.fai \
# /media/40TB/Liuhui/Xsorbifolium/LTR/clustering/Xs_intact_5ltr_clustering.txt \
# Intact_LTR_not_overlap_hit.txt Idmapping.txt


# add large chromosomes or scaffolds to a list
chroms = []
with open(sys.argv[1], 'r') as f:
    for line in f:
        lsp = line.split()
        chroms.append(lsp[0])

# temporary directory
directory = "TEMP"
if not os.path.exists(directory):
    os.makedirs(directory)

# commands used for parallel running
commands = set()
str1 = "python /media/40TB/Liuhui/Xsorbifolium/LTR/bin/Intact_LTR_not_overlap_hit.py "
str2 = " /media/12TB/Xso_genome/genome/genome2.fasta.fai /media/40TB/Liuhui/Xsorbifolium/LTR/clustering/Xs_intact_5ltr_clustering.txt "
str3 = "_intact_LTR_not_overlap_hit.txt"
str4 = "_idmapping.txt"
str5 = None

# processing the large blast output file
with open(sys.argv[2], 'r') as f:
    for line in f:
        qseqid, sseqid, pident, qcovs, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = line.split()
        strand = "+"
        if int(sstart) > int(send):
            sstart, send = send, sstart
            strand = "-"
        if sseqid in chroms:
            if strand == "+":
                str5 = "plus"
                out = directory + "/" + sseqid + "_plus.outfmt6"
                with open(out, 'a') as f:
                    f.write(line)
                commands.add(str1 + out + str2 + directory + "/" + sseqid + "_" + str5 + str3 + " " + directory + "/" + sseqid + "_" + str5 + str4)
            else:
                str5 = "minus"
                out = directory + "/" + sseqid + "_minus.outfmt6"
                with open(out, 'a') as f:
                    f.write(line)
                commands.add(str1 + out + str2 + directory + "/" + sseqid + "_" + str5 + str3 + " " + directory + "/" + sseqid + "_" + str5 + str4)
        else:
            if strand == "+":
                str5 = "plus"
                out = directory + "/" + "unchrom" + "_plus.outfmt6"
                with open(out, 'a') as f:
                    f.write(line)
                commands.add(str1 + out + str2 + directory + "/" + "unchrom" + "_" + str5 + str3 + " " + directory + "/" + "unchrom" + "_" + str5 + str4)
            else:
                str5 = "minus"
                out = directory + "/" + "unchrom" + "_minus.outfmt6"
                with open(out, 'a') as f:
                    f.write(line)
                commands.add(str1 + out + str2 + directory + "/" + "unchrom" + "_" + str5 + str3 + " " + directory + "/" + "unchrom" + "_" + str5 + str4)

# write out the commoands
run_commands = "run_commands.txt"
command_out = open(run_commands, 'w')
for i in commands:
    command_out.write(i + "\n")
command_out.close()

# set the number of cpus
CPU = sys.argv[3]

print "running..."

os.system('ParaFly -c %s -CPU %s' % (run_commands, CPU))

# merging the results
print "combining..."

# Intact_LTR_not_overlap_hit
out1 = open(sys.argv[4], 'w')

# Idmapping
out2 = open(sys.argv[5], 'w')

for c in chroms + ["unchrom"]:
    for s in ["plus", "minus"]:
        temp_str1 = directory + "/" + c + "_" + s + str3
        temp_str2 = directory + "/" + c + "_" + s + str4
        with open(temp_str1, 'r') as f1:
            content1 = f1.read()
            out1.write(content1)
        with open(temp_str2, 'r') as f2:
            content2 = f2.read()
            out2.write(content2)
out1.close()
out2.close()

os.system('rm %s %s.completed' % (run_commands, run_commands))
os.system('rm -r %s' % (directory))
print "finished..."
