import sys

LTR = [0, 0]
LTR_gypsy = [0, 0]
LTR_copia = [0, 0]
LTR_unknown = [0, 0]
with open(sys.argv[1], 'r') as f:
    for line in f:
        chrom, source, feat_type, start, end, score, strand, frame, info = line.rstrip().split("\t")
        d_info = dict([v.split('=') for v in info.split(';')])
        if feat_type == "repeat_region":
            LTR[0] += 1
            LTR[1] += int(end) - int(start) + 1
            if d_info['Classification'] == "LTR/Gypsy":
                LTR_gypsy[0] += 1
                LTR_gypsy[1] += int(end) - int(start) + 1
            elif d_info['Classification'] == "LTR/Copia":
                LTR_copia[0] += 1
                LTR_copia[1] += int(end) - int(start) + 1
            elif d_info['Classification'] == "LTR/unknown":
                LTR_unknown[0] += 1
                LTR_unknown[1] += int(end) - int(start) + 1

print "LTR:", "\t".join(map(str, LTR))
print "LTR/Gypsy:", "\t".join(map(str, LTR_gypsy))
print "LTR/Copia:", "\t".join(map(str, LTR_copia))
print "LTR/unknown:", "\t".join(map(str, LTR_unknown))

