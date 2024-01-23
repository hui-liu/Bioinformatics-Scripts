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
                fas[id].append(line.rstrip().upper())
        for id, seq in fas.iteritems():
            fas[id] = ''.join(seq)
    return fas

def merge_intervals(intervals):
    # https://stackoverflow.com/questions/43600878/merging-overlapping-intervals/43600953
    intervals.sort(key=lambda x: x[0])
    merged = [intervals[0]]
    for current in intervals:
        previous = merged[-1]
        if current[0] <= previous[1]:
            previous[1] = max(previous[1], current[1])
            #merged[-1][1] = max(previous[1], current[1])
        else:
            merged.append(current)
    return merged

repeat = {}
TE = {}
LTR = {}
Copia = {}
Gypsy = {}
DNA = {}
TIR = {}
Helitron = {}
Unclassified = {}
Simple_repeat = {}
Low_complexity = {}

with open(sys.argv[1], 'r') as f:
    for line in f:
        if line[0] == "#": continue
        lsp = line.rstrip().split("\t")
        chr,start,end = lsp[0], lsp[3], lsp[4]
        if int(start) > int(end):
            start, end = end, start
        repeat.setdefault(chr, []).append([int(start), int(end)])
        if "EDTA" in lsp[1]:
            TE.setdefault(chr, []).append([int(start), int(end)])
        if "LTR" in lsp[8]:
            LTR.setdefault(chr, []).append([int(start), int(end)])
        if "Copia" in lsp[8]:
            Copia.setdefault(chr, []).append([int(start), int(end)])
        if "Gypsy" in lsp[8]:
            Gypsy.setdefault(chr, []).append([int(start), int(end)])
        if "DNA" in lsp[8] or "MITE" in lsp[8] or "polinton" in lsp[2]:
            DNA.setdefault(chr, []).append([int(start), int(end)])
        if "TIR" in lsp[2]:
           TIR.setdefault(chr, []).append([int(start), int(end)])
        if "Helitron" in lsp[8]:
           Helitron.setdefault(chr, []).append([int(start), int(end)])
        if "=Unknown;" in lsp[8]:
           Unclassified.setdefault(chr, []).append([int(start), int(end)])
        if "Simple_repeat" in lsp[8]:
           Simple_repeat.setdefault(chr, []).append([int(start), int(end)])
        if "Low_complexity" in lsp[8]:
           Low_complexity.setdefault(chr, []).append([int(start), int(end)])


len_all = []
len_TE = []
len_LTR = []
len_Copia = []
len_Gypsy = []
len_DNA = []
len_TIR = []
len_Helitron = []
len_Unclassified = []
len_Simple_repeat = []
len_Low_complexity = []

for chr in repeat:
    merged_reg = merge_intervals(repeat[chr])
    length = sum([i[1] - i[0] + 1 for i in merged_reg])
    len_all.append(length)

for chr in TE:
    merged_reg = merge_intervals(TE[chr])
    length = sum([i[1] - i[0] + 1 for i in merged_reg])
    len_TE.append(length)

for chr in LTR:
    merged_reg = merge_intervals(LTR[chr])
    length = sum([i[1] - i[0] + 1 for i in merged_reg])
    len_LTR.append(length)

for chr in Copia:
    merged_reg = merge_intervals(Copia[chr])
    length = sum([i[1] - i[0] + 1 for i in merged_reg])
    len_Copia.append(length)

for chr in Gypsy:
    merged_reg = merge_intervals(Gypsy[chr])
    length = sum([i[1] - i[0] + 1 for i in merged_reg])
    len_Gypsy.append(length)

for chr in DNA:
    merged_reg = merge_intervals(DNA[chr])
    length = sum([i[1] - i[0] + 1 for i in merged_reg])
    len_DNA.append(length)

for chr in TIR:
    merged_reg = merge_intervals(TIR[chr])
    length = sum([i[1] - i[0] + 1 for i in merged_reg])
    len_TIR.append(length)

for chr in Helitron:
    merged_reg = merge_intervals(Helitron[chr])
    length = sum([i[1] - i[0] + 1 for i in merged_reg])
    len_Helitron.append(length)

for chr in Unclassified:
    merged_reg = merge_intervals(Unclassified[chr])
    length = sum([i[1] - i[0] + 1 for i in merged_reg])
    len_Unclassified.append(length)

for chr in Simple_repeat:
    merged_reg = merge_intervals(Simple_repeat[chr])
    length = sum([i[1] - i[0] + 1 for i in merged_reg])
    len_Simple_repeat.append(length)

for chr in Low_complexity:
    merged_reg = merge_intervals(Low_complexity[chr])
    length = sum([i[1] - i[0] + 1 for i in merged_reg])
    len_Low_complexity.append(length)

size  = int(sys.argv[2])

print("Total", sum(len_all), round(sum(len_all) * 100.0 / size, 2))
print("TE", sum(len_TE), round(sum(len_TE) * 100.0 / size, 2))
print("LTR", sum(len_LTR), round(sum(len_LTR) * 100.0 / size, 2))
print("Copia", sum(len_Copia), round(sum(len_Copia) * 100.0 / size, 2))
print("Gypsy", sum(len_Gypsy), round(sum(len_Gypsy) * 100.0 / size, 2))
print("DNA", sum(len_DNA), round(sum(len_DNA) * 100.0 / size, 2))
print("TIR", sum(len_TIR), round(sum(len_TIR) * 100.0 / size, 2))
print("Helitron", sum(len_Helitron), round(sum(len_Helitron) * 100.0 / size, 2))
print("Unclassified", sum(len_Unclassified), round(sum(len_Unclassified) * 100.0 / size, 2))
print("Simple_repeat", sum(len_Simple_repeat), round(sum(len_Simple_repeat) * 100.0 / size, 2))
print("Low_complexity", sum(len_Low_complexity), round(sum(len_Low_complexity) * 100.0 / size, 2))

