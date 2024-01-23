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

def get_length(coor_dict):
    size = []
    num = 0
    for chr in coor_dict:
        num += len(coor_dict[chr])
        merged_reg = merge_intervals(coor_dict[chr])
        length = sum([i[1] - i[0] + 1 for i in merged_reg])
        size.append(length)
    return size, num

repeat = {}
TE = {}
LTR = {}
Copia = {}
Gypsy = {}
LTR_other = {}
DNA = {}
TIR = {}
Helitron = {}
Unclassified = {}
TandemRepeat = {}

Mutator_TIR_transposon = {}
CACTA_TIR_transposon = {}
hAT_TIR_transposon = {}
PIF_Harbinger_TIR_transposon = {}
Tc1_Mariner_TIR_transposon = {}
LINE_element = {}
non_LTR_retrotransposon = {} # Non-LTR/Others
all_non_LTR = {}
polinton = {}
Retro = {}

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
            #Retro.setdefault(chr, []).append([int(start), int(end)])
        if "Copia" in lsp[8]:
            Copia.setdefault(chr, []).append([int(start), int(end)])
            Retro.setdefault(chr, []).append([int(start), int(end)])
        if "Gypsy" in lsp[8]:
            Gypsy.setdefault(chr, []).append([int(start), int(end)])
            Retro.setdefault(chr, []).append([int(start), int(end)])
        if "LTR" in lsp[8] and "Copia" not in lsp[8] and "Gypsy" not in lsp[8]:
            LTR_other.setdefault(chr, []).append([int(start), int(end)])
            Retro.setdefault(chr, []).append([int(start), int(end)])
        if "DNA" in lsp[8] or "MITE" in lsp[8] or "polinton" in lsp[2]:
            DNA.setdefault(chr, []).append([int(start), int(end)])
        if "TIR" in lsp[2]:
            TIR.setdefault(chr, []).append([int(start), int(end)])
        if "Helitron" in lsp[8]:
            Helitron.setdefault(chr, []).append([int(start), int(end)])
        if "=Unknown;" in lsp[8]:
            Unclassified.setdefault(chr, []).append([int(start), int(end)])
        if "TandemRepeat" in lsp[8]:
            TandemRepeat.setdefault(chr, []).append([int(start), int(end)])
        if "Mutator_TIR_transposon" == lsp[2]:
            Mutator_TIR_transposon.setdefault(chr, []).append([int(start), int(end)])
        if "CACTA_TIR_transposon" == lsp[2]:
            CACTA_TIR_transposon.setdefault(chr, []).append([int(start), int(end)])
        if "hAT_TIR_transposon" == lsp[2]:
            hAT_TIR_transposon.setdefault(chr, []).append([int(start), int(end)])
        if "PIF_Harbinger_TIR_transposon" == lsp[2]:
            PIF_Harbinger_TIR_transposon.setdefault(chr, []).append([int(start), int(end)])
        if "Tc1_Mariner_TIR_transposon" == lsp[2]:
            Tc1_Mariner_TIR_transposon.setdefault(chr, []).append([int(start), int(end)])
        if "LINE_element" == lsp[2]:
            LINE_element.setdefault(chr, []).append([int(start), int(end)])
            all_non_LTR.setdefault(chr, []).append([int(start), int(end)])
            Retro.setdefault(chr, []).append([int(start), int(end)])
        if "non_LTR_retrotransposon" == lsp[2]:
            non_LTR_retrotransposon.setdefault(chr, []).append([int(start), int(end)])
            all_non_LTR.setdefault(chr, []).append([int(start), int(end)])
            Retro.setdefault(chr, []).append([int(start), int(end)])
        if "polinton" == lsp[2]:
            polinton.setdefault(chr, []).append([int(start), int(end)])

len_all, num_all = get_length(repeat)
len_TE, num_TE = get_length(TE)
len_LTR, num_LTR = get_length(LTR)
len_Copia, num_Copia = get_length(Copia)
len_Gypsy, num_Gypsy = get_length(Gypsy)
len_LTR_other, num_LTR_other = get_length(LTR_other)
len_DNA, num_DNA = get_length(DNA)
len_TIR, num_TIR = get_length(TIR)
len_Helitron, num_Helitron = get_length(Helitron)
len_Unclassified, num_Unclassified = get_length(Unclassified)
len_TandemRepeat, num_TandemRepeat = get_length(TandemRepeat)
len_Mutator_TIR_transposon, num_Mutator_TIR_transposon = get_length(Mutator_TIR_transposon)
len_CACTA_TIR_transposon, num_CACTA_TIR_transposon = get_length(CACTA_TIR_transposon)
len_hAT_TIR_transposon, num_hAT_TIR_transposon = get_length(hAT_TIR_transposon)
len_PIF_Harbinger_TIR_transposon, num_PIF_Harbinger_TIR_transposon = get_length(PIF_Harbinger_TIR_transposon)
len_Tc1_Mariner_TIR_transposon, num_Tc1_Mariner_TIR_transposon = get_length(Tc1_Mariner_TIR_transposon)
len_LINE_element, num_LINE_element = get_length(LINE_element)
len_non_LTR_retrotransposon, num_non_LTR_retrotransposon = get_length(non_LTR_retrotransposon)
len_all_non_LTR, num_all_non_LTR = get_length(all_non_LTR)
len_polinton, num_polinton = get_length(polinton)
len_Retro, num_Retro = get_length(Retro)

size  = int(sys.argv[2])

print("Total", sum(len_all), num_all, round(sum(len_all) * 100.0 / size, 2))
print("TE", sum(len_TE), num_TE, round(sum(len_TE) * 100.0 / size, 2))
print("LTR", sum(len_LTR), num_LTR, round(sum(len_LTR) * 100.0 / size, 2))
print("Copia", sum(len_Copia), num_Copia, round(sum(len_Copia) * 100.0 / size, 2))
print("Gypsy", sum(len_Gypsy), num_Gypsy, round(sum(len_Gypsy) * 100.0 / size, 2))
print("LTR_other", sum(len_LTR_other), num_LTR_other, round(sum(len_LTR_other) * 100.0 / size, 2))

print("DNA", sum(len_DNA), num_DNA, round(sum(len_DNA) * 100.0 / size, 2))

print("TIR", sum(len_TIR), num_TIR, round(sum(len_TIR) * 100.0 / size, 2))
print("CACTA_TIR_transposon", sum(len_CACTA_TIR_transposon),num_CACTA_TIR_transposon, round(sum(len_CACTA_TIR_transposon) * 100.0 / size, 2))
print("Mutator_TIR_transposon", sum(len_Mutator_TIR_transposon),num_Mutator_TIR_transposon, round(sum(len_Mutator_TIR_transposon) * 100.0 / size, 2))
print("hAT_TIR_transposon", sum(len_hAT_TIR_transposon), num_hAT_TIR_transposon, round(sum(len_hAT_TIR_transposon) * 100.0 / size, 2))
print("PIF_Harbinger_TIR_transposon", sum(len_PIF_Harbinger_TIR_transposon),num_PIF_Harbinger_TIR_transposon, round(sum(len_PIF_Harbinger_TIR_transposon) * 100.0 / size, 2))
print("Tc1_Mariner_TIR_transposon", sum(len_Tc1_Mariner_TIR_transposon),num_Tc1_Mariner_TIR_transposon, round(sum(len_Tc1_Mariner_TIR_transposon) * 100.0 / size, 2))
print("polinton", sum(len_polinton), num_polinton, round(sum(len_polinton) * 100.0 / size, 2))

print("non_LTR_retrotransposon", sum(len_non_LTR_retrotransposon), num_non_LTR_retrotransposon, round(sum(len_non_LTR_retrotransposon) * 100.0 / size, 2))
print("LINE_element", sum(len_LINE_element), num_LINE_element, round(sum(len_LINE_element) * 100.0 / size, 2))
print("all_non_LTR",  sum(len_all_non_LTR), num_all_non_LTR, round(sum(len_all_non_LTR) * 100.0 / size, 2))
print("Helitron", sum(len_Helitron), num_Helitron, round(sum(len_Helitron) * 100.0 / size, 2))
print("Unclassified", sum(len_Unclassified), num_Unclassified, round(sum(len_Unclassified) * 100.0 / size, 2))
print("TandemRepeat", sum(len_TandemRepeat), num_TandemRepeat, round(sum(len_TandemRepeat) * 100.0 / size, 2))

print("Retrotransposon", sum(len_Retro), num_Retro, round(sum(len_Retro) * 100.0 / size, 2))

