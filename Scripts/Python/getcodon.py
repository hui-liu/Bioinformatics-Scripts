def locateBin(a, b):
    result = True
    if a < b[0] or b[1] < a:
        result = False
    return result


def getIndex(pos, cds_list):
    for ith, elment in enumerate(cds_list):
        if locateBin(pos, elment):
            return ith

def getcodon(seq, pos, cds_list):
    """
    seq: AGTTATTGATTTTGCTTTTTGCTTGTA
    pos: 10
    cds_list: [[1,8,'+'], [9,20,'+'], [21,27,'+']]
    return: (TTT, 0)
    """
    strand = cds_list[0][2]
    if strand == "+":
        ith_codon = [2, 0, 1]
        index = getIndex(pos, cds_list)
        target_lst = cds_list[0:index] + [[cds_list[index][0], pos]]
        L = 0
        for i in target_lst:
            l = i[1] - i[0] + 1
            L += l
        a = [i for i in range(0, L, 3)]
        b = [a[-1], a[-1]+2]
        return seq[b[0]: b[1]+1], ith_codon[L%3]
    else:
        # [[21, 27, '-'], [9, 20, '-'], [1, 8, '-']]
        cds_list_r = [i[:-1] for i in cds_list[::-1]]
        ith_codon = [0, 2, 1]
        index = getIndex(pos, cds_list_r)
        target_lst = cds_list_r[0:index] + [[pos, cds_list_r[index][1]]]
        print target_lst
        L = 0
        for i in target_lst:
            l = i[1] - i[0] + 1
            L += l
        print L
        a = [i for i in range(0, L, 3)]
        b = [-a[-1]-2, -a[-1]]
        if b[1] == 0:
            return seq[b[0]-1:], ith_codon[L%3]
        else:
            return seq[b[0]-1: b[1]], ith_codon[L%3]

pos = 10
#print getcodon('AGTTATTGATTTTGCTTTTTGCTTGTA', 25, [[1,8,'+'], [9,20,'+'], [21,27,'+']])
print getcodon('AGTTATTGATTTTGCTTTTTGCTTGTA', pos, [[1,8,'-'], [9,20,'-'], [21,27,'-']])
