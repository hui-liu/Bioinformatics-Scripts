import sys
from decimal import Decimal

"""
input:
    qseqid qlen sseqid slen frames pident nident length mismatch gapopen qstart qend sstart send evalue bitscore
output:
    qseqid sseqid q_txon s_taxon evlaue_mant evalue_exp pident percent_match

percent_match is computed by counting the length of the shorter sequence (query or subject) and dividing by the length of *that shorter sequence*
"""


def splitEvalue(string):
    mant = None
    exp = None
    if "0.0" == string:
        mant = string.split(".")[0]
        exp = string.split(".")[1]
    else:
        mant = string.split("e")[0]
        exp = string.split("e")[1]
    return mant, exp

def formatDigit(value):
    if value.split(".")[1] == "0":
        value = value.split(".")[0]
    else:
        value = value
    return value


IN = open(sys.argv[1], 'r')
OUT = open(sys.argv[2], 'w')

for line in IN:
    lsplit = line.split()
    if lsplit[0] == lsplit[2]: continue
    qlen, slen = int(lsplit[1]), int(lsplit[3])
    qstar, qend = [int(lsplit[10]), int(lsplit[11])] if int(lsplit[10]) < int(lsplit[11]) else [int(lsplit[11]), int(lsplit[10])]
    sstar, send = [int(lsplit[12]), int(lsplit[13])] if int(lsplit[12]) < int(lsplit[13]) else [int(lsplit[13]), int(lsplit[12])]
    mant, exp = splitEvalue(lsplit[14])
    pident = str(Decimal(lsplit[5]).quantize(Decimal('0.0')))
    if qlen >= slen:
        match = send - sstar + 1
        percent_match = str(Decimal(str(match / float(slen) * 100)).quantize(Decimal('0.0')))
        OUT.write("\t".join([lsplit[0], lsplit[2], lsplit[0].split("|")[0], lsplit[2].split("|")[0], mant, exp, formatDigit(pident), formatDigit(percent_match)]) + "\n")
    else:
        match = qend - qstar + 1
        percent_match = str(Decimal(str(match / float(qlen) * 100)).quantize(Decimal('0.0')))
        OUT.write("\t".join([lsplit[0], lsplit[2], lsplit[0].split("|")[0], lsplit[2].split("|")[0], mant, exp, formatDigit(pident), formatDigit(percent_match)]) + "\n")

IN.close()
OUT.close()
