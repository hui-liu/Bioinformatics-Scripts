import sys
import gzip
#import numpy as np

def decode(q):
    return ord(q) - 33

n = 0
leng = 0
GC_count = 0
AT_count = 0
Q20 = 0
Q30 = 0
for line in sys.stdin:
    n += 1
    if n % 4 == 2:
        seq = line.strip().upper()
        leng += len(seq)
        gc = seq.count("C") + seq.count("G")
        at = seq.count("A") + seq.count("T")
        GC_count += gc
        AT_count += at
    elif n % 4 == 0:
        values = [decode(x) for x in line.strip()]
        q20 = sum([1 for x in values if x >= 20])
        q30 = sum([1 for x in values if x >= 30])
        Q20 += q20
        Q30 += q30


print("Total number:", int(n/4))
print("Total length (Gb):", round(leng / 1000000000.0, 2))
print("GC (%):", round(GC_count * 1.0 / (GC_count + AT_count) * 100, 2))
print("Q20:", round(Q20 * 1.0 / leng * 100, 2))
print("Q30:", round(Q30 * 1.0 / leng * 100, 2))
