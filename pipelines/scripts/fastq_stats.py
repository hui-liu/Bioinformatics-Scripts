import sys
import gzip

n = 0
leng = 0
GC_count = 0
AT_count = 0
for line in sys.stdin:
    n += 1
    if n % 4 == 2:
        seq = line.strip().upper()
        leng += len(seq)
        gc = seq.count("C") + seq.count("G")
        at = seq.count("A") + seq.count("T")
        GC_count += gc
        AT_count += at

print("Total number:", int(n/4))
print("Total length (Gb):", round(leng / 1000000000.0, 2))
print("GC (%):", round(GC_count * 1.0 / (GC_count + AT_count) * 100, 2))
