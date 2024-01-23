import sys

n = 0
num = 0
outdir = sys.argv[2]
with open(sys.argv[1], 'r') as f:
    for line in f:
        if n % 1000 == 0:
            num += 1
            print outdir + "/chunk" + str(num)
        n += 1
        with open(outdir + "/chunk" + str(num), 'a') as fh:
            fh.write(line)
