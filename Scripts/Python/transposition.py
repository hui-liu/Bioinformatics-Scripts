import sys
"""python transposition.py input output
Row to Column Transposition (转置)"""

in_file = open(sys.argv[1])
out_file = open(sys.argv[2], 'w')

lis = [x.split() for x in in_file]

for x in zip(*lis):
    for y in x:
        print >> out_file, y + '\t',
    print >> out_file # add '\n' throuth 'print'
out_file.close()
in_file.close()
