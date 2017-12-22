import sys
import os
import pandas

# python merge.py dir output

name_list = os.listdir(sys.argv[1])
fram_list = [pandas.read_table("%s/%s" % (sys.argv[1], name)) for name in name_list]
fram = fram_list[0]
genename = fram.columns[0]

for i in range(1,len(fram_list)):
    fram = pandas.merge(fram,fram_list[i], on = genename)

fram.to_csv(sys.argv[2], sep = '\t' , index=False)
