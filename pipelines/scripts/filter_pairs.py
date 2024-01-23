import sys
from collections import OrderedDict
def main(inPairs=sys.stdin, outPairs=sys.stdout):
	d_pairs = OrderedDict()
	for line in inPairs:
		temp = line.rstrip().split()
		pair = tuple(temp[:2])
		if pair in d_pairs or pair[-1::-1] in d_pairs:
			continue
		else:
			d_pairs[pair] = line
	for line in d_pairs.values():
		outPairs.write(line)

if __name__ == '__main__':
	main()
