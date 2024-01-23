
'''
    from Trf import TrfReader
    inTrf = 'trf.out.dat'
    trf_reader = TrfReader(inTrf)
    for sequence in trf_reader.parse():
        for record in sequence.records:
            print sequence.id, record.start, record.period_size, record.copy_number
'''
import sys
import re

class TrfRecord():
	def __init__(self, line):
		self.line = line.split()
		self.start = int(self.line[0])
		self.end = int(self.line[1])
		self.period_size = int(self.line[2])
		self.copy_number = float(self.line[3])
		self.consensus_size = int(self.line[4])
		self.percent_matches = int(self.line[5])
		self.percent_indels = int(self.line[6])
		self.score = int(self.line[7])
		self.A = int(self.line[8])
		self.C = int(self.line[9])
		self.G = int(self.line[10])
		self.T = int(self.line[11])
		self.entropy = float(self.line[12])
		self.motif = self.line[13]
		self.seq = self.line[14]
		self.center = 1.0*(self.start + self.end) / 2
		self.range = (self.start, self.end)
		self.tr_length = self.end - self.start + 1
		self.motif_length = self.period_size

class TrfReader():
	def __init__(self, trf_out):
		if isinstance(trf_out, str):
			self.trf = open(trf_out)
		elif isinstance(trf_out, file):
			self.trf = trf_out
		else:
			raise ValueError('unsupport type %s of %s' % (type(trf_out), trf_out))
	@property
	def records(self):
		for line in self.lines:
			yield TrfRecord(line)
	def parse(self):
		for line in self.trf:
			if line.startswith('Sequence'):
				if 'lines' in dir():
					self.id = sequence_id
					self.lines = lines
					yield self
					lines = []
				else:
					lines = []
				sequence_id = line.rstrip().split()[1]
			elif re.compile(r'\d').match(line):
				line = line.rstrip()
				lines.append(line)
			elif line.startswith('Parameters'):
				self.parameters = line.split(':')[1].strip()
			else:
				continue
		self.id = sequence_id
		self.lines = lines
		yield self

def main(inTrf):
	trf_reader = TrfReader(inTrf)
	for recs in trf_reader.parse():
		#print recs, dir(recs)
		for record in recs.records:
			print recs.id, record.start, record.end, record.start/100000 * 100000, record.period_size, record.consensus_size, record.copy_number, record.tr_length, record.motif

if __name__ == '__main__':
	main(sys.argv[1])


