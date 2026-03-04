import sys

fn = sys.argv[1]

names, seqs, s = [], [], []
with open(fn) as f:
	for line in f:
		if line.startswith(">"):
			names.append(line[1:].strip())
			if s:
				seqs.append("".join(s))
			s = []
			continue
		s.append(line.strip())
seqs.append("".join(s))

aln = list(zip(seqs[0], seqs[1]))

norm = len([x for x in aln if '-' not in x and 'X' not in x])
if norm == 0:
	seqid = 0
else:
	seqid = len([x for x in aln if x[0] == x[1] and x[0] != "X"])/norm

print("{0}\t{1}\t{2:10.6f}".format(names[0], names[1], seqid))
