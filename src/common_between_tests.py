import sys
from collections import defaultdict

comps = sys.argv[1:-1]
output = sys.argv[-1]

comp = defaultdict(lambda: [0 for x in range(len(comps))])

for i, file in enumerate(comps):
    with open(file) as inp:
        for line in inp:
            transcript = line.strip()
            comp[transcript][i] = 1


with open(output, 'w') as out:
    for k, v in comp.items():
        if sum(v) == len(comps):
            out.write(k+"\n")

