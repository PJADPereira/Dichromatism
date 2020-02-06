import sys

files_to_merge = sys.argv[1:-1]
output = sys.argv[-1]

transcripts = []
for file in files_to_merge:
    line_c = 0
    with open(file) as inp:
        for line in inp:
            if line_c == 0:
                line_c += 1
                continue
            else:
                transcripts.append(line.split('\t')[0].strip())

with open(output, 'w') as out:
    out.write("\n".join(list(set(transcripts))))
