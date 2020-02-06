import sys
from collections import defaultdict


first_file = sys.argv[1]
second_file = sys.argv[2]
output_file = sys.argv[3]


first_file_data = defaultdict(list)
second_file_data = defaultdict(list)

with open(first_file) as inp:
    line_count = 0
    for line in inp:
        if line_count == 0:
            header = line
            line_count += 1
            continue

        data = line.strip().split('\t')
        first_file_data[data[0]] = data[1:]

with open(second_file) as inp:
    line_count = 0
    for line in inp:
        if line_count == 0:
            #header = line.strip()
            line_count += 1
            continue

        data = line.strip().split('\t')
        second_file_data[data[0]] = data[1:]

with open(output_file, 'w') as out:
    out.write(header)
    for k, v in second_file_data.items():
        if len(first_file_data[k]) > 0:
            out.write('{}\t{}\n'.format(k, "\t".join(v)))
