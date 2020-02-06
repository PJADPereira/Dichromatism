import sys
conditions = sys.argv[1]
count_matrix = sys.argv[2]
meta_data = sys.argv[3]
output_prefix = sys.argv[4]

#conditions = 'Fchesthousefinch;Mchesthousefinch'
#count_matrix = 'D:/FCUP_Drive/DE/Files/renamed_csv/housefinch_3vs3.isoform.counts.matrix.txt'
#meta_data = 'D:/FCUP_Drive/DE/Files/renamed_csv/housefinch.isoform.meta.txt'
#output_prefix = 'D:/FCUP_Drive/DE/Files/renamed_csv/F_chest_HF_vs_M_chest_HF'

conditions = conditions.split(";")
conditions
data = []
with open(meta_data) as inp:
    for i, line in enumerate(inp):
        if i == 0:
            m_header = line.strip().split(",")
        else:
            data.append(line.strip().split(","))

columns_names = []
with open(output_prefix+"_meta.csv", 'w') as out:
    out.write('{}\n'.format(",".join(m_header)))
    for i in data:
        if "".join(i[1:]) in conditions:
            columns_names.append(i[0])
            out.write('{}\n'.format(",".join(i)))

with open(count_matrix) as inp, open(output_prefix+"_matrix.csv", 'w') as out:
    for i, line in enumerate(inp):
        if i == 0:
            header = line.strip().split(",")

            columns_to_get = [0]+[header.index(x) for x in columns_names]

            out.write('{}\n'.format(
                ",".join([x for i, x in enumerate(header) if i in columns_to_get])))
        else:
            m_data = line.strip().split(',')
            to_write = [x for i, x in enumerate(m_data) if i in columns_to_get]
            out.write('{}\n'.format(",".join(to_write)))
