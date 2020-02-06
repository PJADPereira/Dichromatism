import sys
import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn3_unweighted


contrast_1 = sys.argv[1]
contrast_2 = sys.argv[2]
contrast_3 = sys.argv[3]
labels = tuple(sys.argv[4].split(','))
output = sys.argv[5]

t_contrast_1 = []

with open(contrast_1) as inp:
    for line in inp:
        t_contrast_1.append(line.strip())

t_contrast_2 = []

with open(contrast_2) as inp:
    for line in inp:
        t_contrast_2.append(line.strip())

t_contrast_3 = []

with open(contrast_3) as inp:
    for line in inp:
        t_contrast_3.append(line.strip())


# perform sanity check

assert len(t_contrast_2) == len(list(set(t_contrast_2)))
assert len(t_contrast_1) == len(list(set(t_contrast_1)))
assert len(t_contrast_3) == len(list(set(t_contrast_3)))

# if sanity checks pass convert lists into sets (if they fail there are duplicate transcripts in one of the files which is not expected)

t_contrast_2 = set(t_contrast_2)
t_contrast_1 = set(t_contrast_1)
t_contrast_3 = set(t_contrast_3)

# Plot

venn3_unweighted(subsets=(t_contrast_1, t_contrast_2, t_contrast_3), set_labels=labels, alpha=0.55)
plt.savefig(output, transparent=True, format='pdf')

