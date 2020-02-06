[//]: # ( )

# Differential expression

This document has four main sections:

[RNA-SEQ Quality control](#RNA-SEQ-Quality-control)

[Merging species comparisons](#Merging-species-comparisons)

[Merge DE transcripts inside each contrast](#Merge-DE-transcripts-inside-each-contrast)

[Overlap between contrasts](#Overlap-between-contrasts)

---

&nbsp;

&nbsp;

# RNA-SEQ Quality control

Software:

[Trim Galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)

---

Trim for adapters and low quality bases:

```shell
trim_galore --path_to_cutadapt=*/bin/cutadapt \
            --paired \
            --retain_unpaired \
            --phred33 \
            --output_dir trimmed \
            --length 36 \
            -q 5 \
            --stringency 1 \
            -e 0.1 \
            sampleforde_1.fq \
            sampleforde_2.fq
```

* --path_to_cutadapt: if cutadapt is not in $PATHm, specify the path to its executable.
* --paired: Allows trim galore to work on paired-end files. This means that to pass the validation test, both sequences of a sequence pair are required to have a minimum length specified in --length.
* --retain_unpaired: If only one of the two paired-end reads became too short, the longer read will be written will be retained and written to a unpaired_R1 or unpaired_R2 file.
* --phred33: Use ASCII+33 quality scores as Phred scores for quality trimming.
* --output_dir: All output will be written to the specified folder.
* --length: Discard reads that became shorter than specified length because of either quality or adapter trimming.
* -q: Trim low-quality ends from reads below threshold.
* --stringency: Overlap with adapter sequence required to trim a sequence.
* -e: Maximum allowed error rate.

---

&nbsp;

&nbsp;

# Trinity Count Tables

Software

[Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki)

[DEApp](https://yanli.shinyapps.io/DEApp/)

---

Align reads to reference transcriptome:

```shell
TRINITY_HOME/util/align_and_estimate_abundance.pl \
    --transcript species_specific_transcriptome.fa  \
    --samples_file inputs_table \
    --seqType fq \
    --est_method RSEM \
    --aln_method bowtie2 \
    --SS_lib_type FR \
    --thread_count 40 \
    --prep_reference \
    --trinity_mode \
    --output_dir rsem_outdir

```

* --transcript: reference transcriptome (built in transcriptome.md).
* --samples_file: A tab separated table with the following format:

    |  Condition       |   Repetition            |  R1_reads        |     R2_reads     |
    |:-------:|:-------------:|:----------------:|:----------------:|
    |cond_A   | cond_A_rep1   | A_rep1_left.fq   | A_rep1_right.fq  |
    |cond_A   | cond_A_rep2   | A_rep2_left.fq   | A_rep2_right.fq  |
    |cond_B   | cond_B_rep1   | B_rep1_left.fq   | B_rep2_right.fq  |
    |cond_B   | cond_b_rep2   | B_rep2_left.fq   | B_rep2_right.fq  |
    |         |               |                  |                  |
* --est_method: Abundance estimation method.
* --aln_method: Aligner to use.
* --SS_lib_type: Strand specific library type (R1 - F: R2 - R).
* --thread_count: Number of threads for parallelization.
* --prep_reference: Build index for reference.
* --trinity_mode: Setting --trinity_mode will automatically generate the a gene name to transcript association.
* --output_dir: write all files to output directory.

&nbsp;

&nbsp;

Convert abundance estimate into a count matrix:

```shell
TRINITY_HOME/util/abundance_estimates_to_matrix.pl \
    --est_method RSEM \
    species_specific_transcriptome.fa.gene_trans_map \
    --out_prefix RSEM \
    --name_sample_by_basedir \
     cond_A_rep1/isoforms.results \
     cond_A_rep2/isoforms.results \
     cond_B_rep1/isoforms.results \
     cond_B_rep2/isoforms.results
```

&nbsp;

&nbsp;

At this point you will have the following type of transcript count table:

|ID|cond_A_rep1|cond_A_rep2|cond_A_rep3|cond_B_rep1|cond_B_rep2|cond_B_rep3|cond_C_rep1|cond_C_rep2|cond_C_rep3|cond_D_rep1|cond_D_rep2|cond_D_rep3|
|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
|ENSTGUP00000012686.1|2222|1183|1860|1944|2435|2193|1152|1196|753|1597|1361|1356|
|ENSTGUP00000012628.1|114|48|131|104|127|154|34|56|170|118|85|365|
|ENSTGUP00000006842.1|142|29|169|153|177|169|33|61|1|74|79|20|
|ENSTGUP00000010084.1|5|5|5|7|2|6|2|3|11|8|8|13|
|ENSTGUP00000004756.1|35|15|40|50|80|73|31|33|58|47|88|96|
|ENSTGUP00000011676.1|54|10|56|29|79|86|27|40|1|38|33|12|
| ...|... |...  |...  |...  |...  |...  |...  |...  |...  |...  |...  | ...|

&nbsp;

For which a metadata table can be constructed, e.g:

|Sample|Sex|Patch|Species|
|:---:|:---:|:---:|:---:|
|cond_A_rep1|M|chest|housefinch|
|cond_A_rep2|M|chest|housefinch|
|cond_A_rep3|M|chest|housefinch|
|cond_B_rep1|F|chest|housefinch|
|cond_B_rep2|F|chest|housefinch|
|cond_B_rep3|F|chest|housefinch|
| cond_C_rep1 | M|belly|housefinch|
|...|...|...|...|
| | | | |

&nbsp;

In order to subset counts table to each individual comparison the script [make_inputs.py](./src/make_inputs.py) can be run in the following way:

```shell
python3 make_inputs.py 'Treatment1;Treatment2' path/to/count_table.txt path/to/metadata_table.txt path/to/output
```

* 'Treatment1;Treatment2': String that matches, without any whitespace, the attributes of the samples from the meta-table.
* path/to/count_table.txt: path to the trinity output count table.
* path/to/metadata_table.txt: path to the metadata table.
* path/to/output: path to output files, files will take as prefix the folder name selected.

For example, the following call will pick up all replicates of cond_A and all replicates of cond_B from the count_matrix table and will also subset its metadata and output both files in a folder named MCHFvsFCHF.

```shell
python3 make_inputs.py 'Mchesthousefinch;Fchesthousefinch' counts.matrix.txt meta_data.txt contrast/MCHFvsFCHF
```

Those files are ready to load into DEApp.

&nbsp;

&nbsp;

# Merging species comparisons

To account for potential biases when mapping to different species two tests were performed. First, reads from both species A and B were mapped to species A transcriptome. Second, reads from both species A and B were mapped to species B transcriptome.

In such cases, only transcripts differentially expressed in both tests were considered. Such lists can be obtained by running the [merge_different_species.py](./src/merge_different_species.py) script found in this repository:

```shell
python3 merge_different_species.py SAAvsSBA/all3_overlap.txt SABvsSBB/all3_overlap.txt output/all3overlap_common.txt
```

* SAAvsSBA corresponds to the comparison where both species were mapped to species A transcriptome.
* SABvsSBB corresponds to the comparison where both species were mapped to species B transcriptome.

---

&nbsp;

&nbsp;

# Merge DE transcripts inside each contrast

The total list of differentially expressed transcripts was then compiled for each contrast (Patch, Sex, Species) using the script [merge_tests.py](./src/merge_tests.py):

```shell
python3 merge_tests.py /contrast/test_1/all3_overlap.txt /contrast/test_2/all3_overlap.txt /contrast/test_n/all3_overlap.txt /contrast/DE_transcripts.txt
```

---

&nbsp;

&nbsp;

# Overlap between contrasts

To obtain the overlapping transcripts between each contrast the script [common_between_tests.py](./src/common_betweem_tests.py) can be used:

```shell
python3 common_between_tests.py contrast_1/DE_transcripts.txt contrast_2/DE_transcripts.txt contrast_n/DE_transcripts.txt output/common_between_all.txt
```

Another way to obtain the overlap betweens tests is to draw a Venn Diagram, to do so the [venn.py](./src/venn.py) script can be used. This script requires 3 contrasts.

Note that too run this script two python packages need to be installed, matplotlib and matplotlib_venn. Both can be installed through PIP:

```shell
pip3 install matplotlib
pip3 install matplotlib_venn
```


```shell
python3 venn.py contrast_1/DE_transcripts.txt contrast_2/DE_transcripts.txt contrast_3/DE_transcripts.txt label1,label2,label3 output/venn.pdf
```

* contrast_1/DE_transcripts.txt: replace with path to merged list of differentially expressed transcripts in contrast 1.
* contrast_2/DE_transcripts.txt: replace with path to merged list of differentially expressed transcripts in contrast 2.
* contrast_3/DE_transcripts.txt: replace with path to merged list of differentially expressed transcripts in contrast 3.
* label1,label2,label3: replace with labels to be used in the venn diagram, the order match's the order of the previous arguments.
* output/venn.pdf: replace with path to output pdf.
