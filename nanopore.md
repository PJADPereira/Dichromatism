[//]: # ( )

# Nanopore data analysis

This document has two main sections:

[Trimming of adapter sequences](#Trimming-of-adapter-sequences)

[Mapping to reference genome](#Mapping-to-reference-genome)

---

&nbsp;

&nbsp;

# Trimming of adapter sequences

Software

[Downpore](https://github.com/jteutenberg/downpore)

---
From the basecalled reads, trim adapters using downpore:

```shell
downpore trim -i raw_nanopore_reads.fastq -f ./data/adapters_back.fasta -b ././data/adapters_back.fasta > trimmed_nanopore_reads.fastq
```

* -i: a fastq file containing all reads from a nanopore run.
* -f front-adapters used, adapters are provided as a fasta file in/[downpore_installation_folder]/data/adapters_front.fasta.
* -b back-adapters used, adapters are provided as a fasta file in/[downpore_installation_folder]/data/adapters_back.fasta.

---

&nbsp;

&nbsp;

# Mapping to reference genome

Software:

[NGMLR](https://github.com/philres/ngmlr)

[samtools](https://sourceforge.net/projects/samtools/files/)

---

NGMLR was used to map long reads to the reference genome:

```shell
ngmlr -t 8 -r reference_genome.fasta -q trimmed_nanopore_reads.fastq -o mapped.sam
```

* -t: number of threads working.
* -r: reference genome to map to.
* -q: input reads in fastq format.
* -o: output sam.

Convert sam to bam:

```shell
samtools view -S -b mapped.sam > mapped.bam
```

* -S: input is in sam format.
* -b: output in bam format.

Sort the resulting bam file:

```shell
samtools sort mapped.bam > mapped.sorted.bam
```

Index file bam file:

```shell
samtools index mapped.sorted.bam
```

Subset bam file to capture only the candidate region:

```shell
samtools view -b mapped.sorted.bam "NW_007931177:815000-865000" > candidate.bam
```

Index the new file:

```shell
samtools index candidate.bam.bai
```

---
