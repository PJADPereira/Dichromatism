[//]: # ( )

# Transcriptomics
This document has two main sections:

[RNA-SEQ Quality control](#RNA-SEQ-Quality-control)

[*De novo* transcriptome assemblies](#de-novo-transcriptome-assemblies)

&nbsp;

&nbsp;

# RNA-SEQ Quality control

Software:

[fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

[Rcorrector](https://github.com/mourisl/Rcorrector)

[TranscriptomeAssemblyTools](https://github.com/harvardinformatics/TranscriptomeAssemblyTools)

[Trim Galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)

[cutadapt](https://github.com/marcelm/cutadapt/)

[bowtie2](https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.5.1/)

---

Check quality parameters in fastqc:

```shell
fastqc reads_1.fq.gz reads_2.fq.gz

```

&nbsp;

&nbsp;

Correct RNA-seq reads using a kmer-based error correction method implemented in Rcorrector; For paired reads, two files will be output, with the same name as input but with an added <input_prefix>.cor.<input_suffix>:

```shell
perl run_rcorrector.pl -t 30 -1 reads_1.fq.gz -2 reads_2.fq.gz
```

* -t: threads.
* -1: R1 fastq file.
* -2: R2 fastq file.

&nbsp;

&nbsp;

Remove uncorrectable reads and the corrected tag left from Rcorrector:

```shell
python TranscriptomeAssemblyTools/FilterUncorrectabledPEfastq.py \
        -1 reads_1.cor.fq.gz \
        -2 reads_2.cor.fq.gz \
        -s sample_name
```

* -1: R1 fastq file.
* -2: R2 fastq file.
* -s: sample name for log file.

&nbsp;

&nbsp;

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
            unfix_rm_reads_1.cor.fq \
            unfix_rm_reads_2.cor.fq
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

&nbsp;

&nbsp;

Remove rRNA contaminants, given that we are aligning reads to rRNA, the expected output comes from the paired but unaligned files as no read in the pair aligned to rRNA:

```shell
bowtie2 --nofw --quiet --very-sensitive-local \
        --phred33 \
        -x indexed_rRNA \
        -1 /trimmed/unfix_rm_reads_1.cor.fq.gz \
        -2 /trimmed/unfix_rm_reads_2.cor.fq.gz \
        --threads 40 \
        --met-file bowtie2_metrics.txt \
        --al-conc-gz paired_aligned.fq.gz \
        --un-conc-gz paired_unaligned.fq.gz \
        --al-gz unpaired_aligned.fq.gz \
        --un-gz unpaired_unaligned.fq.gz
```

* --nofw: Don't try to align unpaired reads to the forward reference strand.
* --quiet: Print nothing besides alignments and serious errors.
* --very-sensitive-local: A preset from bowtie2 (equal to using following parameters: -D 20 -R 3 -N 0 -L 20 -i S,1,0.5).
* --phred33: Input quality scores are in ASCII+33.
* -x: index files for the reference sequence.
* -1: R1 fastq file.
* -2: R2 fastq file.
* --threads: number of threads to parallelize the script.
* --met-file: Write bowtie2 metrics to the specified file.
* --al-conc-gz: Write the read pairs that aligned to the reference.
* --un-conc-gz: Write the read pairs that failed to align to the reference.
* --al-gz: Write unpaired reads that aligned to reference.
* --un-gz: Write unpaired reads that failed to align to reference.

&nbsp;

&nbsp;

Re-check quality parameters in fastqc: 

```shell
fastqc paired_unaligned.fq.1.gz paired_unaligned.fq.2.gz
```

---

&nbsp;

&nbsp;

# *De novo* transcriptome assemblies

Software:

[TransAbyss](https://github.com/bcgsc/transabyss)

[Oases](https://github.com/dzerbino/oases/)

[Spades](https://github.com/ablab/spades)

[Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki)

[CD-HIT](https://github.com/weizhongli/cdhit)

[EvidentialGene](http://arthropods.eugenes.org/EvidentialGene/about/EvidentialGene_trassembly_pipe.html)

[BUSCO](https://busco-archive.ezlab.org/v3/)

[blast+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)

[transPS](https://bioinformatics.cs.vt.edu/zhanglab/software/transps/)

---

Assembly with abyss (each kmer size tested was run independently):

```shell
transabyss -k 35 \
           --pe paired_unaligned.fq.1.gz paired_unaligned.fq.2.gz \
           --name kmer35.fa \
           --outdir /*pathto*/abyss_assembly \
           --SS \
           -c 12 \
           --threads 60
```

* -k: kmer size, in this study several kmer sizes were tested (35, 45, 55, 65, 75, 85).
* --pe: paired end files.
* --name: assembly name.
* --outdir: output directory.
* --SS: input reads are strand-specific.
* -c: minimum mean k-mer coverage of a unitig.
* --threads: number of threads for parallel computation.

&nbsp;

&nbsp;

Assembly with Oases:

```shell
python ./oases_pipeline.py \
        -m 33 \
        -M 69 \
        -s 6 \
        -o /*pathto*/oases_assembly \
        -d " -fastq.gz -shortPaired paired_unaligned.fq.1.gz paired_unaligned.fq.2.gz -strand_specific" \
        -p " -ins_length 350 "
```

* -m: minimum kmer size.
* -M: maximum kmer size.
* -s: step increase in kmer size.
* -o: output directory.
* -d: Velveth file descriptors:
  * -fastq.gz: input file format is fastq.gz.
  * -shortPaired: read type.
  * -strand_specific: input reads are strand-specific.
* -p: Oases options to be passed to the command line:
  * -ins_length: expected insert length.

&nbsp;

&nbsp;

Assembly with spades:

```shell
python ./spades.py \
        --rna \
        --ss-fr \
        -1 paired_unaligned.fq.1.gz \
        -2 paired_unaligned.fq.2.gz \
        -o /*pathto*/spades_assembly \
        -t 60
```

* --rna: RNA-seq assembly.
* --ss-fr: strand specific orientation, first read in pair corresponds to forward gene strand.
* -1: R1 fastq file.
* -2: R1 fastq file.
* -o: output folder for the assembly.
* -t: umber of threads for parallel computation.

&nbsp;

&nbsp;

Assembly with trinity:

```shell
Trinity --seqType fq \
        --SS_lib_type FR \
        --max_memory 225G \
        --min_kmer_cov 1 \
        --CPU 64 \
        --left paired_unaligned.fq.1.gz \
        --right paired_unaligned.fq.2.gz \
        --output trinity_assembly
```

* --seqType: Type of input reads.
* --SS_lib_type: Strand specific read orientation.
* --max_memory: suggested max random-access memory (RAM) for Trinity usage.
* --min_kmer_cov: minimum kmer coverage.
* --CPU: number of CPUs to use.
* --left: left reads, equivalent to -1 in other assemblers.
* --right: right reads,equivalent to -2 in other assemblers.
* --output: output folder for the assembly.

&nbsp;

&nbsp;

For each assembler cluster highly sequences (e.g. for oases):

```shell
cd-hit-est -o oases_assembly_0.99.fa -c 0.99 -i oases_assembly.fa -p 1 -d 40 -g 1 -T 10
```

* -o: output filename.
* -c: sequence identity threshold.
* -i: input assembly in fasta format.
* -p: print alignment overlap to .clstr file.
* -d: length of description in .clstr file.
* -g: if one, cluster sequence to the most similar cluster that meets the threshold.
* -T: number of threads for parallel computation.

&nbsp;

&nbsp;

Merge all assemblies after clustering using cat:

```shell
cat oases_assembly_0.99.fa spades_assembly_0.99.fa trinity_assembly_0.99.fa abyss_assembly_0.99.fa > all_assemblies.fasta
```

&nbsp;

&nbsp;

Rename transcripts:

```shell
perl -ane 'if(/\>/){$a++;print ">Locus_$a\n"}else{print;}' all_assemblies.fasta > all_renamed.fasta
```

&nbsp;

&nbsp;

Ensure unique IDs and add prefixes for parameter sets:

```shell
trformat.pl -output all_assemblies.tr -input all_renamed.fasta
```

&nbsp;

&nbsp;

Merge and reduce redundancy across assemblies:

```shell
perl tr2aacds.pl -mrnaseq all_assemblies.tr -NCPU=60 -MAXMEM=200000 1>tr2aacds.log 2>tr2aacds.err
```

* -mrnaseq: input file.
* -NCPU: number of cpus to use.
* -MAXMEM: max random-access memory (RAM) reserved for the program.
* 1> std out (redirect normal console output to file).
* 2> std err (redirect errors to file).

&nbsp;

&nbsp;

Quality assessment of the final transcriptome through BUSCO:

```shell
python run_BUSCO.py \
        -o canary \
        -i all_assemblies_okay.tr \
        -l aves_odb9/ \
        -m transcriptome -f
```

* -o: output folders and files will be labelled with this name.
* -i: Assembled transcriptome.
* -l: Specify location of the BUSCO [lineage](https://busco-archive.ezlab.org/v3/frame_meta.html) data to be used.
* -m: Type of analysis to run.
* -f: Force overwrite if any files are present in output folder.

&nbsp;

&nbsp;

Create a database of zebra finch transcripts:

```shell
makeblastdb -in zebra_finch_proteins.fasta \
        -dbtype prot \
        -out zb_database
```

* -in: input transcriptome.
* -dbtype: type of molecule.
* -out: output database to file.

&nbsp;

&nbsp;

Blast transcriptome to the previously created zebra finch database:

```shell
blastx -db zb_database \
        -query all_assemblies_okay.tr \
        -out blastx_out \
        -outfmt 6 \
        -evalue 0.01 \
        -max_target_seqs 20
```

* -db: database to query to.
* -query: list of transcripts to blast.
* -out: output file.
* -outfmt: format (where 6 is tabular).
* -evalue: [E-value](http://www.metagenomics.wiki/tools/blast/evalue) cutoff for BLAST searches.
* -max_target_seqs: maximum number of sequences.

&nbsp;

&nbsp;

Final filtering and scaffolding guided by the alignment to zebra finch:

```shell
perl transps.pl -t all_assemblies_okay.tr -b blaxtx_out.blastx
```

* -t: transcriptome to scaffold.
* -b: blastx output file.

---
