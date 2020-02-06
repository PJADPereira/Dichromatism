[//]: # ( )

# Whole-genome DNA pools

This document has two main sections:

[Whole-genome sequencing of DNA pools](#Whole-genome-sequencing-of-DNA-pools)

[Genetic mapping using differentiation association and introgression statistics](#genetic-mapping-using-differentiation-association-and-introgression-statistics)

&nbsp;

&nbsp;

# Whole-genome sequencing of DNA pools

Software:

[trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)

[BWA MEM](https://sourceforge.net/projects/bio-bwa/files/)

[Samtools](https://sourceforge.net/projects/samtools/files/)

---

Trimming of adapter sequences:

```shell
trimmomatic-0.39.jar PE input_forward.fq.gz input_reverse.fq.gz R1.fastq.gz R1_unpaired.fastq.gz R2.fastq.gz R2_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:15 SLIDINGWINDOW=4:20 MINLEN:30
```

* PE: Pair-end data, followed by the R1 and R2 fastq input files and R1 paired, R1 unpaired, R2 paired, R2 unpaired output files.
* ILLUMINACLIP: Remove the indicated adapters, arguments are passed in the following order:
    * fastaWithAdapters: Fasta file containing all the adapters, etc.
    * seedMismatches: maximum mismatch count allowed.
    * palindromeClipThreshold: how accurate match's between two 'adapter ligated' reads must be for PE palindrome read alignment. 
    * simpleClipThreshold: how accurate any adapter etc. sequence must be against read.
    * minAdapterLength: Minimum length of adapter to be removed.
    * keepBothReads: Keep forward and reverse reads.
* Leading: Specifies the minimum quality required to keep a base at the beginning of the sequence.
* Trailing: Specifies the minimum quality required to keep a base at the end of the sequence.
* SLIDINGWINDOW: Perform a sliding window trimming. Parameter with two arguments:
    * windowSize: Number of bases to average quality across.
    * requiredQuality: Required quality to pass.
* MINLEN: Minimum length of reads, after trimming, required to be kept.

&nbsp;

&nbsp;

Map all reads in the read-pairs to a reference fasta file:

```shell
bwa mem -t 15 -M reference_genome.fasta R1.fastq.gz R2.fastq.gz > mapped_to_reference.sam
```

* -t: number of threads for bwa mem.
* -M: Mark shorter split hits as secondary.

&nbsp;

&nbsp;

Sort and convert to sam to bam using samtools:

```shell
samtools sort -@10 -O BAM -o output.bam
```

* -@: Set number of sorting and compression threads.
* -o: Write the final sorted output to FILE.
* -O: Write the final output as sam, bam, or cram.

---

&nbsp;

&nbsp;

# Genetic mapping using differentiation association and introgression statistics

Software:

[Samtools](https://sourceforge.net/projects/samtools/files/)

[Popoolation2](https://sourceforge.net/p/popoolation2/wiki/Main/)

[FreeBayes](https://github.com/ekg/freebayes)

[VCFlib](https://github.com/vcflib/vcflib)

---
&nbsp;

&nbsp;

Create mpileup for all pools here sequenced using samtools mpileup:

```shell
samtools mpileup -B -q 20 -Q 20 pool_1.bam pool_2.bam pool_n.bam > all_pools.mpileup

```

* -B: Disable base alignment quality (BAQ) computation.
* -q: Minimum mapping quality for an alignment to be used.
* -Q: Minimum base quality for a base to be considered.

&nbsp;

&nbsp;

Convert mpileup to sync:

```shell
mpileup2sync.jar --input all_pools.mpileup --output all_pools.sync --fastq-type sanger --min-qual 20 --threads 20

```

* --input: input sync file.
* --output: output file of the analysis.
* --fastq-type: determines the type of phred score used illumina = phred 64, sanger = phred 33.
* --min-qual: minimum quality accepted for each base called.
* --threads: number of threads to parallelize the script.

&nbsp;

&nbsp;

Run popoolation2 fst-sliding script:

```shell
perl fst-sliding.pl --input all_pools.sync --output all_pools_w20kb_10kb.fst --suppress-noninformative --min-count 4 --min-coverage 30 --max-coverage 43,47,69,55,51,61,53,60 --min-covered-fraction 0.3 --window-size 20000 --step-size 10000 --pool-size 16:12:12:12:16:39:12:12
```

* --input: input sync file.
* --output: output file of the analysis.
* --suppress-noninformative: suppresses output for windows not containing a SNP. When applied to windows of size one, this option suppresses output for bases that are no SNP.
* --min-count: minimum allele count.
* --min-coverage: minimum coverage accepted at any given position.
* --max-coverage: maximum coverage accepted at any given position.
* --min-covered-fraction: minimum proportion of the window sufficiently covered by the pileup.
* --window-size: size of the window that will scan the genome.
* --step-size: distance to shift the window by each step.
* --pool-size: number of individuals in the pool.
* --fastq-type: determines the type of phred score used illumina = phred 64, sanger = phred 33.
* --min-qual: minimum quality accepted for each base called.
* --threads: number of threads to parallelize the script.

&nbsp;

&nbsp;

Perform a Cochran-Mantel-Haenszel test for repeated tests of independence for every SNP using Popoolation2:

```shell
perl cmh-test.pl --input all_pools.sync --output contrast_1.cmh --min-count 4 --min-coverage 10 --max-coverage 20 --population 1-3,2-4
```

* --input: input sync file.
* --output: output file of the analysis.
* --min-count: minimum allele count.
* --max-coverage: maximum coverage accepted at any given position.
* --population: populations in the sync file to include in the analysis. E.g.: compare populations 1 and 3 contrasted against 3 and 4.

&nbsp;

&nbsp;

SNP variant calling through freebayes:

```shell
./freebayes -F 0.1 -C 3 -H -J \
            --min-mapping-quality 20 \
            --min-base-quality 20 \
            --min-coverage 10 \
            --dont-left-align-indels \
            -f reference_genome.fasta \
            -L bam_list.txt \
            -A file containing ploidy information
            -v variant_call.vcf
 ```

* -F: Require at least this fraction of observations supporting an alternate allele within a single individual in the in order to evaluate the position.
* -C: Require at least this count of observations supporting an alternate allele within a single individual in order to evaluate the position.
* -H: Use a weighted sum of base qualities around an indel, scaled by the distance from the indel.  By default use a minimum BQ in flanking sequence.
* -j:  Use mapping quality of alleles when calculating data likelihoods.
* --min-mapping-quality: Exclude alignments from analysis if they have a mapping quality less than value.
* --min-base-quality: Exclude alleles from analysis if their supporting base quality is less than value.
* --min-coverage: Require at least this coverage to process a site.
* --dont-left-align-indels: Turn off left-alignment of indels.
* -f: Reference sequence.
* -L:  List of bam files to be analysed.
* -A: Read a copy number map from the BED file.
* -v:  Output vcf file.

&nbsp;

&nbsp;

Filter vcf file with vcffilter from vcflib:

```shell
vcffilter -v -f "DP < 33 | DP > 439 | QUAL < 50" variant_call.vcf
```

* -v: inverts the filter
* -f: filter using this parameters:

    * DP: depth / coverage

    * QUAL: Base Quality
