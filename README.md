
# CountMappedReads2 - a tool for NGS read mapping, count reads and calculate normalized abundance

[![license](https://img.shields.io/github/license/mashape/apistatus.svg)](/LICENSE)
[![size](https://img.shields.io/github/size/webcaetano/craft/build/phaser-craft.min.js.svg)]()
[![doi](https://img.shields.io/badge/doi-10.1093%2Fbioinformatics%2Fbtx157-blue.svg?style=flat)](https://doi.org/10.1093/bioinformatics/btx157)

Countmappedreads2 serializes multiple processes including NGS read maping (by bowtie2), conversion to bam (by samtools), count mapped reads (by featureCounts), and coverage/FPKM(RPKM)/TPM calculation (by in-house script).

## requirements
* bowtie2
* samtools
* subread (for featureCounts)
* Ruby (ver >=2.0)

## usage 
```
### CountMappedReads2 ver 2.0.0 (2018-05-02) ###

[description]
CountMappedReads2 - a tool for NGS read mapping, count reads and calculate normalized abundance
CountMappedReads2 serializes multiple processes including NGS read maping (by bowtie2), conversion to bam (by samtools), count mapped reads (by featureCounts),
and coverage/FPKM(RPKM)/TPM calculation (by in-house script).

[usage]
$ CountMappedReads2 <reference fasta> <output dir> {-1 <pe1> -2 <pe2> | -U <up>} [options]

[arguments]
    - reference fasta      -- nucleotide fasta file for read mapping on (e.g. genomes, contigs)
    - output dir           -- output directory (should not exist).

[dependencies]
    - bowtie2
    - samtools
    - subread (for featureCounts)
    - ruby (ver >=2.0)

[options]
  (general)
    -h, --help
    -v, --version
    --gene      [str] gene annotation file -- either of GTF (GTF2), GFF3, or BED (BED6) formatted file.
    --format    [str] to specify format of gene annotation file, either of 'GTF', 'GTF2', 'GFF3', 'BED', or 'BED6'. If not given, format is inferred from extention of gene annotation file.

  (bowtie2)
    -1          [str] (required)      -- bowtie2 '-1' option (for paired-end read 1)
    -2          [str] (required)      -- bowtie2 '-2' option (for paired-end read 2)
    -U          [str] (required)      -- bowtie2 '-U' option (for unpaired reads)
    --score-min [str]                 -- bowtie2 '--score-min' option under '--end-to-end' mapping (default: 'L,0,-0.6')

  (featureCounts)
    --feature   [str] (default: gene)         -- specify feature name for count (some type written in 3rd column in GTF/GFF format). If '--format' is BED(6), all lines are included for count, regardless of '--feature'.
    --attr      [str] (default: <feature>_id) -- speficy attribute to count by (written in 9th column in GTF/GFF format). If '--format' is BED(6), 4th column is used to count by, regardless of '--attr'.

  (computation)
    --threads   [int] (default: 1)    -- number of threads for run of bowtie2, samtools, and featureCounts
    --mem       [str] (default: 1G)   -- maximum memory per thread for 'samtools sort'; suffix K/M/G is recognized (e.g., 800M)

[output files]
count/sequence.result.tsv                                -- tab separated file reporting  sequence abundance (BED6 and 7: count, 8: coverage, 9: FPKM(RPKM), 10: TPM)
count/<feature>.result.tsv (only when '--gene' is given) -- tab separated file reporting <feature> abundance (BED6 and 7: count, 8: coverage, 9: FPKM(RPKM), 10: TPM)

[note]
CountMappedReads2 does not perform quality control of reads. It is recommended to do QC and use QCed read as input.
If '--gene' option is given, CountMappedReads2 returns count and normalized abundance both for input sequences (e.g., genomes, contigs) and for elements (e.g., genes)
In this tool, FPKM(RPKM) and TPM are defined as described in http://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/

The '--score-min' option of bowtie2 controls a minimum acceptable level of mapping quality.
For example, with the 'end-to-end' mode, '--score-min L,0,-0.6' means an aligned region shows at least 90% identity, assuming that a read quality is high enough and a gap is not open.
Similarity, 'L,0,-0.3' means at least 95% identity, 'L,0,-0.9' means at least 85% identity, and 'L,0,-1.2' means at least 80% identity under the same assumption.
For deteils of the '--score-min' option and the alignment score, see bowtie2 manual page (http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml).
```

