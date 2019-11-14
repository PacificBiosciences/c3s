<h1 align="center">c3s</h1>
<p align="center">Consensus of CCS reads</p>

***

Generates one consensus sequence of all input CCS reads,
using the Partial Order Alignment implementation [spoa](https://github.com/rvaser/spoa).

## Formats
Allowed input formats are
 - PacBio CCS `.bam`
 - PacBio dataset `.xml`
 - FASTQ (plain) `.fastq` or `.fq`
 - FASTQ (gzip or bgzip) `.fastq.gz` or `.fq.gz`

 Allowed output formats are
 - FASTA (plain) `.fasta` or `.fa`
 - FASTA (bgzip) `.fasta.gz` or `.fa.gz`


## How to use
Example invocation

```
c3s input.ccs.bam out.fa
```
