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

## DISCLAIMER
THIS WEBSITE AND CONTENT AND ALL SITE-RELATED SERVICES, INCLUDING ANY DATA, ARE PROVIDED "AS IS," WITH ALL FAULTS, WITH NO REPRESENTATIONS OR WARRANTIES OF ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTIES OF MERCHANTABILITY, SATISFACTORY QUALITY, NON-INFRINGEMENT OR FITNESS FOR A PARTICULAR PURPOSE. YOU ASSUME TOTAL RESPONSIBILITY AND RISK FOR YOUR USE OF THIS SITE, ALL SITE-RELATED SERVICES, AND ANY THIRD PARTY WEBSITES OR APPLICATIONS. NO ORAL OR WRITTEN INFORMATION OR ADVICE SHALL CREATE A WARRANTY OF ANY KIND. ANY REFERENCES TO SPECIFIC PRODUCTS OR SERVICES ON THE WEBSITES DO NOT CONSTITUTE OR IMPLY A RECOMMENDATION OR ENDORSEMENT BY PACIFIC BIOSCIENCES.
