#!/bin/bash

module load bioinfo-tools
module load samtools

samtools view -bh \
    /proj/a2009002/webexport/opendata/HiSeqX_CEPH/CEP-1-7/03-BAM/CEP-1-7.clean.dedup.recal.bam \
    20 > CEP-1-7.clean.dedup.recal.chr20.bam
